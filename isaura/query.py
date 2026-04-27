import gc
import time

import numpy as np
import pandas as pd

from isaura.logging import logger
from isaura.utils import cpu_cnt, mem_gb_lim, rss_mb


def query(conn, header, wanted, file_glob, columns="*", tmpdir="/tmp", preserve_order=False):
  """Run a batched DuckDB query and return all results as a single DataFrame.

  Convenience wrapper around query_batched that concatenates all batches.

  Args:
      conn: DuckDB connection pre-configured with S3 credentials.
      header: Column name to filter on (e.g. "input").
      wanted: List of values to look up.
      file_glob: S3 URI glob or list of S3 URIs pointing to Parquet files.
      columns: SQL column expression (default "*" for all columns).
      tmpdir: Directory for DuckDB to spill temporary data to disk.
      preserve_order: If True, results are returned in the same order as wanted.

  Returns:
      DataFrame of all matching rows, or an empty DataFrame if none found.
  """
  chunks = list(
    query_batched(conn, header, wanted, file_glob, columns=columns, tmpdir=tmpdir, preserve_order=preserve_order)
  )
  if not chunks:
    return pd.DataFrame()
  return pd.concat(chunks, ignore_index=True)


def _configure_conn(conn, tmpdir):
  """Set DuckDB memory limit, temp directory, and thread count based on available resources."""
  mem_lim = mem_gb_lim()
  threads = cpu_cnt(ratio=0.8)
  try:
    conn.execute(f"SET memory_limit='{mem_lim}GB'")
    conn.execute(f"SET temp_directory='{tmpdir}'")
    conn.execute("PRAGMA enable_object_cache")
    conn.execute(f"SET threads TO {threads}")
  except Exception:
    pass
  return mem_lim, threads


def _src_expr(file_glob):
  """Build a DuckDB read_parquet() SQL expression from a glob string or list of URIs."""
  if isinstance(file_glob, list):
    escaped = ", ".join(f"'{u}'" for u in file_glob)
    return f"read_parquet([{escaped}])", f"{len(file_glob)} files"
  return f"read_parquet('{file_glob}')", f"glob={file_glob}"


def query_batched(
  conn, header, wanted, file_glob,
  batch_size=10000, columns="*", tmpdir="/tmp", preserve_order=False,
):
  """Query Parquet files on S3 via DuckDB and yield results as DataFrames in batches.

  Registers the wanted list as a temporary DuckDB table and runs a single SQL
  query with an IN filter. Memory limit and thread count are set automatically
  based on available system resources.

  Args:
      conn: DuckDB connection pre-configured with S3 credentials.
      header: Column name to filter on (e.g. "input").
      wanted: List of values to look up.
      file_glob: S3 URI glob or list of S3 URIs pointing to Parquet files.
      batch_size: Number of rows per yielded DataFrame.
      columns: SQL column expression (default "*" for all columns).
      tmpdir: Directory for DuckDB to spill temporary data to disk.
      preserve_order: If True, results are returned in the same order as wanted.

  Yields:
      DataFrames of matching rows in batch_size chunks.
  """
  if not wanted:
    return

  mem_lim, threads = _configure_conn(conn, tmpdir)
  wanted_list = list(wanted)
  src, src_desc = _src_expr(file_glob)

  logger.info(
    f"[query_batched] inputs={len(wanted_list)} "
    f"batch_size={batch_size} mem={mem_lim}GB threads={threads} src={src_desc}"
  )

  if preserve_order:
    sql = f"""
      WITH p AS (
        SELECT {columns} FROM {src}
        WHERE {header} IN (SELECT {header} FROM __wanted_inputs)
      )
      SELECT p.* FROM p
      JOIN __wanted_inputs w ON p.{header} = w.{header}
      ORDER BY w.__o
    """
    wdf = pd.DataFrame({
      header: wanted_list,
      "__o": np.arange(len(wanted_list), dtype=np.int64),
    })
  else:
    sql = f"""
      SELECT {columns} FROM {src}
      WHERE {header} IN (SELECT {header} FROM __wanted_inputs)
    """
    wdf = pd.DataFrame({header: wanted_list})

  conn.register("__wanted_inputs", wdf)
  total_rows = 0
  n_batches = 0
  t0 = time.perf_counter()

  try:
    for batch in conn.execute(sql).fetch_record_batch(batch_size):
      n_batches += 1
      df = batch.to_pandas(split_blocks=True, self_destruct=True)
      total_rows += len(df)
      logger.debug(f"[query_batched] batch#{n_batches} rows={len(df)} total={total_rows} rss={rss_mb():.0f}MB")
      yield df
  finally:
    conn.unregister("__wanted_inputs")

  dt = time.perf_counter() - t0
  logger.info(f"[query_batched] done rows={total_rows} batches={n_batches} elapsed={dt:.2f}s rss={rss_mb():.0f}MB")


def chunked_query_batched(
  conn, header, wanted, file_glob,
  slice_size=2000, batch_size=2000,
  columns="*", tmpdir="/tmp", preserve_order=False,
):
  """Query wide Parquet files by slicing the wanted list into smaller chunks.

  Used for models with 100+ output columns where loading all wanted inputs in
  a single DuckDB query would exhaust memory. Each slice of slice_size inputs
  is queried independently with reduced memory and thread limits, and results
  are yielded incrementally.

  Args:
      conn: DuckDB connection pre-configured with S3 credentials.
      header: Column name to filter on (e.g. "input").
      wanted: List of values to look up.
      file_glob: S3 URI glob or list of S3 URIs pointing to Parquet files.
      slice_size: Number of wanted inputs per DuckDB query.
      batch_size: Number of rows per yielded DataFrame within each slice.
      columns: SQL column expression (default "*" for all columns).
      tmpdir: Directory for DuckDB to spill temporary data to disk.
      preserve_order: If True, results within each slice are returned in wanted order.

  Yields:
      DataFrames of matching rows in batch_size chunks.
  """
  if not wanted:
    return

  wanted_list = list(wanted)
  n_total = len(wanted_list)
  n_slices = (n_total + slice_size - 1) // slice_size
  src, src_desc = _src_expr(file_glob)

  mem_lim = mem_gb_lim(ratio=0.4, floor_gb=1)
  wide_threads = max(1, min(4, cpu_cnt(ratio=0.3)))

  try:
    conn.execute(f"SET memory_limit='{mem_lim}GB'")
    conn.execute(f"SET temp_directory='{tmpdir}'")
    conn.execute("PRAGMA enable_object_cache")
    conn.execute(f"SET threads TO {wide_threads}")
    conn.execute("SET preserve_insertion_order=false")
  except Exception:
    pass

  logger.info(
    f"[chunked_query] inputs={n_total} slice={slice_size} "
    f"slices={n_slices} mem={mem_lim}GB threads={wide_threads} src={src_desc}"
  )

  sql = f"""
    SELECT {columns} FROM {src}
    WHERE {header} IN (SELECT {header} FROM __slice)
  """

  total_rows = 0
  t0 = time.perf_counter()

  for si in range(n_slices):
    lo = si * slice_size
    hi = min(lo + slice_size, n_total)
    chunk_wanted = wanted_list[lo:hi]
    wdf = pd.DataFrame({header: chunk_wanted})

    conn.register("__slice", wdf)
    slice_rows = {}
    try:
      for batch in conn.execute(sql).fetch_record_batch(batch_size):
        df = batch.to_pandas(split_blocks=True, self_destruct=True)
        if preserve_order:
          for _, row in df.iterrows():
            k = row.get(header)
            if k is not None:
              slice_rows[str(k).strip()] = row
        else:
          total_rows += len(df)
          yield df
        del df
    finally:
      conn.unregister("__slice")

    if preserve_order and slice_rows:
      ordered = []
      for k in chunk_wanted:
        row = slice_rows.get(str(k).strip())
        if row is not None:
          ordered.append(row)
      if ordered:
        out = pd.DataFrame(ordered)
        total_rows += len(out)
        yield out
        del out
      del ordered

    del wdf, chunk_wanted, slice_rows
    gc.collect()

    logger.info(
      f"[chunked_query] slice {si+1}/{n_slices} "
      f"rows={total_rows} rss={rss_mb():.0f}MB"
    )

  try:
    conn.execute("SET preserve_insertion_order=true")
  except Exception:
    pass

  dt = time.perf_counter() - t0
  logger.info(f"[chunked_query] done rows={total_rows} slices={n_slices} elapsed={dt:.2f}s rss={rss_mb():.0f}MB")
