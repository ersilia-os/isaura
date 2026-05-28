import os

import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from isaura.const import INPUT_C, STREAM_DENSE_BATCH_ROWS, STREAM_DENSE_FILE_RATIO
from isaura.logging import logger
from isaura.utils import cleanup_temp_dir, make_temp, rss_mb


def _yield_filtered_chunks(filtered, header, batch_size, wanted_set, progress_cb):
  """Slice a filtered PyArrow Table into DataFrames and track resolved molecules.

  Yields each slice as a DataFrame and removes matched molecules from
  wanted_set so the caller knows which inputs are still outstanding.

  Args:
      filtered: PyArrow Table already filtered to matching rows.
      header: Column name identifying each molecule (e.g. "input").
      batch_size: Maximum rows per yielded DataFrame.
      wanted_set: Mutable set of molecule strings still to be found. Modified in place.
      progress_cb: Optional callable receiving the row count of each batch.

  Yields:
      Tuple of (DataFrame, n_newly_resolved).
  """
  remaining_delta = 0
  for start in range(0, filtered.num_rows, batch_size):
    chunk = filtered.slice(start, batch_size).to_pandas(split_blocks=True, self_destruct=True)
    matched = set(chunk[header].astype(str).str.strip())
    n = len(matched & wanted_set)
    wanted_set -= matched
    remaining_delta += n
    if progress_cb:
      progress_cb(len(chunk))
    yield chunk, n
    del chunk


def stream_parquet_filtered(
  store,
  bucket,
  prefix,
  wanted,
  header="input",
  batch_size=10000,
  dense_file_ratio=STREAM_DENSE_FILE_RATIO,
  progress=None,
):
  """Stream matching rows from Parquet chunk files on MinIO without loading all files into memory.

  Downloads one chunk file at a time, filters it for the wanted molecules,
  and yields matching rows as DataFrames. Stops early once all wanted molecules
  have been found. Two read strategies are used per file:

  - **Dense mode**: when most rows in a file are expected to match (ratio >= dense_file_ratio),
    reads the file in large batches and filters each batch in memory.
  - **Row-group mode**: reads only the input column first to check for matches,
    then fetches the full row group only if matches exist (saves I/O on sparse queries).

  Args:
      store: MinioStore used to download chunk files.
      bucket: Source bucket.
      prefix: Key prefix under which chunk files are listed.
      wanted: List or set of molecule input strings to look up.
      header: Column name identifying each molecule (default "input").
      batch_size: Maximum rows per yielded DataFrame.
      dense_file_ratio: Fraction of a file's rows that must be wanted to use dense mode.
      progress: Optional ReadProgress instance for live progress updates.

  Yields:
      DataFrames of matching rows, in file order.
  """
  wanted_set = set(wanted) if not isinstance(wanted, set) else wanted
  remaining = len(wanted_set)
  wanted_arr = pa.array(list(wanted_set))
  tmpdir = make_temp("isaura_stream_")
  keys = []
  total_rows_yielded = 0
  n_chunks_yielded = 0
  try:
    keys = sorted(
      (
        obj["Key"]
        for obj in store.list_keys(bucket, prefix)
        if obj["Key"].endswith(".parquet") and "/chunk_" in obj["Key"]
      )
    )

    logger.debug(f"[stream] starting: {len(keys)} files, {remaining} wanted, batch_size={batch_size}")
    if progress is not None:
      progress.update(
        stage="starting",
        files_done=0,
        files_total=len(keys),
        found_rows=0,
        emitted_rows=0,
        unresolved=remaining,
      )

    for ki, key in enumerate(keys):
      if remaining <= 0:
        break
      if remaining < len(wanted_arr) * 0.5:
        wanted_arr = pa.array(list(wanted_set))
      if progress is not None:
        progress.update(stage="scanning file", files_done=ki, files_total=len(keys), unresolved=remaining)

      local = os.path.join(tmpdir, f"s_{ki}.parquet")
      try:
        store.download_file(bucket, key, local)
      except Exception as e:
        logger.warning(f"[stream] skip {key}: {e}")
        continue

      try:
        pf = pq.ParquetFile(local)
        n_rg = pf.metadata.num_row_groups
        n_rows = pf.metadata.num_rows
        use_dense = (
          bool(dense_file_ratio) and n_rows > 0 and remaining >= max(1, int(n_rows * dense_file_ratio))
        )

        schema_names = set(pf.schema_arrow.names)
        parquet_col = next((c for c in [header] + INPUT_C if c in schema_names), header)

        if use_dense:
          dense_batch_rows = max(batch_size, STREAM_DENSE_BATCH_ROWS)
          for batch in pf.iter_batches(batch_size=dense_batch_rows):
            if remaining <= 0:
              break
            table = pa.Table.from_batches([batch])
            mask = pc.is_in(table.column(parquet_col), wanted_arr)
            if pc.any(mask).as_py() is not True:
              continue
            filtered = table.filter(mask)
            if filtered.num_rows == 0:
              continue
            for start in range(0, filtered.num_rows, batch_size):
              chunk = filtered.slice(start, batch_size).to_pandas(split_blocks=True, self_destruct=True)
              if parquet_col != header:
                chunk = chunk.rename(columns={parquet_col: header})
              matched = set(chunk[header].astype(str).str.strip())
              remaining -= len(matched & wanted_set)
              wanted_set -= matched
              n_chunks_yielded += 1
              total_rows_yielded += len(chunk)
              if progress is not None:
                progress.update(
                  stage=f"yielding matches {ki + 1}/{len(keys)}",
                  files_done=ki + 1,
                  files_total=len(keys),
                  found_rows=total_rows_yielded,
                  emitted_rows=total_rows_yielded,
                  unresolved=remaining,
                )
              yield chunk
              del chunk
        else:
          for rg_idx in range(n_rg):
            if remaining <= 0:
              break
            key_col = pf.read_row_group(rg_idx, columns=[parquet_col]).column(parquet_col)
            mask = pc.is_in(key_col, wanted_arr)
            if pc.any(mask).as_py() is not True:
              del key_col, mask
              continue
            filtered = pf.read_row_group(rg_idx).filter(mask)
            del key_col, mask
            if filtered.num_rows == 0:
              del filtered
              continue
            for start in range(0, filtered.num_rows, batch_size):
              chunk = filtered.slice(start, batch_size).to_pandas(split_blocks=True, self_destruct=True)
              if parquet_col != header:
                chunk = chunk.rename(columns={parquet_col: header})
              matched = set(chunk[header].astype(str).str.strip())
              remaining -= len(matched & wanted_set)
              wanted_set -= matched
              n_chunks_yielded += 1
              total_rows_yielded += len(chunk)
              if progress is not None:
                progress.update(
                  stage=f"yielding matches {ki + 1}/{len(keys)}",
                  files_done=ki + 1,
                  files_total=len(keys),
                  found_rows=total_rows_yielded,
                  emitted_rows=total_rows_yielded,
                  unresolved=remaining,
                )
              yield chunk
              del chunk
            del filtered

        logger.debug(
          f"[stream] file {ki + 1}/{len(keys)} done key={key.split('/')[-1]} "
          f"yielded={total_rows_yielded} remaining={remaining} rss={rss_mb():.0f}MB"
        )
        if progress is not None:
          progress.update(stage="advancing", files_done=ki + 1, files_total=len(keys), unresolved=remaining)
        del pf
      except Exception as e:
        logger.warning(f"[stream] error reading {key}: {e}")
      finally:
        try:
          os.remove(local)
        except Exception:
          pass
  finally:
    logger.debug(
      f"[stream] finished: files={len(keys)} yielded={total_rows_yielded} "
      f"chunks={n_chunks_yielded} remaining={remaining} rss={rss_mb():.0f}MB"
    )
    cleanup_temp_dir(tmpdir)


def stream_parquet_filtered_ordered(
  store,
  bucket,
  prefix,
  wanted,
  header="input",
  batch_size=10000,
  progress=None,
):
  """Stream matching rows from MinIO in the same order as the wanted input list.

  Wraps stream_parquet_filtered and buffers resolved rows until they can be
  emitted in the original wanted order. Molecules not found in the store are
  emitted as None-filled placeholder rows so the output always has exactly
  len(wanted) rows.

  Args:
      store: MinioStore used to download chunk files.
      bucket: Source bucket.
      prefix: Key prefix under which chunk files are listed.
      wanted: Ordered list of molecule input strings.
      header: Column name identifying each molecule (default "input").
      batch_size: Maximum rows per yielded DataFrame.
      progress: Optional ReadProgress instance for live progress updates.

  Yields:
      DataFrames of rows aligned to the wanted order, with None rows for misses.
  """
  wanted_list = [str(v).strip() for v in wanted if str(v).strip()]
  if not wanted_list:
    return

  unresolved = set(wanted_list)
  resolved = {}
  schema_cols = None
  next_idx = 0
  emitted = 0

  def missing_row(key):
    base = {header: key}
    if schema_cols:
      for col in schema_cols:
        if col != header:
          base[col] = None
    return base

  def make_frame(rows):
    if not rows:
      return pd.DataFrame(columns=schema_cols or [header])
    if schema_cols:
      return pd.DataFrame.from_records(rows, columns=schema_cols)
    return pd.DataFrame.from_records(rows)

  def flush_ready():
    nonlocal next_idx, emitted
    rows = []
    while next_idx < len(wanted_list):
      key = wanted_list[next_idx]
      row = resolved.get(key)
      if row is None:
        if key in unresolved:
          break
        rows.append(missing_row(key))
      else:
        rows.append(dict(row))
      next_idx += 1
      emitted += 1
      if len(rows) >= batch_size:
        yield make_frame(rows)
        rows = []
    if rows:
      yield make_frame(rows)

  source = stream_parquet_filtered(
    store,
    bucket,
    prefix,
    set(unresolved),
    header=header,
    batch_size=batch_size,
    progress=progress,
  )
  try:
    for chunk in source:
      if chunk is None or chunk.empty:
        continue
      if schema_cols is None:
        schema_cols = list(chunk.columns)
      for row in chunk.to_dict("records"):
        key = str(row.get(header) or "").strip()
        if key and key not in resolved:
          resolved[key] = row
          unresolved.discard(key)
      for ready in flush_ready():
        yield ready
      if not unresolved:
        break
  finally:
    close = getattr(source, "close", None)
    if close is not None:
      close()

  unresolved_before_fill = len(unresolved)
  unresolved.clear()
  if next_idx < len(wanted_list):
    tail = []
    for key in wanted_list[next_idx:]:
      row = resolved.get(key)
      tail.append(dict(row) if row is not None else missing_row(key))
      if len(tail) >= batch_size:
        yield make_frame(tail)
        tail = []
    if tail:
      yield make_frame(tail)
    if unresolved_before_fill:
      logger.warning(f"[stream-ordered] missing={unresolved_before_fill} rows; emitted blank placeholders")
