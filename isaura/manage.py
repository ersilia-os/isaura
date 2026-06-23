import csv, datetime, json, os, shutil, sys, time, uuid, pyarrow.parquet as pq, pandas as pd
from collections import defaultdict, Counter
from contextlib import AbstractContextManager
from concurrent.futures import ThreadPoolExecutor, as_completed
from rich.progress import (
  Progress,
  SpinnerColumn,
  TextColumn,
  BarColumn,
  TaskProgressColumn,
  TimeElapsedColumn,
  TimeRemainingColumn,
)
from isaura.base import _BaseTransfer, BloomIndex, TrancheState, MinioStore, DuckDBMinio, _S3RangeFile, upload_catalog_json
from isaura.helpers import (
  ACCESS_FILE,
  BLOOM_FILENAME,
  CATALOG_FILE,
  CHECKPOINT_EVERY,
  INDEX_FILE,
  STREAM_PARQUET_THRESHOLD,
  DEFAULT_BUCKET_NAME as PUB,
  DEFAULT_PRIVATE_BUCKET_NAME as PRI,
  INPUT_C,
  MINIO_ENDPOINT_CLOUD,
  MIN_NNS_RESULT_SIZE,
  MINIO_CLOUD_AK as mcak,
  MINIO_CLOUD_SK as mcsk,
  MINIO_PRIV_CLOUD_AK as mcpak,
  MINIO_PRIV_CLOUD_SK as mcpsk,
  ReadProgress,
  StreamingCsvSink,
  WRITE_INPUT_CHUNK_ROWS,
  bind_temp_dirs,
  release_temp_dirs,
  chunk_row_limit,
  logger,
  console,
  output_dimension_from_metadata,
  rss_mb,
  fetch_schema_from_github,
  get_acc_key,
  get_apprx,
  get_base,
  get_coll,
  get_idx_key,
  get_pref,
  group_inputs,
  hive_prefix,
  list_parquet_keys,
  make_temp,
  query,
  query_batched,
  split_csv,
  stream_parquet_filtered,
  stream_parquet_filtered_ordered,
  track_write_progress,
  write_access_file,
)
from isaura.const import (
  WIDE_READ_SLICE,
  WIDE_REORDER_BUCKET_ROWS,
  WIDE_REORDER_MAX_OPEN_BUCKETS,
  WIDE_REORDER_MIN_FREE_MB,
)
from isaura.query import chunked_query_batched


class IsauraChecker(AbstractContextManager):
  """Lightweight context manager for checking and registering molecule inputs.

  Wraps a BloomIndex to provide a clean interface for deduplication checks
  without setting up a full writer. Useful when you only need to know whether
  inputs are already stored before running an expensive model, or to register
  new inputs incrementally without writing any data.

  Args:
      bucket: MinIO bucket to operate on.
      model_id: Ersilia model identifier.
      model_version: Model version string.
      store: Optional pre-built MinioStore. Creates one from env vars if absent.
      bloom_filename: Object key suffix for the Bloom filter file.
      checkpoint_every: Number of additions between automatic persist() calls.
  """

  def __init__(
    self,
    bucket,
    model_id,
    model_version,
    store=None,
    bloom_filename=BLOOM_FILENAME,
    checkpoint_every=CHECKPOINT_EVERY,
  ):
    self.bucket = bucket
    self.base_prefix = get_base(model_id, model_version)
    self.store = store or MinioStore()
    self.local_dir = make_temp("isaura_checker_")
    bind_temp_dirs(self, self.local_dir)
    self.checkpoint_every = checkpoint_every
    self.bi = BloomIndex(
      self.store, self.bucket, self.base_prefix, self.local_dir, bloom_filename=bloom_filename
    )
    self._added = 0

  def seen(self, v):
    """Return True if v has been previously stored (Bloom filter check)."""
    return self.bi.seen(v)

  def seen_many(self, vs):
    """Check multiple values at once. Returns dict mapping each value to [bool_seen, bucket]."""
    return self.bi.seen_many(vs)

  def register(self, v, rc=None):
    """Add v to the Bloom filter and trigger a checkpoint persist if needed."""
    self.bi.register(v, rc=rc)
    self._added += 1
    if self._added >= self.checkpoint_every:
      self.bi.persist()
      self._added = 0

  def close(self):
    """Persist any pending Bloom filter additions and release temporary directories."""
    try:
      if self._added:
        self.bi.persist()
        self._added = 0
    finally:
      release_temp_dirs(self)

  def __enter__(self):
    return self

  def __exit__(self, et, ev, tb):
    self.close()


class IsauraWriter:
  """Reads a CSV of model outputs and stores them as chunked Parquet files on MinIO.

  Deduplicates against the existing Bloom filter before writing, so re-running
  on a dataset that is already partially stored will only add the new rows.
  Writes are buffered in memory and flushed to MinIO in chunks of max_rows.
  On completion, uploads updated access metadata and catalog.json.

  Args:
      input_csv: Path to the CSV file containing model outputs.
      model_id: Ersilia model identifier.
      model_version: Model version string.
      access: Access level for the stored data ("public" or "private").
      bucket: Target MinIO bucket. Defaults to the public bucket.
      endpoint: MinIO endpoint URL. Defaults to local MinIO.
      access_key: MinIO access key.
      secrete: MinIO secret key.
  """

  def __init__(
    self,
    input_csv,
    model_id,
    model_version,
    bucket=None,
    endpoint=None,
    access_key=None,
    secrete=None,
  ):
    self.input_csv = input_csv
    self.model_id = model_id
    self.model_version = model_version
    self.bucket = bucket or PUB
    self.base_prefix = get_base(self.model_id, model_version)
    self.access_key = get_acc_key(self.base_prefix)
    self.model_meta = fetch_schema_from_github(self.model_id)
    self.output_dimension = output_dimension_from_metadata(self.model_meta)
    self.max_rows = int(chunk_row_limit(self.output_dimension))
    self.store = MinioStore(endpoint=endpoint, access=access_key, secret=secrete)
    self.store.require_bucket(self.bucket)
    self.tmpdir = make_temp("isaura_writter_")
    self.access = self._read_bucket_access()
    bind_temp_dirs(self, self.tmpdir)
    self.metadata_path = os.path.join(self.tmpdir, ACCESS_FILE)
    self.bi = BloomIndex(
      self.store,
      self.bucket,
      self.base_prefix,
      self.tmpdir,
      bloom_filename=os.getenv("BLOOM_FILENAME", BLOOM_FILENAME),
    )
    self.chunk_state = TrancheState(
      self.store,
      self.bucket,
      self.base_prefix,
      self.tmpdir,
      self.max_rows,
      output_dimension=self.output_dimension,
    )
    self.buffers = []
    self.buf_rows = 0
    self.schema_cols = None
    logger.debug(
      f"writer init: bucket={self.bucket} base={self.base_prefix} csv={self.input_csv} output_dim={self.output_dimension} chunk_limit={self.max_rows}"
    )

  def _read_bucket_access(self):
    """Resolve access level: hardcoded for canonical buckets, read from access.json for custom ones."""
    if self.bucket == PUB:
      return "public"
    if self.bucket == PRI:
      return "private"
    local = os.path.join(self.tmpdir, f"bucket_access_{uuid.uuid4().hex}.json")
    try:
      self.store.download_file(self.bucket, ACCESS_FILE, local)
      with open(local, "r", encoding="utf-8") as f:
        return json.load(f).get("access", "public")
    except Exception:
      return "public"

  def _load_metadata(self):
    """Download and parse the existing access metadata file, or return None if absent."""
    local = os.path.join(self.tmpdir, f"{uuid.uuid4().hex}.json")
    try:
      self.store.download_file(self.bucket, self.access_key, local)
      with open(local, "r", encoding="utf-8") as f:
        return json.load(f)
    except Exception:
      return None
    finally:
      try:
        os.remove(local)
      except Exception:
        pass

  def _upload_metadata(self, inputs):
    """Merge newly written inputs into the access metadata file and upload to MinIO."""
    if not self.metadata_path or not self.bucket:
      return
    try:
      existed = self._load_metadata()
      write_access_file(existed, inputs, self.access, self.metadata_path)
      self.store.upload_file(self.metadata_path, self.bucket, f"{self.base_prefix}/{ACCESS_FILE}")
      logger.debug(f"{ACCESS_FILE} -> minio://{self.bucket}/{self.base_prefix}/{ACCESS_FILE}")
    except Exception as e:
      logger.warning(f"metadata upload failed: {e}")

  def _set_schema(self, row):
    """Infer and store column schema from the first row seen."""
    if self.schema_cols is None:
      self.schema_cols = list(row.keys())
      logger.debug(f"writer schema: {self.schema_cols[:10]}")

  def _flush_if_needed(self):
    """Flush the in-memory buffer to MinIO if it has reached max_rows."""
    if self.buf_rows < self.max_rows or not self.buffers:
      return
    self.chunk_state.flush(self.buffers, self.schema_cols)
    self.buffers = []
    self.buf_rows = 0

  def write(self, df=None, show_progress: bool = True):
    """Write model outputs to MinIO, skipping rows already in the store.

    Reads from self.input_csv by default. Pass a DataFrame via df to write
    from memory instead. Deduplicates against the Bloom filter, buffers new
    rows, flushes to Parquet chunks, and finalises the Bloom filter, access
    metadata, and catalog.json when done.

    Args:
        df: Optional DataFrame to write instead of reading from input_csv.
        show_progress: Whether to display a Rich progress bar.
    """
    total = dupes = 0
    new = []
    old_entries = len(self.bi.index) if self.bi.index is not None else 0
    if df is not None:
      total_rows = len(df)
      chunk_iter = (
        df.iloc[start : start + WRITE_INPUT_CHUNK_ROWS] for start in range(0, len(df), WRITE_INPUT_CHUNK_ROWS)
      )
    else:
      try:
        with open(self.input_csv, "rb") as _f:
          total_rows = sum(1 for _ in _f) - 1  # subtract header
      except Exception:
        total_rows = None
      chunk_iter = pd.read_csv(self.input_csv, chunksize=WRITE_INPUT_CHUNK_ROWS)
    progress = None
    task_id = None
    try:
      if show_progress:
        progress = Progress(
          SpinnerColumn(),
          TextColumn("[progress.description]{task.description}"),
          BarColumn(),
          TaskProgressColumn(text_format="[bold bright_cyan]{task.percentage:>3.0f}%[/]"),
          TextColumn("{task.completed}" + ("/{task.total}" if total_rows else "")),
          TimeElapsedColumn(),
          TimeRemainingColumn(),
          console=logger.console,
          transient=False,
        )
        progress.__enter__()
        task_id = progress.add_task(
          f"Writing [bold]{self.model_id}[/bold] → [bold]{self.bucket}[/bold]",
          total=total_rows,
        )
      for chunk in chunk_iter:
        if chunk is None or chunk.empty:
          continue
        if self.schema_cols is None:
          self.schema_cols = list(chunk.columns)
          logger.debug(f"writer schema: {self.schema_cols[:10]}")
        if "input" in chunk.columns:
          inputs = chunk["input"].astype(str).str.strip()
        elif "smiles" in chunk.columns:
          inputs = chunk["smiles"].astype(str).str.strip()
        else:
          inputs = pd.Series([""] * len(chunk), index=chunk.index, dtype="object")
        non_blank = inputs != ""
        if hasattr(self.bi, "sbf"):
          sbf = self.bi.sbf
          not_seen = non_blank & ~inputs.map(lambda v: v in sbf)
        else:
          not_seen = non_blank & ~inputs.map(self.bi.seen)
        new_df = chunk.loc[not_seen]
        added = len(new_df)
        total += added
        dupes += len(chunk) - added
        if added:
          added_inputs = inputs.loc[not_seen].tolist()
          new.extend(added_inputs)
          self.buffers.extend(new_df.to_dict("records"))
          self.buf_rows += added
          for smi in added_inputs:
            self.bi.register(smi, rc=(1, 1))
          self._flush_if_needed()
        if progress is not None and task_id is not None:
          progress.advance(task_id, len(chunk))
    finally:
      if progress is not None:
        progress.__exit__(None, None, None)
    if self.buffers:
      self.chunk_state.flush(self.buffers, self.schema_cols)
      self.buffers = []
      self.buf_rows = 0
    self.bi.persist()
    if new:
      self._upload_metadata(new)
    entries = len(self.bi.index) if self.bi.index is not None else 0
    chunks = len(self.chunk_state._list_chunks())
    upload_catalog_json(self.store, self.bucket, self.base_prefix, entries, chunks, self.tmpdir)
    actual_new = entries - old_entries
    actual_dupes = (total + dupes) - actual_new
    dupes_str = f", [dim]{actual_dupes} duplicates skipped[/dim]" if actual_dupes else ""
    console.print(f"[green]✓[/green] [bold]{self.model_id}/{self.model_version}[/bold] → [bold]{self.bucket}[/bold]: [bold]{actual_new}[/bold] molecules written{dupes_str}")

  def close(self):
    """Release temporary directories created during construction."""
    release_temp_dirs(self)

  def __enter__(self):
    return self

  def __exit__(self, et, ev, tb):
    self.close()


class IsauraReader:
  """Retrieves precalculated model outputs from MinIO for a list of molecule inputs.

  Supports two lookup modes:
  - **Exact**: checks the Bloom filter to confirm all inputs are stored, then
    queries Parquet files via DuckDB or streaming.
  - **Approximate**: replaces inputs with their nearest-neighbor equivalents
    from the Milvus NNS service, then retrieves those stored outputs.

  Automatically switches between DuckDB and incremental Parquet streaming
  based on query size, and between normal and wide-table paths based on
  output dimension.

  Args:
      model_id: Ersilia model identifier.
      model_version: Model version string.
      input_csv: Path to a CSV file with an "input" or "smiles" column.
      approximate: If True, use nearest-neighbor search instead of exact lookup.
      bucket: Source MinIO bucket. Defaults to the public bucket.
      endpoint: MinIO endpoint URL.
      access_key: MinIO access key.
      secrete: MinIO secret key.
  """

  def __init__(
    self,
    model_id,
    model_version,
    input_csv,
    approximate,
    bucket=None,
    endpoint=None,
    access_key=None,
    secrete=None,
  ):
    self.model_id = model_id
    self.model_version = model_version
    self.approximate = approximate
    self.input_csv = input_csv
    self.bucket = bucket or PUB
    self.endpoint = endpoint
    self.base = get_base(self.model_id, model_version)
    self.pref = get_pref(self.model_id, self.model_version)
    self.collection = get_coll(self.model_id, self.model_version)
    self.index_key = get_idx_key(self.base)
    self.tmpdir = make_temp("isaura_reader_")
    bind_temp_dirs(self, self.tmpdir)
    self.store = MinioStore(endpoint=self.endpoint, access=access_key, secret=secrete)
    self.store.require_bucket(self.bucket)
    self.duck = DuckDBMinio(endpoint=self.endpoint, access=access_key, secret=secrete)
    self.bi = BloomIndex(
      self.store,
      self.bucket,
      self.base,
      self.tmpdir,
      bloom_filename=os.getenv("BLOOM_FILENAME", BLOOM_FILENAME),
      load_index=False,
    )
    try:
      self.model_meta = fetch_schema_from_github(self.model_id)
    except Exception as e:
      logger.warning(f"[read] metadata fetch failed for {self.model_id}: {e}")
      self.model_meta = {}
    self.output_dimension = output_dimension_from_metadata(self.model_meta)
    self.wide_model = self.output_dimension is not None and self.output_dimension >= 100
    logger.debug(
      f"reader init bucket={self.bucket} base={self.base} csv={self.input_csv} output_dim={self.output_dimension} wide={self.wide_model}"
    )

  def close(self):
    release_temp_dirs(self)

  def __enter__(self):
    return self

  def __exit__(self, et, ev, tb):
    self.close()

  def _wanted(self, df=None):
    """Parse molecule inputs from input_csv or a DataFrame.

    Returns:
        Tuple of (list_of_input_strings, column_header_name).
    """
    if df is None and hasattr(self, "_cached_wanted"):
      return self._cached_wanted
    wanted = []
    if df is not None:
      logger.debug(f"[read:wanted] parsing inputs from dataframe rows={len(df)}")
      columns = list(df.columns)
      header = self._resolve_input_header(columns, source="input dataframe")
      for v in df[header].tolist():
        v = ("" if v is None else str(v)).strip()
        if v:
          wanted.append(v)
    else:
      logger.debug(f"[read:wanted] parsing inputs from csv={self.input_csv}")
      with open(self.input_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        header = self._resolve_input_header(reader.fieldnames or [], source=self.input_csv)
        for row in reader:
          v = (row.get(header) or "").strip()
          if v:
            wanted.append(v)
    logger.debug(f"[read:wanted] collected inputs={len(wanted)} header={header}")
    if df is None:
      self._cached_wanted = (wanted, header)
    return (wanted, header)

  def _resolve_input_header(self, columns, source):
    """Return the molecule column among `columns`, or exit with a clear error.

    The input must have a column named one of INPUT_C ("input" or "smiles").
    A wrong/missing header is the usual cause of a "0 found" result, so fail
    loudly here naming the expected vs actual headers instead of silently
    collecting no inputs.
    """
    header = next((c for c in INPUT_C if c in columns), None)
    if header is None:
      logger.error(
        f"No molecule column found in {source}. Expected a header named one of "
        f"{list(INPUT_C)}, but found: {list(columns) or 'no columns'}. "
        f"Rename the column to 'input' or 'smiles' and try again."
      )
      sys.exit(1)
    return header

  def _prepare_read(self, df=None):
    """Validate inputs and resolve them to queryable molecule strings.

    In exact mode: checks the Bloom filter and exits if any inputs are missing.
    In approximate mode: loads the NNS index and replaces inputs with their
    nearest stored neighbors.

    Returns:
        Tuple of (start_time, resolved_wanted_list, header_column_name).
    """
    if hasattr(self, "_cached_prepare"):
      return self._cached_prepare
    index = None
    t0 = time.time()
    wanted, header = self._wanted(df=df)
    if self.approximate:
      logger.debug(f"[read:prepare] approximate mode — loading index for collection={self.collection}")
      index = self._load_index()
      if index and len(index) < MIN_NNS_RESULT_SIZE:
        logger.error(
          f"Minimum precalculation size for enabling nearest neighbor search is {MIN_NNS_RESULT_SIZE}, found {len(index)}. Aborting the Ops!"
        )
        sys.exit(1)
      st = time.perf_counter()
      wanted = get_apprx(wanted, self.collection)
      et = time.perf_counter()
      logger.debug(f"Approximate inputs are retrieved {len(wanted)} in {et - st:.2f} seconds!")
      group_inputs(wanted, index, bloom=self.bi)
    else:
      logger.debug(f"[read:prepare] exact mode — checking bloom filter for {len(wanted)} inputs")
      bloom_t0 = time.perf_counter()
      missing = [v for v in wanted if not self.bi.seen(v)]
      bloom_dt = time.perf_counter() - bloom_t0
      logger.debug(f"[read:prepare] bloom check done in {bloom_dt:.3f}s missing={len(missing)}")
      if missing:
        logger.error(
          f"inputs not indexed: {missing[:5]}{('...' if len(missing) > 5 else '')} total_missing={len(missing)}"
        )
        sys.exit(1)
    prep_dt = time.time() - t0
    logger.debug(f"[read:prepare] ready inputs={len(wanted)} header={header} prep_time={prep_dt:.2f}s")
    self._cached_prepare = (t0, wanted, header)
    return self._cached_prepare

  def _detect_parquet_col(self, header):
    """Return the molecule column name actually used in stored Parquet files.

    Downloads the first chunk file and inspects its schema. Falls back to
    header if no known INPUT_C column is found.
    """
    prefix = hive_prefix(self.base) + "/"
    try:
      keys = [
        obj["Key"]
        for obj in self.store.list_keys(self.bucket, prefix)
        if obj["Key"].endswith(".parquet") and "/chunk_" in obj["Key"]
      ]
      if not keys:
        return header
      local = os.path.join(self.tmpdir, f"schema_probe_{uuid.uuid4().hex}.parquet")
      try:
        self.store.download_file(self.bucket, keys[0], local)
        schema_names = set(pq.ParquetFile(local).schema_arrow.names)
        return next((c for c in [header] + INPUT_C if c in schema_names), header)
      finally:
        try:
          os.remove(local)
        except Exception:
          pass
    except Exception:
      return header

  def _make_read_source(self, wanted, header, batch_size=10000, ordered=True):
    """Choose a read strategy and return (mode, payload).

    Three strategies based on query size and output width:
    - "wide": Parquet streaming for models with 100+ output columns.
    - "stream": Parquet streaming for large queries (>= STREAM_PARQUET_THRESHOLD).
    - "duckdb": DuckDB batched query for smaller queries.

    Returns:
        Tuple of (mode_string, source_or_prefix).
    """
    n = len(wanted)

    def _rename(batches, src, dst):
      for chunk in batches:
        if src in chunk.columns:
          chunk = chunk.rename(columns={src: dst})
        yield chunk

    if self.wide_model:
      # Wide models stream chunk files; read_batched/read handle ordering
      # (ordered → concat+reorder; unordered pull → straight passthrough).
      prefix = hive_prefix(self.base) + "/"
      logger.debug(f"[read] wide streaming n={n} ordered={ordered} bucket={self.bucket}")
      return "wide", prefix
    if n >= STREAM_PARQUET_THRESHOLD:
      prefix = hive_prefix(self.base) + "/"
      logger.debug(f"[read] streaming n={n} bucket={self.bucket}")
      return "stream", prefix
    files = list_parquet_keys(self.store, self.bucket, self.base)
    parquet_col = self._detect_parquet_col(header)
    logger.debug(f"[read] duckdb n={n} bucket={self.bucket}")
    raw = query_batched(
      self.duck.con, parquet_col, wanted, files, batch_size=batch_size, tmpdir=self.tmpdir, preserve_order=ordered
    )
    if parquet_col == header:
      return "duckdb", raw
    return "duckdb", _rename(raw, parquet_col, header)

  def _load_index(self):
    """Download and parse the JSON index for this model. Returns [] on failure."""
    local = os.path.join(self.tmpdir, f"{uuid.uuid4().hex}.json")
    downloaded = False
    try:
      try:
        self.store.download_file(self.bucket, self.index_key, local)
        downloaded = True
      except Exception as e:
        logger.error(
          f"Exception occurred when downloading the index (bucket={self.bucket}, key={self.index_key}, local={local}): {e}"
        )
        return []
      if not downloaded or not os.path.exists(local):
        return []
      try:
        with open(local, "r", encoding="utf-8") as f:
          return json.load(f)
      except Exception as e:
        logger.error(f"Exception occurred when reading/parsing index JSON (local={local}): {e}")
        return []
    finally:
      try:
        if os.path.exists(local):
          os.remove(local)
      except Exception as e:
        logger.error(f"Exception occurred when removing temp index file (local={local}): {e}")

  def _reorder(self, wanted, header, result_df):
    """Reorder result_df rows to match the order of wanted inputs.

    Rows missing from the result are included as None-filled placeholders
    so the output always has exactly len(wanted) rows.

    Args:
        wanted: Original ordered list of input values.
        header: Column name used to match rows.
        result_df: Unordered DataFrame of retrieved rows.

    Returns:
        DataFrame aligned to the wanted order.
    """
    import numpy as np

    keys = result_df[header].astype(str).str.strip()
    lookup = {}
    for i, k in enumerate(keys):
      if k not in lookup:
        lookup[k] = i
    n = len(wanted)
    idx = np.empty(n, dtype=np.intp)
    missing_mask = np.zeros(n, dtype=bool)
    for i, inp in enumerate(wanted):
      k = str(inp).strip()
      pos = lookup.get(k, -1)
      if pos >= 0:
        idx[i] = pos
      else:
        idx[i] = 0
        missing_mask[i] = True
    matched = int(n - missing_mask.sum())
    out = result_df.iloc[idx].reset_index(drop=True)
    if missing_mask.any():
      out.loc[missing_mask] = None
      out.loc[missing_mask, header] = [str(wanted[i]).strip() for i in range(n) if missing_mask[i]]
    logger.info(f"[reorder] wanted={n} matched={matched} missing={n - matched}")
    return out

  def _reorder_batched(self, wanted, header, result_df, batch_size=10000):
    """Reorder result_df to match wanted, then yield it in batch_size chunks."""
    ordered = self._reorder(wanted, header, result_df)
    for start in range(0, len(ordered), batch_size):
      yield ordered.iloc[start : start + batch_size].copy()

  def _reorder_external(self, wanted, header, source_iter, batch_size=10000):
    """Reorder streamed rows to wanted order with bounded memory (external sort).

    Consumes ``source_iter`` — DataFrames of found rows in arbitrary store order
    (e.g. from stream_parquet_filtered) — and yields DataFrames in exact
    ``wanted`` order, batch_size rows at a time. Rows are scattered to on-disk
    "bucket" Parquet files keyed by input position (Pass A), then gathered one
    bucket at a time in order (Pass B), so peak RAM is ~one bucket rather than
    the full result set. This is the memory-safe replacement for
    concat + _reorder on wide models.

    Guarantees: output has exactly len(wanted) rows; misses appear as
    None-filled rows in their input position; duplicate inputs are emitted once
    per occurrence; pairing is matched by ``header`` value, with a guard that
    raises if any row would land against the wrong input (a fatal data error).
    """
    import numpy as np
    import pyarrow as pa

    n = len(wanted)
    if n == 0:
      return

    # input value -> list of positions (preserves order and duplicates)
    wanted_norm = [str(v).strip() for v in wanted]
    positions = {}
    for i, k in enumerate(wanted_norm):
      positions.setdefault(k, []).append(i)

    # Widen the bucket span so the number of open writers/files stays under the
    # cap, then size buckets so one fits comfortably in RAM during the gather.
    span = max(WIDE_REORDER_BUCKET_ROWS, (n + WIDE_REORDER_MAX_OPEN_BUCKETS - 1) // WIDE_REORDER_MAX_OPEN_BUCKETS)
    n_buckets = (n + span - 1) // span

    spill = os.path.join(self.tmpdir, f"reorder_{uuid.uuid4().hex}")
    os.makedirs(spill, exist_ok=True)
    try:
      free_mb = shutil.disk_usage(spill).free // (1024 * 1024)
      if free_mb < WIDE_REORDER_MIN_FREE_MB:
        logger.warning(
          f"[reorder] low spill space: {free_mb}MB free at {spill} (< {WIDE_REORDER_MIN_FREE_MB}MB)"
        )
    except Exception:
      pass

    writers = {}        # bucket idx -> pq.ParquetWriter
    bucket_paths = {}   # bucket idx -> path
    locked_schema = None
    scattered = 0
    try:
      # --- Pass A: scatter found rows to per-position-range bucket files ---
      try:
        for chunk in source_iter:
          if chunk is None or len(chunk) == 0:
            continue
          keys = chunk[header].astype(str).str.strip()
          # Fan each row out to every position its input occupies (duplicates).
          pos_lists = keys.map(positions.get)
          sub = chunk.assign(__pos=pos_lists.values)
          sub = sub[sub["__pos"].notna()]
          if sub.empty:
            continue
          sub = sub.explode("__pos")
          sub["__order"] = sub["__pos"].astype(np.int64)
          sub = sub.drop(columns="__pos")
          scattered += len(sub)
          if locked_schema is None:
            tbl = pa.Table.from_pandas(sub, preserve_index=False)
            locked_schema = tbl.schema
          else:
            tbl = pa.Table.from_pandas(sub, schema=locked_schema, preserve_index=False)
          buckets = sub["__order"].to_numpy() // span
          for b in np.unique(buckets):
            b = int(b)
            part = tbl.filter(pa.array(buckets == b))
            w = writers.get(b)
            if w is None:
              bp = os.path.join(spill, f"bucket_{b}.parquet")
              bucket_paths[b] = bp
              w = pq.ParquetWriter(bp, locked_schema)
              writers[b] = w
            w.write_table(part)
      finally:
        for w in writers.values():
          try:
            w.close()
          except Exception:
            pass

      logger.debug(
        f"[reorder] scattered rows={scattered} buckets={len(bucket_paths)} span={span} "
        f"n={n} rss={rss_mb():.0f}MB"
      )

      # --- Pass B: gather one bucket at a time, in input order ---
      cols = [c for c in (locked_schema.names if locked_schema is not None else [header]) if c != "__order"]
      for b in range(n_buckets):
        lo = b * span
        hi = min((b + 1) * span, n)
        idxrange = list(range(lo, hi))
        want_slice = wanted_norm[lo:hi]
        bp = bucket_paths.get(b)
        if bp is not None:
          bdf = pq.read_table(bp).to_pandas()
          if "__order" in bdf.columns:
            bdf = bdf.drop_duplicates(subset="__order", keep="first").set_index("__order")
          aligned = bdf.reindex(idxrange)
        else:
          aligned = pd.DataFrame(index=idxrange)
        aligned = aligned.reindex(columns=cols)
        # Fatal pairing guard: any present row must match its position's input.
        present = aligned[header].notna().to_numpy()
        if present.any():
          stored = aligned[header].astype(str).str.strip().to_numpy()
          mism = present & (stored != np.asarray(want_slice, dtype=object))
          if mism.any():
            j = int(np.argmax(mism))
            raise ValueError(
              f"[reorder] FATAL input/output mismatch at position {lo + j}: stored input "
              f"{stored[j]!r} does not match requested {want_slice[j]!r}"
            )
        # Misses keep blank descriptors; set header for every position.
        aligned[header] = want_slice
        aligned = aligned.reset_index(drop=True)
        for start in range(0, len(aligned), batch_size):
          yield aligned.iloc[start : start + batch_size].copy()
    finally:
      shutil.rmtree(spill, ignore_errors=True)

  def read(self, output_csv=None, df=None):
    """Retrieve stored outputs and return them as a DataFrame (or write to CSV).

    Validates inputs, picks the right read strategy, streams results, and
    optionally writes them to output_csv. Returns an empty DataFrame if
    output_csv is set (results are on disk) or if no inputs were found.

    Args:
        output_csv: Optional path to write results to. If set, returns empty DataFrame.
        df: Optional DataFrame of inputs to use instead of input_csv.

    Returns:
        DataFrame of results, or empty DataFrame if written to output_csv.
    """
    t0, wanted, header = self._prepare_read(df=df)
    if not wanted:
      logger.debug("[read] no inputs — returning empty frame")
      return pd.DataFrame()
    try:
      total_rows = 0
      n_chunks = 0
      mode, payload = self._make_read_source(wanted, header, ordered=True)
      with ReadProgress(total_inputs=len(wanted), console=logger.console, description=f"Reading [bold]{self.model_id}[/bold] → [bold]{self.bucket}[/bold]") as progress:
        if mode == "wide":
          # Memory-bounded external sort yields rows already in input order.
          source = self._reorder_external(
            wanted,
            header,
            stream_parquet_filtered(
              self.store,
              self.bucket,
              payload,
              wanted,
              header=header,
              batch_size=10000,
              progress=progress,
            ),
            10000,
          )
        elif mode == "stream":
          source = stream_parquet_filtered_ordered(
            self.store,
            self.bucket,
            payload,
            wanted,
            header=header,
            progress=progress,
          )
        else:
          source = payload
        parts = []
        for chunk in source:
          n_chunks += 1
          parts.append(chunk)
          total_rows += len(chunk)
          progress.update(
            stage="scanning", emitted_rows=total_rows, unresolved=max(0, len(wanted) - total_rows)
          )
        result = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        del parts
      if output_csv and not result.empty:
        logger.debug(f"[read] writing csv to {output_csv}")
        with StreamingCsvSink(output_csv) as sink:
          sink.write_df(result)
        out = pd.DataFrame()
      elif output_csv:
        out = pd.DataFrame()
      else:
        out = result
      elapsed = time.time() - t0
      logger.debug(f"Query successfully fetched for a given inputs in {elapsed:.2f} sec")
    except Exception as e:
      logger.error(f"[read] query failed: {e}")
      sys.exit(1)
    rate = total_rows / elapsed if elapsed > 0 and total_rows else 0.0
    logger.debug(
      f"[read] done model={self.pref} inputs={len(wanted)} matched={total_rows} "
      f"chunks={n_chunks} elapsed={time.time() - t0:.2f}s rate={rate:.1f}/s rss={rss_mb():.0f}MB"
    )
    return out

  def read_batched(self, batch_size=10000, output_csv=None, df=None, ordered=True):
    """Retrieve stored outputs as a generator of DataFrames, optionally writing to CSV.

    Memory-efficient alternative to read() for large result sets. Yields one
    DataFrame per batch. If output_csv is set, also streams results to disk.

    Args:
        batch_size: Number of rows per yielded DataFrame.
        output_csv: Optional path to write all results to as they are streamed.
        df: Optional DataFrame of inputs to use instead of input_csv.
        ordered: When True (default), rows are emitted in the same order as the
            wanted inputs, with None-filled placeholders for misses. When False,
            rows are streamed straight from the parquet files without buffering
            or reordering — bounded memory, found rows only. Callers that don't
            care about row order (e.g. pull, which re-indexes by molecule) should
            pass ordered=False to avoid materializing the full result set, which
            OOMs for wide models (output dimension >= 100) over large inputs.

    Yields:
        DataFrames of matching rows in batch_size chunks.
    """
    t0, wanted, header = self._prepare_read(df=df)
    if not wanted:
      logger.info("[read_batched] no inputs — returning")
      return
    total_rows = 0
    n_chunks = 0
    sink = None
    source = None
    try:
      mode, payload = self._make_read_source(wanted, header, batch_size=batch_size, ordered=ordered)
      with ReadProgress(
        total_inputs=len(wanted),
        console=logger.console,
        description=f"Reading [bold]{self.model_id}[/bold] → [bold]{self.bucket}[/bold]",
        writing_description=(f"Writing [bold]{self.model_id}[/bold] → [bold]{output_csv}[/bold]" if output_csv else None),
      ) as progress:
        if mode == "wide":
          if ordered:
            # Memory-bounded external sort: scatter found rows to on-disk
            # buckets keyed by input position, then gather in order. Replaces
            # the concat + full-frame _reorder that OOMs on wide models.
            source = self._reorder_external(
              wanted,
              header,
              stream_parquet_filtered(
                self.store,
                self.bucket,
                payload,
                wanted,
                header=header,
                batch_size=batch_size,
                progress=progress,
              ),
              batch_size,
            )
          else:
            source = stream_parquet_filtered(
              self.store,
              self.bucket,
              payload,
              wanted,
              header=header,
              batch_size=batch_size,
              progress=progress,
            )
        elif mode == "stream":
          if ordered:
            source = stream_parquet_filtered_ordered(
              self.store,
              self.bucket,
              payload,
              wanted,
              header=header,
              batch_size=batch_size,
              progress=progress,
            )
          else:
            source = stream_parquet_filtered(
              self.store,
              self.bucket,
              payload,
              wanted,
              header=header,
              batch_size=batch_size,
              progress=progress,
            )
        else:
          source = payload
        if output_csv:
          logger.debug(f"[read_batched] streaming csv to {output_csv}")
          sink = StreamingCsvSink(output_csv)
          sink.__enter__()
        for chunk in source:
          n_chunks += 1
          if sink is not None:
            sink.write_df(chunk)
          total_rows += len(chunk)
          progress.update(
            stage="emitting", emitted_rows=total_rows, unresolved=max(0, len(wanted) - total_rows)
          )
          yield chunk
    except Exception as e:
      logger.error(f"[read_batched] query failed: {e}")
      sys.exit(1)
    finally:
      if sink is not None:
        sink.__exit__(None, None, None)
        logger.debug(f"[read_batched] csv closed rows={total_rows} path={output_csv}")
      close = getattr(source, "close", None)
      if close is not None:
        close()
    elapsed = time.time() - t0
    rate = total_rows / elapsed if elapsed > 0 and total_rows else 0.0
    logger.debug(
      f"[read_batched] done model={self.pref} inputs={len(wanted)} matched={total_rows} "
      f"chunks={n_chunks} elapsed={elapsed:.2f}s rate={rate:.1f}/s rss={rss_mb():.0f}MB"
    )


class IsauraCopy(_BaseTransfer):
  """Copies a model's data from a staging bucket to the public and/or private buckets.

  Routing is determined by the access metadata: each input is sent to
  isaura-public or isaura-private based on its "access" field.
  If output_dir is set, dumps all objects to the local filesystem instead.
  """

  def copy(self):
    """Run the copy operation. Loads access metadata and routes rows to the right buckets."""
    logger.debug(f"[copy] starting model={self.model_id} v={self.model_version} bucket={self.bucket}")
    meta_local, meta = self._load_metadata()
    return self._copy(meta_local, meta)


class IsauraMover(_BaseTransfer):
  """Moves a model's data between buckets by copying then deleting the source."""

  def move(self):
    """Copy data to the destination buckets, then delete all source objects."""
    t0 = time.time()
    logger.info(f"[move] starting model={self.model_id} v={self.model_version} bucket={self.bucket}")
    meta_local, meta = self._load_metadata()
    self._copy(meta_local, meta)
    n = self._delete()
    logger.success(f"[move] done wiped={n} elapsed={time.time() - t0:.1f}s")


class IsauraRemover(_BaseTransfer):
  """Deletes a model's data from a bucket, or all objects in a project bucket.

  If project_name is set, deletes the entire project bucket without needing
  a model_id or model_version. Otherwise delegates to _BaseTransfer._delete().
  """

  def __init__(self, model_id=None, model_version=None, bucket=None, project_name=None):
    self._project_name = project_name
    if project_name:
      self._store = MinioStore()
    else:
      super().__init__(model_id, model_version, bucket)

  def remove(self):
    """Delete all objects for this model (or entire project bucket if project_name is set)."""
    if self._project_name:
      store = self._store
      n = store.delete_prefix(self._project_name, "")
      logger.info(f"removed all objects in project={self._project_name} count={n}")
    else:
      n = self._delete()
      logger.info(f"removed objects={n}")


class IsauraMolRemover:
  """Removes specific molecules from a model's stored data in a bucket.

  Downloads each Parquet chunk, filters out the requested rows, and re-uploads
  the modified chunks. Rebuilds the Bloom filter and JSON index from scratch,
  and updates the per-model access.json and catalog.json.

  Args:
      model_id: Ersilia model identifier.
      model_version: Model version string.
      bucket: Target MinIO bucket.
      input_csv: Path to a CSV file with an "input" or "smiles" column.
  """

  def __init__(self, model_id, model_version, bucket, input_csv):
    self.model_id = model_id
    self.model_version = model_version
    self.bucket = bucket
    self.input_csv = input_csv
    self.base_prefix = get_base(model_id, model_version)
    self.store = MinioStore()
    self.store.require_bucket(bucket)
    self.tmpdir = make_temp("isaura_mol_remover_")
    bind_temp_dirs(self, self.tmpdir)

  def _read_inputs(self):
    """Parse the input CSV and return a list of molecule strings."""
    inputs = []
    with open(self.input_csv, newline="", encoding="utf-8") as f:
      for row in csv.DictReader(f):
        v = (row.get("input") or row.get("smiles") or "").strip()
        if v:
          inputs.append(v)
    return inputs

  def _list_chunk_keys(self):
    """Return a sorted list of Parquet chunk keys for this model."""
    pref = hive_prefix(self.base_prefix) + "/"
    keys = [
      obj["Key"] for obj in self.store.list_keys(self.bucket, pref)
      if obj["Key"].endswith(".parquet") and "/chunk_" in obj["Key"]
    ]
    def _idx(key):
      try:
        return int(os.path.basename(key).split("_")[1].split(".")[0])
      except Exception:
        return 0
    return sorted(keys, key=_idx)

  def _rewrite_chunk(self, chunk_key, to_remove):
    """Download a chunk, filter out to_remove rows, re-upload (or delete if empty).

    Returns:
        Tuple of (rows_before, rows_after).
    """
    import pyarrow as pa
    local = os.path.join(self.tmpdir, f"chunk_dl_{uuid.uuid4().hex}.parquet")
    local_out = os.path.join(self.tmpdir, f"chunk_out_{uuid.uuid4().hex}.parquet")
    try:
      self.store.download_file(self.bucket, chunk_key, local)
      pf = pq.ParquetFile(local)
      schema_names = pf.schema_arrow.names
      input_col = next((c for c in ["input", "smiles"] if c in schema_names), schema_names[0])
      rows_before = pf.metadata.num_rows
      df = pf.read().to_pandas()
      mask = ~df[input_col].astype(str).str.strip().isin(to_remove)
      filtered = df.loc[mask].reset_index(drop=True)
      rows_after = len(filtered)
      if rows_after == 0:
        self.store.client.delete_object(Bucket=self.bucket, Key=chunk_key)
        logger.debug(f"[remove] deleted empty chunk {chunk_key}")
      else:
        table = pa.Table.from_pandas(filtered, preserve_index=False)
        pq.write_table(table, local_out)
        self.store.upload_file(local_out, self.bucket, chunk_key)
        logger.debug(f"[remove] rewrote chunk {chunk_key}: {rows_before}→{rows_after} rows")
      return (rows_before, rows_after)
    except Exception as e:
      logger.error(f"[remove] chunk {chunk_key}: {e}")
      return (0, 0)
    finally:
      for f in [local, local_out]:
        try:
          os.remove(f)
        except Exception:
          pass

  def _update_access_json(self, to_remove):
    """Remove deleted molecules from the per-model access.json and re-upload."""
    acc_key = get_acc_key(self.base_prefix)
    local = os.path.join(self.tmpdir, f"access_upd_{uuid.uuid4().hex}.json")
    try:
      self.store.download_file(self.bucket, acc_key, local)
      with open(local, "r", encoding="utf-8") as f:
        data = json.load(f)
      if isinstance(data, list):
        filtered = [d for d in data if (d.get("input") or "").strip() not in to_remove]
        with open(local, "w", encoding="utf-8") as f:
          json.dump(filtered, f, indent=2)
        self.store.upload_file(local, self.bucket, acc_key)
        logger.debug(f"[remove] access.json: {len(data)}→{len(filtered)}")
    except Exception as e:
      logger.debug(f"[remove] access.json update skipped: {e}")
    finally:
      try:
        os.remove(local)
      except Exception:
        pass

  def remove(self, show_progress=True):
    """Remove the specified molecules from the store.

    Returns:
        Tuple of (n_removed, n_not_found).
    """
    to_remove_list = self._read_inputs()
    if not to_remove_list:
      logger.error("No valid molecules found in input file.")
      sys.exit(1)
    to_remove = set(to_remove_list)

    bi = BloomIndex(self.store, self.bucket, self.base_prefix, self.tmpdir)
    if not bi.index:
      logger.error(f"No index found for {self.model_id}/{self.model_version} in {self.bucket}.")
      sys.exit(1)

    actually_present = to_remove & set(bi.index.keys())
    n_not_found = len(to_remove) - len(actually_present)

    if not actually_present:
      return (0, len(to_remove))

    chunk_keys = self._list_chunk_keys()
    if not chunk_keys:
      logger.error("No chunk files found.")
      sys.exit(1)

    total_removed = 0
    progress = None
    task_id = None
    try:
      if show_progress:
        progress = Progress(
          SpinnerColumn(),
          TextColumn("[progress.description]{task.description}"),
          BarColumn(),
          TextColumn("{task.completed}/{task.total}"),
          TimeElapsedColumn(),
          console=logger.console,
          transient=False,
        )
        progress.__enter__()
        task_id = progress.add_task(
          f"Removing from [bold]{self.model_id}[/bold] → [bold]{self.bucket}[/bold]",
          total=len(chunk_keys),
        )
      for key in chunk_keys:
        rows_before, rows_after = self._rewrite_chunk(key, actually_present)
        total_removed += rows_before - rows_after
        if progress is not None and task_id is not None:
          progress.advance(task_id, 1)
    finally:
      if progress is not None:
        progress.__exit__(None, None, None)

    for k in actually_present:
      bi.index.pop(k, None)
    from pybloom_live import ScalableBloomFilter
    bi.sbf = ScalableBloomFilter(
      mode=ScalableBloomFilter.SMALL_SET_GROWTH,
      initial_capacity=max(1000000, len(bi.index)),
      error_rate=1e-14,
    )
    for k in bi.index:
      bi.sbf.add(k)
    bi._added = 1
    bi.persist()

    self._update_access_json(actually_present)

    remaining_chunks = self._list_chunk_keys()
    upload_catalog_json(
      self.store, self.bucket, self.base_prefix, len(bi.index), len(remaining_chunks), self.tmpdir
    )

    return (total_removed, n_not_found)

  def close(self):
    """Release temporary directories."""
    release_temp_dirs(self)

  def __enter__(self):
    return self

  def __exit__(self, et, ev, tb):
    self.close()


class IsauraPull(_BaseTransfer):
  """Downloads precalculated outputs from the cloud MinIO and stores them locally.

  Uses IsauraReader pointed at the cloud endpoint to retrieve the data for
  the requested inputs, then writes them to the local MinIO via _pull_batched.
  """

  def __init__(self, model_id, model_version, bucket, input_csv, output_dir=None):
    super().__init__(model_id, model_version, bucket, output_dir)
    self.model_id = model_id
    self.model_version = model_version
    self.bucket = bucket
    self.input_csv = input_csv

  def pull(self):
    """Read outputs from cloud MinIO for the given inputs and store them locally."""
    t0 = time.time()
    input_name = os.path.basename(self.input_csv) if self.input_csv else ""
    console.print(f"Pulling [bold]{self.model_id}[/bold] ({self.model_version}) from [bold]{self.bucket}[/bold]")
    with IsauraReader(
      model_id=self.model_id,
      model_version=self.model_version,
      input_csv=self.input_csv,
      approximate=False,
      bucket=self.bucket,
      endpoint=MINIO_ENDPOINT_CLOUD,
      access_key=mcak,
      secrete=mcsk,
    ) as r:
      with console.status("Loading inputs...", spinner="dots"):
        wanted, _ = r._wanted()
      console.print(f"[green]✓[/green] {len(wanted):,} inputs loaded from {input_name}")
      with console.status("Checking cloud index...", spinner="dots"):
        r._prepare_read()
      console.print(f"[green]✓[/green] All {len(wanted):,} found in cloud index")
      fetch_t0 = time.time()
      out = self._pull_batched(r.read_batched(ordered=False))
      console.print(f"[green]✓[/green] {out[0]:,} rows fetched and stored locally ({time.time() - fetch_t0:.1f}s)")
    console.print(f"[green]✓[/green] [bold]{self.model_id}/{self.model_version}[/bold] pulled [bold]{out[0]:,}[/bold] rows from [bold]{self.bucket}[/bold] in {time.time() - t0:.1f}s")
    logger.debug(f"pulled objects={out}")


class IsauraPush:
  """Uploads a model's local Parquet data to the cloud MinIO in parallel.

  Uploads all chunk files from both the public and private local buckets to
  their cloud counterparts (isaura-public and isaura-private). Bloom filter
  and index files are merged rather than blindly overwritten, so cloud entries
  already present are preserved.
  """

  def __init__(self, model_id, model_version, bucket):
    self.model_id = model_id
    self.model_version = model_version
    self.bucket = bucket

  def _bucket_exists(self, store, bucket):
    """Return True if the given bucket exists and is accessible."""
    try:
      store.client.head_bucket(Bucket=bucket)
      return True
    except Exception:
      return False

  def _relay_file(self, local_store, src_bucket, cloud_store, dst_bucket, key, tmpdir):
    """Download a file from local MinIO and upload it to cloud MinIO.

    Returns:
        Tuple of (key, None) on success, or (None, error_string) on failure.
    """
    local = os.path.join(tmpdir, uuid.uuid4().hex)
    try:
      local_store.download_file(src_bucket, key, local)
      cloud_store.upload_file(local, dst_bucket, key)
      return (key, None)
    except Exception as e:
      return (None, f"[push] failed {key}: {e}")
    finally:
      try:
        os.remove(local)
      except Exception:
        pass

  def _merge_bloom_index(self, local_store, src_bucket, cloud_store, dst_bucket, tmpdir):
    """Load cloud bloom+index, merge local entries in, persist back to cloud. Returns (entries, chunks)."""
    base = get_base(self.model_id, self.model_version)
    cloud_bi = BloomIndex(cloud_store, dst_bucket, base, os.path.join(tmpdir, "cloud"))
    local_bi = BloomIndex(local_store, src_bucket, base, os.path.join(tmpdir, "local"))
    if local_bi.index:
      for k, v in local_bi.index.items():
        if cloud_bi.index is not None and k not in cloud_bi.index:
          cloud_bi.index[k] = v
    if local_bi.index:
      for k in local_bi.index:
        cloud_bi.sbf.add(k)
    cloud_bi._added = 1
    cloud_bi.persist()
    merged_entries = len(cloud_bi.index or {})
    # count cloud chunks from listing
    pref = hive_prefix(base) + "/"
    cloud_chunks = sum(
      1 for obj in cloud_store.list_keys(dst_bucket, pref)
      if obj["Key"].endswith(".parquet") and "/chunk_" in obj["Key"]
    )
    logger.debug(f"[push] merged bloom+index: entries={merged_entries} chunks={cloud_chunks}")
    return (merged_entries, cloud_chunks)

  def push(self):
    """Upload all local Parquet chunks to cloud MinIO and merge the Bloom filter and index."""
    local_store = MinioStore()
    prefix = get_pref(self.model_id, self.model_version) + "/"
    skip_suffixes = (BLOOM_FILENAME, INDEX_FILE, CATALOG_FILE)
    has_data = False
    for access, src_bucket, mck, mcs in [
      ("public", PUB, mcak, mcsk),
      ("private", PRI, mcpak, mcpsk),
    ]:
      if not self._bucket_exists(local_store, src_bucket):
        logger.warning(f"[push] local bucket {src_bucket} does not exist. Skipping {access}.")
        continue
      keys = [
        obj["Key"] for obj in local_store.list_keys(src_bucket, prefix)
        if not obj["Key"].endswith(skip_suffixes)
      ]
      if not keys:
        logger.debug(f"[push] no data for {self.model_id} in {src_bucket}. Skipping {access}.")
        continue
      has_data = True
      dst_bucket = f"isaura-{access}"
      cloud_store = MinioStore(endpoint=MINIO_ENDPOINT_CLOUD, access=mck, secret=mcs)
      cloud_store.ensure_bucket(dst_bucket)
      tmpdir = make_temp("isaura_push_")
      logger.debug(
        f"[push] uploading {len(keys)} objects for {self.model_id}/{self.model_version} "
        f"({access}) to cloud"
      )
      try:
        done, errors = 0, []
        with logger.console.status(
          f"Pushing {len(keys)} objects to cloud ({access})...", spinner="dots"
        ):
          with ThreadPoolExecutor(max_workers=4) as pool:
            futs = [
              pool.submit(self._relay_file, local_store, src_bucket, cloud_store, dst_bucket, k, tmpdir)
              for k in keys
            ]
            for fut in as_completed(futs):
              key, err = fut.result()
              if key is not None:
                done += 1
              if err is not None:
                errors.append(err)
        for err in errors:
          logger.warning(err)
        logger.debug(f"[push] {access} done: {done}/{len(keys)} objects uploaded to cloud")
        with logger.console.status("Merging bloom filter and index...", spinner="dots"):
          merged_entries, cloud_chunks = self._merge_bloom_index(
            local_store, src_bucket, cloud_store, dst_bucket, tmpdir
          )
        base = get_base(self.model_id, self.model_version)
        upload_catalog_json(cloud_store, dst_bucket, base, merged_entries, cloud_chunks, tmpdir)
      finally:
        try:
          shutil.rmtree(tmpdir)
        except Exception:
          pass
    if not has_data:
      logger.error("No data found in any default bucket for a given model! Aborting push.")
      sys.exit(1)
    console.print(f"[green]✓[/green] [bold]{self.model_id}/{self.model_version}[/bold] pushed to cloud")


class IsauraInspect:
  """Inspects models, buckets, and inputs in MinIO without modifying any data.

  Supports both local and cloud MinIO. Provides methods for listing models,
  checking which inputs are available, counting rows, and reading Parquet
  metadata (columns, row counts) via HTTP range requests to avoid full downloads.

  Args:
      model_id: If set, scopes all operations to this model.
      model_version: If set, further scopes to this version.
      cloud: If True, connects to cloud MinIO using cloud credentials.
      project_name: If set, operates on a custom bucket instead of the defaults.
      access: Which buckets to inspect — "public", "private", or "both".
      endpoint: Custom MinIO endpoint. Ignored if cloud=True.
      parquet_block_size: Block size for _S3RangeFile when reading Parquet metadata.
      max_small_json_bytes: Size limit for JSON files fetched inline (larger files are skipped).
      heavy_index: If True, loads the full JSON index even in cloud mode.
  """

  def __init__(
    self,
    model_id=None,
    model_version=None,
    cloud=False,
    project_name=None,
    access="both",
    endpoint=None,
    parquet_block_size=262144,
    max_small_json_bytes=20000000,
    heavy_index=False,
  ):
    self.model_id, self.model_version = (model_id, model_version)
    self.cloud, self.project_name, self.access = (cloud, project_name, access)
    self.endpoint = MINIO_ENDPOINT_CLOUD if cloud else endpoint
    self.parquet_block_size = int(parquet_block_size)
    self.max_small_json_bytes = int(max_small_json_bytes)
    self.heavy_index = bool(heavy_index)
    self._cache = {}
    self._rowcount_cache = {}

  def buckets(self):
    """Return the list of bucket names to operate on based on access and project_name."""
    if self.project_name:
      return [self.project_name]
    if self.access == "public":
      return [PUB]
    if self.access == "private":
      return [PRI]
    return [PUB, PRI]

  def _creds(self, bucket):
    """Return MinioStore constructor kwargs for the given bucket (local or cloud credentials)."""
    if not self.cloud:
      return {"endpoint": self.endpoint}
    if bucket == PUB:
      return {"endpoint": MINIO_ENDPOINT_CLOUD, "access": mcak, "secret": mcsk}
    if bucket == PRI:
      return {"endpoint": MINIO_ENDPOINT_CLOUD, "access": mcpak, "secret": mcpsk}
    return {"endpoint": MINIO_ENDPOINT_CLOUD}

  def _clients(self, bucket):
    k = (bucket, self.cloud, self.endpoint)
    if k not in self._cache:
      c = self._creds(bucket)
      self._cache[k] = (MinioStore(**c),)
    return self._cache[k]

  def _client(self, bucket):
    return self._clients(bucket)[0].client

  def _paginator(self, bucket):
    return self._client(bucket).get_paginator("list_objects_v2")

  def list_prefixes(self, bucket, prefix=""):
    """Yield common prefixes (virtual directories) under prefix, one level deep."""
    p = self._paginator(bucket)
    for page in p.paginate(Bucket=bucket, Prefix=prefix, Delimiter="/"):
      for cp in page.get("CommonPrefixes", []):
        yield cp["Prefix"]

  def list_objects(self, bucket, prefix=""):
    """Yield all object keys under prefix (paginated, recursive)."""
    p = self._paginator(bucket)
    for page in p.paginate(Bucket=bucket, Prefix=prefix):
      for obj in page.get("Contents", []):
        yield obj["Key"]

  def iter_object_meta(self, bucket, prefix=""):
    """Yield metadata dicts (Key, Size, LastModified, ETag, StorageClass) for all objects under prefix."""
    p = self._paginator(bucket)
    for page in p.paginate(Bucket=bucket, Prefix=prefix):
      for obj in page.get("Contents", []):
        yield {
          "Key": obj.get("Key", ""),
          "Size": obj.get("Size", 0),
          "LastModified": obj.get("LastModified"),
          "ETag": obj.get("ETag"),
          "StorageClass": obj.get("StorageClass"),
        }

  def get_json(self, bucket, key, max_bytes=None):
    """Fetch and parse a JSON object from MinIO. Returns None if too large or on error."""
    client = self._client(bucket)
    lim = self.max_small_json_bytes if max_bytes is None else int(max_bytes)
    try:
      try:
        h = client.head_object(Bucket=bucket, Key=key)
        if h and int(h.get("ContentLength", 0) or 0) > lim:
          return None
      except Exception:
        pass
      o = client.get_object(Bucket=bucket, Key=key)
      body = o["Body"].read()
      try:
        o["Body"].close()
      except Exception:
        pass
      return json.loads(body.decode("utf-8"))
    except Exception:
      return None

  def iter_models(self, bucket, prefix_filter=""):
    """Yield (model_id, model_version) tuples for all models found in the bucket."""
    if self.model_id:
      m = self.model_id
      if self.model_version:
        yield (m, self.model_version)
        return
      for v_pref in self.list_prefixes(bucket, f"{m}/"):
        mv = v_pref[len(f"{m}/") :].strip("/")
        if mv:
          yield (m, mv)
      return
    for m_pref in self.list_prefixes(bucket, ""):
      mid = m_pref.strip("/")
      if prefix_filter and (not mid.startswith(prefix_filter)):
        continue
      for v_pref in self.list_prefixes(bucket, m_pref):
        mv = v_pref[len(m_pref) :].strip("/")
        if mv:
          yield (mid, mv)

  def load_index(self, bucket, model_id, model_version):
    """Load the JSON index for a model. Returns {} in cloud mode unless heavy_index=True."""
    base = get_base(model_id, model_version)
    key = get_idx_key(base)
    if self.heavy_index or not self.cloud:
      return self.get_json(bucket, key, max_bytes=10**18) or {}
    return {}

  def load_metadata(self, bucket, model_id, model_version):
    """Load the access metadata file for a model. Returns None if absent."""
    base = get_base(model_id, model_version)
    return self.get_json(bucket, f"{base}/{ACCESS_FILE}")

  def _indices_union(self, force=False):
    union, owner = ({}, {})
    for b in self.buckets():
      for mid, mv in self.iter_models(b):
        if force:
          idx = self.get_json(b, get_idx_key(get_base(mid, mv)), max_bytes=10**18) or {}
        else:
          idx = self.load_index(b, mid, mv)
        for smi in idx.keys():
          if smi not in union:
            union[smi] = True
            owner[smi] = b
    return (union, owner)

  def list_available(self, output_csv=None):
    """Return a DataFrame of all stored inputs with a single column: input."""
    _, owner = self._indices_union(force=self.cloud)
    df = pd.DataFrame([{"input": smi} for smi in owner])
    if output_csv:
      df.to_csv(output_csv, index=False)
    return df

  def inspect_inputs(self, input_csv, output_csv=None):
    """Check which inputs from input_csv are available in the store.

    Returns a DataFrame with a single column: input.
    """
    _, owner = self._indices_union(force=self.cloud)
    with open(input_csv, newline="", encoding="utf-8") as f:
      wanted = [(r.get("input") or "").strip() for r in csv.DictReader(f) if (r.get("input") or "").strip()]
    df = pd.DataFrame([{"input": s} for s in wanted if s in owner])
    if output_csv:
      df.to_csv(output_csv, index=False)
    return df

  def find_chunk_keys(self, bucket, model_id, model_version):
    """Yield (key, size_bytes) for every Parquet chunk file for this model."""
    pref = get_pref(model_id, model_version)
    for obj in self.iter_object_meta(bucket, pref):
      k = obj.get("Key") or ""
      if "/chunk_" in k and k.endswith(".parquet"):
        yield (k, int(obj.get("Size") or 0))

  def parquet_num_rows(self, bucket, parquet_key, size=None):
    """Read the row count from a Parquet file's footer without downloading the full file.

    Uses _S3RangeFile for range-based HTTP access. Returns None on failure.
    Results are cached for the lifetime of this instance.
    """
    if pq is None or not parquet_key:
      return None
    ck = (bucket, parquet_key)
    if ck in self._rowcount_cache:
      return self._rowcount_cache[ck]
    client = self._client(bucket)
    try:
      if size is None:
        h = client.head_object(Bucket=bucket, Key=parquet_key)
        size = int(h.get("ContentLength", 0) or 0)
      size = int(size or 0)
      if size <= 0:
        return None
      f = _S3RangeFile(client, bucket, parquet_key, size=size, block_size=self.parquet_block_size)
      pf = pq.ParquetFile(f)
      md = pf.metadata
      n = int(md.num_rows) if md is not None else None
      try:
        f.close()
      except Exception:
        pass
      if n is not None:
        self._rowcount_cache[ck] = n
      return n
    except Exception:
      return None

  def entries_from_chunks(self, bucket, model_id, model_version):
    """Return (total_row_count, chunk_count) for a model.

    Tries catalog.json first for a fast answer. Falls back to reading row
    counts from each Parquet footer in parallel if catalog.json is absent.
    """
    # try catalog.json first
    base = get_base(model_id, model_version)
    cat = self.get_json(bucket, f"{base}/{CATALOG_FILE}")
    if cat and "entries" in cat:
      return (int(cat["entries"]), int(cat.get("chunks", 0)))
    # fallback: read parquet footers
    chunks = list(self.find_chunk_keys(bucket, model_id, model_version))
    if not chunks:
      return (0, 0)
    with ThreadPoolExecutor(max_workers=min(8, len(chunks))) as pool:
      futs = {pool.submit(self.parquet_num_rows, bucket, k, sz): k for k, sz in chunks}
      total = 0
      for fut in as_completed(futs):
        n = fut.result()
        if n is not None:
          total += int(n)
    return (total, len(chunks))

  def inspect_models(self, bucket, prefix_filter=""):
    """Return a list of model info dicts for all models in the bucket.

    In cloud mode without heavy_index, uses the fast listing path
    (_inspect_models_from_listing). Otherwise loads the full JSON index.
    """
    if self.cloud and (not self.heavy_index):
      return self._inspect_models_from_listing(bucket, prefix_filter)
    out = []
    for mid, mv in self.iter_models(bucket, prefix_filter=prefix_filter):
      out.append({
        "bucket": bucket, "model_id": mid, "model_version": mv,
        "model": f"{mid}/{mv}", "entries": len(self.load_index(bucket, mid, mv)), "chunks": None,
      })
    return out

  def _inspect_models_from_listing(self, bucket, prefix_filter=""):
    """Single LIST for discovery. catalog.json for fast counts, footer reads as fallback."""
    stats = defaultdict(lambda: {"chunks": [], "bytes": 0, "has_catalog": False, "catalog_key": None})
    for obj in self.iter_object_meta(bucket, prefix_filter):
      k = obj.get("Key", "")
      parts = k.split("/")
      if len(parts) < 2:
        continue
      mid, mv = parts[0], parts[1]
      if k.endswith(f"/{CATALOG_FILE}"):
        stats[(mid, mv)]["has_catalog"] = True
        stats[(mid, mv)]["catalog_key"] = k
        continue
      if "/chunk_" not in k or not k.endswith(".parquet"):
        continue
      sz = int(obj.get("Size") or 0)
      s = stats[(mid, mv)]
      s["chunks"].append((k, sz))
      s["bytes"] += sz

    row_counts = {}

    # fast path: fetch catalog.json for models that have it
    catalog_models = [(mid, mv, s["catalog_key"]) for (mid, mv), s in stats.items() if s["has_catalog"]]
    if catalog_models:
      def _fetch_cat(key):
        return self.get_json(bucket, key)
      with ThreadPoolExecutor(max_workers=min(16, len(catalog_models))) as pool:
        futs = {pool.submit(_fetch_cat, ck): (mid, mv) for mid, mv, ck in catalog_models}
        for fut in as_completed(futs):
          cat = fut.result()
          if cat and "entries" in cat:
            row_counts[futs[fut]] = int(cat["entries"])

    # slow fallback: footer reads for models without catalog.json
    fallback_tasks = []
    for (mid, mv), s in stats.items():
      if (mid, mv) in row_counts:
        continue
      for k, sz in s["chunks"]:
        fallback_tasks.append(((mid, mv), k, sz))
    if fallback_tasks:
      fallback_counts = defaultdict(int)
      with ThreadPoolExecutor(max_workers=min(16, len(fallback_tasks))) as pool:
        futs = {
          pool.submit(self.parquet_num_rows, bucket, k, sz): (mid, mv)
          for (mid, mv), k, sz in fallback_tasks
        }
        for fut in as_completed(futs):
          n = fut.result()
          if n is not None:
            fallback_counts[futs[fut]] += int(n)
      row_counts.update(fallback_counts)

    out = []
    for (mid, mv), s in sorted(stats.items()):
      b = s["bytes"]
      if b >= 1 << 30:
        size = f"[bold red]{b / (1 << 30):.1f} GB[/]"
      elif b >= 100 * (1 << 20):
        size = f"[yellow]{b / (1 << 20):.1f} MB[/]"
      elif b >= 1 << 20:
        size = f"[green]{b / (1 << 20):.1f} MB[/]"
      elif b >= 1 << 10:
        size = f"[dim]{b / (1 << 10):.1f} KB[/]"
      else:
        size = f"[dim]{b} B[/]"
      rows = row_counts.get((mid, mv), 0)
      out.append({
        "bucket": bucket, "model_id": mid, "model_version": mv,
        "model": f"{mid}/{mv}", "entries": f"{rows:,}", "size": size, "chunks": len(s["chunks"]),
      })
    return out

  def find_any_chunk_key(self, bucket, model_id, model_version):
    """Return the key of the first Parquet chunk found for a model, or None."""
    pref = get_pref(model_id, model_version)
    for k in self.list_objects(bucket, pref):
      if "/chunk_" in k and k.endswith(".parquet"):
        return k
    return None

  def parquet_columns(self, bucket, parquet_key):
    """Return the column names of a Parquet file by reading its footer via range request."""
    if not parquet_key or pq is None:
      return None
    client = self._client(bucket)
    try:
      h = client.head_object(Bucket=bucket, Key=parquet_key)
      size = int(h.get("ContentLength", 0) or 0)
      if size <= 0:
        return None
      f = _S3RangeFile(client, bucket, parquet_key, size=size, block_size=self.parquet_block_size)
      pf = pq.ParquetFile(f)
      cols = list(pf.schema_arrow.names)
      try:
        f.close()
      except Exception:
        pass
      return cols
    except Exception:
      return None

  def duckdb_columns(self, bucket, parquet_key):
    return self.parquet_columns(bucket, parquet_key)


class IsauraStat:
  """Computes storage statistics across all models in one or more buckets.

  Produces a JSON report with per-model molecule counts, storage sizes, column
  counts, and GitHub metadata. In cloud mode uses a fast parallel path
  (_compute_cloud) that makes a single LIST per bucket and fetches counts from
  catalog.json. In local mode reads the full JSON index for exact counts.

  Args:
      project_name: If set, scopes stats to a custom bucket.
      access: Which buckets to include — "public", "private", or "both".
      cloud: If True, connects to cloud MinIO.
      endpoint: Custom MinIO endpoint.
      include_columns: If True, reads the number of output columns from each model's first chunk.
      include_column_names: If True, also includes the full list of column names.
      schema_version: Version tag included in the output JSON.
      producer: Producer label included in the output JSON.
      max_workers: Max parallel threads for fetching model data.
  """

  _MAX_WORKERS = 8

  def __init__(
    self,
    project_name=None,
    access="both",
    cloud=False,
    endpoint=None,
    include_columns=True,
    include_column_names=False,
    schema_version="1",
    producer="isaura stats",
    max_workers=None,
  ):
    self.insp = IsauraInspect(cloud=cloud, project_name=project_name, access=access, endpoint=endpoint)
    self.include_columns = bool(include_columns)
    self.include_column_names = bool(include_column_names)
    self.schema_version = str(schema_version)
    self.producer = str(producer)
    self.max_workers = int(max_workers) if max_workers is not None else self._MAX_WORKERS

  def _percentile(self, sorted_vals, p):
    """Return the p-th percentile of a pre-sorted list, or None if empty."""
    if not sorted_vals:
      return None
    n = len(sorted_vals)
    i = int(round(p / 100.0 * (n - 1)))
    i = max(0, min(i, n - 1))
    return sorted_vals[i]

  def _use_chunks(self):
    """Return True if the fast cloud path (chunk-based counting) should be used."""
    return self.insp.cloud and (not self.insp.heavy_index)

  def _collect_objects(self, bucket, model_id, model_version):
    """List all objects for a model and return (total_bytes, total_gb, chunks, first_chunk_key)."""
    pref = get_pref(model_id, model_version)
    total_bytes = 0
    chunks = []
    first_chunk = None
    for obj in self.insp.iter_object_meta(bucket, pref):
      total_bytes += int(obj.get("Size") or 0)
      k = obj.get("Key") or ""
      if "/chunk_" in k and k.endswith(".parquet"):
        sz = int(obj.get("Size") or 0)
        chunks.append((k, sz))
        if first_chunk is None:
          first_chunk = k
    total_gb = round(total_bytes / 1024**3, 6)
    return (total_bytes, total_gb, chunks, first_chunk)

  def _count_rows_from_chunks(self, bucket, chunks):
    """Sum row counts across a list of (key, size) chunk tuples in parallel."""
    total = 0
    if not chunks:
      return total
    with ThreadPoolExecutor(max_workers=min(8, len(chunks))) as pool:
      futs = {pool.submit(self.insp.parquet_num_rows, bucket, k, sz): k for k, sz in chunks}
      for fut in as_completed(futs):
        n = fut.result()
        if n is not None:
          total += int(n)
    return total

  def _process_model(self, bucket, mid, mv, use_chunks):
    """Build a stats row dict for a single model. Used in local mode."""
    total_bytes, total_gb, chunks, first_chunk = self._collect_objects(bucket, mid, mv)
    if use_chunks:
      # try catalog.json first
      base = get_base(mid, mv)
      cat = self.insp.get_json(bucket, f"{base}/{CATALOG_FILE}")
      if cat and "entries" in cat:
        mol_count = int(cat["entries"])
      else:
        mol_count = self._count_rows_from_chunks(bucket, chunks)
      idx = None
    else:
      idx = self.insp.load_index(bucket, mid, mv)
      mol_count = len(idx)
    meta_out = fetch_schema_from_github(mid)
    cols = None
    ncols = None
    if self.include_columns and first_chunk:
      cols = self.insp.parquet_columns(bucket, first_chunk)
      ncols = len(cols) if cols else None
    row = {
      "bucket": bucket,
      "model_id": mid,
      "model_version": mv,
      "model_key": f"{bucket}:{mid}:{mv}",
      "model": f"{mid}/{mv}",
      "molecules": mol_count,
      "total_bytes": total_bytes,
      "total_gb": total_gb,
      "n_columns": ncols,
      "metadata": meta_out,
    }
    if self.include_column_names:
      row["columns"] = cols
    return (row, idx)

  def _compute_cloud(self, buckets):
    """Fast cloud path: single LIST, catalog.json for counts, parallel GitHub + column fetches."""
    generated_at = datetime.datetime.now(datetime.timezone.utc).isoformat()
    # single LIST per bucket — gather everything
    model_data = defaultdict(lambda: {"bucket": None, "total_bytes": 0, "chunks": [], "first_chunk": None,
                                       "catalog_key": None})
    for b in buckets:
      for obj in self.insp.iter_object_meta(b, ""):
        k = obj.get("Key", "")
        parts = k.split("/")
        if len(parts) < 2:
          continue
        mid, mv = parts[0], parts[1]
        key = (b, mid, mv)
        d = model_data[key]
        d["bucket"] = b
        sz = int(obj.get("Size") or 0)
        d["total_bytes"] += sz
        if k.endswith(f"/{CATALOG_FILE}"):
          d["catalog_key"] = k
        elif "/chunk_" in k and k.endswith(".parquet"):
          d["chunks"].append((k, sz))
          if d["first_chunk"] is None:
            d["first_chunk"] = k

    tasks = sorted(model_data.keys())
    if not tasks:
      return self._empty_result(generated_at, buckets)

    # parallel: catalog.json + GitHub metadata + parquet columns
    def _enrich(key):
      b, mid, mv = key
      d = model_data[key]
      # row count from catalog.json or fallback
      mol_count = 0
      if d["catalog_key"]:
        cat = self.insp.get_json(b, d["catalog_key"])
        if cat and "entries" in cat:
          mol_count = int(cat["entries"])
      if not mol_count and d["chunks"]:
        mol_count = self._count_rows_from_chunks(b, d["chunks"])
      # GitHub metadata
      meta_out = fetch_schema_from_github(mid)
      # columns
      cols, ncols = None, None
      if self.include_columns and d["first_chunk"]:
        cols = self.insp.parquet_columns(b, d["first_chunk"])
        ncols = len(cols) if cols else None
      total_bytes = d["total_bytes"]
      row = {
        "bucket": b, "model_id": mid, "model_version": mv,
        "model_key": f"{b}:{mid}:{mv}", "model": f"{mid}/{mv}",
        "molecules": mol_count, "total_bytes": total_bytes,
        "total_gb": round(total_bytes / 1024**3, 6),
        "n_columns": ncols, "metadata": meta_out,
      }
      if self.include_column_names:
        row["columns"] = cols
      return row

    models = []
    models_total_by_bucket = defaultdict(int)
    with ThreadPoolExecutor(max_workers=min(self.max_workers, len(tasks))) as pool:
      futs = {pool.submit(_enrich, k): k for k in tasks}
      for fut in as_completed(futs):
        row = fut.result()
        models.append(row)
        models_total_by_bucket[row["bucket"]] += 1
    models.sort(key=lambda r: (r["bucket"], r["model_id"], r["model_version"]))
    molecules_total = sum(m["molecules"] for m in models)
    return {
      "schema_version": self.schema_version, "producer": self.producer,
      "generated_at_utc": generated_at, "buckets": buckets,
      "models_total": len(models), "models_total_by_bucket": dict(models_total_by_bucket),
      "models": models,
      "models_per_molecule": {
        "molecules_total_unique": molecules_total,
        "min": None, "max": None, "p50": None, "p95": None, "histogram": [],
      },
    }

  def _empty_result(self, generated_at, buckets):
    """Return a zero-filled stats result dict for when no models are found."""
    return {
      "schema_version": self.schema_version, "producer": self.producer,
      "generated_at_utc": generated_at, "buckets": buckets,
      "models_total": 0, "models_total_by_bucket": {},
      "models": [],
      "models_per_molecule": {
        "molecules_total_unique": 0,
        "min": None, "max": None, "p50": None, "p95": None, "histogram": [],
      },
    }

  def compute(self):
    """Compute and return the full stats report as a dict.

    In cloud mode delegates to _compute_cloud (fast, parallel, catalog.json-based).
    In local mode iterates models sequentially and reads the full JSON index.

    Returns:
        Dict with schema_version, generated_at_utc, models list, and aggregate stats.
    """
    generated_at = datetime.datetime.now(datetime.timezone.utc).isoformat()
    buckets = self.insp.buckets()
    use_chunks = self._use_chunks()
    if use_chunks:
      return self._compute_cloud(buckets)
    tasks = []
    for b in buckets:
      for mid, mv in self.insp.iter_models(b):
        tasks.append((b, mid, mv))
    models = []
    models_per_molecule = Counter()
    models_total_by_bucket = defaultdict(int)
    for b, mid, mv in tasks:
      row, idx = self._process_model(b, mid, mv, False)
      if idx is not None:
        for smi in idx.keys():
          models_per_molecule[smi] += 1
      models.append(row)
      models_total_by_bucket[b] += 1
    counts = sorted(models_per_molecule.values()) if models_per_molecule else []
    hist = Counter(counts)
    hist_out = [{"models": k, "molecules": v} for k, v in sorted(hist.items(), key=lambda x: x[0])]
    molecules_total = sum((m["molecules"] for m in models))
    return {
      "schema_version": self.schema_version,
      "producer": self.producer,
      "generated_at_utc": generated_at,
      "buckets": buckets,
      "models_total": len(models),
      "models_total_by_bucket": dict(models_total_by_bucket),
      "models": models,
      "models_per_molecule": {
        "molecules_total_unique": len(models_per_molecule) if models_per_molecule else molecules_total,
        "min": min(counts) if counts else None,
        "max": max(counts) if counts else None,
        "p50": self._percentile(counts, 50),
        "p95": self._percentile(counts, 95),
        "histogram": hist_out,
      },
    }

  def write_json(self, output_path):
    """Compute stats and write the result to output_path as formatted JSON.

    Uses an atomic write (tmp file + rename) to avoid leaving a partial file.

    Returns:
        The output_path that was written.
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    tmp = output_path + ".tmp"
    data = self.compute()
    with open(tmp, "w", encoding="utf-8") as f:
      json.dump(data, f, indent=2, ensure_ascii=False)
    os.replace(tmp, output_path)
    return output_path
