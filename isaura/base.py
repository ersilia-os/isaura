import boto3, duckdb, gc, json, os, pandas as pd, pickle, pyarrow as pa, pyarrow.parquet as pq, requests, sqlite3, sys, time, uuid
from isaura.const import INPUT_C as _INPUT_C, INDEX_SQLITE_FILE, INDEX_FORMAT, INDEX_GRANULARITY
from botocore.config import Config
from boto3.s3.transfer import TransferConfig
from pybloom_live import ScalableBloomFilter
from isaura.helpers import (
  MINIO_ENDPOINT,
  MINIO_LOCAL_AK,
  MINIO_LOCAL_SK,
  MAX_ROWS,
  CHECKPOINT_EVERY,
  BLOOM_FILENAME,
  ACCESS_FILE,
  CATALOG_FILE,
  INDEX_FILE,
  STREAM_PARQUET_THRESHOLD,
  DEFAULT_BUCKET_NAME as pub_bucket,
  DEFAULT_PRIVATE_BUCKET_NAME as priv_bucket,
  MINIO_MULTIPART_THRESHOLD,
  MINIO_MULTIPART_CHUNKSIZE,
  MINIO_MAX_CONCURRENCY,
  MINIO_MAX_POOL_CONNECTIONS,
  WIDE_TARGET_FILE_BYTES,
  WIDE_TARGET_ROWGROUP_BYTES,
  logger,
  console,
  rss_mb,
  is_wide,
  chunk_row_limit,
  chunk_write_batch_rows,
  parquet_writer_kwargs,
  build_typed_array,
  resolve_write_types,
  fetch_schema_from_github,
  get_acc_key,
  get_base,
  get_coll,
  get_idx_key,
  output_dimension_from_metadata,
  get_pref,
  hive_prefix,
  list_parquet_keys,
  make_temp,
  post_apprx,
  query,
  query_batched,
  bind_temp_dirs,
  release_temp_dirs,
  stream_parquet_filtered,
  cpu_cnt,
)


class _S3RangeFile:
  """Seekable file-like object backed by S3/MinIO HTTP range requests.

  Lets libraries like PyArrow read Parquet footer metadata from MinIO
  without downloading the entire file. Fetches fixed-size blocks on
  demand and caches them for the lifetime of the object.

  Args:
      client: boto3 S3 client.
      bucket: Bucket name.
      key: Object key.
      size: Total object size in bytes.
      block_size: Number of bytes per cached block (default 256 KB).
  """

  def __init__(self, client, bucket, key, size, block_size=262144):
    self._client = client
    self._bucket = bucket
    self._key = key
    self._size = int(size or 0)
    self._block_size = int(block_size)
    self._pos = 0
    self._cache = {}
    self.closed = False

  def readable(self):
    """Return True — this stream supports reading."""
    return True

  def seekable(self):
    """Return True — this stream supports seeking."""
    return True

  def tell(self):
    """Return the current byte offset."""
    return self._pos

  def seek(self, offset, whence=0):
    """Move the current position.

    Args:
        offset: Byte offset relative to whence.
        whence: 0 = start, 1 = current position, 2 = end.

    Returns:
        The new absolute byte position.
    """
    if whence == 0:
      new_pos = int(offset)
    elif whence == 1:
      new_pos = int(self._pos + offset)
    elif whence == 2:
      new_pos = int(self._size + offset)
    else:
      raise ValueError("invalid whence")
    if new_pos < 0:
      new_pos = 0
    if new_pos > self._size:
      new_pos = self._size
    self._pos = new_pos
    return self._pos

  def _get_block(self, block_idx):
    """Fetch and cache a single block by index, using an S3 range request."""
    if block_idx in self._cache:
      return self._cache[block_idx]
    start = block_idx * self._block_size
    if start >= self._size:
      data = b""
      self._cache[block_idx] = data
      return data
    end = min(self._size - 1, start + self._block_size - 1)
    resp = self._client.get_object(Bucket=self._bucket, Key=self._key, Range=f"bytes={start}-{end}")
    body = resp["Body"].read()
    try:
      resp["Body"].close()
    except Exception:
      pass
    self._cache[block_idx] = body
    return body

  def read(self, n=-1):
    """Read up to n bytes from the current position.

    Args:
        n: Number of bytes to read. Reads to end-of-file if -1 or None.

    Returns:
        Bytes read. Returns b"" when at end-of-file or if closed.
    """
    if self.closed:
      return b""
    if n is None or n < 0:
      n = self._size - self._pos
    n = int(n)
    if n <= 0 or self._pos >= self._size:
      return b""
    end_pos = min(self._size, self._pos + n)
    out = bytearray()
    while self._pos < end_pos:
      block_idx = self._pos // self._block_size
      block_off = self._pos % self._block_size
      block = self._get_block(block_idx)
      if not block:
        break
      take = min(len(block) - block_off, end_pos - self._pos)
      if take <= 0:
        break
      out += block[block_off : block_off + take]
      self._pos += take
    return bytes(out)

  def close(self):
    """Release the block cache and mark the stream as closed."""
    self.closed = True
    self._cache.clear()


class MinioStore:
  """boto3 S3 client wrapper for a MinIO endpoint.

  Checks that the MinIO server is reachable on construction and exits
  if it is not. Provides upload, download, listing, and delete helpers
  used throughout the rest of the codebase.

  Args:
      endpoint: MinIO HTTP endpoint URL. Defaults to MINIO_ENDPOINT env var.
      access: Access key. Defaults to MINIO_LOCAL_AK env var.
      secret: Secret key. Defaults to MINIO_LOCAL_SK env var.
      multipart_threshold: File size (bytes) above which uploads are split into chunks.
      multipart_chunksize: Size of each upload chunk in bytes.
      max_concurrency: Number of chunks to upload in parallel.
      use_threads: Whether to use threads for transfers.
  """

  def __init__(
    self,
    endpoint=None,
    access=None,
    secret=None,
    multipart_threshold=None,
    multipart_chunksize=None,
    max_concurrency=None,
    use_threads=True,
  ):
    self.endpoint = endpoint or MINIO_ENDPOINT
    self.access = access or MINIO_LOCAL_AK
    self.secret = secret or MINIO_LOCAL_SK
    self.client = boto3.client(
      "s3",
      endpoint_url=self.endpoint,
      aws_access_key_id=self.access,
      aws_secret_access_key=self.secret,
      region_name="us-east-1",
      config=Config(
        signature_version="s3v4",
        s3={"addressing_style": "path"},
        max_pool_connections=MINIO_MAX_POOL_CONNECTIONS,
        tcp_keepalive=True,
        retries={"mode": "adaptive", "max_attempts": 5},
      ),
    )
    self.transfer_config = TransferConfig(
      multipart_threshold=multipart_threshold or MINIO_MULTIPART_THRESHOLD,
      multipart_chunksize=multipart_chunksize or MINIO_MULTIPART_CHUNKSIZE,
      max_concurrency=max_concurrency or MINIO_MAX_CONCURRENCY,
      use_threads=use_threads,
    )
    if not self.ping(self.client):
      sys.exit(1)

  def ping(self, store):
    """Return True if the MinIO health endpoint responds with 200, False otherwise."""
    try:
      url = f"{self.endpoint.rstrip('/')}/minio/health/ready"
      resp = requests.get(url, timeout=3)
      if resp.status_code == 200:
        logger.debug(f"MinIO server healthy at {url}")
        return True
      logger.error(
        f"MinIO health check failed. Maybe the minio containers not running[{resp.status_code}]: {resp.text.strip()}"
      )
    except requests.exceptions.RequestException as e:
      logger.error(f"MinIO health request error. Maybe the minio containers not running: {e}")
    except Exception as e:
      logger.error(f"Unexpected error during MinIO health check: {e}")
    return False

  def credential_check(self) -> bool:
    """Return True if the current credentials are accepted by the server."""
    try:
      self.client.list_buckets()
      return True
    except Exception:
      return False

  def ensure_bucket(self, bucket):
    """Create the bucket if it does not already exist. Exits on failure."""
    try:
      self.client.head_bucket(Bucket=bucket)
    except Exception:
      try:
        self.client.create_bucket(Bucket=bucket)
      except Exception as e:
        logger.error(e)
        sys.exit(1)

  def require_bucket(self, bucket):
    """Exit with a helpful error if the bucket does not exist."""
    try:
      self.client.head_bucket(Bucket=bucket)
    except Exception:
      logger.error(
        f"Project [bold]{bucket}[/bold] does not exist. "
        f"Run [bold]isaura info[/bold] to see available projects or [bold]isaura create -pn {bucket} --access <public|private>[/bold] to create it."
      )
      sys.exit(1)

  def download_file(self, bucket, key, local):
    """Download an object from MinIO to a local path. Raises on missing key."""
    try:
      self.client.download_file(bucket, key, local, Config=self.transfer_config)
    except Exception as e:
      logger.debug(f"The file {key} is not found in bucket {bucket}. Details -> {e}")
      raise

  def upload_file(self, local, bucket, key, extra_args=None):
    """Upload a local file to MinIO at the given bucket and key."""
    self.client.upload_file(local, bucket, key, ExtraArgs=extra_args or {}, Config=self.transfer_config)

  def list_keys(self, bucket, prefix):
    """Yield all S3 object metadata dicts under the given prefix (paginated)."""
    p = self.client.get_paginator("list_objects_v2")
    for page in p.paginate(Bucket=bucket, Prefix=prefix):
      for obj in page.get("Contents", []):
        yield obj

  def delete_prefix(self, bucket, prefix):
    """Delete all objects under a prefix in batches of 1000. Returns the count deleted."""
    batch = []
    deleted = 0
    for obj in self.list_keys(bucket, prefix):
      batch.append({"Key": obj["Key"]})
      if len(batch) == 1000:
        self.client.delete_objects(Bucket=bucket, Delete={"Objects": batch})
        deleted += len(batch)
        batch = []
    if batch:
      self.client.delete_objects(Bucket=bucket, Delete={"Objects": batch})
      deleted += len(batch)
    return deleted


class DuckDBMinio:
  """Singleton DuckDB connection for querying Parquet files directly on MinIO.

  DuckDB runs entirely in RAM (no database file) and uses its httpfs extension
  to read Parquet over S3. On first construction the connection is configured
  with the MinIO endpoint and credentials so that subsequent SQL queries can
  reference s3:// paths without any extra setup.

  Re-construction with identical settings is a no-op, so this can be safely
  instantiated multiple times across the codebase.

  Args:
      endpoint: MinIO HTTP endpoint URL.
      access: Access key.
      secret: Secret key.
      threads: Number of DuckDB query threads. Defaults to all available CPUs.
  """

  _instance = None

  def __new__(cls, *a, **kw):
    if cls._instance is None:
      cls._instance = super().__new__(cls)
    return cls._instance

  def __init__(self, endpoint=None, access=None, secret=None, threads=None):
    endpoint = endpoint or MINIO_ENDPOINT
    access = access or MINIO_LOCAL_AK
    secret = secret or MINIO_LOCAL_SK
    threads = threads or max(2, os.cpu_count() or 2)
    cfg = (endpoint, access, secret, threads)
    if getattr(self, "_cfg", None) == cfg:
      return
    if getattr(self, "con", None):
      self.con.close()
    self.endpoint, self.access, self.secret = (endpoint, access, secret)
    self.con = duckdb.connect(":memory:")
    self.con.execute("PRAGMA enable_object_cache")
    self.con.execute(f"PRAGMA threads={cpu_cnt(ratio=1)};")
    self.con.execute("INSTALL httpfs; LOAD httpfs;")
    self.con.execute("SET s3_access_key_id=?", [access])
    self.con.execute("SET s3_secret_access_key=?", [secret])
    self.con.execute("SET s3_endpoint=?", [endpoint.replace("http://", "").replace("https://", "")])
    self.con.execute("SET s3_region='us-east-1';")
    self.con.execute("SET s3_use_ssl=?", [endpoint.startswith("https://")])
    self.con.execute("SET s3_url_style='path'")
    self._cfg = cfg

  def close(self):
    """Close the DuckDB connection and release its memory."""
    if getattr(self, "con", None):
      self.con.close()
      self.con = None


class BloomIndex:
  """Bloom filter and JSON index for the Parquet store, persisted on MinIO.

  The Bloom filter answers "have we already stored this molecule?" without
  scanning any Parquet files, preventing duplicate writes cheaply. The JSON
  index maps each molecule SMILES to its (row, chunk) coordinates so that
  exact lookups can open just the one relevant chunk file instead of all of them.

  Both files are downloaded from MinIO on construction (or created fresh if
  absent) and re-uploaded on every call to persist(). An automatic checkpoint
  upload is triggered every CHECKPOINT_EVERY additions.

  Args:
      store: MinioStore used for uploads and downloads.
      bucket: Bucket where the Bloom filter and index files live.
      base_prefix: Key prefix for this model (e.g. "eos1234/v1").
      local_dir: Local directory used to cache the downloaded files.
      bloom_filename: Object key suffix for the Bloom filter file.
      error_rate: Acceptable false-positive rate for the Bloom filter (default ~0).
      initial_capacity: Expected number of entries at creation time.
      load_index: If False, skip loading the JSON index. Useful during bulk
          writes where only duplicate-checking is needed, not coordinate lookup.
  """

  def __init__(
    self,
    store,
    bucket,
    base_prefix,
    local_dir,
    bloom_filename=BLOOM_FILENAME,
    error_rate=1e-14,
    initial_capacity=1000000,
    load_index=True,
  ):
    self.store = store
    self.bucket = bucket
    self.base = base_prefix.strip("/")
    self.bloom_key = f"{self.base}/{bloom_filename}"
    self.index_key = f"{self.base}/{INDEX_FILE}"
    self.local_bloom = os.path.join(local_dir, bloom_filename)
    self.local_index = os.path.join(local_dir, "index.json")
    os.makedirs(local_dir, exist_ok=True)
    try:
      self.store.download_file(self.bucket, self.bloom_key, self.local_bloom)
      with open(self.local_bloom, "rb") as f:
        self.sbf = pickle.load(f)
      logger.debug(f"loaded bloom {self.bucket}/{self.bloom_key}")
    except Exception:
      self.sbf = ScalableBloomFilter(
        mode=ScalableBloomFilter.SMALL_SET_GROWTH, initial_capacity=initial_capacity, error_rate=error_rate
      )
      logger.debug("created new bloom")
    if load_index:
      try:
        self.store.download_file(self.bucket, self.index_key, self.local_index)
        with open(self.local_index, "r", encoding="utf-8") as f:
          self.index = json.load(f)
        logger.debug(f"loaded index {self.bucket}/{self.index_key} entries={len(self.index)}")
      except Exception:
        self.index = {}
        logger.debug("created new index")
    else:
      self.index = None
      logger.debug("skipped index load")
    self._added = 0

  def seen(self, v):
    """Return True if v is in the Bloom filter (molecule was previously stored)."""
    return v in self.sbf

  def seen_many(self, vs):
    """Check multiple values against the Bloom filter at once.

    Returns:
        Dict mapping each value to [bool_seen, bucket_name].
    """
    sbf = self.sbf
    b = self.bucket
    return {v: [v in sbf, b] for v in vs}

  def rc(self, v):
    """Return the (row, chunk) coordinates for v from the JSON index, or None if absent."""
    return self.index.get(v) if self.index is not None else None

  def register(self, v, rc=None):
    """Add v to the Bloom filter and optionally record its (row, chunk) coordinates.

    Triggers a checkpoint persist to MinIO every CHECKPOINT_EVERY additions.

    Args:
        v: Molecule SMILES string to register.
        rc: Optional (row, chunk) tuple to store in the JSON index.
    """
    self.sbf.add(v)
    if rc is not None and self.index is not None and (v not in self.index):
      self.index[v] = list(rc)
    self._added += 1
    if self._added >= CHECKPOINT_EVERY:
      self.persist()

  def persist(self, retries=3):
    """Save the Bloom filter and JSON index to MinIO with atomic local writes.

    Each file is written to a `.tmp` sibling first, then renamed to avoid
    leaving a corrupt file if the process is interrupted. Retries the upload
    up to `retries` times before raising.

    Args:
        retries: Number of upload attempts before raising RuntimeError.

    Raises:
        RuntimeError: If the upload fails after all retries.
    """
    tb = f"{self.local_bloom}.tmp"
    with open(tb, "wb") as f:
      pickle.dump(self.sbf, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tb, self.local_bloom)
    for attempt in range(1, retries + 1):
      try:
        self.store.upload_file(self.local_bloom, self.bucket, self.bloom_key)
        break
      except Exception as e:
        if attempt < retries:
          logger.warning(f"bloom upload attempt {attempt}/{retries} failed: {e}")
          time.sleep(1)
        else:
          raise RuntimeError(
            f"bloom upload failed after {retries} attempts — aborting to prevent duplicates: {e}"
          ) from e
    if self.index is not None:
      ti = f"{self.local_index}.tmp"
      with open(ti, "w", encoding="utf-8") as f:
        json.dump(self.index, f, separators=(",", ":"), ensure_ascii=False)
      os.replace(ti, self.local_index)
      for attempt in range(1, retries + 1):
        try:
          self.store.upload_file(self.local_index, self.bucket, self.index_key)
          break
        except Exception as e:
          if attempt < retries:
            logger.warning(f"index upload attempt {attempt}/{retries} failed: {e}")
            time.sleep(1)
          else:
            raise RuntimeError(
              f"index upload failed after {retries} attempts: {e}"
            ) from e
    self._added = 0


class LocationSink:
  """On-disk capture of {key -> (chunk, row_group)} as wide chunks are written (Workstream 4).

  Backed by a local SQLite `loc` table with NO index yet — step 2 adds the unique index
  and `meta` marker and uploads it as index.sqlite. Rows are inserted in batches so peak
  RAM stays ~one batch, never the full ~1.35M-key set (the ~300 MB the JSON location dict
  cost — the OOM this workstream removes)."""

  def __init__(self, path, batch=50000):
    self.path = path
    self.batch = batch
    self.count = 0
    self._buf = []
    self.conn = sqlite3.connect(path)
    self.conn.execute("PRAGMA journal_mode=OFF")
    self.conn.execute("PRAGMA synchronous=OFF")
    self.conn.execute("CREATE TABLE IF NOT EXISTS loc (key TEXT, chunk INTEGER, rg INTEGER)")

  def add(self, keys, chunk, rgs):
    """Buffer (key, chunk, rg) triples; keys and rgs are equal-length sequences."""
    self._buf.extend((k, chunk, rg) for k, rg in zip(keys, rgs))
    self.count += len(keys)
    if len(self._buf) >= self.batch:
      self._drain()

  def _drain(self):
    if not self._buf:
      return
    self.conn.executemany("INSERT INTO loc(key, chunk, rg) VALUES (?, ?, ?)", self._buf)
    self.conn.commit()
    self._buf = []

  def finalize_index(self):
    """Turn the raw loc table into the queryable index: a unique key index for fast
    batched lookups + a `meta` marker table (format/granularity/indexed_chunks/entries).
    Commits and closes; the file at self.path is then ready to upload. Returns entries."""
    self._drain()
    self.conn.execute("CREATE UNIQUE INDEX IF NOT EXISTS idx_loc_key ON loc(key)")
    (entries,) = self.conn.execute("SELECT COUNT(*) FROM loc").fetchone()
    chunks = [r[0] for r in self.conn.execute("SELECT DISTINCT chunk FROM loc ORDER BY chunk")]
    self.conn.execute("CREATE TABLE IF NOT EXISTS meta (k TEXT PRIMARY KEY, v TEXT)")
    self.conn.executemany(
      "INSERT OR REPLACE INTO meta(k, v) VALUES (?, ?)",
      [
        ("format", INDEX_FORMAT),
        ("granularity", INDEX_GRANULARITY),
        ("indexed_chunks", json.dumps(chunks)),
        ("entries", str(entries)),
      ],
    )
    self.conn.commit()
    self.conn.close()
    return entries

  def close(self):
    self._drain()
    self.conn.close()


def build_location_index(loc_sink, store, bucket, base, retries=3):
  """Workstream 4 step 2: finalize a LocationSink into index.sqlite and upload it to
  {base}/index.sqlite (a single-object PUT is atomic). No-op when nothing was captured —
  we do NOT overwrite an existing index with an empty one. Returns entries indexed.

  Covers only the keys written in THIS session; an incremental merge that adds to existing
  chunks yields a partial index (correct via the read-side full-scan fallback). Full-rewrite
  paths (precalc, pull, W5) write everything in one session, so the index is complete."""
  if loc_sink is None:
    return 0
  if loc_sink.count == 0:
    loc_sink.close()
    return 0
  entries = loc_sink.finalize_index()
  key = f"{base.strip('/')}/{INDEX_SQLITE_FILE}"
  for attempt in range(1, retries + 1):
    try:
      store.upload_file(loc_sink.path, bucket, key)
      break
    except Exception as e:
      if attempt < retries:
        logger.warning(f"index.sqlite upload attempt {attempt}/{retries} failed: {e}")
        time.sleep(1)
      else:
        raise RuntimeError(f"index.sqlite upload failed after {retries} attempts: {e}") from e
  logger.debug(f"[loc] uploaded index.sqlite entries={entries} -> {bucket}/{key}")
  return entries


def read_index_keys(store, bucket, base, tmpdir):
  """Return the key list from {base}/index.sqlite (a loc-v1 index), or None if there is no
  usable index. Materializes the key list — for admin/merge paths (push, inspect), NOT the
  hot read path (which uses the bounded batched lookup instead)."""
  local = os.path.join(tmpdir, f"idxkeys_{uuid.uuid4().hex}.sqlite")
  try:
    store.download_file(bucket, f"{base.strip('/')}/{INDEX_SQLITE_FILE}", local)
  except Exception:
    return None
  try:
    conn = sqlite3.connect(f"file:{local}?mode=ro", uri=True)
    try:
      fmt = conn.execute("SELECT v FROM meta WHERE k='format'").fetchone()
      if not fmt or fmt[0] != INDEX_FORMAT:
        return None
      return [r[0] for r in conn.execute("SELECT key FROM loc")]
    finally:
      conn.close()
  except Exception:
    return None
  finally:
    try:
      os.remove(local)
    except Exception:
      pass


def scan_chunk_locations(store, bucket, base, tmpdir, key_cols=_INPUT_C):
  """Rebuild-from-data counterpart to write-time capture (Option B): read the input/smiles
  column back from every chunk_*.parquet under {base}/data/, one row group at a time (no
  output decode), into a fresh LocationSink recording exact (key, chunk, row_group). Used to
  rebuild the index after a delete rewrites chunks, and as the offline backfill for legacy
  models. Bounded RAM — one chunk file + one string column at a time. Returns the sink."""
  sink = LocationSink(os.path.join(tmpdir, f"rescan_{uuid.uuid4().hex}.sqlite"))
  prefix = f"{hive_prefix(base)}/"
  keys = sorted(
    k["Key"]
    for k in store.list_keys(bucket, prefix)
    if k["Key"].endswith(".parquet") and "/chunk_" in k["Key"]
  )
  for okey in keys:
    try:
      cidx = int(os.path.basename(okey).split("_")[1].split(".")[0])
    except Exception:
      continue
    local = os.path.join(tmpdir, f"rescan_{uuid.uuid4().hex}.parquet")
    try:
      store.download_file(bucket, okey, local)
      pf = pq.ParquetFile(local)
      kc = next((c for c in key_cols if c in pf.schema_arrow.names), None)
      if kc is None:
        continue
      for rg in range(pf.metadata.num_row_groups):
        col = pf.read_row_group(rg, columns=[kc]).column(kc).to_pylist()
        sink.add(col, cidx, [rg] * len(col))
    except Exception as e:
      logger.warning(f"[rescan] skip {okey}: {e}")
    finally:
      try:
        os.remove(local)
      except Exception:
        pass
  return sink


class ChunkState:
  """Manages chunked Parquet writes to MinIO for a single model.

  Rows are written into sequentially numbered chunk files
  (chunk_1.parquet, chunk_2.parquet, ...) stored under a Hive-style prefix
  on MinIO. The last chunk file stays "open" (accepting appends) until it
  reaches max_rows, at which point a new chunk is started.

  Each flush downloads the open chunk if it exists, appends the new rows,
  and re-uploads it. Full chunks are uploaded once and never touched again.

  Args:
      store: MinioStore used for uploads and downloads.
      bucket: Target bucket.
      base_prefix: Key prefix for this model (e.g. "eos1234/v1").
      tmpdir: Local directory for temporary files during writes.
      max_rows: Maximum number of rows per chunk file.
      output_dimension: Number of output columns. Used to tune Parquet
          writer settings (compression, row group size) for wide tables.
  """

  def __init__(self, store, bucket, base_prefix, tmpdir, max_rows, output_dimension=None):
    self.store = store
    self.bucket = bucket
    self.base = base_prefix.strip("/")
    self.tmpdir = tmpdir
    self.max_rows = max_rows
    self.wide = is_wide(output_dimension)
    self.write_batch_rows = chunk_write_batch_rows(output_dimension)
    self.pq_ctor_kw, self.pq_wt_kw = parquet_writer_kwargs(output_dimension)
    # {column_name: arrow type} from run_columns.csv; set by the writer once the
    # schema is known. Empty => fall back to per-column type inference.
    self.column_types = {}
    # W4 location capture (wide-only): the writer attaches a LocationSink and rows'
    # (chunk, row_group) are recorded as they are written. key_col is resolved from
    # the schema on first capture. None => no capture (narrow / not wired).
    self.loc_sink = None
    self.key_col = None
    self.state = {}
    # Wide byte-sized accumulation (Workstream 3): file/row-group sizes are
    # derived from the first typed table (_calibrate); rows stream into one
    # ParquetWriter kept open on local disk across flushes and uploaded at
    # file_rows, with the trailing partial flushed by finalize_chunks().
    self._ensured = False
    self._calibrated = False
    self._open_writer = None
    self._open_path = None
    self._open_schema = None
    self._open_idx = None
    self._open_rows = 0
    # Rows accumulated across flushes but not yet cut into a row group; coalesced
    # to whole row groups (write_batch_rows) so the byte-sized row-group target is
    # honored regardless of how small each incoming flush is.
    self._pending = []
    self._pending_rows = 0
    self._pending_schema = None

  def _rows_to_frame(self, rows, schema_cols):
    """Convert a list of row dicts to a DataFrame with exactly schema_cols columns."""
    df = pd.DataFrame.from_records(rows)
    for col in schema_cols:
      if col not in df.columns:
        df[col] = None
    return df[schema_cols]

  def _normalize_scalar(self, value):
    """Return None for NaN/None values, otherwise return the value unchanged."""
    if value is None:
      return None
    try:
      if pd.isna(value):
        return None
    except Exception:
      pass
    return value

  def _stringify_scalar(self, value):
    """Convert a value to a string, encoding bytes and serialising collections to JSON."""
    if value is None:
      return None
    if isinstance(value, bytes):
      return value.decode("utf-8", errors="replace")
    if isinstance(value, (list, dict)):
      return json.dumps(value, ensure_ascii=False)
    if isinstance(value, (tuple, set)):
      return json.dumps(list(value), ensure_ascii=False)
    return str(value)

  def _build_array(self, values, col=None):
    """Build a PyArrow array from a list of Python values.

    When the column has a declared type (from run_columns.csv via the writer),
    build that exact type and hard-fail on bad data. Otherwise fall back to native
    type inference, then to string for heterogeneous/non-serialisable columns.
    """
    target = self.column_types.get(col) if col is not None else None
    if target is not None:
      return build_typed_array(values, target)
    normalized = [self._normalize_scalar(v) for v in values]
    try:
      return pa.array(normalized, from_pandas=True)
    except Exception:
      return pa.array([self._stringify_scalar(v) for v in normalized], type=pa.string())

  def _frame_to_table(self, df, schema_cols):
    """Convert a DataFrame slice to a PyArrow Table with schema_cols columns."""
    df = self._ensure_cols(df, schema_cols)
    arrays = [self._build_array(df[col].tolist(), col) for col in schema_cols]
    return pa.Table.from_arrays(arrays, names=schema_cols)

  def _rows_to_table(self, rows, schema_cols):
    """Convert a list of row dicts to a PyArrow Table with schema_cols columns."""
    cols = {col: [row.get(col) if isinstance(row, dict) else None for row in rows] for col in schema_cols}
    arrays = [self._build_array(cols[col], col) for col in schema_cols]
    return pa.Table.from_arrays(arrays, names=schema_cols)

  def _ensure_cols(self, df, schema_cols):
    """Add any missing schema columns as None and return df restricted to schema_cols."""
    for col in schema_cols:
      if col not in df.columns:
        df[col] = None
    return df[schema_cols]

  def _iter_tables(self, rows, schema_cols):
    """Yield PyArrow Tables from a row list in write_batch_rows-sized slices."""
    for start in range(0, len(rows), self.write_batch_rows):
      yield self._rows_to_table(rows[start : start + self.write_batch_rows], schema_cols)

  def _iter_tables_df(self, df, schema_cols):
    """Yield PyArrow Tables from a DataFrame in write_batch_rows-sized slices."""
    df = self._ensure_cols(df, schema_cols)
    for start in range(0, len(df), self.write_batch_rows):
      frame = df.iloc[start : start + self.write_batch_rows]
      yield self._frame_to_table(frame, schema_cols)

  def _rows_in_remote(self, key):
    """Download a remote Parquet file and return its row count. Returns 0 on failure."""
    local = os.path.join(self.tmpdir, f"inspect_{uuid.uuid4().hex}.parquet")
    try:
      self.store.download_file(self.bucket, key, local)
      try:
        return pq.ParquetFile(local).metadata.num_rows
      except:
        return len(pd.read_parquet(local))
    except:
      return 0
    finally:
      try:
        os.remove(local)
      except:
        pass

  def ensure(self):
    """Initialise chunk state by inspecting existing chunks on MinIO.

    Sets self.state["data"] with the next chunk index to write, whether
    there is currently an open (not-yet-full) chunk, and its row count.
    Must be called before any flush.
    """
    keys, t = (self._list_chunks(), "data")
    if not keys:
      self.state[t] = {"next": 1, "open": None, "rows": 0}
      logger.debug("chunk new: next=1")
      return
    last = keys[-1]
    idx = self._chunk_idx(last)
    self.state[t] = {"next": idx + 1, "open": None, "rows": 0}
    logger.debug(f"chunk next: {idx + 1}")

  def _list_chunks(self):
    """Return a sorted list of Parquet chunk keys for this model from MinIO."""
    pref = hive_prefix(self.base) + "/"
    keys = []
    for obj in self.store.list_keys(self.bucket, pref):
      k = obj["Key"]
      if k.endswith(".parquet") and "/chunk_" in k:
        keys.append(k)
    return sorted(keys, key=self._chunk_idx)

  def _chunk_idx(self, key):
    """Parse and return the integer index from a chunk filename (e.g. chunk_3.parquet → 3)."""
    base = os.path.basename(key)
    name, _ = os.path.splitext(base)
    try:
      return int(name.split("_")[1])
    except:
      return 1

  def _write_chunk(self, rows, schema_cols, idx, mode="new", existing_local=None):
    """Write a list of row dicts to chunk_{idx}.parquet and upload to MinIO.

    Args:
        rows: Row dicts to write.
        schema_cols: Column order for the Parquet schema.
        idx: Chunk index (determines the filename).
        mode: "new" to create a fresh file, "append" to merge into existing_local.
        existing_local: Path to the already-downloaded open chunk (append mode only).

    Returns:
        The MinIO object key of the written chunk.
    """
    os_key = f"{hive_prefix(self.base)}/chunk_{idx}.parquet"
    local = existing_local or os.path.join(self.tmpdir, f"chunk_{uuid.uuid4().hex}.parquet")
    ctor_kw, wt_kw = self.pq_ctor_kw, self.pq_wt_kw
    if mode == "append" and existing_local:
      tmp = local + ".merged"
      old_pf = pq.ParquetFile(existing_local)
      schema = old_pf.schema_arrow
      with pq.ParquetWriter(tmp, schema, **ctor_kw) as writer:
        for batch in old_pf.iter_batches(batch_size=self.write_batch_rows):
          writer.write_batch(batch)
        for table in self._iter_tables(rows, schema_cols):
          writer.write_table(table.cast(schema), **wt_kw)
      os.replace(tmp, local)
    else:
      tmp = local + ".tmp"
      writer = None
      try:
        for table in self._iter_tables(rows, schema_cols):
          if writer is None:
            writer = pq.ParquetWriter(tmp, table.schema, **ctor_kw)
          writer.write_table(table, **wt_kw)
      finally:
        if writer is not None:
          writer.close()
      os.replace(tmp, local)
    self.store.upload_file(local, self.bucket, os_key)
    if not existing_local:
      try:
        os.remove(local)
      except:
        pass
    return os_key

  def _write_chunk_df(self, df, schema_cols, idx, mode="new", existing_local=None):
    """Write a DataFrame to chunk_{idx}.parquet and upload to MinIO.

    Same as _write_chunk but accepts a DataFrame instead of a list of row dicts.

    Args:
        df: DataFrame to write.
        schema_cols: Column order for the Parquet schema.
        idx: Chunk index (determines the filename).
        mode: "new" to create a fresh file, "append" to merge into existing_local.
        existing_local: Path to the already-downloaded open chunk (append mode only).

    Returns:
        The MinIO object key of the written chunk.
    """
    os_key = f"{hive_prefix(self.base)}/chunk_{idx}.parquet"
    local = existing_local or os.path.join(self.tmpdir, f"chunk_{uuid.uuid4().hex}.parquet")
    ctor_kw, wt_kw = self.pq_ctor_kw, self.pq_wt_kw
    if mode == "append" and existing_local:
      tmp = local + ".merged"
      old_pf = pq.ParquetFile(existing_local)
      schema = old_pf.schema_arrow
      with pq.ParquetWriter(tmp, schema, **ctor_kw) as writer:
        for batch in old_pf.iter_batches(batch_size=self.write_batch_rows):
          writer.write_batch(batch)
        for table in self._iter_tables_df(df, schema_cols):
          writer.write_table(table.cast(schema), **wt_kw)
      os.replace(tmp, local)
    else:
      tmp = local + ".tmp"
      writer = None
      try:
        for table in self._iter_tables_df(df, schema_cols):
          if writer is None:
            writer = pq.ParquetWriter(tmp, table.schema, **ctor_kw)
          writer.write_table(table, **wt_kw)
      finally:
        if writer is not None:
          writer.close()
      os.replace(tmp, local)
    self.store.upload_file(local, self.bucket, os_key)
    if not existing_local:
      try:
        os.remove(local)
      except:
        pass
    return os_key

  def _ensure_once(self):
    """Initialise chunk state from MinIO on the first wide flush only.

    Unlike the narrow path (which re-lists on every flush), the wide path keeps
    a writer open across flushes, so it must inspect the existing chunks exactly
    once to learn the next chunk index.
    """
    if not self._ensured:
      self.ensure()
      self._ensured = True

  def _calibrate(self, table):
    """Size the chunk file (self.max_rows) and row group from the table's measured
    bytes/row, so both are bounded in bytes regardless of width/dtype."""
    self._calibrated = True
    n = table.num_rows
    if n <= 0:
      return
    bytes_per_row = max(1.0, table.nbytes / n)
    file_rows = max(1, int(WIDE_TARGET_FILE_BYTES / bytes_per_row))
    rg_rows = max(1, int(WIDE_TARGET_ROWGROUP_BYTES / bytes_per_row))
    rg_rows = min(rg_rows, file_rows)
    file_rows = max(rg_rows, (file_rows // rg_rows) * rg_rows)  # whole row groups per file
    self.max_rows = file_rows
    self.write_batch_rows = rg_rows
    self.pq_wt_kw["row_group_size"] = rg_rows
    logger.debug(
      f"[chunk] calibrated bytes/row={bytes_per_row:.0f} file_rows={file_rows} rg_rows={rg_rows}"
    )

  def _open_new_writer(self, schema):
    """Open a fresh ParquetWriter on local disk for the next chunk index."""
    self._open_idx = self.state["data"]["next"]
    self._open_path = os.path.join(self.tmpdir, f"chunk_{uuid.uuid4().hex}.parquet")
    self._open_schema = schema
    self._open_writer = pq.ParquetWriter(self._open_path, schema, **self.pq_ctor_kw)
    self._open_rows = 0

  def _close_and_upload_open(self):
    """Close the open chunk writer, upload it to MinIO, and advance the index."""
    if self._open_writer is None:
      return
    self._open_writer.close()
    self._open_writer = None
    os_key = f"{hive_prefix(self.base)}/chunk_{self._open_idx}.parquet"
    self.store.upload_file(self._open_path, self.bucket, os_key)
    logger.debug(f"[chunk] uploaded chunk idx={self._open_idx} rows={self._open_rows}")
    try:
      os.remove(self._open_path)
    except Exception:
      pass
    self.state["data"]["next"] = self._open_idx + 1
    self._open_path = None
    self._open_idx = None
    self._open_rows = 0

  def _flush_accum(self, tables):
    """Accumulate typed tables and cut them into whole row groups (write_batch_rows) as
    they fill, streaming into the open chunk and rolling at file_rows. This makes the
    byte-sized row-group target real even when each incoming flush is small; pyarrow
    never merges across write_table calls, so we coalesce here instead. The trailing
    partial row group + chunk are flushed by finalize_chunks(); peak RAM ~= one row group."""
    for table in tables:
      if table.num_rows == 0:
        continue
      if not self._calibrated:
        self._calibrate(table)
      if self._pending_schema is None:
        self._pending_schema = table.schema
      self._pending.append(table.cast(self._pending_schema))
      self._pending_rows += table.num_rows
      if self._pending_rows >= self.write_batch_rows:
        self._drain_pending(final=False)

  def _drain_pending(self, final):
    """Write accumulated rows into the open chunk in whole row groups. Mid-run writes
    only complete row groups (write_batch_rows); the trailing partial is written when final."""
    if not self._pending_rows:
      return
    rg = self.write_batch_rows
    combined = self._pending[0] if len(self._pending) == 1 else pa.concat_tables(self._pending)
    n_write = combined.num_rows if final else (combined.num_rows // rg) * rg
    if n_write <= 0:
      return
    self._write_to_open(combined.slice(0, n_write), final=final)
    remainder = combined.slice(n_write)
    self._pending = [remainder.combine_chunks()] if remainder.num_rows else []
    self._pending_rows = remainder.num_rows

  def _write_to_open(self, table, final=False):
    """Append a table to the open chunk, rolling to a new chunk at file_rows; row groups
    are sized by pq_wt_kw['row_group_size'] (= write_batch_rows)."""
    n = table.num_rows
    offset = 0
    while offset < n:
      if self._open_writer is None:
        self._open_new_writer(table.schema)
      take = min(self.max_rows - self._open_rows, n - offset)
      sub = table.slice(offset, take)
      self._open_writer.write_table(sub.cast(self._open_schema), **self.pq_wt_kw)
      if self.loc_sink is not None:
        self._capture_locs(sub, self._open_idx, self._open_rows, final)
      self._open_rows += take
      offset += take
      if self._open_rows >= self.max_rows:
        self._close_and_upload_open()

  def _capture_locs(self, table, chunk_idx, start_row, final):
    """Record (key, chunk_idx, row_group) for each row in `table` into the location sink.

    Rows occupy [start_row, start_row + n) within chunk_idx; row_group = row // rg_rows.
    The W3 write path only ever cuts whole row groups (write_batch_rows) except the final
    trailing partial, so mid-run writes MUST be row-group-aligned for the ordinal to be
    exact. Assert it: a mislocated key is not caught by the read-side full-scan fallback
    (that only fires on keys fully absent), so it would return a silent blank row."""
    rg = self.write_batch_rows
    n = table.num_rows
    if not final and (start_row % rg or n % rg):
      raise RuntimeError(
        f"[loc] unaligned wide write (start_row={start_row} n={n} rg={rg}); "
        "row-group index would mislocate keys — aborting"
      )
    if self.key_col is None:
      self.key_col = next((c for c in _INPUT_C if c in table.schema.names), None)
    if self.key_col is None:
      return
    keys = table.column(self.key_col).to_pylist()
    rgs = [(start_row + i) // rg for i in range(n)]
    self.loc_sink.add(keys, chunk_idx, rgs)

  def finalize_chunks(self):
    """Flush the trailing partial row group and open chunk. Call once after all flushes.

    No-op for narrow models (which upload every chunk inside flush/flush_df).
    """
    if not self.wide:
      return
    if self._pending_rows:
      self._drain_pending(final=True)
    if self._open_writer is not None and self._open_rows > 0:
      self._close_and_upload_open()

  def flush(self, rows, schema_cols):
    """Write a list of row dicts to MinIO as one or more Parquet chunk files.

    Appends to the currently open chunk if space remains, then opens new
    chunks as needed until all rows are written. Each completed chunk is
    uploaded to MinIO immediately.

    Args:
        rows: List of row dicts to write.
        schema_cols: Ordered column names for the Parquet schema.
    """
    if self.wide:
      self._ensure_once()
      if rows:
        self._flush_accum(self._iter_tables(rows, schema_cols))
      return
    if not rows:
      return self.ensure()
    self.ensure()
    st = self.state["data"]
    remaining = len(rows)
    start = 0
    logger.debug(f"flush: chunk rows={remaining} cols={len(schema_cols or [])} limit={self.max_rows}")
    if st["open"]:
      tmp = os.path.join(self.tmpdir, f"open_{uuid.uuid4().hex}.parquet")
      try:
        self.store.download_file(self.bucket, st["open"], tmp)
        space = self.max_rows - st["rows"]
        take = min(space, remaining)
        if take > 0:
          part = rows[start : start + take]
          self._write_chunk(part, schema_cols, st["next"], mode="append", existing_local=tmp)
          st["rows"] += take
          remaining -= take
          start += take
          logger.debug(f"flush: appended chunk idx={st['next']} +{take} -> {st['rows']}")
        if st["rows"] >= self.max_rows:
          st["next"] += 1
          st["open"] = None
          st["rows"] = 0
          logger.debug(f"flush: closed chunk next={st['next']}")
      finally:
        try:
          os.remove(tmp)
        except:
          pass
    while remaining > 0:
      take = min(self.max_rows, remaining)
      part = rows[start : start + take]
      os_key = self._write_chunk(part, schema_cols, st["next"], mode="new")
      if take < self.max_rows:
        st["open"] = os_key
        st["rows"] = take
        logger.debug(f"flush: new open chunk idx={st['next']} rows={take}")
      else:
        st["next"] += 1
        st["open"] = None
        st["rows"] = 0
        logger.debug(f"flush: full chunk idx={st['next'] - 1} rows={take}")
      remaining -= take
      start += take

  def flush_df(self, df, schema_cols):
    """Write a DataFrame to MinIO as one or more Parquet chunk files.

    Same as flush but accepts a DataFrame instead of a list of row dicts.

    Args:
        df: DataFrame to write.
        schema_cols: Ordered column names for the Parquet schema.
    """
    if self.wide:
      self._ensure_once()
      if df is not None and not df.empty:
        self._flush_accum(self._iter_tables_df(df, schema_cols))
      return
    if df is None or df.empty:
      return self.ensure()
    self.ensure()
    st = self.state["data"]
    remaining = len(df)
    start = 0
    logger.debug(f"flush_df: chunk rows={remaining} cols={len(schema_cols or [])} limit={self.max_rows}")
    if st["open"]:
      tmp = os.path.join(self.tmpdir, f"open_{uuid.uuid4().hex}.parquet")
      try:
        self.store.download_file(self.bucket, st["open"], tmp)
        space = self.max_rows - st["rows"]
        take = min(space, remaining)
        if take > 0:
          part = df.iloc[start : start + take]
          self._write_chunk_df(part, schema_cols, st["next"], mode="append", existing_local=tmp)
          st["rows"] += take
          remaining -= take
          start += take
          logger.debug(f"flush_df: appended chunk idx={st['next']} +{take} -> {st['rows']}")
        if st["rows"] >= self.max_rows:
          st["next"] += 1
          st["open"] = None
          st["rows"] = 0
      finally:
        try:
          os.remove(tmp)
        except:
          pass
    while remaining > 0:
      take = min(self.max_rows, remaining)
      part = df.iloc[start : start + take]
      os_key = self._write_chunk_df(part, schema_cols, st["next"], mode="new")
      if take < self.max_rows:
        st["open"] = os_key
        st["rows"] = take
      else:
        st["next"] += 1
        st["open"] = None
        st["rows"] = 0
      remaining -= take
      start += take


TrancheState = ChunkState


class _SinkWriter:
  """High-level writer that deduplicates and buffers rows before writing to MinIO.

  Combines a BloomIndex (to skip already-stored molecules) with a ChunkState
  (to manage chunked Parquet writes). Incoming DataFrames are checked against
  the Bloom filter, new rows are buffered in memory, and the buffer is flushed
  to MinIO once it reaches max_rows.

  Args:
      store: MinioStore used for uploads and downloads.
      bucket: Target bucket.
      model_id: Ersilia model identifier (e.g. "eos1234").
      model_version: Model version string.
      tmpdir: Local directory for temporary files.
      max_rows: Buffer size before an automatic flush. Defaults to MAX_ROWS.
      output_dimension: Number of output columns, used to tune Parquet settings.
  """

  def __init__(self, store, bucket, model_id, model_version, tmpdir, max_rows=None, output_dimension=None):
    self.store = store
    self.bucket = bucket
    self.model_id = model_id
    self.model_version = model_version
    self.base = get_base(model_id, model_version)
    self.tmpdir = tmpdir
    self.max_rows = int(max_rows or MAX_ROWS)
    self.store.ensure_bucket(self.bucket)
    # W4: wide models track molecules via index.sqlite + bloom, not the JSON index — skip
    # loading/writing it (saves ~300 MB RAM on wide writes/pulls). Narrow keeps the JSON index.
    self.bi = BloomIndex(self.store, self.bucket, self.base, tmpdir, load_index=not is_wide(output_dimension))
    self.chunk_state = ChunkState(
      self.store, self.bucket, self.base, tmpdir, self.max_rows, output_dimension=output_dimension
    )
    # W4: capture real (chunk, row_group) locations during wide writes (step 2 builds
    # index.sqlite from this). Wide-only; narrow path untouched.
    if self.chunk_state.wide:
      self.chunk_state.loc_sink = LocationSink(os.path.join(tmpdir, "loc_index.sqlite"))
    self.buffers = []
    self.buf_rows = 0
    self.schema_cols = None
    self.last_added_inputs = []

  def add_rows(self, df):
    """Filter out duplicates and blank inputs, buffer new rows, and auto-flush if full.

    Args:
        df: DataFrame with an "input" or "smiles" column identifying each molecule.

    Returns:
        Number of new rows added (after deduplication).
    """
    try:
      if df is None or df.empty:
        logger.debug("[sink] empty batch")
        self.last_added_inputs = []
        return 0
      n_in = len(df)
      if self.schema_cols is None:
        self.schema_cols = list(df.columns)
        # Resolve declared column types from run_columns.csv (enforces the contract
        # and self-heals string-typed numeric data on pull write-back).
        self.chunk_state.column_types = resolve_write_types(self.model_id, self.schema_cols)
      if "input" in df.columns:
        inputs = df["input"].astype(str).str.strip()
      elif "smiles" in df.columns:
        inputs = df["smiles"].astype(str).str.strip()
      else:
        inputs = pd.Series([""] * len(df), index=df.index, dtype="object")
      non_blank = inputs != ""
      sbf = self.bi.sbf
      not_seen = non_blank & ~inputs.map(lambda v: v in sbf)
      skipped_blank = int((~non_blank).sum())
      skipped_seen = int(non_blank.sum() - not_seen.sum())
      new_df = df.loc[not_seen]
      added = len(new_df)
      self.last_added_inputs = []
      if added == 0:
        logger.debug(f"[sink] in={n_in} added=0 seen={skipped_seen} blank={skipped_blank}")
        return 0
      self.last_added_inputs = inputs.loc[not_seen].tolist()
      for smi in inputs.loc[not_seen]:
        self.bi.register(smi, rc=(1, 1))
      self.buffers.append(new_df)
      self.buf_rows += added
      logger.debug(
        f"[sink] in={n_in} added={added} seen={skipped_seen} blank={skipped_blank} buf={self.buf_rows}/{self.max_rows}"
      )
      if self.buf_rows >= self.max_rows:
        self._flush()
      return added
    except RuntimeError:
      raise  # column-contract violation must surface, never be silently swallowed
    except Exception as e:
      logger.error(f"[sink] error: {e}")
      return 0

  def _flush(self):
    """Concatenate buffered DataFrames and write them to MinIO via ChunkState."""
    parts = self.buffers
    if not parts:
      return
    merged = pd.concat(parts, ignore_index=True)
    logger.debug(f"[sink] flushing {len(merged)} rows")
    self.chunk_state.flush_df(merged, self.schema_cols)
    self.buffers = []
    self.buf_rows = 0
    del merged
    gc.collect()

  def finalize(self, metadata_local=None, schema_cols=None):
    """Flush any remaining buffered rows, persist the Bloom filter and index, and upload catalog.json.

    Should be called once after all add_rows() calls are complete.

    Args:
        metadata_local: Optional local path to an access metadata file to upload alongside the data.
        schema_cols: Column schema to use if none was inferred during add_rows().
    """
    schema_cols = schema_cols or self.schema_cols
    if self.buffers:
      merged = pd.concat(self.buffers, ignore_index=True)
      self.chunk_state.flush_df(merged, schema_cols)
      self.buffers = []
      self.buf_rows = 0
      del merged
    # Upload the trailing partial wide chunk held open across flushes (no-op for narrow).
    self.chunk_state.finalize_chunks()
    # W4 step 2: build + upload index.sqlite from the captured locations (wide-only).
    build_location_index(self.chunk_state.loc_sink, self.store, self.bucket, self.base)
    gc.collect()
    self.bi.persist()
    if metadata_local:
      try:
        self.store.upload_file(metadata_local, self.bucket, f"{self.base}/{ACCESS_FILE}")
      except Exception:
        pass
    entries = len(self.bi.index) if self.bi.index is not None else len(self.bi.sbf)
    chunks = len(self.chunk_state._list_chunks())
    upload_catalog_json(self.store, self.bucket, self.base, entries, chunks, self.tmpdir)


def upload_catalog_json(store, bucket, base_prefix, entries, chunks, tmpdir):
  """Write a catalog.json with entry and chunk counts to MinIO.

  catalog.json is a small summary file used by the catalog and stats commands
  to report model sizes without scanning all Parquet files.

  Args:
      store: MinioStore used for the upload.
      bucket: Target bucket.
      base_prefix: Key prefix for this model.
      entries: Total number of stored molecule entries.
      chunks: Total number of Parquet chunk files.
      tmpdir: Local directory to write the temporary JSON file.
  """
  cat = {"entries": int(entries), "chunks": int(chunks)}
  local = os.path.join(tmpdir, f"catalog_{uuid.uuid4().hex}.json")
  try:
    with open(local, "w", encoding="utf-8") as f:
      json.dump(cat, f, separators=(",", ":"))
    key = f"{base_prefix.strip('/')}/{CATALOG_FILE}"
    store.upload_file(local, bucket, key)
    logger.debug(f"[catalog] {CATALOG_FILE} -> {bucket}/{key} entries={entries} chunks={chunks}")
  except Exception as e:
    logger.warning(f"catalog.json upload failed: {e}")
  finally:
    try:
      os.remove(local)
    except Exception:
      pass


class _BaseTransfer:
  """Base class for all isaura manager operations (read, write, copy, move, push, pull).

  Sets up the shared infrastructure each operation needs: a MinioStore, a
  DuckDB connection, a BloomIndex for the target model, and two temporary
  directories (one for general use, one for the SinkWriter). Implements
  context-manager protocol so resources are cleaned up with `with` blocks.

  Subclasses in manage.py implement specific operations by calling the
  helper methods here (select_rows, _pull, _copy, _dump, etc.).

  Args:
      model_id: Ersilia model identifier (e.g. "eos1234").
      model_version: Model version string.
      bucket: MinIO bucket to operate on.
      output_dir: Local directory for dump operations (copy to local filesystem).
  """

  def __init__(self, model_id, model_version, bucket, output_dir=None):
    self.model_id = model_id
    self.model_version = model_version
    self.bucket = bucket
    self.output_dir = output_dir
    self.base = get_pref(model_id, model_version)
    self.tranches = get_base(model_id, model_version)
    self.collection = get_coll(self.model_id, self.model_version)
    self.store = MinioStore()
    self.store.require_bucket(self.bucket)
    self.tmpdir = make_temp("isaura_xfer_")
    self.tmpdir_sinkw = make_temp("isaura_sinkw_")
    bind_temp_dirs(self, self.tmpdir, self.tmpdir_sinkw)
    self.duck = DuckDBMinio(endpoint=self.store.endpoint, access=self.store.access, secret=self.store.secret)
    self.bi = BloomIndex(
      self.store,
      self.bucket,
      self.tranches,
      self.tmpdir,
      bloom_filename=os.getenv("BLOOM_FILENAME", BLOOM_FILENAME),
    )
    self._model_meta = None
    self._output_dimension = None

  def close(self):
    """Release temporary directories created during construction."""
    release_temp_dirs(self)

  def __enter__(self):
    return self

  def __exit__(self, et, ev, tb):
    self.close()

  def _get_model_meta(self):
    """Fetch and cache model schema metadata from GitHub. Returns {} on failure."""
    if self._model_meta is None:
      try:
        self._model_meta = fetch_schema_from_github(self.model_id)
      except Exception as e:
        logger.warning(f"[schema] metadata fetch failed for {self.model_id}: {e}")
        self._model_meta = {}
    return self._model_meta

  def _get_output_dimension(self):
    """Return the number of output columns for this model, cached after first call."""
    if self._output_dimension is None:
      self._output_dimension = output_dimension_from_metadata(self._get_model_meta())
    return self._output_dimension

  def _chunk_row_limit(self):
    """Return the max rows per chunk file, adjusted for output width."""
    return chunk_row_limit(self._get_output_dimension())

  def _download_if_exists(self, key, local):
    """Download key to local path. Exits the process if the download fails."""
    try:
      self.store.download_file(self.bucket, key, local)
      return True
    except Exception as e:
      logger.error(f"Either model id or version or project name are maybe specified wrongly. Results -> {e}")
      sys.exit(1)
      return False

  def _load_metadata(self):
    """Download and parse the access metadata file for this model.

    Returns:
        Tuple of (local_path, metadata_dict), or None on failure.
    """
    local = os.path.join(self.tmpdir, ACCESS_FILE)
    acc = get_acc_key(self.tranches)
    try:
      if self._download_if_exists(acc, local):
        with open(local, "r", encoding="utf-8") as f:
          return (local, json.load(f))
    except Exception as e:
      logger.error(f"{e}. Possible causes -> 1. No project 2. Incorrect model id or version")

  def _load_index(self):
    """Download and parse index.json for this model.

    Returns:
        Tuple of (index_dict, local_path).

    Raises:
        RuntimeError: If index.json is not found on MinIO.
    """
    key = get_idx_key(self.tranches)
    local = os.path.join(self.tmpdir, "index_src.json")
    if not self._download_if_exists(key, local):
      raise RuntimeError("index.json not found")
    with open(local, "r", encoding="utf-8") as f:
      return (json.load(f), local)

  def select_rows(self, wanted, input_col="input"):
    """Return a DataFrame of all rows matching the wanted input values.

    Queries all Parquet chunk files via DuckDB, preserving input order.

    Args:
        wanted: List of molecule SMILES (or other input values) to look up.
        input_col: Column name to filter on (default "input").

    Returns:
        DataFrame of matching rows, ordered by wanted.
    """
    files = list_parquet_keys(self.store, self.bucket, self.tranches)
    return query(self.duck.con, input_col, wanted, files, preserve_order=True)

  def select_rows_batched(self, wanted, input_col="input", batch_size=10000):
    """Yield DataFrames of matching rows in batches, switching strategy based on query size.

    For large queries (>= STREAM_PARQUET_THRESHOLD inputs) uses incremental
    Parquet streaming to avoid loading all chunk files into DuckDB at once.
    For smaller queries uses DuckDB directly.

    Args:
        wanted: List of molecule SMILES to look up.
        input_col: Column name to filter on (default "input").
        batch_size: Number of rows per yielded DataFrame.

    Yields:
        DataFrames of matching rows.
    """
    n = len(wanted) if hasattr(wanted, "__len__") else None
    logger.debug(
      f"[select_rows_batched] starting n={n} input_col={input_col} batch_size={batch_size} bucket={self.bucket}"
    )
    source = None
    try:
      if n is not None and n >= STREAM_PARQUET_THRESHOLD:
        prefix = hive_prefix(self.tranches) + "/"
        logger.debug(f"[select_rows_batched] streaming mode n={n} bucket={self.bucket} prefix={prefix}")
        source = stream_parquet_filtered(
          self.store, self.bucket, prefix, wanted, header=input_col, batch_size=batch_size
        )
      else:
        files = list_parquet_keys(self.store, self.bucket, self.tranches)
        logger.debug(
          f"[select_rows_batched] duckdb mode n={n} bucket={self.bucket} files={(len(files) if isinstance(files, list) else files)}"
        )
        source = query_batched(
          self.duck.con, input_col, wanted, files, batch_size=batch_size, tmpdir=self.tmpdir
        )
      for chunk in source:
        yield chunk
    finally:
      close = getattr(source, "close", None)
      if close is not None:
        close()

  def _delete(self):
    """Delete all objects for this model from the bucket. Returns the count deleted."""
    return self.store.delete_prefix(self.bucket, self.tranches)

  def _pull(self, df, index):
    """Write a single DataFrame to the local bucket and register inputs with the NNS service.

    Args:
        df: DataFrame of model outputs with an "input" column.
        index: Unused legacy parameter (kept for API compatibility).

    Returns:
        Tuple of (rows_added, 0).
    """
    t0 = time.time()
    wt = _SinkWriter(
      self.store,
      self.bucket,
      self.model_id,
      self.model_version,
      self.tmpdir,
      max_rows=self._chunk_row_limit(),
      output_dimension=self._get_output_dimension(),
    )
    tp, tu = (0, 0)
    apprx_inputs = []
    _tp = wt.add_rows(df)
    if _tp != 0:
      mol_col = next((c for c in _INPUT_C if c in df.columns), df.columns[0])
      apprx_inputs.extend(df[mol_col].astype(str).tolist())
    tp += _tp
    if wt:
      wt.finalize(schema_cols=list(df.keys()))
    if apprx_inputs:
      apprx_df = pd.DataFrame({"input": apprx_inputs})
      post_apprx(apprx_df, self.collection)
      del apprx_df
    logger.debug(f"[pull] done rows={tp} elapsed={time.time() - t0:.1f}s rss={rss_mb():.0f}MB")
    return (tp, tu)

  def _pull_batched(self, chunks):
    """Write an iterable of DataFrames to the local bucket and register with the NNS service.

    Args:
        chunks: Iterable of DataFrames, each with an "input" column.

    Returns:
        Tuple of (total_rows_added, 0).
    """
    t0 = time.time()
    wt = _SinkWriter(
      self.store,
      self.bucket,
      self.model_id,
      self.model_version,
      self.tmpdir,
      max_rows=self._chunk_row_limit(),
      output_dimension=self._get_output_dimension(),
    )
    tp = 0
    apprx_inputs = []
    schema_cols = None
    mol_col = None
    for chunk in chunks:
      if chunk is None or chunk.empty:
        continue
      if schema_cols is None:
        schema_cols = list(chunk.columns)
        mol_col = next((c for c in _INPUT_C if c in chunk.columns), schema_cols[0])
      added = wt.add_rows(chunk)
      tp += added
      if added:
        apprx_inputs.extend(chunk[mol_col].astype(str).tolist())
    with console.status("Storing to local MinIO...", spinner="dots"):
      wt.finalize(schema_cols=schema_cols)
    actual = len(wt.bi.index) if wt.bi.index is not None else len(wt.bi.sbf)
    if apprx_inputs:
      try:
        post_apprx(pd.DataFrame({"input": apprx_inputs}), self.collection)
      except Exception as e:
        logger.warning(f"[pull] post_apprx failed: {e}")
    logger.debug(
      f"[pull] done rows={actual} chunk_limit={self._chunk_row_limit()} output_dim={self._get_output_dimension()} elapsed={time.time() - t0:.1f}s rss={rss_mb():.0f}MB"
    )
    return (actual, 0)

  def _dump(self):
    """Download every object in the bucket to self.output_dir, preserving key paths."""
    os.makedirs(self.output_dir, exist_ok=True)
    paginator = self.store.client.get_paginator("list_objects_v2")
    pages = paginator.paginate(Bucket=self.bucket)
    found = False
    for page in pages:
      for obj in page.get("Contents", []):
        found = True
        key = obj["Key"]
        local_path = os.path.join(self.output_dir, key)
        os.makedirs(os.path.dirname(local_path), exist_ok=True)
        logger.info(f"dumping from minio://{self.bucket}/{key} -> {local_path}")
        self.store.download_file(self.bucket, key, local_path)
    if not found:
      logger.info(f"bucket empty: {self.bucket}")

  def _copy(self, meta_local, meta_list):
    """Route rows from the source bucket to the public and/or private buckets based on access metadata.

    Reads all rows for the wanted inputs, splits them by access level ("private" →
    isaura-private, everything else → isaura-public), and writes each subset via
    a _SinkWriter. Also registers all copied inputs with the NNS service.

    If self.output_dir is set, falls back to _dump() instead (local filesystem copy).

    Args:
        meta_local: Local path to the access metadata file to upload alongside the data.
        meta_list: List of dicts with "input" and "access" keys defining routing.

    Returns:
        Tuple of (rows_written_to_private, rows_written_to_public).
    """
    t0 = time.time()
    if self.output_dir is not None:
      self._dump()
      return
    route = {}
    for d in meta_list:
      inp = (d.get("input") or "").strip()
      if not inp:
        continue
      acc = (d.get("access") or "").lower()
      route[inp] = "priv" if acc == "private" else "pub"
    all_wanted = list(route.keys())
    logger.info(f"[copy] inputs={len(all_wanted)} rss={rss_mb():.0f}MB")
    has_priv = any((v == "priv" for v in route.values()))
    has_pub = any((v == "pub" for v in route.values()))
    w_priv = (
      _SinkWriter(
        self.store,
        priv_bucket,
        self.model_id,
        self.model_version,
        self.tmpdir_sinkw,
        max_rows=self._chunk_row_limit(),
        output_dimension=self._get_output_dimension(),
      )
      if has_priv
      else None
    )
    w_pub = (
      _SinkWriter(
        self.store,
        pub_bucket,
        self.model_id,
        self.model_version,
        self.tmpdir_sinkw,
        max_rows=self._chunk_row_limit(),
        output_dimension=self._get_output_dimension(),
      )
      if has_pub
      else None
    )
    tp, tu, n_batches = (0, 0, 0)
    schema_cols = None
    for df in self.select_rows_batched(all_wanted):
      if df is None or df.empty:
        continue
      n_batches += 1
      if schema_cols is None:
        schema_cols = list(df.columns)
      inputs = df["input"].astype(str).str.strip()
      routes = inputs.map(route.get)
      priv_mask = routes == "priv"
      pub_mask = routes == "pub"
      if w_priv and priv_mask.any():
        tp += w_priv.add_rows(df.loc[priv_mask])
      if w_pub and pub_mask.any():
        tu += w_pub.add_rows(df.loc[pub_mask])
      del df
      if n_batches % 20 == 0:
        gc.collect()
        logger.debug(
          f"[copy] batch={n_batches} priv={tp} pub={tu} elapsed={time.time() - t0:.1f}s rss={rss_mb():.0f}MB"
        )
    if w_priv:
      w_priv.finalize(metadata_local=meta_local, schema_cols=schema_cols)
    if w_pub:
      w_pub.finalize(metadata_local=meta_local, schema_cols=schema_cols)
    total = tp + tu
    if total > 0:
      try:
        apprx_df = pd.DataFrame({"input": all_wanted[:total]})
        post_apprx(apprx_df, self.collection)
        del apprx_df
      except Exception as e:
        logger.warning(f"[copy] post_apprx failed: {e}")
    del route
    gc.collect()
    elapsed = time.time() - t0
    rate = total / elapsed if elapsed > 0 else 0
    logger.debug(
      f"[copy] done priv={tp} pub={tu} total={total} elapsed={elapsed:.1f}s rate={rate:.0f}/s rss={rss_mb():.0f}MB"
    )
    return (tp, tu)
