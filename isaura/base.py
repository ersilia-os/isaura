import boto3, duckdb, gc, json, os, pandas as pd, pickle, pyarrow as pa, pyarrow.parquet as pq, requests, sys, time, uuid
from botocore.config import Config
from boto3.s3.transfer import TransferConfig
from collections import defaultdict

from pybloom_live import ScalableBloomFilter

from isaura.helpers import (
  MINIO_ENDPOINT,
  MINIO_LOCAL_AK,
  MINIO_LOCAL_SK,
  MAX_ROWS_PER_FILE,
  IMMUTABLE_CHUNK_COLS_THRESHOLD,
  CHECKPOINT_EVERY,
  BLOOM_FILENAME,
  ACCESS_FILE,
  INDEX_FILE,
  STREAM_PARQUET_THRESHOLD,
  DEFAULT_BUCKET_NAME as pub_bucket,
  DEFAULT_PRIVATE_BUCKET_NAME as priv_bucket,
  logger,
  rss_mb,
  get_acc_key,
  get_base,
  get_coll,
  get_idx_key,
  get_files_glob,
  get_pref,
  group_inputs,
  hive_prefix,
  make_temp,
  post_apprx,
  query,
  query_batched,
  stream_parquet_filtered,
  cpu_cnt,
)


class _S3RangeFile:
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
    return True

  def seekable(self):
    return True

  def tell(self):
    return self._pos

  def seek(self, offset, whence=0):
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
    self.closed = True
    self._cache.clear()


class MinioStore:
  def __init__(
    self,
    endpoint=None,
    access=None,
    secret=None,
    multipart_threshold=5 * 1024 * 1024,
    multipart_chunksize=5 * 1024 * 1024,
    max_concurrency=2,
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
      config=Config(signature_version="s3v4", s3={"addressing_style": "path"}),
    )
    self.transfer_config = TransferConfig(
      multipart_threshold=multipart_threshold,
      multipart_chunksize=multipart_chunksize,
      max_concurrency=max_concurrency,
      use_threads=use_threads,
    )
    if not self.ping(self.client):
      sys.exit(1)

  def ping(self, store):
    try:
      url = f"{self.endpoint.rstrip('/')}/minio/health/ready"
      resp = requests.get(url, timeout=3)

      if resp.status_code == 200:
        logger.info(f"MinIO server healthy at {url}")
        return True
      logger.error(
        f"MinIO health check failed. Maybe the minio containers not running[{resp.status_code}]: {resp.text.strip()}"
      )
    except requests.exceptions.RequestException as e:
      logger.error(f"MinIO health request error. Maybe the minio containers not running: {e}")
    except Exception as e:
      logger.error(f"Unexpected error during MinIO health check: {e}")
    return False

  def ensure_bucket(self, bucket):
    try:
      self.client.head_bucket(Bucket=bucket)
    except Exception as e:
      logger.error(e)
      try:
        self.client.create_bucket(Bucket=bucket)
      except Exception as e:
        logger.error(e)
        sys.exit(1)

  def download_file(self, bucket, key, local):
    try:
      self.client.download_file(bucket, key, local, Config=self.transfer_config)
    except Exception as e:
      logger.warning(f"The file {key} is not found in bucket {bucket}. Details -> {e}")
      raise

  def upload_file(self, local, bucket, key, extra_args=None):
    self.client.upload_file(local, bucket, key, ExtraArgs=extra_args or {}, Config=self.transfer_config)

  def list_keys(self, bucket, prefix):
    p = self.client.get_paginator("list_objects_v2")
    for page in p.paginate(Bucket=bucket, Prefix=prefix):
      for obj in page.get("Contents", []):
        yield obj

  def delete_prefix(self, bucket, prefix):
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


class DuckDBMinio:  # Singleton design here
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

    self.endpoint, self.access, self.secret = endpoint, access, secret
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
    if getattr(self, "con", None):
      self.con.close()
      self.con = None


class BloomIndex:
  def __init__(
    self,
    store,
    bucket,
    base_prefix,
    local_dir,
    bloom_filename=BLOOM_FILENAME,
    error_rate=1e-14,
    initial_capacity=1_000_000,
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
      logger.info(f"loaded bloom {self.bucket}/{self.bloom_key}")
    except Exception:
      self.sbf = ScalableBloomFilter(
        mode=ScalableBloomFilter.SMALL_SET_GROWTH,
        initial_capacity=initial_capacity,
        error_rate=error_rate,
      )
      logger.info("created new bloom")
    if load_index:
      try:
        self.store.download_file(self.bucket, self.index_key, self.local_index)
        with open(self.local_index, "r", encoding="utf-8") as f:
          self.index = json.load(f)
        logger.info(f"loaded index {self.bucket}/{self.index_key} entries={len(self.index)}")
      except Exception:
        self.index = {}
        logger.info("created new index")
    else:
      self.index = None
      logger.info("skipped index load")
    self._added = 0

  def seen(self, v):
    return v in self.sbf

  def seen_many(self, vs):
    sbf = self.sbf
    b = self.bucket
    return {v: [v in sbf, b] for v in vs}

  def rc(self, v):
    return self.index.get(v) if self.index is not None else None

  def register(self, v, rc=None):
    self.sbf.add(v)
    if rc is not None and self.index is not None and v not in self.index:
      self.index[v] = list(rc)
    self._added += 1
    if self._added >= CHECKPOINT_EVERY:
      self.persist()

  def persist(self):
    tb = f"{self.local_bloom}.tmp"
    with open(tb, "wb") as f:
      pickle.dump(self.sbf, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tb, self.local_bloom)
    try:
      self.store.upload_file(self.local_bloom, self.bucket, self.bloom_key)
    except Exception as e:
      logger.warning(f"bloom upload failed: {e}")
    if self.index is not None:
      ti = f"{self.local_index}.tmp"
      with open(ti, "w", encoding="utf-8") as f:
        json.dump(self.index, f, separators=(",", ":"), ensure_ascii=False)
      os.replace(ti, self.local_index)
      try:
        self.store.upload_file(self.local_index, self.bucket, self.index_key)
      except Exception as e:
        logger.warning(f"index upload failed: {e}")
    self._added = 0


class TrancheState:
  def __init__(self, store, bucket, base_prefix, tmpdir, max_rows):
    self.store = store
    self.bucket = bucket
    self.base = base_prefix.strip("/")
    self.tmpdir = tmpdir
    self.max_rows = max_rows
    self.write_batch_rows = 2048
    self.state = {}

  def _use_immutable_chunks(self, schema_cols):
    return bool(schema_cols) and len(schema_cols) >= IMMUTABLE_CHUNK_COLS_THRESHOLD

  def _rows_to_frame(self, rows, schema_cols):
    df = pd.DataFrame.from_records(rows)
    for col in schema_cols:
      if col not in df.columns:
        df[col] = None
    return df[schema_cols]

  def _ensure_cols(self, df, schema_cols):
    for col in schema_cols:
      if col not in df.columns:
        df[col] = None
    return df[schema_cols]

  def _iter_tables(self, rows, schema_cols):
    for start in range(0, len(rows), self.write_batch_rows):
      frame = self._rows_to_frame(rows[start : start + self.write_batch_rows], schema_cols)
      yield pa.Table.from_pandas(frame, preserve_index=False)

  def _iter_tables_df(self, df, schema_cols):
    df = self._ensure_cols(df, schema_cols)
    for start in range(0, len(df), self.write_batch_rows):
      frame = df.iloc[start : start + self.write_batch_rows]
      yield pa.Table.from_pandas(frame, preserve_index=False)

  def _rows_in_remote(self, key):
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
    keys, t = self._list_chunks(), "data"
    if not keys:
      self.state[t] = {"next": 1, "open": None, "rows": 0}
      logger.info(f"tranche new: next=1")
      return

    last = keys[-1]
    idx = self._chunk_idx(last)
    self.state[t] = {"next": idx + 1, "open": None, "rows": 0}
    logger.info(f"tranche next: {idx + 1}")

  def _list_chunks(self):
    pref = hive_prefix(self.base) + "/"
    keys = []
    for obj in self.store.list_keys(self.bucket, pref):
      k = obj["Key"]
      if k.endswith(".parquet") and "/chunk_" in k:
        keys.append(k)
    return sorted(keys, key=self._chunk_idx)

  def _chunk_idx(self, key):
    base = os.path.basename(key)
    name, _ = os.path.splitext(base)
    try:
      return int(name.split("_")[1])
    except:
      return 1

  def _write_chunk(self, rows, schema_cols, idx, mode="new", existing_local=None):
    os_key = f"{hive_prefix(self.base)}/chunk_{idx}.parquet"
    local = existing_local or os.path.join(self.tmpdir, f"chunk_{uuid.uuid4().hex}.parquet")

    if mode == "append" and existing_local:
      tmp = local + ".merged"
      old_pf = pq.ParquetFile(existing_local)
      schema = old_pf.schema_arrow
      with pq.ParquetWriter(tmp, schema) as writer:
        for batch in old_pf.iter_batches(batch_size=self.write_batch_rows):
          writer.write_batch(batch)
        for table in self._iter_tables(rows, schema_cols):
          writer.write_table(table.cast(schema))
      os.replace(tmp, local)
    else:
      tmp = local + ".tmp"
      writer = None
      try:
        for table in self._iter_tables(rows, schema_cols):
          if writer is None:
            writer = pq.ParquetWriter(tmp, table.schema)
          writer.write_table(table)
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
    os_key = f"{hive_prefix(self.base)}/chunk_{idx}.parquet"
    local = existing_local or os.path.join(self.tmpdir, f"chunk_{uuid.uuid4().hex}.parquet")

    if mode == "append" and existing_local:
      tmp = local + ".merged"
      old_pf = pq.ParquetFile(existing_local)
      schema = old_pf.schema_arrow
      with pq.ParquetWriter(tmp, schema) as writer:
        for batch in old_pf.iter_batches(batch_size=self.write_batch_rows):
          writer.write_batch(batch)
        for table in self._iter_tables_df(df, schema_cols):
          writer.write_table(table.cast(schema))
      os.replace(tmp, local)
    else:
      tmp = local + ".tmp"
      writer = None
      try:
        for table in self._iter_tables_df(df, schema_cols):
          if writer is None:
            writer = pq.ParquetWriter(tmp, table.schema)
          writer.write_table(table)
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

  def flush(self, rows, schema_cols):
    if not rows:
      return self.ensure()
    self.ensure()
    st = self.state["data"]
    remaining = len(rows)
    start = 0
    immutable = self._use_immutable_chunks(schema_cols)

    logger.info(
      f"flush: tranche rows={remaining} cols={len(schema_cols or [])} mode={'immutable' if immutable else 'append'}"
    )

    if immutable:
      while remaining > 0:
        take = min(self.max_rows, remaining)
        part = rows[start : start + take]
        os_key = self._write_chunk(part, schema_cols, st["next"], mode="new")
        logger.info(f"flush: immutable tranche idx={st['next']} rows={take}")
        st["next"] += 1
        st["open"] = None
        st["rows"] = 0
        remaining -= take
        start += take
      return

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
          logger.info(f"flush: appended tranche idx={st['next']} +{take} -> {st['rows']}")

        if st["rows"] >= self.max_rows:
          st["next"] += 1
          st["open"] = None
          st["rows"] = 0
          logger.info(f"flush: closed tranche next={st['next']}")
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
        logger.info(f"flush: new open tranche idx={st['next']} rows={take}")
      else:
        st["next"] += 1
        st["open"] = None
        st["rows"] = 0
        logger.info(f"flush: full tranche idx={st['next'] - 1} rows={take}")

      remaining -= take
      start += take

  def flush_df(self, df, schema_cols):
    if df is None or df.empty:
      return self.ensure()
    self.ensure()
    st = self.state["data"]
    remaining = len(df)
    start = 0
    immutable = self._use_immutable_chunks(schema_cols)

    logger.info(
      f"flush_df: tranche rows={remaining} cols={len(schema_cols or [])} mode={'immutable' if immutable else 'append'}"
    )

    if immutable:
      while remaining > 0:
        take = min(self.max_rows, remaining)
        part = df.iloc[start : start + take]
        os_key = self._write_chunk_df(part, schema_cols, st["next"], mode="new")
        logger.info(f"flush_df: immutable tranche idx={st['next']} rows={take}")
        st["next"] += 1
        st["open"] = None
        st["rows"] = 0
        remaining -= take
        start += take
      return

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
          logger.info(f"flush_df: appended idx={st['next']} +{take} -> {st['rows']}")

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


class _SinkWriter:
  def __init__(self, store, bucket, model_id, model_version, tmpdir):
    self.store = store
    self.bucket = bucket
    self.model_id = model_id
    self.model_version = model_version
    self.base = get_base(model_id, model_version)
    self.tmpdir = tmpdir
    self.max_rows = MAX_ROWS_PER_FILE
    self.store.ensure_bucket(self.bucket)
    self.bi = BloomIndex(self.store, self.bucket, self.base, tmpdir)
    self.tranche = TrancheState(self.store, self.bucket, self.base, tmpdir, self.max_rows)
    self.buffers = defaultdict(list)
    self.buf_rows = defaultdict(int)
    self.schema_cols = None

  def add_rows(self, r, c, df):
    try:
      if df is None or df.empty:
        logger.debug(f"[sink] ({r},{c}) empty batch")
        return 0

      n_in = len(df)
      if self.schema_cols is None:
        self.schema_cols = list(df.columns)

      inputs = df["input"].astype(str).str.strip()
      non_blank = inputs != ""
      sbf = self.bi.sbf
      not_seen = non_blank & ~inputs.map(lambda v: v in sbf)

      skipped_blank = int((~non_blank).sum())
      skipped_seen = int(non_blank.sum() - not_seen.sum())
      new_df = df.loc[not_seen]
      added = len(new_df)

      if added == 0:
        logger.info(
          f"[sink] ({r},{c}) in={n_in} added=0 seen={skipped_seen} blank={skipped_blank}"
        )
        return 0

      for smi in inputs.loc[not_seen]:
        self.bi.register(smi, rc=(r, c))

      self.buffers[(r, c)].append(new_df)
      self.buf_rows[(r, c)] += added

      logger.info(
        f"[sink] ({r},{c}) in={n_in} added={added} seen={skipped_seen} "
        f"blank={skipped_blank} buf={self.buf_rows[(r, c)]}/{self.max_rows}"
      )

      if self.buf_rows[(r, c)] >= self.max_rows:
        self._flush_key(r, c)

      return added

    except Exception as e:
      logger.error(f"[sink] ({r},{c}) error: {e}")
      return 0

  def _flush_key(self, r, c):
    parts = self.buffers[(r, c)]
    if not parts:
      return
    merged = pd.concat(parts, ignore_index=True)
    logger.info(f"[sink] ({r},{c}) flushing {len(merged)} rows")
    self.tranche.flush_df(merged, self.schema_cols)
    self.buffers[(r, c)] = []
    self.buf_rows[(r, c)] = 0
    del merged
    gc.collect()

  def finalize(self, metadata_local=None, schema_cols=None):
    schema_cols = schema_cols or self.schema_cols
    for k in list(self.buffers.keys()):
      if self.buffers[k]:
        merged = pd.concat(self.buffers[k], ignore_index=True)
        self.tranche.flush_df(merged, schema_cols)
        self.buffers[k] = []
        self.buf_rows[k] = 0
        del merged
    gc.collect()
    self.bi.persist()
    if metadata_local:
      try:
        self.store.upload_file(metadata_local, self.bucket, f"{self.base}/{ACCESS_FILE}")
      except Exception:
        pass


class _BaseTransfer:
  def __init__(self, model_id, model_version, bucket, output_dir=None):
    self.model_id = model_id
    self.model_version = model_version
    self.bucket = bucket
    self.output_dir = output_dir
    self.base = get_pref(model_id, model_version)
    self.tranches = get_base(model_id, model_version)
    self.collection = get_coll(self.model_id, self.model_version)
    self.store = MinioStore()
    self.tmpdir = make_temp("isaura_xfer_")
    self.tmpdir_sinkw = make_temp("isaura_sinkw_")
    self.duck = DuckDBMinio(
      endpoint=self.store.endpoint,
      access=self.store.access,
      secret=self.store.secret,
    )
    self.bi = BloomIndex(
      self.store,
      self.bucket,
      self.tranches,
      self.tmpdir,
      bloom_filename=os.getenv("BLOOM_FILENAME", BLOOM_FILENAME),
    )

  def _download_if_exists(self, key, local):
    try:
      self.store.download_file(self.bucket, key, local)
      return True
    except Exception as e:
      logger.error(f"Either model id or version or project name are maybe specified wrongly. Results -> {e}")
      sys.exit(1)
      return False

  def _load_metadata(self):
    local = os.path.join(self.tmpdir, ACCESS_FILE)
    acc = get_acc_key(self.tranches)
    try:
      if self._download_if_exists(acc, local):
        with open(local, "r", encoding="utf-8") as f:
          return local, json.load(f)
    except Exception as e:
      logger.error(f"{e}. Possible causes -> 1. No project 2. Incorrect model id or version")

  def _load_index(self):
    key = get_idx_key(self.tranches)
    local = os.path.join(self.tmpdir, "index_src.json")
    if not self._download_if_exists(key, local):
      raise RuntimeError("index.json not found")
    with open(local, "r", encoding="utf-8") as f:
      return json.load(f), local

  def select_rows(self, wanted, input_col="input"):
    return query(self.duck.con, input_col, wanted, get_files_glob(self.bucket, self.tranches))

  def select_rows_batched(self, wanted, input_col="input", batch_size=10_000):
    n = len(wanted) if hasattr(wanted, '__len__') else None
    if n is not None and n >= STREAM_PARQUET_THRESHOLD:
      prefix = hive_prefix(self.tranches) + "/"
      logger.info(f"[select] streaming mode n={n} bucket={self.bucket} prefix={prefix}")
      yield from stream_parquet_filtered(
        self.store, self.bucket, prefix, wanted, header=input_col, batch_size=batch_size
      )
    else:
      files = get_files_glob(self.bucket, self.tranches)
      for chunk in query_batched(self.duck.con, input_col, wanted, files, batch_size=batch_size, tmpdir=self.tmpdir):
        yield chunk

  def _delete(self):
    return self.store.delete_prefix(self.bucket, self.tranches)

  def _pull(self, df, index):
    t0 = time.time()
    input = df["input"].astype(str).to_list()
    gp = group_inputs(input, index, bloom=self.bi, force=True)
    wt = _SinkWriter(self.store, self.bucket, self.model_id, self.model_version, self.tmpdir)
    tp, tu = 0, 0
    apprx_inputs = []
    if gp is None:
      return 0, 0
    for (r, c), want in gp.items():
      df_ = df[df["input"].isin(list(want))]
      _tp = wt.add_rows(r, c, df_)
      if _tp != 0:
        apprx_inputs.extend(df_["input"].astype(str).tolist())
      tp += _tp
      del df_
    if wt:
      wt.finalize(schema_cols=list(df.keys()))
    if apprx_inputs:
      apprx_df = pd.DataFrame({"input": apprx_inputs})
      post_apprx(apprx_df, self.collection)
      del apprx_df
    logger.info(f"[pull] done rows={tp} elapsed={time.time()-t0:.1f}s rss={rss_mb():.0f}MB")
    return tp, tu

  def _dump(self):
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

    index, _ = self._load_index()
    rc_lookup = {}
    for inp in all_wanted:
      if inp in index:
        r, c = index[inp]
        rc_lookup[inp] = (int(r), int(c))
    del index
    gc.collect()

    has_priv = any(v == "priv" for v in route.values())
    has_pub = any(v == "pub" for v in route.values())

    w_priv = (
      _SinkWriter(self.store, priv_bucket, self.model_id, self.model_version, self.tmpdir_sinkw)
      if has_priv
      else None
    )
    w_pub = (
      _SinkWriter(self.store, pub_bucket, self.model_id, self.model_version, self.tmpdir_sinkw)
      if has_pub
      else None
    )

    tp, tu, n_batches = 0, 0, 0
    schema_cols = None

    for df in self.select_rows_batched(all_wanted):
      if df is None or df.empty:
        continue
      n_batches += 1
      if schema_cols is None:
        schema_cols = list(df.columns)

      inputs = df["input"].astype(str).str.strip()
      rcs = inputs.map(rc_lookup.get)
      routes = inputs.map(route.get)
      valid_rc = rcs.notna()

      rc_r = rcs.apply(lambda v: v[0] if v is not None else None)
      rc_c = rcs.apply(lambda v: v[1] if v is not None else None)

      priv_mask = (routes == "priv") & valid_rc
      pub_mask = (routes == "pub") & valid_rc

      if w_priv and priv_mask.any():
        for (r, c), grp in df.loc[priv_mask].groupby([rc_r.loc[priv_mask], rc_c.loc[priv_mask]]):
          tp += w_priv.add_rows(int(r), int(c), grp)

      if w_pub and pub_mask.any():
        for (r, c), grp in df.loc[pub_mask].groupby([rc_r.loc[pub_mask], rc_c.loc[pub_mask]]):
          tu += w_pub.add_rows(int(r), int(c), grp)

      del df
      if n_batches % 20 == 0:
        gc.collect()
        logger.info(
          f"[copy] batch={n_batches} priv={tp} pub={tu} "
          f"elapsed={time.time()-t0:.1f}s rss={rss_mb():.0f}MB"
        )

    if w_priv:
      w_priv.finalize(metadata_local=meta_local, schema_cols=schema_cols)
    if w_pub:
      w_pub.finalize(metadata_local=meta_local, schema_cols=schema_cols)

    total = tp + tu
    if total > 0:
      try:
        apprx_df = pd.DataFrame({"input": list(rc_lookup.keys())[:total]})
        post_apprx(apprx_df, self.collection)
        del apprx_df
      except Exception as e:
        logger.warning(f"[copy] post_apprx failed: {e}")

    del rc_lookup, route
    gc.collect()
    elapsed = time.time() - t0
    rate = total / elapsed if elapsed > 0 else 0
    logger.success(
      f"[copy] done priv={tp} pub={tu} total={total} "
      f"elapsed={elapsed:.1f}s rate={rate:.0f}/s rss={rss_mb():.0f}MB"
    )
    return tp, tu
