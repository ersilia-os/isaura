import csv, datetime, json, os, shutil, sys, time, uuid, pyarrow.parquet as pq, pandas as pd
from collections import defaultdict, Counter
from contextlib import AbstractContextManager
from concurrent.futures import ThreadPoolExecutor, as_completed
from isaura.base import _BaseTransfer, BloomIndex, TrancheState, MinioStore, DuckDBMinio, _S3RangeFile
from isaura.helpers import (
  ACCESS_FILE,
  BLOOM_FILENAME,
  CHECKPOINT_EVERY,
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
  StreamingCsvSink,
  chunk_row_limit,
  logger,
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
  query_batched,
  split_csv,
  stream_parquet_filtered,
  stream_parquet_filtered_ordered,
  track_write_progress,
  write_access_file,
)


class IsauraChecker(AbstractContextManager):
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
    self.checkpoint_every = checkpoint_every
    self.bi = BloomIndex(
      self.store, self.bucket, self.base_prefix, self.local_dir, bloom_filename=bloom_filename
    )
    self._added = 0

  def seen(self, v):
    return self.bi.seen(v)

  def seen_many(self, vs):
    return self.bi.seen_many(vs)

  def register(self, v, rc=None):
    self.bi.register(v, rc=rc)
    self._added += 1
    if self._added >= self.checkpoint_every:
      self.bi.persist()
      self._added = 0

  def close(self):
    if self._added:
      self.bi.persist()
      self._added = 0

  def __enter__(self):
    return self

  def __exit__(self, et, ev, tb):
    self.close()


class IsauraWriter:
  def __init__(
    self,
    input_csv,
    model_id,
    model_version,
    access="public",
    bucket=None,
    endpoint=None,
    access_key=None,
    secrete=None,
  ):
    self.input_csv = input_csv
    self.model_id = model_id
    self.model_version = model_version
    self.access = access
    self.bucket = bucket or PUB
    self.base_prefix = get_base(self.model_id, model_version)
    self.access_key = get_acc_key(self.base_prefix)
    self.model_meta = fetch_schema_from_github(self.model_id)
    self.output_dimension = output_dimension_from_metadata(self.model_meta)
    self.max_rows = int(chunk_row_limit(self.output_dimension))
    self.store = MinioStore(endpoint=endpoint, access=access_key, secret=secrete)
    self.store.ensure_bucket(self.bucket)
    self.tmpdir = make_temp("isaura_writter_")
    self.metadata_path = os.path.join(self.tmpdir, ACCESS_FILE)
    self.bi = BloomIndex(
      self.store,
      self.bucket,
      self.base_prefix,
      self.tmpdir,
      bloom_filename=os.getenv("BLOOM_FILENAME", BLOOM_FILENAME),
    )
    self.chunk_state = TrancheState(self.store, self.bucket, self.base_prefix, self.tmpdir, self.max_rows)
    self.buffers = []
    self.buf_rows = 0
    self.schema_cols = None
    logger.info(
      f"writer init: bucket={self.bucket} base={self.base_prefix} csv={self.input_csv} output_dim={self.output_dimension} chunk_limit={self.max_rows}"
    )

  def _load_metadata(self):
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
    if not self.metadata_path or not self.bucket:
      return
    try:
      existed = self._load_metadata()
      write_access_file(existed, inputs, self.access, self.metadata_path)
      self.store.upload_file(self.metadata_path, self.bucket, f"{self.base_prefix}/{ACCESS_FILE}")
      logger.info(f"{ACCESS_FILE} -> minio://{self.bucket}/{self.base_prefix}/{ACCESS_FILE}")
    except Exception as e:
      logger.warning(f"metadata upload failed: {e}")

  def _set_schema(self, row):
    if self.schema_cols is None:
      self.schema_cols = list(row.keys())
      logger.info(f"writer schema: {self.schema_cols[:10]}")

  def _flush_if_needed(self):
    if self.buf_rows < self.max_rows or not self.buffers:
      return
    self.chunk_state.flush(self.buffers, self.schema_cols)
    self.buffers = []
    self.buf_rows = 0

  def write(self, df=None, show_progress: bool = True):
    total = dupes = 0
    new = []
    rows_iter, total_rows = (None, None)
    if df is not None:
      rows_list = df.to_dict("records")
      rows_iter = rows_list
      total_rows = len(rows_list)
    else:
      f = open(self.input_csv, newline="", encoding="utf-8")
      reader = csv.DictReader(f)
      rows_iter = reader
      total_rows = None
    try:
      if show_progress:
        rows_iter = track_write_progress(
          rows_iter, total=total_rows, description="Writing rows", console=logger.console
        )
      for row in rows_iter:
        self._set_schema(row)
        smi = (row.get("input") or row.get("smiles") or "").strip()
        if not smi:
          continue
        if not self.bi.seen(smi):
          new.append(smi)
          self.buffers.append(dict(row))
          self.buf_rows += 1
          self.bi.register(smi, rc=(1, 1))
          total += 1
          self._flush_if_needed()
          if self.bi._added >= int(CHECKPOINT_EVERY):
            self.bi.persist()
        else:
          dupes += 1
          continue
    finally:
      if df is None and "f" in locals():
        try:
          f.close()
        except Exception:
          pass
    if self.buffers:
      self.chunk_state.flush(self.buffers, self.schema_cols)
      self.buffers = []
      self.buf_rows = 0
    self.bi.persist()
    if new:
      self._upload_metadata(new)
    logger.info(f"write done: new={total} dupes={dupes}")

  def close(self):
    try:
      shutil.rmtree(self.tmpdir)
    except Exception:
      pass

  def __enter__(self):
    return self

  def __exit__(self, et, ev, tb):
    self.close()


class IsauraReader:
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
    self.store = MinioStore(endpoint=self.endpoint, access=access_key, secret=secrete)
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
    logger.info(
      f"reader init bucket={self.bucket} base={self.base} csv={self.input_csv} output_dim={self.output_dimension} wide={self.wide_model}"
    )

  def _wanted(self, df=None):
    wanted, header_set = ([], set())
    rows = df.to_dict("records") if df is not None else None
    if rows is None:
      logger.debug(f"[read:wanted] parsing inputs from csv={self.input_csv}")
      with open(self.input_csv, newline="", encoding="utf-8") as f:
        rows = csv.DictReader(f)
        for row in rows:
          h = INPUT_C[0] if row.get(INPUT_C[0]) else INPUT_C[1]
          v = (row.get(h) or "").strip()
          if v:
            wanted.append(v)
          if h:
            header_set.add(h)
    else:
      logger.debug(f"[read:wanted] parsing inputs from dataframe rows={len(rows)}")
      for row in rows:
        h = INPUT_C[0] if row.get(INPUT_C[0]) else INPUT_C[1]
        v = (row.get(h) or "").strip()
        if v:
          wanted.append(v)
        if h:
          header_set.add(h)
    header = list(header_set)[0] if header_set else "smiles"
    logger.info(f"[read:wanted] collected inputs={len(wanted)} header={header}")
    return (wanted, header)

  def _prepare_read(self, df=None):
    index = None
    t0 = time.time()
    wanted, header = self._wanted(df=df)
    if self.approximate:
      logger.info(f"[read:prepare] approximate mode — loading index for collection={self.collection}")
      index = self._load_index()
      if index and len(index) < MIN_NNS_RESULT_SIZE:
        logger.error(
          f"Minimum precalculation size for enabling nearest neighbor search is {MIN_NNS_RESULT_SIZE}, found {len(index)}. Aborting the Ops!"
        )
        sys.exit(1)
      st = time.perf_counter()
      wanted = get_apprx(wanted, self.collection)
      et = time.perf_counter()
      logger.info(f"Approximate inputs are retrieved {len(wanted)} in {et - st:.2f} seconds!")
      group_inputs(wanted, index, bloom=self.bi)
    else:
      logger.info(f"[read:prepare] exact mode — checking bloom filter for {len(wanted)} inputs")
      bloom_t0 = time.perf_counter()
      missing = [v for v in wanted if not self.bi.seen(v)]
      bloom_dt = time.perf_counter() - bloom_t0
      logger.info(f"[read:prepare] bloom check done in {bloom_dt:.3f}s missing={len(missing)}")
      if missing:
        logger.error(
          f"inputs not indexed: {missing[:5]}{('...' if len(missing) > 5 else '')} total_missing={len(missing)}"
        )
        sys.exit(1)
    prep_dt = time.time() - t0
    logger.info(f"[read:prepare] ready inputs={len(wanted)} header={header} prep_time={prep_dt:.2f}s")
    return (t0, wanted, header)

  def _load_index(self):
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

  def read(self, output_csv=None, df=None):
    t0, wanted, header = self._prepare_read(df=df)
    if not wanted:
      logger.info("[read] no inputs to query — returning empty frame")
      return pd.DataFrame()
    try:
      if self.wide_model:
        prefix = hive_prefix(self.base) + "/"
        logger.info(
          f"[read] wide-model ordered streaming mode n={len(wanted)} bucket={self.bucket} prefix={prefix} output_dim={self.output_dimension}"
        )
        source = stream_parquet_filtered_ordered(self.store, self.bucket, prefix, wanted, header=header)
      elif len(wanted) >= STREAM_PARQUET_THRESHOLD:
        prefix = hive_prefix(self.base) + "/"
        logger.info(f"[read] streaming mode n={len(wanted)} bucket={self.bucket} prefix={prefix}")
        source = stream_parquet_filtered(self.store, self.bucket, prefix, wanted, header=header)
      else:
        files = list_parquet_keys(self.store, self.bucket, self.base)
        logger.info(
          f"[read] duckdb mode n={len(wanted)} bucket={self.bucket} files={(len(files) if isinstance(files, list) else files)}"
        )
        source = query_batched(self.duck.con, header, wanted, files, tmpdir=self.tmpdir, preserve_order=True)
      total_rows = 0
      n_chunks = 0
      if output_csv:
        logger.info(f"[read] streaming csv output to {output_csv}")
        collect = []
        with StreamingCsvSink(output_csv) as sink:
          for chunk in source:
            n_chunks += 1
            sink.write_df(chunk)
            collect.append(chunk)
            total_rows += len(chunk)
            logger.debug(
              f"[read] chunk #{n_chunks} rows={len(chunk)} total_rows={total_rows} cols={len(chunk.columns)} rss={rss_mb():.0f}MB"
            )
        out = pd.concat(collect, ignore_index=True) if collect else pd.DataFrame()
        del collect
        logger.info(f"[read] wrote csv chunks={n_chunks} rows={total_rows} path={output_csv}")
      else:
        parts = []
        for chunk in source:
          n_chunks += 1
          parts.append(chunk)
          total_rows += len(chunk)
          logger.debug(
            f"[read] chunk #{n_chunks} rows={len(chunk)} total_rows={total_rows} cols={len(chunk.columns)} rss={rss_mb():.0f}MB"
          )
        out = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        del parts
      elapsed = time.time() - t0
      logger.success(f"Query successfully fetched for a given inputs in {elapsed:.2f} sec")
    except Exception as e:
      logger.error(f"[read] query failed: {e}")
      sys.exit(1)
    rate = total_rows / elapsed if elapsed > 0 and total_rows else 0.0
    te = time.time() - t0
    logger.info(
      f"[read] done model+version={self.pref} bucket={self.bucket} inputs={len(wanted)} matched={total_rows} chunks={n_chunks} (query elapsed)={elapsed:.2f}s rate={rate:.1f}/s (total time elapsed)={te:.2f}s rss={rss_mb():.0f}MB"
    )
    return out

  def read_batched(self, batch_size=10000, output_csv=None, df=None):
    t0, wanted, header = self._prepare_read(df=df)
    if not wanted:
      logger.info("[read_batched] no inputs to query — returning")
      return
    total_rows = 0
    n_chunks = 0
    sink = None
    try:
      if self.wide_model:
        prefix = hive_prefix(self.base) + "/"
        logger.info(
          f"[read_batched] wide-model ordered streaming mode n={len(wanted)} batch_size={batch_size} bucket={self.bucket} prefix={prefix} output_dim={self.output_dimension}"
        )
        source = stream_parquet_filtered_ordered(
          self.store, self.bucket, prefix, wanted, header=header, batch_size=batch_size
        )
      elif len(wanted) >= STREAM_PARQUET_THRESHOLD:
        prefix = hive_prefix(self.base) + "/"
        logger.info(
          f"[read_batched] streaming mode n={len(wanted)} batch_size={batch_size} bucket={self.bucket} prefix={prefix}"
        )
        source = stream_parquet_filtered(
          self.store, self.bucket, prefix, wanted, header=header, batch_size=batch_size
        )
      else:
        files = list_parquet_keys(self.store, self.bucket, self.base)
        logger.info(
          f"[read_batched] duckdb mode n={len(wanted)} batch_size={batch_size} bucket={self.bucket} files={(len(files) if isinstance(files, list) else files)}"
        )
        source = query_batched(
          self.duck.con, header, wanted, files, batch_size=batch_size, tmpdir=self.tmpdir, preserve_order=True
        )
      if output_csv:
        logger.info(f"[read_batched] streaming csv output to {output_csv}")
        sink = StreamingCsvSink(output_csv)
        sink.__enter__()
      for chunk in source:
        n_chunks += 1
        if sink is not None:
          sink.write_df(chunk)
        total_rows += len(chunk)
        logger.debug(
          f"[read_batched] chunk #{n_chunks} rows={len(chunk)} total_rows={total_rows} cols={len(chunk.columns)} rss={rss_mb():.0f}MB"
        )
        yield chunk
    except Exception as e:
      logger.error(f"[read_batched] query failed: {e}")
      sys.exit(1)
    finally:
      if sink is not None:
        sink.__exit__(None, None, None)
        logger.info(f"[read_batched] csv sink closed rows={total_rows} path={output_csv}")
    elapsed = time.time() - t0
    rate = total_rows / elapsed if elapsed > 0 and total_rows else 0.0
    logger.success(
      f"[read_batched] done model+version={self.pref} bucket={self.bucket} inputs={len(wanted)} matched={total_rows} chunks={n_chunks} elapsed={elapsed:.2f}s rate={rate:.1f}/s rss={rss_mb():.0f}MB"
    )


class IsauraCopy(_BaseTransfer):
  def copy(self):
    logger.info(f"[copy] starting model={self.model_id} v={self.model_version} bucket={self.bucket}")
    meta_local, meta = self._load_metadata()
    return self._copy(meta_local, meta)


class IsauraMover(_BaseTransfer):
  def move(self):
    t0 = time.time()
    logger.info(f"[move] starting model={self.model_id} v={self.model_version} bucket={self.bucket}")
    meta_local, meta = self._load_metadata()
    self._copy(meta_local, meta)
    n = self._delete()
    logger.success(f"[move] done wiped={n} elapsed={time.time() - t0:.1f}s")


class IsauraRemover(_BaseTransfer):
  def __init__(self, model_id=None, model_version=None, bucket=None, project_name=None):
    self._project_name = project_name
    if project_name:
      self._store = MinioStore()
    else:
      super().__init__(model_id, model_version, bucket)

  def remove(self):
    if self._project_name:
      store = self._store
      n = store.delete_prefix(self._project_name, "")
      logger.info(f"removed all objects in project={self._project_name} count={n}")
    else:
      n = self._delete()
      logger.info(f"removed objects={n}")


class IsauraPull(_BaseTransfer):
  def __init__(self, model_id, model_version, bucket, input_csv, output_dir=None):
    super().__init__(model_id, model_version, bucket, output_dir)
    self.model_id = model_id
    self.model_version = model_version
    self.bucket = bucket
    self.input_csv = input_csv

  def pull(self):
    logger.info(
      f"Pulling precalculation from the cloud for model={self.model_id}, version={self.model_version}, bucket={self.bucket}"
    )
    r = IsauraReader(
      model_id=self.model_id,
      model_version=self.model_version,
      input_csv=self.input_csv,
      approximate=False,
      bucket=self.bucket,
      endpoint=MINIO_ENDPOINT_CLOUD,
      access_key=mcak,
      secrete=mcsk,
    )
    out = self._pull_batched(r.read_batched())
    logger.info(f"pulled objects={out}")


class IsauraPush:
  def __init__(self, model_id, model_version, bucket):
    self.model_id = model_id
    self.model_version = model_version
    self.bucket = bucket

  def push(self):
    insp = IsauraInspect(model_id=self.model_id, model_version=self.model_version, access="both", cloud=False)
    df = insp.list_available()
    if df.empty:
      logger.error("No data found in any default bucket for a given model! Aborting push.")
      sys.exit(1)
    files = split_csv(df)
    if not files:
      logger.error("No data found in any default bucket! Aborting push.")
      sys.exit(1)
    file1 = files[0]
    file2 = files[1] if len(files) > 1 else None
    if not file2:
      logger.warning("Private bucket has no data! Skipping push for it.")
    for access, file, mck, mcs in [("public", file1, mcak, mcsk), ("private", file2, mcpak, mcpsk)]:
      if not file:
        continue
      r = IsauraReader(
        model_id=self.model_id,
        model_version=self.model_version,
        input_csv=file,
        approximate=False,
        bucket=f"isaura-{access}",
      )
      with IsauraWriter(
        input_csv=file,
        model_id=self.model_id,
        model_version=self.model_version,
        bucket=f"isaura-{access}",
        access=None if access == "public" else "private",
        endpoint=MINIO_ENDPOINT_CLOUD,
        access_key=mck,
        secrete=mcs,
      ) as w:
        for chunk in r.read_batched():
          w.write(df=chunk, show_progress=False)


class IsauraInspect:
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
    if self.project_name:
      return [self.project_name]
    if self.access == "public":
      return [PUB]
    if self.access == "private":
      return [PRI]
    return [PUB, PRI]

  def _creds(self, bucket):
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
    p = self._paginator(bucket)
    for page in p.paginate(Bucket=bucket, Prefix=prefix, Delimiter="/"):
      for cp in page.get("CommonPrefixes", []):
        yield cp["Prefix"]

  def list_objects(self, bucket, prefix=""):
    p = self._paginator(bucket)
    for page in p.paginate(Bucket=bucket, Prefix=prefix):
      for obj in page.get("Contents", []):
        yield obj["Key"]

  def iter_object_meta(self, bucket, prefix=""):
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
    base = get_base(model_id, model_version)
    key = get_idx_key(base)
    if self.heavy_index or not self.cloud:
      return self.get_json(bucket, key, max_bytes=10**18) or {}
    return {}

  def load_metadata(self, bucket, model_id, model_version):
    base = get_base(model_id, model_version)
    return self.get_json(bucket, f"{base}/{ACCESS_FILE}")

  def _indices_union(self):
    union, owner = ({}, {})
    for b in self.buckets():
      for mid, mv in self.iter_models(b):
        idx = self.load_index(b, mid, mv)
        for smi in idx.keys():
          if smi not in union:
            union[smi] = True
            owner[smi] = b
    return (union, owner)

  def list_available(self, output_csv=None):
    _, owner = self._indices_union()
    df = pd.DataFrame([{"input": smi, "bucket": b} for smi, b in owner.items()])
    if output_csv:
      df.to_csv(output_csv, index=False)
    return df

  def inspect_inputs(self, input_csv, output_csv=None):
    _, owner = self._indices_union()
    with open(input_csv, newline="", encoding="utf-8") as f:
      wanted = [(r.get("input") or "").strip() for r in csv.DictReader(f) if (r.get("input") or "").strip()]
    df = pd.DataFrame([{"input": s, "available": s in owner, "bucket": owner.get(s, "")} for s in wanted])
    if output_csv:
      df.to_csv(output_csv, index=False)
    return df

  def find_chunk_keys(self, bucket, model_id, model_version):
    pref = get_pref(model_id, model_version)
    for obj in self.iter_object_meta(bucket, pref):
      k = obj.get("Key") or ""
      if "/chunk_" in k and k.endswith(".parquet"):
        yield (k, int(obj.get("Size") or 0))

  def parquet_num_rows(self, bucket, parquet_key, size=None):
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
    total = 0
    n_chunks = 0
    for k, sz in self.find_chunk_keys(bucket, model_id, model_version):
      n_chunks += 1
      n = self.parquet_num_rows(bucket, k, size=sz)
      if n is not None:
        total += int(n)
    return (total, n_chunks)

  def inspect_models(self, bucket, prefix_filter=""):
    out = []
    for mid, mv in self.iter_models(bucket, prefix_filter=prefix_filter):
      if self.cloud and (not self.heavy_index):
        entries, chunks = self.entries_from_chunks(bucket, mid, mv)
        out.append({
          "bucket": bucket,
          "model_id": mid,
          "model_version": mv,
          "model": f"{mid}/{mv}",
          "entries": int(entries),
          "chunks": int(chunks),
        })
      else:
        out.append({
          "bucket": bucket,
          "model_id": mid,
          "model_version": mv,
          "model": f"{mid}/{mv}",
          "entries": len(self.load_index(bucket, mid, mv)),
          "chunks": None,
        })
    return out

  def find_any_chunk_key(self, bucket, model_id, model_version):
    pref = get_pref(model_id, model_version)
    for k in self.list_objects(bucket, pref):
      if "/chunk_" in k and k.endswith(".parquet"):
        return k
    return None

  def parquet_columns(self, bucket, parquet_key):
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
    if not sorted_vals:
      return None
    n = len(sorted_vals)
    i = int(round(p / 100.0 * (n - 1)))
    i = max(0, min(i, n - 1))
    return sorted_vals[i]

  def _use_chunks(self):
    return self.insp.cloud and (not self.insp.heavy_index)

  def _collect_objects(self, bucket, model_id, model_version):
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
    total = 0
    for k, sz in chunks:
      n = self.insp.parquet_num_rows(bucket, k, size=sz)
      if n is not None:
        total += int(n)
    return total

  def _process_model(self, bucket, mid, mv, use_chunks):
    total_bytes, total_gb, chunks, first_chunk = self._collect_objects(bucket, mid, mv)
    if use_chunks:
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

  def compute(self):
    generated_at = datetime.datetime.now(datetime.timezone.utc).isoformat()
    buckets = self.insp.buckets()
    use_chunks = self._use_chunks()
    tasks = []
    for b in buckets:
      for mid, mv in self.insp.iter_models(b):
        tasks.append((b, mid, mv))
    models = []
    models_per_molecule = Counter()
    models_total_by_bucket = defaultdict(int)
    if use_chunks and len(tasks) > 1:
      results = [None] * len(tasks)
      with ThreadPoolExecutor(max_workers=min(self.max_workers, len(tasks))) as pool:
        futures = {}
        for i, (b, mid, mv) in enumerate(tasks):
          fut = pool.submit(self._process_model, b, mid, mv, True)
          futures[fut] = i
        for fut in as_completed(futures):
          results[futures[fut]] = fut.result()
      for row, _ in results:
        models.append(row)
        models_total_by_bucket[row["bucket"]] += 1
    else:
      for b, mid, mv in tasks:
        row, idx = self._process_model(b, mid, mv, use_chunks)
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
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    tmp = output_path + ".tmp"
    data = self.compute()
    with open(tmp, "w", encoding="utf-8") as f:
      json.dump(data, f, indent=2, ensure_ascii=False)
    os.replace(tmp, output_path)
    return output_path
