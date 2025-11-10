import boto3, duckdb, json, os, pandas as pd, pickle, pyarrow as pa, pyarrow.parquet as pq, requests, sys, uuid

from isaura.helpers import (
  MINIO_ENDPOINT,
  MINIO_LOCAL_AK,
  MINIO_LOCAL_SK,
  MAX_ROWS_PER_FILE,
  CHECKPOINT_EVERY,
  BLOOM_FILENAME,
  ACCESS_FILE,
  INDEX_FILE,
  logger,
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
)

from botocore.config import Config
from boto3.s3.transfer import TransferConfig
from collections import defaultdict


from pybloom_live import ScalableBloomFilter


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
      logger.warning(f"The file {key} is not found in bucket {bucket}. Operation Abort!. Details -> {e}")

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
    self.con.execute(f"PRAGMA threads={threads};")
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
    error_rate=0.001,
    initial_capacity=1_000_000,
  ):
    self.store = store
    self.bucket = bucket
    self.base = base_prefix.strip("/")
    self.bloom_key = f"{self.base}/{BLOOM_FILENAME}"
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
    try:
      self.store.download_file(self.bucket, self.index_key, self.local_index)
      with open(self.local_index, "r", encoding="utf-8") as f:
        self.index = json.load(f)
      logger.info(f"loaded index {self.bucket}/{self.index_key} entries={len(self.index)}")
    except Exception:
      self.index = {}
      logger.info("created new index")
    self._added = 0

  def seen(self, v):
    return v in self.sbf

  def rc(self, v):
    return self.index.get(v)

  def register(self, v, rc=None):
    self.sbf.add(v)
    if rc is not None and v not in self.index:
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
    self.st = {"next": 1, "open": None, "rows": 0}
    self._w = None
    self._schema = None

  def _pref(self):
    return hive_prefix(self.base) + "/"

  def _ok(self, i):
    return f"{self._pref()}chunk_{i}.parquet"

  def _lp(self, i):
    return os.path.join(self.tmpdir, f"open_{i}.parquet")

  def _ls(self):
    ks = []
    for o in self.store.list_keys(self.bucket, self._pref()):
      k = o["Key"]
      if k.endswith(".parquet") and "/chunk_" in k:
        ks.append(k)
    ks.sort(key=self._idx)
    return ks

  def _idx(self, k):
    b = os.path.basename(k)
    n, _ = os.path.splitext(b)
    try: return int(n.split("_")[1])
    except: return 1

  def _esc(self, cols):
    if cols is None: return []
    return list(cols)

  def _ens(self, cols):
    if self._schema is None:
      self._schema = pa.schema([pa.field(c, pa.string()) for c in cols])

  def _wopen(self, i, cols):
    if self._w is None:
      self._ens(cols)
      self._w = pq.ParquetWriter(self._lp(i), self._schema, compression="zstd", use_dictionary=True)

  def _append(self, df, cols):
    cols = self._esc(cols)
    for c in cols:
      if c not in df.columns:
        df[c] = None
    df = df[cols]
    t = pa.Table.from_pandas(df, schema=self._schema, preserve_index=False)
    self._w.write_table(t)
    self.st["rows"] += t.num_rows

  def _closeup(self, i):
    if self._w:
      self._w.close()
      self._w = None
      local = self._lp(i)
      self.store.upload_file(local, self.bucket, self._ok(i))
      try: os.remove(local)
      except: pass

  def ensure(self):
    ks = self._ls()
    if not ks:
      self.st = {"next": 1, "open": self._lp(1), "rows": 0}
      return
    n = self._idx(ks[-1]) + 1
    self.st = {"next": n, "open": self._lp(n), "rows": 0}

  def flush(self, rows, schema_cols):
    if not rows: return
    self.ensure()
    df = pd.DataFrame(rows)
    cols = self._esc(schema_cols)
    i = self.st["next"]
    if not self.st["open"]:
      self.st["open"] = self._lp(i)
    r = len(df); s = 0
    while r > 0:
      space = self.max_rows - self.st["rows"]
      take = min(space, r)
      if take > 0:
        part = df.iloc[s:s+take]
        self._wopen(i, cols)
        self._append(part, cols)
        r -= take; s += take
      if self.st["rows"] >= self.max_rows:
        self._closeup(i)
        self.st["next"] += 1
        i = self.st["next"]
        self.st["open"] = self._lp(i)
        self.st["rows"] = 0
  

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

  def add_rows(self, r, c, df):
    if df.empty:
      return 0
    new_rows = []
    for _, row in df.iterrows():
      smi = row["input"]
      if smi and not self.bi.seen(smi):
        new_rows.append(dict(row))
        self.bi.register(smi, rc=(r, c))
    if not new_rows:
      return 0
    self.buffers[(r, c)].extend(new_rows)
    if len(self.buffers[(r, c)]) >= self.max_rows:
      self.tranche.flush(r, c, self.buffers[(r, c)], list(df.keys()))
      self.buffers[(r, c)].clear()
    return len(new_rows)

  def finalize(self, metadata_local=None, schema_cols=None):
    for k in list(self.buffers.keys()):
      if self.buffers[k]:
        r, c = k
        self.tranche.flush(r, c, self.buffers[k], schema_cols)
        self.buffers[k].clear()
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
    self.duck = DuckDBMinio(
      endpoint=self.store.endpoint,
      access=self.store.access,
      secret=self.store.secret,
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

  def select_rows(self, r, c, wanted, input_col="input"):
    return query(self.duck.con, input_col, wanted, get_files_glob(self.bucket, self.tranches))

  def _delete(self):
    return self.store.delete_prefix(self.bucket, self.tranches)

  def _pull(self, df, index):
    input = df["input"].astype(str).to_list()
    gp = group_inputs(input, index, force=True)
    wt = _SinkWriter(self.store, self.bucket, self.model_id, self.model_version, self.tmpdir)
    tp, tu, dfs = 0, 0, []
    if gp is None:
      return 0, 0
    for (r, c), want in gp.items():
      df_ = df[df["input"].isin(list(want))]
      _tp = wt.add_rows(r, c, df)
      if _tp != 0:
        dfs.append(df_)
      tp += _tp
    if wt:
      wt.finalize(schema_cols=list(df.keys()))
    if dfs:
      dfs = pd.concat(dfs, ignore_index=True)
      post_apprx(dfs, self.collection)
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
    if self.output_dir is not None:
      self._dump()
      return
    priv = [
      (d.get("input") or "").strip()
      for d in meta_list
      if (d.get("access") or "").lower() == "private" and (d.get("input") or "").strip()
    ]
    pub = [
      (d.get("input") or "").strip()
      for d in meta_list
      if (d.get("access") or "").lower() == "public" and (d.get("input") or "").strip()
    ]
    index, _ = self._load_index()
    gp = group_inputs(priv, index) if priv else {}
    gu = group_inputs(pub, index) if pub else {}
    w_priv = (
      _SinkWriter(self.store, "isaura-private", self.model_id, self.model_version, self.tmpdir)
      if gp
      else None
    )
    w_pub = (
      _SinkWriter(self.store, "isaura-public", self.model_id, self.model_version, self.tmpdir) if gu else None
    )
    tp, tu, dfs = 0, 0, []
    for (r, c), want in gp.items():
      df = self.select_rows(r, c, want)
      _tp = w_priv.add_rows(r, c, df)
      if _tp != 0:
        dfs.append(df)
      tp += _tp
    for (r, c), want in gu.items():
      df = self.select_rows(r, c, want)
      _tu = w_pub.add_rows(r, c, df)
      if _tu != 0:
        dfs.append(df)
      tu += _tu
    if w_priv:
      w_priv.finalize(metadata_local=meta_local, schema_cols=list(df.keys()))
    if w_pub:
      w_pub.finalize(metadata_local=meta_local, schema_cols=list(df.keys()))
    if dfs:
      dfs = pd.concat(dfs, ignore_index=True)
      post_apprx(dfs, self.collection)
    return tp, tu
