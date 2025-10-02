import boto3, duckdb, json, os, pickle, tempfile, uuid

import pandas as pd

from isaura.helpers import (
  MINIO_ENDPOINT,
  MINIO_ACCESS_KEY,
  MINIO_SECRET_KEY,
  STORE_DIRECTORY,
  MAX_ROWS_PER_FILE,
  CHECKPOINT_EVERY,
  BLOOM_FILENAME,
  logger,
)

from botocore.config import Config
from collections import defaultdict

from pybloom_live import ScalableBloomFilter


class _BloomIndex:
  def __init__(self, s3, bucket, base_prefix, local_dir):
    self.s3 = s3
    self.bucket = bucket
    self.base = base_prefix.strip("/")
    self.bloom_key = f"{self.base}/bloom.pkl"
    self.index_key = f"{self.base}/index.json"
    self.local_bloom = os.path.join(local_dir, BLOOM_FILENAME)
    self.local_index = os.path.join(local_dir, "index.json")
    self._added = 0
    os.makedirs(local_dir, exist_ok=True)
    try:
      self.s3.download_file(self.bucket, self.bloom_key, self.local_bloom)
      with open(self.local_bloom, "rb") as f:
        self.sbf = pickle.load(f)
      logger.info(f"downloaded bloom: s3://{self.bucket}/{self.bloom_key}")
    except Exception:
      self.sbf = ScalableBloomFilter(
        mode=ScalableBloomFilter.SMALL_SET_GROWTH,
        initial_capacity=1_000_000,
        error_rate=0.001,
      )
      logger.info(f"created new bloom for {self.bucket}/{self.base}")
    try:
      self.s3.download_file(self.bucket, self.index_key, self.local_index)
      with open(self.local_index, "r", encoding="utf-8") as f:
        self.index = json.load(f)
      logger.info(
        f"downloaded index: s3://{self.bucket}/{self.index_key} entries={len(self.index)}"
      )
    except Exception:
      self.index = {}
      logger.info(f"starting empty index for {self.bucket}/{self.base}")

  def seen(self, v):
    return v in self.sbf

  def register(self, v, rc=None):
    self.sbf.add(v)
    if rc is not None and v not in self.index:
      self.index[v] = list(rc)
    self._added += 1
    if self._added >= CHECKPOINT_EVERY:
      self.persist()

  def persist(self):
    tmp = f"{self.local_bloom}.tmp"
    with open(tmp, "wb") as f:
      pickle.dump(self.sbf, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, self.local_bloom)
    try:
      self.s3.upload_file(self.local_bloom, self.bucket, self.bloom_key)
    except Exception as e:
      logger.warning(f"failed to upload bloom: {e}")
    tmp = f"{self.local_index}.tmp"
    with open(tmp, "w", encoding="utf-8") as f:
      json.dump(self.index, f, separators=(",", ":"), ensure_ascii=False)
    os.replace(tmp, self.local_index)
    try:
      self.s3.upload_file(self.local_index, self.bucket, self.index_key)
    except Exception as e:
      logger.warning(f"failed to upload index: {e}")
    self._added = 0


class _SinkWriter:
  def __init__(self, s3, bucket, model_id, model_version, tmpdir):
    self.s3 = s3
    self.bucket = bucket
    self.model_id = model_id
    self.model_version = model_version
    self.base = f"{model_id}/{model_version}/tranches"
    self.tmpdir = tmpdir
    self.max_rows = MAX_ROWS_PER_FILE
    self.state = {}
    self.buffers = defaultdict(list)
    self.bi = _BloomIndex(self.s3, self.bucket, self.base, tmpdir)

  def _tranche_prefix(self, r, c):
    return f"{self.base}/tranche_{r}_{c}"

  def _list_chunks(self, r, c):
    prefix = self._tranche_prefix(r, c) + "/"
    paginator = self.s3.get_paginator("list_objects_v2")
    keys = []
    for page in paginator.paginate(Bucket=self.bucket, Prefix=prefix):
      for o in page.get("Contents", []):
        k = o["Key"]
        if k.endswith(".parquet") and "/chunk_" in k:
          keys.append(k)
    return sorted(keys)

  def _chunk_idx(self, key):
    b = os.path.basename(key)
    n = os.path.splitext(b)[0]
    try:
      return int(n.split("_")[1])
    except:
      return 1

  def _num_rows_remote(self, key):
    local = os.path.join(self.tmpdir, f"inspect_{uuid.uuid4().hex}.parquet")
    try:
      self.s3.download_file(self.bucket, key, local)
      df = pd.read_parquet(local)
      return len(df)
    except:
      return 0
    finally:
      try:
        os.remove(local)
      except:
        pass

  def _ensure_state(self, r, c):
    t = (r, c)
    if t in self.state:
      return
    keys = self._list_chunks(r, c)
    if not keys:
      self.state[t] = {"next": 1, "open": None, "rows": 0}
      return
    last = keys[-1]
    n = self._num_rows_remote(last)
    idx = self._chunk_idx(last)
    if n < self.max_rows:
      self.state[t] = {"next": idx, "open": last, "rows": n}
    else:
      self.state[t] = {"next": idx + 1, "open": None, "rows": 0}

  def _write_chunk(self, df, r, c, idx, mode="new", existing_local=None):
    prefix = self._tranche_prefix(r, c)
    os_key = f"{prefix}/chunk_{idx}.parquet"
    local = existing_local or os.path.join(
      self.tmpdir, f"chunk_{uuid.uuid4().hex}.parquet"
    )
    if mode == "append" and existing_local:
      old = pd.read_parquet(existing_local)
      df = pd.concat([old, df], ignore_index=True)
    df.to_parquet(local, index=False)
    self.s3.upload_file(local, self.bucket, os_key)
    if not existing_local:
      try:
        os.remove(local)
      except:
        pass
    return os_key, len(df)

  def _flush(self, r, c):
    t = (r, c)
    buf = self.buffers.get(t, [])
    if not buf:
      return
    self._ensure_state(r, c)
    st = self.state[t]
    df = pd.DataFrame(buf)
    rem = len(df)
    i = 0
    if st["open"]:
      tmp = os.path.join(self.tmpdir, f"open_{uuid.uuid4().hex}.parquet")
      try:
        self.s3.download_file(self.bucket, st["open"], tmp)
        space = self.max_rows - st["rows"]
        take = min(space, rem)
        if take > 0:
          part = df.iloc[i : i + take]
          self._write_chunk(part, r, c, st["next"], mode="append", existing_local=tmp)
          st["rows"] += take
          rem -= take
          i += take
        if st["rows"] >= self.max_rows:
          st["next"] += 1
          st["open"] = None
          st["rows"] = 0
      finally:
        try:
          os.remove(tmp)
        except:
          pass
    while rem > 0:
      take = min(self.max_rows, rem)
      part = df.iloc[i : i + take]
      os_key, _ = self._write_chunk(part, r, c, st["next"], mode="new")
      if take < self.max_rows:
        st["open"] = os_key
        st["rows"] = take
      else:
        st["next"] += 1
        st["open"] = None
        st["rows"] = 0
      rem -= take
      i += take
    self.buffers[t].clear()

  def add_rows(self, r, c, df):
    if df.empty:
      return 0
    new_rows = []
    for _, row in df.iterrows():
      smi = row["input"]
      if smi and not self.bi.seen(smi):
        new_rows.append(dict(row))
        self.bi.register(smi, (r, c))
    if not new_rows:
      return 0
    self.buffers[(r, c)].extend(new_rows)
    if len(self.buffers[(r, c)]) >= self.max_rows:
      self._flush(r, c)
    return len(new_rows)

  def finalize(self, metadata_local=None):
    for r, c in list(self.buffers.keys()):
      if self.buffers[(r, c)]:
        self._flush(r, c)
    self.bi.persist()
    if metadata_local:
      meta_tranches_key = f"{self.base}/metadata.json"
      try:
        self.s3.upload_file(metadata_local, self.bucket, meta_tranches_key)
      except Exception as e:
        logger.warning(f"failed to upload metadata to {self.bucket}: {e}")


class _BaseTransfer:
  def __init__(self, model_id, model_version, source_bucket):
    self.model_id = model_id
    self.model_version = model_version
    self.source_bucket = source_bucket
    self.base = f"{model_id}/{model_version}"
    self.tranches = f"{self.base}/tranches"
    self.s3 = boto3.client(
      "s3",
      endpoint_url=MINIO_ENDPOINT,
      aws_access_key_id=MINIO_ACCESS_KEY,
      aws_secret_access_key=MINIO_SECRET_KEY,
      config=Config(signature_version="s3v4", s3={"addressing_style": "path"}),
    )
    self.tmpdir = tempfile.mkdtemp(prefix="isaura_xfer_", dir=STORE_DIRECTORY)
    self.con = duckdb.connect(database=":memory:")
    self.con.execute("INSTALL httpfs; LOAD httpfs;")
    ep = MINIO_ENDPOINT.replace("http://", "").replace("https://", "")
    use_ssl = not MINIO_ENDPOINT.startswith("http://")
    self.con.execute("SET s3_access_key_id=?", [MINIO_ACCESS_KEY])
    self.con.execute("SET s3_secret_access_key=?", [MINIO_SECRET_KEY])
    self.con.execute("SET s3_endpoint=?", [ep])
    self.con.execute("SET s3_region='us-east-1'")
    self.con.execute("SET s3_use_ssl=?", [use_ssl])
    self.con.execute("SET s3_url_style='path'")

  def _download_if_exists(self, key, local):
    try:
      self.s3.download_file(self.source_bucket, key, local)
      return True
    except Exception:
      return False

  def _load_metadata(self):
    m1 = f"{self.base}/metadata.json"
    m2 = f"{self.tranches}/metadata.json"
    local = os.path.join(self.tmpdir, "metadata.json")
    if self._download_if_exists(m1, local) or self._download_if_exists(m2, local):
      with open(local, "r", encoding="utf-8") as f:
        data = json.load(f)
      return local, data
    raise RuntimeError(
      f"metadata.json not found at s3://{self.source_bucket}/{m1} or {m2}"
    )

  def _load_source_index(self):
    key = f"{self.tranches}/index.json"
    local = os.path.join(self.tmpdir, "index_src.json")
    if not self._download_if_exists(key, local):
      raise RuntimeError(f"index.json not found at s3://{self.source_bucket}/{key}")
    with open(local, "r", encoding="utf-8") as f:
      return json.load(f)

  def _group_by_tranche(self, inputs, index_dict):
    missing = [s for s in inputs if s not in index_dict]
    if missing:
      raise RuntimeError(
        f"inputs not indexed: {missing[:5]}{'...' if len(missing) > 5 else ''} total_missing={len(missing)}"
      )
    g = defaultdict(set)
    for s in inputs:
      r, c = index_dict[s]
      g[(int(r), int(c))].add(s)
    return g

  def _list_chunks(self, r, c):
    prefix = f"{self.tranches}/tranche_{r}_{c}/"
    paginator = self.s3.get_paginator("list_objects_v2")
    keys = []
    for page in paginator.paginate(Bucket=self.source_bucket, Prefix=prefix):
      for o in page.get("Contents", []):
        k = o["Key"]
        if k.endswith(".parquet") and "/chunk_" in k:
          keys.append(k)
    return sorted(keys)

  def _select_rows(self, r, c, wanted_set, input_col="input"):
    if not wanted_set:
      return pd.DataFrame()
    files_glob = (
      f"s3://{self.source_bucket}/{self.tranches}/tranche_{r}_{c}/chunk_*.parquet"
    )
    dfw = pd.DataFrame({input_col: list(wanted_set)})
    view = f"wanted_{r}_{c}_{uuid.uuid4().hex[:6]}"
    self.con.register(view, dfw)
    q = f"SELECT t.* FROM read_parquet('{files_glob}') t INNER JOIN {view} w ON t.{input_col}=w.{input_col}"
    out = self.con.execute(q).fetchdf()
    self.con.unregister(view)
    return out

  def _delete_tranches_tree(self):
    prefix = f"{self.tranches}/"
    paginator = self.s3.get_paginator("list_objects_v2")
    pages = paginator.paginate(Bucket=self.source_bucket, Prefix=prefix)
    batch = []
    deleted = 0
    for page in pages:
      for obj in page.get("Contents", []):
        batch.append({"Key": obj["Key"]})
        if len(batch) == 1000:
          self.s3.delete_objects(Bucket=self.source_bucket, Delete={"Objects": batch})
          deleted += len(batch)
          batch = []
    if batch:
      self.s3.delete_objects(Bucket=self.source_bucket, Delete={"Objects": batch})
      deleted += len(batch)
    logger.info(
      f"deleted objects under s3://{self.source_bucket}/{prefix} count={deleted}"
    )

  def _copy_to_buckets(self, meta_local, meta_list):
    priv = [
      (d.get("input") or "").strip()
      for d in meta_list
      if (d.get("access_level") or "").lower() == "private"
      and (d.get("input") or "").strip()
    ]
    pub = [
      (d.get("input") or "").strip()
      for d in meta_list
      if (d.get("access_level") or "").lower() == "public"
      and (d.get("input") or "").strip()
    ]
    logger.info(f"metadata: private={len(priv)} public={len(pub)}")
    src_index = self._load_source_index()
    logger.info(f"source index entries={len(src_index)}")
    priv_g = self._group_by_tranche(priv, src_index)
    pub_g = self._group_by_tranche(pub, src_index)
    dest_priv = (
      _SinkWriter(
        self.s3, "isaura-private", self.model_id, self.model_version, self.tmpdir
      )
      if priv_g
      else None
    )
    dest_pub = (
      _SinkWriter(
        self.s3, "isaura-public", self.model_id, self.model_version, self.tmpdir
      )
      if pub_g
      else None
    )
    total_priv = 0
    total_pub = 0
    for (r, c), want in priv_g.items():
      if not self._list_chunks(r, c):
        continue
      df = self._select_rows(r, c, want)
      total_priv += dest_priv.add_rows(r, c, df)
      logger.info(f"private tranche=({r},{c}) copied={len(df)}")
    for (r, c), want in pub_g.items():
      if not self._list_chunks(r, c):
        continue
      df = self._select_rows(r, c, want)
      total_pub += dest_pub.add_rows(r, c, df)
      logger.info(f"public tranche=({r},{c}) copied={len(df)}")
    if dest_priv:
      dest_priv.finalize(metadata_local=meta_local)
    if dest_pub:
      dest_pub.finalize(metadata_local=meta_local)
    logger.info(f"copy totals: private_new={total_priv} public_new={total_pub}")
    return total_priv, total_pub


class IsauraCopy(_BaseTransfer):
  def copy(self):
    meta_local, meta = self._load_metadata()
    return self._copy_to_buckets(meta_local, meta)


class IsauraRemover(_BaseTransfer):
  def remove(self):
    logger.info(
      f"removing tranches for {self.model_id}/{self.model_version} in {self.source_bucket}"
    )
    self._delete_tranches_tree()


class IsauraMover(_BaseTransfer):
  def move(self):
    meta_local, meta = self._load_metadata()
    self._copy_to_buckets(meta_local, meta)
    logger.info(
      f"wiping source tranches tree for {self.model_id}/{self.model_version} in {self.source_bucket}"
    )
    self._delete_tranches_tree()
