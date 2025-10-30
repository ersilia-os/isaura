import csv, json, os, shutil, sys, time, uuid

import pandas as pd
from isaura.base import _BaseTransfer, BloomIndex, TrancheState, MinioStore, DuckDBMinio
from isaura.helpers import (
  ACCESS_FILE,
  BLOOM_FILENAME,
  CHECKPOINT_EVERY,
  DEFAULT_BUCKET_NAME as PUB,
  DEFAULT_PRIVATE_BUCKET_NAME as PRI,
  INPUT_C,
  MAX_ROWS,
  MINIO_ENDPOINT_CLOUD,
  MINIO_CLOUD_AK as mcak,
  MINIO_CLOUD_SK as mcsk,
  MINIO_PRIV_CLOUD_AK as mcpak,
  MINIO_PRIV_CLOUD_SK as mcpsk,
  logger,
  get_acc_key,
  get_apprx,
  get_base,
  get_coll,
  get_files_glob,
  get_idx_key,
  get_pref,
  group_inputs,
  make_temp,
  query,
  split_csv,
  spinner,
  tranche_coordinates,
  write_access_file,
)


from collections import defaultdict
from contextlib import AbstractContextManager


class IsauraChecker(AbstractContextManager):
  def __init__(
    self,
    bucket,
    base_prefix,
    store=None,
    store_directory=".",
    bloom_filename=BLOOM_FILENAME,
    checkpoint_every=CHECKPOINT_EVERY,
  ):
    self.bucket = bucket
    self.base_prefix = base_prefix
    self.store = store or MinioStore()
    self.local_dir = store_directory
    self.checkpoint_every = checkpoint_every
    self.bi = BloomIndex(
      self.store,
      self.bucket,
      self.base_prefix,
      self.local_dir,
      bloom_filename=bloom_filename,
    )
    self._added = 0

  def seen(self, v):
    return self.bi.seen(v)

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
    self.max_rows = int(MAX_ROWS)
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
    self.tranche = TrancheState(self.store, self.bucket, self.base_prefix, self.tmpdir, self.max_rows)
    self.buffers = defaultdict(list)
    self.schema_cols = None
    logger.info(f"writer init: bucket={self.bucket} base={self.base_prefix} csv={self.input_csv}")

  def _load_metadata(self):
    local = os.path.join(self.tmpdir, f"{uuid.uuid4().hex}.json")
    try:
      self.store.download_file(self.bucket, self.access_key, local)
      with open(local, "r", encoding="utf-8") as f:
        return json.load(f)
    except Exception as e:
      return None
    finally:
      try:
        os.remove(local)
      except:
        pass

  def _upload_metadata(self, inputs):
    if not self.metadata_path or not self.bucket:
      return
    try:
      existed = self._load_metadata()
      write_access_file(existed, inputs, self.access, self.metadata_path)
      self.store.upload_file(self.metadata_path, self.bucket, f"{self.base_prefix}/{ACCESS_FILE}")
      logger.info(f"{ACCESS_FILE} -> s3://{self.bucket}/{self.base_prefix}/{ACCESS_FILE}")
    except Exception as e:
      logger.warning(f"metadata upload failed: {e}")

  def _set_schema(self, row):
    if self.schema_cols is None:
      self.schema_cols = list(row.keys())
      logger.info(f"writer schema: {self.schema_cols[:10]}")

  def _flush_if_needed(self, r, c):
    buf = self.buffers[(r, c)]
    if len(buf) >= self.max_rows:
      self.tranche.flush(buf, self.schema_cols)
      self.buffers[(r, c)].clear()

  def write(self, df=None):
    total = dupes = 0
    new = []
    rows = df.to_dict("records") if df is not None else None
    if rows is None:
      with open(self.input_csv, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    for row in rows:
      self._set_schema(row)
      smi = (row.get("input") or row.get("smiles") or "").strip()
      if not smi:
        continue
      if self.bi.seen(smi):
        dupes += 1
        continue
      try:
        r, c, _, _ = tranche_coordinates(smi)
        if r < 1 or c < 1:
          raise ValueError("Tranche coordinates must be >= 1")
      except Exception:
        logger.warning("invalid SMILES or tranche mapping; skipped")
        continue
      new.append(smi)
      self.buffers[(r, c)].append(dict(row))
      self.bi.register(smi, rc=(r, c))
      total += 1
      self._flush_if_needed(r, c)
      if self.bi._added >= int(CHECKPOINT_EVERY):
        self.bi.persist()
    for (r, c), buf in list(self.buffers.items()):
      if buf:
        self.tranche.flush(buf, self.schema_cols)
        self.buffers[(r, c)].clear()
    self.bi.persist()
    if new:
      self._upload_metadata(new)
    logger.info(f"write done: new={total} dupes={dupes}")

  def close(self):
    try:
      shutil.rmtree(self.tmpdir)
    except:
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
    logger.info(f"reader init bucket={self.bucket} base={self.base} csv={self.input_csv}")

  def _load_index(self):
    local = os.path.join(self.tmpdir, f"{uuid.uuid4().hex}.json")
    self.store.download_file(self.bucket, self.index_key, local)
    try:
      with open(local, "r", encoding="utf-8") as f:
        return json.load(f)
    finally:
      try:
        os.remove(local)
      except:
        pass

  def read(self, output_csv=None, df=None):
    t0, wanted, header_set = time.time(), [], set()
    rows = df.to_dict("records") if df is not None else None
    if rows is None:
      with open(self.input_csv, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    for row in rows:
      h = INPUT_C[0] if row.get(INPUT_C[0]) else INPUT_C[1]
      v = (row.get(h) or "").strip()
      if v:
        wanted.append(v)
      if h:
        header_set.add(h)
    if self.approximate:
      st = time.perf_counter()
      wanted = get_apprx(wanted, self.collection)
      et = time.perf_counter()
      logger.info(f"Approximate inputs are retrieved {len(wanted)} in {et - st:.2f} seconds!")
    header = list(header_set)[0] if header_set else "smiles"
    index = self._load_index()
    group_inputs(wanted, index)
    if not wanted:
      return pd.DataFrame()
    try:
      files = get_files_glob(self.bucket, self.base)
      out = spinner("Fetching queries. Please wait!", query, self.duck.con, header, wanted, files)
    except Exception as e:
      logger.error(e)
      sys.exit(1)
    if output_csv:
      out.to_csv(output_csv, index=False)
      logger.info(f"wrote csv rows={len(out)} path={output_csv}")
    elapsed = time.time() - t0
    rate = (len(out) / elapsed) if elapsed > 0 and len(out) else 0.0
    logger.info(
      f"read done model+version={self.pref}"
      f"bucket={self.bucket} inputs={len(wanted)} matched={len(out)} "
      f"elapsed={elapsed:.2f}s rate={rate:.1f}/s"
    )
    return out


class IsauraInspect:
  def __init__(self, model_id, model_version, cloud, project_name=None, access="both"):
    self.mid, self.mv, self.proj, self.acc, self.cloud = model_id, model_version, project_name, access, cloud
    self.base = get_base(self.mid, self.mv)
    self.idx_key = get_idx_key(self.base)
    endpoint = MINIO_ENDPOINT_CLOUD if cloud else None
    self.s = MinioStore(endpoint=endpoint)
    logger.info(f"inspect init model={self.mid} version={self.mv} project={self.proj} access={self.acc}")

  def _buckets(self):
    return (
      [self.proj]
      if self.proj
      else ([PUB] if self.acc == "public" else [PRI] if self.acc == "private" else [PUB, PRI])
    )

  def _idx(self, b):
    try:
      if self.cloud and b == PUB:
        self.s = MinioStore(endpoint=MINIO_ENDPOINT_CLOUD, access=mcak, secret=mcsk)
      if self.cloud and b == PRI:
        self.s = MinioStore(endpoint=MINIO_ENDPOINT_CLOUD, access=mcpak, secret=mcpsk)
      o = self.s.client.get_object(Bucket=b, Key=self.idx_key)
      d = json.loads(o["Body"].read().decode("utf-8"))
      logger.info(f"loaded index: bucket={b} entries={len(d)}")
      return d
    except Exception as e:
      logger.warning(e)
      return None

  def _union(self):
    if self.proj:
      try:
        d = self._idx(self.proj)
      except Exception as e:
        raise RuntimeError(f"project bucket not available or missing index.json: {self.proj} ({e})")
      return {s: self.proj for s in d.keys()}
    own = {}
    for b in self._buckets():
      try:
        for s in self._idx(b).keys():
          own[s] = b
      except Exception as e:
        logger.info(f"no index in {b}: {e}")
    return own

  def _candidate_buckets(self):
    if self.proj:
      return [self.proj]
    if self.acc == "public":
      return [PUB]
    if self.acc == "private":
      return [PRI]
    return [PUB, PRI]

  def _load_indices_union(self):
    buckets = self._candidate_buckets()
    union = {}
    owners = {}
    for b in buckets:
      idx = self._idx(b)
      if idx:
        for smi in idx.keys():
          if smi not in union:
            union[smi] = True
            owners[smi] = b
    return union, owners

  def list_available(self, output_csv=None):
    _, owner = self._load_indices_union()
    rows = [{"input": smi, "bucket": b} for smi, b in owner.items()]
    df = pd.DataFrame(rows)
    if output_csv:
      df.to_csv(output_csv, index=False)
      logger.info(f"inspect list wrote={len(df)} path={output_csv}")
    return df

  def inspect_inputs(self, input_csv, output_csv=None):
    own = self._union()
    with open(input_csv, newline="", encoding="utf-8") as f:
      wanted = [(r.get("input") or "").strip() for r in csv.DictReader(f) if (r.get("input") or "").strip()]
    logger.info(f"parsed inputs csv={input_csv} count={len(wanted)}")
    df = pd.DataFrame([{"input": s, "available": s in own, "bucket": own.get(s, "")} for s in wanted])
    if output_csv:
      df.to_csv(output_csv, index=False)
      logger.info(f"inspect inputs wrote={len(df)} path={output_csv}")
    return df

  def inspect_models(self, project_name, prefix_filter=""):
    if self.cloud and project_name == PUB:
      self.s = MinioStore(endpoint=MINIO_ENDPOINT_CLOUD, access=mcak, secret=mcsk)
    if self.cloud and project_name == PRI:
      self.s = MinioStore(endpoint=MINIO_ENDPOINT_CLOUD, access=mcpak, secret=mcpsk)
    c = self.s.client
    rows = []
    p = c.get_paginator("list_objects_v2")

    def list_prefixes(pref):
      try:
        for page in p.paginate(Bucket=project_name, Prefix=pref, Delimiter="/"):
          for cp in page.get("CommonPrefixes", []):
            yield cp["Prefix"]
      except Exception as e:
        logger.error(e)

    def count_chunks(base):
      tr = set()
      ch = 0
      for page in p.paginate(Bucket=project_name, Prefix=base):
        for obj in page.get("Contents", []):
          k = obj["Key"]
          if "/chunk_" in k and k.endswith(".parquet"):
            ch += 1
            i = k.find("row=")
            if i != -1:
              tr.add(k[i:].split("/")[0])
      return len(tr), ch

    for m_pref in list_prefixes(""):
      model = m_pref.strip("/")
      if prefix_filter and not model.startswith(prefix_filter):
        continue
      for v_pref in list_prefixes(m_pref):
        ver = v_pref[len(m_pref) :].strip("/")
        try:
          obj = c.get_object(Bucket=project_name, Key=get_idx_key(get_base(model, ver)))
          idx = json.loads(obj["Body"].read().decode("utf-8"))
          entries = len(idx)
        except Exception as e:
          logger.error(e)
          entries = 0
        tr, ch = count_chunks(get_pref(model, ver))
        rows.append({"model": get_pref(model, ver), "entries": entries, "rows": tr, "chunks": ch})
    return rows


class IsauraCopy(_BaseTransfer):
  def copy(self):
    meta_local, meta = self._load_metadata()
    return self._copy(meta_local, meta)


class IsauraMover(_BaseTransfer):
  def move(self):
    meta_local, meta = self._load_metadata()
    self._copy(meta_local, meta)
    n = self._delete()
    logger.info(f"move wiped objects={n}")


class IsauraRemover(_BaseTransfer):
  def remove(self):
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
    out = self._pull(r.read(), self._load_index())
    logger.info(f"pulled objects={out}")


class IsauraPush:
  def __init__(self, model_id, model_version, bucket):
    self.model_id = model_id
    self.model_version = model_version
    self.bucket = bucket

  def push(self):
    insp = IsauraInspect(
      model_id=self.model_id,
      model_version=self.model_version,
      access="both",
      cloud=False,
    )
    df = insp.list_available()
    files = split_csv(df)

    if not files:
      logger.error("No data found in any default bucket! Aborting pull.")
      sys.exit(1)

    file1 = files[0]
    file2 = files[1] if len(files) > 1 else None

    if not file2:
      logger.warning("Private bucket has no data! Skipping pull for it.")

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
      df = r.read()
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
        w.write(df=df)
