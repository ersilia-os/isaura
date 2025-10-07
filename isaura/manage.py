import csv, json, os, sys, tempfile, time, uuid

import pandas as pd
from isaura.base import _BaseTransfer, BloomIndex, TrancheState, MinioStore, DuckDBMinio
from isaura.helpers import (
  DEFAULT_BUCKET_NAME as PUB,
  DEFAULT_PRIVATE_BUCKET_NAME as PRI,
  STORE_DIRECTORY,
  CHECKPOINT_EVERY,
  BLOOM_FILENAME,
  MAX_ROWS,
  logger,
  tranche_coordinates,
  write_access_file,
  progress,
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
    acess_level=None,
    bucket=None,
    allowed_inputs=None,
    store=None,
  ):
    self.input_csv = input_csv
    self.model_id = model_id
    self.model_version = model_version
    self.access_level = acess_level
    self.bucket = bucket or PUB
    self.base_prefix = f"{model_id}/{model_version}/tranches"
    self.allowed_inputs = set(allowed_inputs) if allowed_inputs else None
    self.store = store or MinioStore()
    self.store.ensure_bucket(self.bucket)
    self.tmpdir = tempfile.mkdtemp(prefix="isaura_", dir=STORE_DIRECTORY)
    self.metadata = os.path.join(self.tmpdir, "access.json")
    self.bi = BloomIndex(
      store=self.store,
      bucket=self.bucket,
      base_prefix=self.base_prefix,
      local_dir=self.tmpdir,
      bloom_filename=BLOOM_FILENAME,
    )
    self.tranche = TrancheState(self.store, self.bucket, self.base_prefix, self.tmpdir, MAX_ROWS)
    self.buffers = defaultdict(list)
    self.schema_cols = None
    logger.info(f"writer init: bucket={self.bucket} base={self.base_prefix} csv={self.input_csv}")

  def _upload_metadata(self):
    if self.bucket not in (PUB, PRI):
      if self.access_level is None:
        logger.warning(
          f"You provide a custom project name but access type is not provided! Please provide access type as either of [public, private, both]. System exiting."
        )
        sys.exit(0)
    if self.bucket in (PUB, PRI):
      return
    try:
      write_access_file(self.input_csv, self.access_level, self.metadata)
      self.store.upload_file(self.metadata, self.bucket, f"{self.base_prefix}/access.json")
      logger.info(f"{self.metadata} -> s3://{self.bucket}/{self.base_prefix}/access.json")
    except Exception as e:
      logger.warning(f"metadata upload failed: {e}")

  def _set_schema(self, row):
    if self.schema_cols is None:
      self.schema_cols = list(row.keys())
      logger.info(f"writer schema: {self.schema_cols}")

  def _flush_if_needed(self, r, c):
    buf = self.buffers[(r, c)]
    if len(buf) >= MAX_ROWS:
      self.tranche.flush(r, c, buf, self.schema_cols)
      self.buffers[(r, c)].clear()

  def write(self):
    self._upload_metadata()
    total = dupes = 0
    with open(self.input_csv, newline="", encoding="utf-8") as f:
      total_rows = sum(1 for _ in csv.DictReader(f))
    desc = f"Writing {self.model_id}/{self.model_version} from {os.path.basename(self.input_csv)}"
    with progress(desc=desc, total_bytes=total_rows, transient=True) as (pg, tid):
      with open(self.input_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
          self._set_schema(row)
          smi = (row.get("input") or "").strip()
          if not smi:
            pg.advance(tid, 1)
            continue
          if self.allowed_inputs is not None and smi not in self.allowed_inputs:
            pg.advance(tid, 1)
            continue
          if self.bi.seen(smi):
            dupes += 1
            pg.advance(tid, 1)
            continue
          try:
            r, c, _, _ = tranche_coordinates(smi)
          except Exception:
            logger.warning("invalid SMILES skipped")
            pg.advance(tid, 1)
            continue
          self.buffers[(r, c)].append(dict(row))
          self.bi.register(smi, rc=(r, c))
          total += 1
          self._flush_if_needed(r, c)
          if self.bi._added >= CHECKPOINT_EVERY:
            self.bi.persist()
          pg.advance(tid, 1)
    for (r, c), buf in list(self.buffers.items()):
      if buf:
        self.tranche.flush(r, c, buf, self.schema_cols)
        self.buffers[(r, c)].clear()
    self.bi.persist()
    self.close()
    logger.info(f"write done: new={total} dupes={dupes}")

  def close(self):
    try:
      os.rmdir(self.tmpdir)
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
    bucket=None,
    store=None,
    endpoint=None,
    access=None,
    secret=None,
  ):
    self.model_id = model_id
    self.model_version = model_version
    self.input_csv = input_csv
    self.bucket = bucket or PUB
    self.base_prefix = f"{model_id}/{model_version}/tranches"
    self.index_key = f"{self.base_prefix}/index.json"
    self.store = store or MinioStore(endpoint=endpoint, access=access, secret=secret)
    self.tmpdir = tempfile.mkdtemp(prefix="isaura_reader_", dir=STORE_DIRECTORY)
    self.duck = DuckDBMinio(
      endpoint=self.store.endpoint,
      access=self.store.access,
      secret=self.store.secret,
    )
    logger.info(f"reader init bucket={self.bucket} base={self.base_prefix} csv={self.input_csv}")

  def _tranche_prefix(self, r, c):
    return f"{self.base_prefix}/tranche_{r}_{c}"

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

  def read(self, input_col="input", output_csv=None):
    t0 = time.time()
    wanted = []
    with open(self.input_csv, newline="", encoding="utf-8") as f:
      for row in csv.DictReader(f):
        v = (row.get(input_col) or "").strip()
        if v:
          wanted.append(v)
    if not wanted:
      return pd.DataFrame()
    index = self._load_index()
    missing = [s for s in wanted if s not in index]
    if missing:
      raise RuntimeError(
        f"inputs not indexed: {missing[:5]}{'...' if len(missing) > 5 else ''} total_missing={len(missing)}"
      )
    groups = defaultdict(set)
    for s in wanted:
      r, c = index[s]
      groups[(int(r), int(c))].add(s)
    results = []
    order_map = {s: i for i, s in enumerate(wanted)}
    for (r, c), smi_set in groups.items():
      files_glob = f"s3://{self.bucket}/{self._tranche_prefix(r, c)}/chunk_*.parquet"
      dfw = pd.DataFrame({input_col: list(smi_set)})
      view = f"w_{r}_{c}_{uuid.uuid4().hex[:6]}"
      self.duck.con.register(view, dfw)
      q = f"SELECT t.* FROM read_parquet('{files_glob}') t INNER JOIN {view} w ON t.{input_col}=w.{input_col}"
      part = self.duck.con.execute(q).fetchdf()
      self.duck.con.unregister(view)
      if not part.empty:
        part["__o"] = part[input_col].map(order_map)
        results.append(part)
    out = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
    if not out.empty and "__o" in out.columns:
      out = out.sort_values("__o").drop(columns="__o").reset_index(drop=True)
    if output_csv:
      out.to_csv(output_csv, index=False)
    elapsed = time.time() - t0
    rate = (len(out) / elapsed) if elapsed > 0 and len(out) else 0.0
    logger.info(f"read matched={len(out)} elapsed={elapsed:.2f}s rate={rate:.1f}/s")
    return out


class IsauraInspect:
  def __init__(self, model_id, model_version, project_name=None, access="both"):
    self.mid, self.mv, self.proj, self.acc = model_id, model_version, project_name, access
    self.base = f"{self.mid}/{self.mv}/tranches"
    self.idx_key = f"{self.base}/index.json"
    self.s = MinioStore()
    logger.info(
      f"inspect init model={self.mid} version={self.mv} project={self.proj} access={self.acc}"
    )

  def _buckets(self):
    return (
      [self.proj]
      if self.proj
      else ([PUB] if self.acc == "public" else [PRI] if self.acc == "private" else [PUB, PRI])
    )

  def _idx(self, b):
    o = self.s.client.get_object(Bucket=b, Key=self.idx_key)
    d = json.loads(o["Body"].read().decode("utf-8"))
    logger.info(f"loaded index: bucket={b} entries={len(d)}")
    return d

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

  def inspect_inputs(self, input_csv, output_csv=None):
    own = self._union()
    with open(input_csv, newline="", encoding="utf-8") as f:
      wanted = [
        (r.get("input") or "").strip() for r in csv.DictReader(f) if (r.get("input") or "").strip()
      ]
    logger.info(f"parsed inputs csv={input_csv} count={len(wanted)}")
    df = pd.DataFrame([
      {"input": s, "available": s in own, "bucket": own.get(s, "")} for s in wanted
    ])
    if output_csv:
      df.to_csv(output_csv, index=False)
      logger.info(f"inspect inputs wrote={len(df)} path={output_csv}")
    return df

  def inspect_models(self, project_name, prefix_filter=""):
    c = self.s.client
    rows = []
    p = c.get_paginator("list_objects_v2")

    def list_prefixes(pref):
      for page in p.paginate(Bucket=project_name, Prefix=pref, Delimiter="/"):
        for cp in page.get("CommonPrefixes", []):
          yield cp["Prefix"]

    def count_chunks(base):
      tr = set()
      ch = 0
      for page in p.paginate(Bucket=project_name, Prefix=f"{base}/tranches/"):
        for obj in page.get("Contents", []):
          k = obj["Key"]
          if "/chunk_" in k and k.endswith(".parquet"):
            ch += 1
            i = k.find("tranche_")
            if i != -1:
              tr.add(k[i:].split("/")[0])
      return len(tr), ch

    for m_pref in list_prefixes(""):
      model = m_pref.strip("/")
      if prefix_filter and not model.startswith(prefix_filter):
        continue
      for v_pref in list_prefixes(m_pref):
        ver = v_pref[len(m_pref) :].strip("/")
        idx_key = f"{model}/{ver}/tranches/index.json"
        try:
          obj = c.get_object(Bucket=project_name, Key=idx_key)
          idx = json.loads(obj["Body"].read().decode("utf-8"))
          entries = len(idx)
        except Exception:
          entries = 0
        tr, ch = count_chunks(f"{model}/{ver}")
        rows.append({"model": f"{model}/{ver}", "entries": entries, "tranches": tr, "chunks": ch})
    return rows


class IsauraCopy(_BaseTransfer):
  def copy(self):
    meta_local, meta = self._load_metadata()
    return self._copy_to_buckets(meta_local, meta)


class IsauraMover(_BaseTransfer):
  def move(self):
    self.copy()
    n = self._delete_tranches_tree()
    logger.info(f"move wiped objects={n}")


class IsauraRemover(_BaseTransfer):
  def remove(self):
    n = self._delete_tranches_tree()
    logger.info(f"removed objects={n}")
