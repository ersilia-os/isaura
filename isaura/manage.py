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
  INPUT_C,
  logger,
  tranche_coordinates,
  write_access_file,
  progress,
  get_apprx,
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
    acess_level="public",
    bucket=None,
    allowed_inputs=None,
    metadata_path=None,
    store=None,
    use_hive=True,
  ):
    self.input_csv = input_csv
    self.model_id = model_id
    self.model_version = model_version
    self.access_level = acess_level
    self.bucket = bucket or os.getenv("DEFAULT_BUCKET_NAME", "isaura-public")
    self.base_prefix = f"{model_id}/{model_version}/tranches"
    self.metadata_path = metadata_path
    self.allowed_inputs = set(allowed_inputs) if allowed_inputs else None
    self.max_rows = int(os.getenv("MAX_ROWS_PER_FILE", str(MAX_ROWS)))
    self.use_hive = use_hive
    self.store = store or MinioStore()
    self.store.ensure_bucket(self.bucket)
    self.tmpdir = tempfile.mkdtemp(prefix="isaura_", dir=os.getenv("STORE_DIRECTORY", STORE_DIRECTORY))
    self.bi = BloomIndex(
      self.store,
      self.bucket,
      self.base_prefix,
      self.tmpdir,
      bloom_filename=os.getenv("BLOOM_FILENAME", BLOOM_FILENAME),
    )
    self.tranche = TrancheState(
      self.store, self.bucket, self.base_prefix, self.tmpdir, self.max_rows, use_hive=self.use_hive
    )
    self.buffers = defaultdict(list)
    self.schema_cols = None
    logger.info(
      f"writer init: bucket={self.bucket} base={self.base_prefix} hive={self.use_hive} csv={self.input_csv}"
    )

  def _upload_metadata(self):
    if not self.metadata_path:
      return
    try:
      self.store.upload_file(self.metadata_path, self.bucket, f"{self.base_prefix}/metadata.json")
      logger.info(f"metadata.json -> s3://{self.bucket}/{self.base_prefix}/metadata.json")
    except Exception as e:
      logger.warning(f"metadata upload failed: {e}")

  def _set_schema(self, row):
    if self.schema_cols is None:
      self.schema_cols = list(row.keys())
      logger.info(f"writer schema: {self.schema_cols[:10]}")

  def _flush_if_needed(self, r, c):
    buf = self.buffers[(r, c)]
    if len(buf) >= self.max_rows:
      self.tranche.flush(r, c, buf, self.schema_cols)
      self.buffers[(r, c)].clear()

  def write(self):
    self._upload_metadata()
    total = dupes = 0
    with open(self.input_csv, newline="", encoding="utf-8") as f:
      reader = csv.DictReader(f)
      for row in reader:
        self._set_schema(row)
        smi = (row.get("input") or row.get("smiles")).strip()
        if not smi:
          continue
        if self.allowed_inputs is not None and smi not in self.allowed_inputs:
          continue
        if self.bi.seen(smi):
          dupes += 1
          continue
        try:
          r, c, _, _ = tranche_coordinates(smi)
        except Exception:
          logger.warning("invalid SMILES skipped")
          continue
        self.buffers[(r, c)].append(dict(row))
        self.bi.register(smi, rc=(r, c))
        total += 1
        self._flush_if_needed(r, c)
        if self.bi._added >= int(os.getenv("CHECKPOINT_EVERY", str(CHECKPOINT_EVERY))):
          self.bi.persist()
    for (r, c), buf in list(self.buffers.items()):
      if buf:
        self.tranche.flush(r, c, buf, self.schema_cols)
        self.buffers[(r, c)].clear()
    self.bi.persist()
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
    approximate,
    bucket=None,
    store=None,
    endpoint=None,
    access=None,
    secret=None,
  ):
    self.model_id = model_id
    self.approximate = approximate
    self.model_version = model_version
    self.input_csv = input_csv
    self.bucket = bucket or PUB
    self.base = f"{model_id}/{model_version}/tranches"
    self.index_key = f"{self.base}/index.json"
    self.store = store or MinioStore(endpoint=endpoint, access=access, secret=secret)
    self.tmpdir = tempfile.mkdtemp(prefix="isaura_reader_", dir=STORE_DIRECTORY)
    self.duck = DuckDBMinio(endpoint=self.store.endpoint, access=self.store.access, secret=self.store.secret)
    logger.info(f"reader init bucket={self.bucket} base={self.base} csv={self.input_csv}")

  def _hive_prefix(self, r, c):
    return f"{self.base}/row={r}/col={c}"

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

  def _group_inputs(self, wanted, index):
    miss = [s for s in wanted if s not in index]
    if miss:
      raise RuntimeError(
        f"inputs not indexed: {miss[:5]}{'...' if len(miss) > 5 else ''} total_missing={len(miss)}"
      )
    g = defaultdict(set)
    for s in wanted:
      r, c = index[s]
      g[(int(r), int(c))].add(s)
    return g

  def _sizes_for_groups(self, groups):
    sizes = {}
    total = 0
    for r, c in groups.keys():
      pref = self._hive_prefix(r, c) + "/"
      s = 0
      for obj in self.store.list_keys(self.bucket, pref):
        k = obj["Key"]
        if k.endswith(".parquet") and "/chunk_" in k:
          s += int(obj.get("Size", 0))
      sizes[(r, c)] = s
      total += s
    return sizes, total

  def read(self, output_csv=None, transient_progress=True):
    t0, wanted, header = time.time(), [], set()
    with open(self.input_csv, newline="", encoding="utf-8") as f:
      for row in csv.DictReader(f):
        h = INPUT_C[0] if row.get(INPUT_C[0]) else INPUT_C[1]
        v = (row.get(h)).strip()
        if v:
          wanted.append(v)
        if h not in header:
          header.add(h)
    if self.approximate:
      st = time.perf_counter()
      wanted = get_apprx(wanted)
      et = time.perf_counter()
      logger.info(f"Approximate inputs are retrieved {len(wanted)} in {et - st:.2f} seconds!")

    header = list(header)[0]
    if not wanted:
      return pd.DataFrame()
    index = self._load_index()
    groups = self._group_inputs(wanted, index)
    sizes, total_bytes = self._sizes_for_groups(groups)
    order_map = {s: i for i, s in enumerate(wanted)}
    results = []
    desc = f"Fetching hive partitions {self.model_id}/{self.model_version} ({len(wanted)} inputs)"
    with progress(desc, total_bytes or 0, transient=transient_progress) as (prog, task_id):
      for (r, c), _ in groups.items():
        files = f"s3://{self.bucket}/{self._hive_prefix(r, c)}/chunk_*.parquet"
        q = f"SELECT * FROM read_parquet('{files}', hive_partitioning=1)"
        part = self.duck.con.execute(q).fetchdf()
        if not part.empty:
          part["__o"] = part[header].map(order_map)
          results.append(part)
        if total_bytes:
          prog.update(task_id, advance=sizes.get((r, c), 0))
    out = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
    if not out.empty and "__o" in out.columns:
      out = out.sort_values("__o").drop(columns="__o").reset_index(drop=True)
    if output_csv:
      out = (
        out[out[header].isin(wanted)]
        .assign(__o=lambda d: d[header].map({s: i for i, s in enumerate(wanted)}))
        .sort_values("__o")
        .drop(columns=["__o", "row", "col"], errors="ignore")
        .reset_index(drop=True)
      )
      out.to_csv(output_csv, index=False)
      logger.info(f"wrote csv rows={len(out)} path={output_csv}")
    elapsed = time.time() - t0
    rate = (len(out) / elapsed) if elapsed > 0 and len(out) else 0.0
    logger.info(
      f"read done model={self.model_id} version={self.model_version} "
      f"bucket={self.bucket} inputs={len(wanted)} matched={len(out)} "
      f"elapsed={elapsed:.2f}s rate={rate:.1f}/s"
    )
    return out


class IsauraInspect:
  def __init__(self, model_id, model_version, project_name=None, access="both"):
    self.mid, self.mv, self.proj, self.acc = model_id, model_version, project_name, access
    self.base = f"{self.mid}/{self.mv}/tranches"
    self.idx_key = f"{self.base}/index.json"
    self.s = MinioStore()
    logger.info(f"inspect init model={self.mid} version={self.mv} project={self.proj} access={self.acc}")

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
      wanted = [(r.get("input") or "").strip() for r in csv.DictReader(f) if (r.get("input") or "").strip()]
    logger.info(f"parsed inputs csv={input_csv} count={len(wanted)}")
    df = pd.DataFrame([{"input": s, "available": s in own, "bucket": own.get(s, "")} for s in wanted])
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
