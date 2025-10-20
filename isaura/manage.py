import csv, json, os, shutil, time, uuid

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
  logger,
  filter_out,
  get_apprx,
  get_base,
  get_desc,
  get_files,
  get_idx_key,
  get_pref,
  group_inputs,
  make_temp,
  progress,
  query,
  write_access_file,
  tranche_coordinates,
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
  ):
    self.input_csv = input_csv
    self.model_id = model_id
    self.model_version = model_version
    self.access = access
    self.bucket = bucket or PUB
    self.base_prefix = get_base(self.model_id, model_version)
    self.max_rows = int(MAX_ROWS)
    self.store = MinioStore()
    self.store.ensure_bucket(self.bucket)
    self.tmpdir = make_temp("isaura_")
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

  def _upload_metadata(self):
    if not self.metadata_path or not self.bucket:
      return
    try:
      write_access_file(self.input_csv, self.access, self.metadata_path)
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
        if self.bi._added >= int(CHECKPOINT_EVERY):
          self.bi.persist()
    for (r, c), buf in list(self.buffers.items()):
      if buf:
        self.tranche.flush(r, c, buf, self.schema_cols)
        self.buffers[(r, c)].clear()
    self.bi.persist()
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
  ):
    self.model_id = model_id
    self.model_version = model_version
    self.approximate = approximate
    self.input_csv = input_csv
    self.bucket = bucket or PUB
    self.base = get_base(self.model_id, model_version)
    self.index_key = get_idx_key(self.base)
    self.pref = get_pref(self.model_id, self.model_version)
    self.store = MinioStore()
    self.tmpdir = make_temp("isaura_reader_")
    self.duck = DuckDBMinio()
    logger.info(f"reader init bucket={self.bucket} base={self.base} csv={self.input_csv}")

  def _hive_prefix(self, r, c):
    return 

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

  def _sizes_for_groups(self, groups):
    sizes, total = {}, 0
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

  def read(self, output_csv=None):
    t0, wanted, header, objc, results = time.time(), [], set(), "__o", []
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
      wanted = get_apprx(wanted, self.pref)
      et = time.perf_counter()
      logger.info(f"Approximate inputs are retrieved {len(wanted)} in {et - st:.2f} seconds!")

    header = list(header)[0]
    if not wanted:
      return pd.DataFrame()
    index = self._load_index()
    groups = group_inputs(wanted, index)
    sizes, total_bytes = self._sizes_for_groups(groups)
    order_map = {s: i for i, s in enumerate(wanted)}
    with progress(get_desc(self.pref, wanted), total_bytes or 0, transient=True) as (prog, task_id):
      for (r, c), _ in groups.items():
        part = self.duck.con.execute(query(get_files(self.bucket, r, c, self.base))).fetchdf()
        if not part.empty:
          part[objc] = part[header].map(order_map)
          results.append(part)
        if total_bytes:
          prog.update(task_id, advance=sizes.get((r, c), 0))
    out = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
    if not out.empty and objc in out.columns:
      out = out.sort_values(objc).drop(columns=objc).reset_index(drop=True)
    if output_csv:
      filter_out(out, objc, wanted, header).to_csv(output_csv, index=False)
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
  def __init__(self, model_id, model_version, project_name=None, access="both"):
    self.mid, self.mv, self.proj, self.acc = model_id, model_version, project_name, access
    self.base = get_base(self.mid, self.mv)
    self.idx_key = get_idx_key(self.base)
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
      for page in p.paginate(Bucket=project_name, Prefix=base):
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
        try:
          obj = c.get_object(Bucket=project_name, Key=get_idx_key(self.base))
          idx = json.loads(obj["Body"].read().decode("utf-8"))
          entries = len(idx)
        except Exception:
          entries = 0
        tr, ch = count_chunks(get_pref(model, ver))
        rows.append({"model": get_pref(model, ver), "entries": entries, "tranches": tr, "chunks": ch})
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


class IsauraPull(_BaseTransfer):
  def remove(self):
    n = self._delete_tranches_tree()
    logger.info(f"removed objects={n}")
