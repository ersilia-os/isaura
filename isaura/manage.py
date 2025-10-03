import csv, json, os, tempfile, time, uuid

import pandas as pd
from isaura.base import _BaseTransfer, BloomIndex, TrancheState, S3Store, DuckDBS3
from isaura.helpers import (
  DEFAULT_BUCKET_NAME,
  STORE_DIRECTORY,
  MAX_ROWS_PER_FILE,
  CHECKPOINT_EVERY,
  BLOOM_FILENAME,
  logger,
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
    self.store = store or S3Store()
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


class IsauraWriter(AbstractContextManager):
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
  ):
    self.input_csv = input_csv
    self.model_id = model_id
    self.model_version = model_version
    self.access_level = acess_level
    self.bucket = bucket or DEFAULT_BUCKET_NAME
    self.base_prefix = f"{model_id}/{model_version}/tranches"
    self.metadata_path = metadata_path
    self.allowed_inputs = set(allowed_inputs) if allowed_inputs else None
    self.max_rows = MAX_ROWS_PER_FILE
    self.tmpdir = tempfile.mkdtemp(prefix="isaura_", dir=STORE_DIRECTORY)
    self.store = store or S3Store()
    self.store.ensure_bucket(self.bucket)
    self.checker = IsauraChecker(
      self.bucket, self.base_prefix, store=self.store, store_directory=STORE_DIRECTORY
    )
    self.bi = self.checker.bi
    self.tranche = TrancheState(
      self.store, self.bucket, self.base_prefix, self.tmpdir, self.max_rows
    )
    self.buffers = defaultdict(list)
    logger.info(f"writer init bucket={self.bucket} base={self.base_prefix} csv={self.input_csv}")

  def _upload_metadata(self):
    if not self.metadata_path:
      return
    try:
      self.store.upload_file(self.metadata_path, self.bucket, f"{self.base_prefix}/metadata.json")
      logger.info("uploaded metadata.json")
    except Exception as e:
      logger.warning(f"metadata upload failed: {e}")

  def _flush_if_needed(self, r, c):
    if len(self.buffers[(r, c)]) >= self.max_rows:
      self.tranche.flush(r, c, self.buffers[(r, c)])
      self.buffers[(r, c)].clear()

  def write(self):
    self._upload_metadata()
    total = 0
    dupes = 0
    with open(self.input_csv, newline="", encoding="utf-8") as f:
      for row in csv.DictReader(f):
        smi = (row.get("input") or "").strip()
        if not smi:
          continue
        if self.allowed_inputs is not None and smi not in self.allowed_inputs:
          continue
        if self.checker.seen(smi):
          dupes += 1
          continue
        try:
          r, c = tranche_coordinates(smi)
        except Exception:
          continue
        self.buffers[(r, c)].append(dict(row))
        self.checker.register(smi, rc=(r, c))
        total += 1
        self._flush_if_needed(r, c)
        if self.bi._added >= CHECKPOINT_EVERY:
          self.bi.persist()
    for r, c in list(self.buffers.keys()):
      if self.buffers[(r, c)]:
        self.tranche.flush(r, c, self.buffers[(r, c)])
        self.buffers[(r, c)].clear()
    self.bi.persist()
    logger.info(f"write done new={total} dupes={dupes}")

  def close(self):
    try:
      os.rmdir(self.tmpdir)
    except:
      pass
    self.checker.close()

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
    region=None,
  ):
    self.model_id = model_id
    self.model_version = model_version
    self.input_csv = input_csv
    self.bucket = bucket or DEFAULT_BUCKET_NAME
    self.base_prefix = f"{model_id}/{model_version}/tranches"
    self.index_key = f"{self.base_prefix}/index.json"
    self.store = store or S3Store(endpoint=endpoint, access=access, secret=secret, region=region)
    self.tmpdir = tempfile.mkdtemp(prefix="isaura_reader_", dir=STORE_DIRECTORY)
    self.duck = DuckDBS3(
      endpoint=self.store.endpoint,
      access=self.store.access,
      secret=self.store.secret,
      region=self.store.region,
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


class IsauraCleaner:
  def __init__(self, model_id, model_version, store=None):
    self.model_id = model_id
    self.model_version = model_version
    self.store = store or S3Store()

  def clean(self):
    total = {}
    for b in ["isaura-public", "isaura-private"]:
      base = f"{self.model_id}/{self.model_version}/tranches/"
      n = self.store.delete_prefix(b, base)
      total[b] = n
    return total
