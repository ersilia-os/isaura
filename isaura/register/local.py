import boto3, csv, duckdb, json, os, pickle, tempfile, time, uuid

import pandas as pd
import pyarrow.parquet as pq

from isaura.helpers import (
  MINIO_ENDPOINT,
  MINIO_ACCESS_KEY,
  MINIO_SECRET_KEY,
  STORE_DIRECTORY,
  MAX_ROWS_PER_FILE,
  CHECKPOINT_EVERY,
  BLOOM_FILENAME,
  DEFAULT_BUCKET_NAME,
  tranche_coordinates,
  logger,
  make_fetching_progress,
)

from botocore.config import Config
from collections import defaultdict
from contextlib import AbstractContextManager

from pybloom_live import ScalableBloomFilter


class IsauraChecker(AbstractContextManager):
  def __init__(
    self,
    bucket,
    base_prefix,
    store_directory=".",
    bloom_filename="bloom.pkl",
    checkpoint_every=50000,
    error_rate=0.001,
    initial_capacity=1000000,
    s3=None,
  ):
    self.bucket = bucket
    self.base_prefix = base_prefix.strip("/")
    self.bloom_key = f"{self.base_prefix}/bloom.pkl"
    self.local_dir = store_directory
    self.local_bloom = os.path.join(store_directory, bloom_filename)
    self.checkpoint_every = checkpoint_every
    self._added_since_save = 0
    self.s3 = s3
    if not os.path.isdir(self.local_dir):
      os.makedirs(self.local_dir, exist_ok=True)
    try:
      self.s3.head_object(Bucket=self.bucket, Key=self.bloom_key)
      self.s3.download_file(self.bucket, self.bloom_key, self.local_bloom)
      logger.info(
        f"downloaded bloom from s3: s3://{self.bucket}/{self.bloom_key} -> {self.local_bloom}"
      )
    except Exception as e:
      logger.info(f"no remote bloom found or failed to download: {e}")
    try:
      with open(self.local_bloom, "rb") as f:
        self.sbf = pickle.load(f)
      logger.info(f"loaded bloom locally: {self.local_bloom}")
    except FileNotFoundError:
      self.sbf = ScalableBloomFilter(
        mode=ScalableBloomFilter.SMALL_SET_GROWTH,
        initial_capacity=initial_capacity,
        error_rate=error_rate,
      )
      logger.info(
        f"created new bloom: initial_capacity={initial_capacity} error_rate={error_rate}"
      )

  def _persist(self):
    tmp = f"{self.local_bloom}.tmp"
    with open(tmp, "wb") as f:
      pickle.dump(self.sbf, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, self.local_bloom)
    logger.info(f"saved bloom locally: {self.local_bloom}")
    try:
      self.s3.upload_file(self.local_bloom, self.bucket, self.bloom_key)
      logger.info(f"uploaded bloom to s3: s3://{self.bucket}/{self.bloom_key}")
    except Exception as e:
      logger.error(f"failed to upload bloom to s3: {e}")
    self._added_since_save = 0

  def seen(self, v):
    return v in self.sbf

  def register(self, v):
    self.sbf.add(v)
    self._added_since_save += 1
    if self.checkpoint_every and self._added_since_save >= self.checkpoint_every:
      self._persist()

  def close(self):
    if self._added_since_save:
      self._persist()

  def __enter__(self):
    return self

  def __exit__(self, exc_type, exc, tb):
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
  ):
    self.access_level = acess_level
    self.input_csv = input_csv
    self.model_id = model_id
    self.model_version = model_version
    self.bucket = bucket if bucket is not None else DEFAULT_BUCKET_NAME
    self.base_prefix = f"{model_id}/{model_version}/tranches"
    self.index_key = f"{self.base_prefix}/index.json"
    self.meta_key = f"{model_id}/{model_version}/metadata.json"
    self.meta_tranches_key = f"{self.base_prefix}/metadata.json"
    self.local_index = os.path.join(STORE_DIRECTORY, "index.json")
    self.metadata_path = metadata_path
    self._index_dirty = 0
    self.max_rows = int(MAX_ROWS_PER_FILE)
    self.allowed_inputs = set(allowed_inputs) if allowed_inputs else None
    self.tmpdir = tempfile.mkdtemp(prefix="isaura_", dir=STORE_DIRECTORY)
    self.s3 = boto3.client(
      "s3",
      endpoint_url=MINIO_ENDPOINT,
      aws_access_key_id=MINIO_ACCESS_KEY,
      aws_secret_access_key=MINIO_SECRET_KEY,
      config=Config(signature_version="s3v4", s3={"addressing_style": "path"}),
    )
    self._ensure_bucket()
    self.checker = IsauraChecker(
      bucket=self.bucket,
      base_prefix=self.base_prefix,
      store_directory=STORE_DIRECTORY,
      bloom_filename=BLOOM_FILENAME,
      checkpoint_every=CHECKPOINT_EVERY,
      s3=self.s3,
    )
    self.buffers = defaultdict(list)
    self.tranche_state = {}
    self._load_index()
    logger.info(
      f"writer init: bucket={self.bucket} prefix={self.base_prefix} input_csv={self.input_csv} allowed_inputs={len(self.allowed_inputs) if self.allowed_inputs else 'ALL'}"
    )

  def _ensure_bucket(self):
    try:
      self.s3.head_bucket(Bucket=self.bucket)
      logger.info(f"bucket exists: {self.bucket}")
    except Exception:
      self.s3.create_bucket(Bucket=self.bucket)
      logger.info(f"created bucket: {self.bucket}")

  def _load_index(self):
    try:
      self.s3.head_object(Bucket=self.bucket, Key=self.index_key)
      self.s3.download_file(self.bucket, self.index_key, self.local_index)
      with open(self.local_index, "r", encoding="utf-8") as f:
        self.index = json.load(f)
      logger.info(
        f"downloaded index: s3://{self.bucket}/{self.index_key} -> {self.local_index} entries={len(self.index)}"
      )
    except Exception:
      self.index = {}
      logger.info("starting new index")

  def _persist_index(self, force=False):
    if not force and self._index_dirty < CHECKPOINT_EVERY:
      return
    tmp = f"{self.local_index}.tmp"
    with open(tmp, "w", encoding="utf-8") as f:
      json.dump(self.index, f, separators=(",", ":"), ensure_ascii=False)
    os.replace(tmp, self.local_index)
    try:
      self.s3.upload_file(self.local_index, self.bucket, self.index_key)
      logger.info(
        f"uploaded index to s3: s3://{self.bucket}/{self.index_key} entries={len(self.index)}"
      )
    except Exception as e:
      logger.error(f"failed to upload index: {e}")
    self._index_dirty = 0

  def _maybe_upload_metadata(self):
    if not self.metadata_path:
      return
    try:
      self.s3.upload_file(self.metadata_path, self.bucket, self.meta_tranches_key)
      logger.info(
        f"uploaded metadata.json to tranches root: s3://{self.bucket}/{self.meta_tranches_key}"
      )
    except Exception as e:
      logger.error(f"failed to upload metadata.json to tranches: {e}")

  def _tranche_prefix(self, row, col):
    return f"{self.base_prefix}/tranche_{row}_{col}"

  def _list_chunks(self, row, col):
    prefix = self._tranche_prefix(row, col) + "/"
    paginator = self.s3.get_paginator("list_objects_v2")
    page_iter = paginator.paginate(Bucket=self.bucket, Prefix=prefix)
    keys = []
    for page in page_iter:
      for obj in page.get("Contents", []):
        k = obj["Key"]
        if k.endswith(".parquet") and "/chunk_" in k:
          keys.append(k)
    logger.info(f"list chunks: tranche={prefix} count={len(keys)}")
    return sorted(keys)

  def _chunk_index_from_key(self, key):
    base = os.path.basename(key)
    name, _ = os.path.splitext(base)
    try:
      return int(name.split("_")[1])
    except:
      return 0

  def _num_rows_in_parquet(self, local_path):
    try:
      pf = pq.ParquetFile(local_path)
      return pf.metadata.num_rows
    except:
      try:
        return len(pd.read_parquet(local_path))
      except:
        return 0

  def _ensure_tranche_state(self, row, col):
    tkey = (row, col)
    if tkey in self.tranche_state:
      return
    keys = self._list_chunks(row, col)
    if not keys:
      self.tranche_state[tkey] = {
        "next_chunk_idx": 1,
        "open_chunk_key": None,
        "open_chunk_rows": 0,
      }
      logger.info(f"tranche state new: ({row},{col}) -> next=1")
      return
    last = keys[-1]
    local = os.path.join(self.tmpdir, f"inspect_{uuid.uuid4().hex}.parquet")
    n = 0
    try:
      self.s3.download_file(self.bucket, last, local)
      n = self._num_rows_in_parquet(local)
      logger.info(f"inspect last chunk: {last} rows={n}")
    except Exception as e:
      logger.error(f"failed to inspect last chunk: {e}")
    try:
      os.remove(local)
    except:
      pass
    idx = self._chunk_index_from_key(last)
    if n < self.max_rows:
      self.tranche_state[tkey] = {
        "next_chunk_idx": idx,
        "open_chunk_key": last,
        "open_chunk_rows": n,
      }
      logger.info(f"tranche state open: ({row},{col}) -> idx={idx} rows={n}")
    else:
      self.tranche_state[tkey] = {
        "next_chunk_idx": idx + 1,
        "open_chunk_key": None,
        "open_chunk_rows": 0,
      }
      logger.info(f"tranche state rotate: ({row},{col}) -> next={idx + 1}")

  def _write_chunk(self, df, row, col, chunk_idx, mode="new", existing_local=None):
    prefix = self._tranche_prefix(row, col)
    os_key = f"{prefix}/chunk_{chunk_idx}.parquet"
    local = existing_local or os.path.join(
      self.tmpdir, f"chunk_{uuid.uuid4().hex}.parquet"
    )
    if mode == "append" and existing_local:
      old = pd.read_parquet(existing_local)
      df = pd.concat([old, df], ignore_index=True)
    df.to_parquet(local, index=False)
    try:
      self.s3.upload_file(local, self.bucket, os_key)
      logger.info(f"uploaded: rows={len(df)} s3://{self.bucket}/{os_key}")
    except Exception as e:
      logger.error(f"upload failed: {e}")
      raise
    if not existing_local:
      try:
        os.remove(local)
      except:
        pass
    return os_key, len(df)

  def _flush_tranche_buffer(self, row, col):
    tkey = (row, col)
    buf = self.buffers.get(tkey, [])
    if not buf:
      return
    self._ensure_tranche_state(row, col)
    state = self.tranche_state[tkey]
    df = pd.DataFrame(buf)
    remaining = len(df)
    start = 0
    logger.info(f"flush tranche buffer: ({row},{col}) rows={remaining}")
    if state["open_chunk_key"]:
      tmp_local = os.path.join(self.tmpdir, f"open_{uuid.uuid4().hex}.parquet")
      try:
        self.s3.download_file(self.bucket, state["open_chunk_key"], tmp_local)
        space = self.max_rows - state["open_chunk_rows"]
        take = min(space, remaining)
        if take > 0:
          part = df.iloc[start : start + take]
          _, _ = self._write_chunk(
            part,
            row,
            col,
            state["next_chunk_idx"],
            mode="append",
            existing_local=tmp_local,
          )
          state["open_chunk_rows"] += take
          remaining -= take
          start += take
          logger.info(
            f"appended to open chunk: ({row},{col}) idx={state['next_chunk_idx']} added={take} total={state['open_chunk_rows']}"
          )
        if state["open_chunk_rows"] >= self.max_rows:
          state["next_chunk_idx"] += 1
          state["open_chunk_key"] = None
          state["open_chunk_rows"] = 0
          logger.info(f"open chunk closed, next idx={state['next_chunk_idx']}")
      finally:
        try:
          os.remove(tmp_local)
        except:
          pass
    while remaining > 0:
      take = min(self.max_rows, remaining)
      part = df.iloc[start : start + take]
      os_key, _ = self._write_chunk(part, row, col, state["next_chunk_idx"], mode="new")
      if take < self.max_rows:
        state["open_chunk_key"] = os_key
        state["open_chunk_rows"] = take
        logger.info(
          f"created new open chunk: ({row},{col}) idx={state['next_chunk_idx']} rows={take}"
        )
      else:
        state["next_chunk_idx"] += 1
        state["open_chunk_key"] = None
        state["open_chunk_rows"] = 0
        logger.info(
          f"created full chunk: ({row},{col}) idx={state['next_chunk_idx'] - 1} rows={take}"
        )
      remaining -= take
      start += take
    self.buffers[tkey].clear()

  def write(self):
    total = 0
    skipped = 0
    if self.metadata_path:
      self._maybe_upload_metadata()
    with open(self.input_csv, newline="", encoding="utf-8") as f:
      reader = csv.DictReader(f)
      for row in reader:
        smi = (row.get("input") or "").strip()
        if not smi:
          continue
        if self.allowed_inputs is not None and smi not in self.allowed_inputs:
          continue
        if self.checker.seen(smi):
          skipped += 1
          continue
        try:
          r, c, mw, lp = tranche_coordinates(smi)
        except Exception as e:
          logger.warning(f"invalid smiles skipped: {e}")
          continue
        self.buffers[(r, c)].append(dict(row))
        if smi not in getattr(self, "index", {}):
          self.index[smi] = [r, c]
          self._index_dirty += 1
        self.checker.register(smi)
        total += 1
        if len(self.buffers[(r, c)]) >= self.max_rows:
          self._flush_tranche_buffer(r, c)
        if self._index_dirty >= CHECKPOINT_EVERY:
          self._persist_index()
    for r, c in list(self.buffers.keys()):
      if self.buffers[(r, c)]:
        self._flush_tranche_buffer(r, c)
    self._persist_index(force=True)
    logger.info(f"write done: new={total} skipped={skipped}")

  def close(self):
    self.checker.close()
    try:
      os.rmdir(self.tmpdir)
    except:
      pass

  def __enter__(self):
    return self

  def __exit__(self, exc_type, exc, tb):
    self.close()


class IsauraReader:
  def __init__(
    self,
    model_id,
    model_version,
    input_csv,
    bucket=None,
    minio_endpoint=None,
    minio_access_key=None,
    minio_secret_key=None,
  ):
    self.model_id = model_id
    self.model_version = model_version
    self.input_csv = input_csv
    self.bucket = bucket or DEFAULT_BUCKET_NAME
    endpoint = minio_endpoint or MINIO_ENDPOINT
    access = minio_access_key or MINIO_ACCESS_KEY
    secret = minio_secret_key or MINIO_SECRET_KEY
    self.base_prefix = f"{model_id}/{model_version}/tranches"
    self.index_key = f"{self.base_prefix}/index.json"
    self.s3 = boto3.client(
      "s3",
      endpoint_url=endpoint,
      aws_access_key_id=access,
      aws_secret_access_key=secret,
      config=Config(signature_version="s3v4", s3={"addressing_style": "path"}),
    )
    self.tmpdir = tempfile.mkdtemp(prefix="isaura_reader_", dir=STORE_DIRECTORY)
    self.con = duckdb.connect(database=":memory:")
    self.con.execute("INSTALL httpfs; LOAD httpfs;")
    ep = endpoint.replace("http://", "").replace("https://", "")
    use_ssl = not endpoint.startswith("http://")
    self.con.execute("SET s3_access_key_id=?", [access])
    self.con.execute("SET s3_secret_access_key=?", [secret])
    self.con.execute("SET s3_endpoint=?", [ep])
    self.con.execute("SET s3_region='us-east-1'")
    self.con.execute("SET s3_use_ssl=?", [use_ssl])
    self.con.execute("SET s3_url_style='path'")
    logger.info(
      f"IsauraReader init bucket={self.bucket} prefix={self.base_prefix} input_csv={self.input_csv}"
    )

  def _tranche_prefix(self, row, col):
    return f"{self.base_prefix}/tranche_{row}_{col}"

  def _list_chunks_with_sizes(self, row, col):
    prefix = self._tranche_prefix(row, col) + "/"
    paginator = self.s3.get_paginator("list_objects_v2")
    keys = []
    for page in paginator.paginate(Bucket=self.bucket, Prefix=prefix):
      for obj in page.get("Contents", []):
        k = obj["Key"]
        if k.endswith(".parquet") and "/chunk_" in k:
          keys.append((k, int(obj.get("Size", 0))))
    logger.info(
      f"found {len(keys)} chunks for tranche=({row},{col}) total_size={sum(sz for _, sz in keys)} bytes"
    )
    return sorted(keys, key=lambda x: x[0])

  def _load_index(self):
    local = os.path.join(self.tmpdir, f"{uuid.uuid4().hex}.json")
    self.s3.download_file(self.bucket, self.index_key, local)
    try:
      with open(local, "r", encoding="utf-8") as f:
        return json.load(f)
    finally:
      try:
        os.remove(local)
      except:
        pass

  def read(self, input_col="input", output_file=None, transient_progress=True):
    t0 = time.time()
    wanted = []
    with open(self.input_csv, newline="", encoding="utf-8") as f:
      for row in csv.DictReader(f):
        v = (row.get(input_col) or "").strip()
        if v:
          wanted.append(v)
    if not wanted:
      logger.info(f"no inputs found in csv={self.input_csv}")
      return pd.DataFrame()
    logger.info(f"parsed {len(wanted)} inputs from csv={self.input_csv}")
    index = self._load_index()
    logger.info(f"loaded index entries={len(index)}")
    missing = [s for s in wanted if s not in index]
    if missing:
      raise RuntimeError(
        f"inputs not indexed: {missing[:5]}{'...' if len(missing) > 5 else ''} total_missing={len(missing)}"
      )
    groups = {}
    for smi in wanted:
      rc = index[smi]
      groups.setdefault(tuple(rc), set()).add(smi)
    tranche_sets = []
    total_bytes = 0
    for (r, c), smi_set in groups.items():
      keys = self._list_chunks_with_sizes(r, c)
      tranche_sets.append((r, c, smi_set, keys))
      total_bytes += sum(sz for _, sz in keys)
    logger.info(
      f"total expected download size={total_bytes} bytes across {sum(len(keys) for _, _, _, keys in tranche_sets)} chunks"
    )
    results = []
    order_map = {smi: i for i, smi in enumerate(wanted)}
    desc = f"Fetching {self.model_id}/{self.model_version} ({len(wanted)} inputs)"
    with make_fetching_progress(transient=transient_progress) as progress:
      task_id = progress.add_task("download", total=total_bytes or 0, desc=desc)
      progressed = 0
      for r, c, smi_set, keys in tranche_sets:
        if not keys:
          logger.info(f"no chunks for tranche=({r},{c})")
          continue
        logger.info(
          f"querying tranche=({r},{c}) keys={len(keys)} inputs={len(smi_set)}"
        )
        files_glob = f"s3://{self.bucket}/{self._tranche_prefix(r, c)}/chunk_*.parquet"
        df_wanted = pd.DataFrame({input_col: list(smi_set)})
        view_name = f"wanted_{r}_{c}_{uuid.uuid4().hex[:6]}"
        self.con.register(view_name, df_wanted)
        q = f"SELECT t.* FROM read_parquet('{files_glob}') t INNER JOIN {view_name} w ON t.{input_col}=w.{input_col}"
        part = self.con.execute(q).fetchdf()
        self.con.unregister(view_name)
        if not part.empty:
          part["__order"] = part[input_col].map(order_map)
          results.append(part)
          logger.info(f"matched {len(part)} rows in tranche=({r},{c})")
        progressed += sum(sz for _, sz in keys)
        if total_bytes:
          progress.update(task_id, completed=min(progressed, total_bytes))
    out = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
    if not out.empty and "__order" in out.columns:
      out = out.sort_values("__order").drop(columns="__order").reset_index(drop=True)
    logger.info(f"collected total matched rows={len(out)}")
    if output_file:
      out.to_csv(output_file, index=False)
      logger.info(f"wrote result csv path={output_file} rows={len(out)}")
    elapsed = time.time() - t0
    rate = (len(out) / elapsed) if elapsed > 0 and len(out) else 0.0
    logger.info(
      f"read done model={self.model_id} version={self.model_version} bucket={self.bucket} inputs={len(wanted)} matched={len(out)} elapsed={elapsed:.2f}s rate={rate:.1f} rows/s"
    )
    return out
