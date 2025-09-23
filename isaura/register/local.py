import boto3, csv, os, pickle, tempfile, uuid

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
        "downloaded bloom from s3: s3://%s/%s -> %s",
        self.bucket,
        self.bloom_key,
        self.local_bloom,
      )
    except Exception as e:
      logger.info("no remote bloom found or failed to download: %s", e)
    try:
      with open(self.local_bloom, "rb") as f:
        self.sbf = pickle.load(f)
      logger.info("loaded bloom locally: %s", self.local_bloom)
    except FileNotFoundError:
      self.sbf = ScalableBloomFilter(
        mode=ScalableBloomFilter.SMALL_SET_GROWTH,
        initial_capacity=initial_capacity,
        error_rate=error_rate,
      )
      logger.info(
        "created new bloom: initial_capacity=%d error_rate=%s",
        initial_capacity,
        error_rate,
      )

  def _persist(self):
    tmp = f"{self.local_bloom}.tmp"
    with open(tmp, "wb") as f:
      pickle.dump(self.sbf, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, self.local_bloom)
    logger.info("saved bloom locally: %s", self.local_bloom)
    try:
      self.s3.upload_file(self.local_bloom, self.bucket, self.bloom_key)
      logger.info("uploaded bloom to s3: s3://%s/%s", self.bucket, self.bloom_key)
    except Exception as e:
      logger.error("failed to upload bloom to s3: %s", e)
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
    acess_level,
    bucket=None,
  ):
    self.access_level = acess_level
    self.input_csv = input_csv
    self.model_id = model_id
    self.model_version = model_version
    self.bucket = bucket or DEFAULT_BUCKET_NAME
    self.base_prefix = f"{model_id}/{model_version}/tranches"
    self.max_rows = int(MAX_ROWS_PER_FILE)
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
      bucket=bucket,
      base_prefix=self.base_prefix,
      store_directory=STORE_DIRECTORY,
      bloom_filename=BLOOM_FILENAME,
      checkpoint_every=CHECKPOINT_EVERY,
      s3=self.s3,
    )
    self.buffers = defaultdict(list)
    self.tranche_state = {}
    logger.info(
      "writer init: bucket=%s prefix=%s input_csv=%s",
      self.bucket,
      self.base_prefix,
      self.input_csv,
    )

  def _ensure_bucket(self):
    try:
      self.s3.head_bucket(Bucket=self.bucket)
      logger.info("bucket exists: %s", self.bucket)
    except Exception:
      self.s3.create_bucket(Bucket=self.bucket)
      logger.info("created bucket: %s", self.bucket)

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
    logger.info("list chunks: tranche=%s count=%d", prefix, len(keys))
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
      logger.info("tranche state new: (%d,%d) -> next=1", row, col)
      return
    last = keys[-1]
    local = os.path.join(self.tmpdir, f"inspect_{uuid.uuid4().hex}.parquet")
    n = 0
    try:
      self.s3.download_file(self.bucket, last, local)
      n = self._num_rows_in_parquet(local)
      logger.info("inspect last chunk: %s rows=%d", last, n)
    except Exception as e:
      logger.error("failed to inspect last chunk: %s", e)
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
      logger.info("tranche state open: (%d,%d) -> idx=%d rows=%d", row, col, idx, n)
    else:
      self.tranche_state[tkey] = {
        "next_chunk_idx": idx + 1,
        "open_chunk_key": None,
        "open_chunk_rows": 0,
      }
      logger.info("tranche state rotate: (%d,%d) -> next=%d", row, col, idx + 1)

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
      logger.info("uploaded: rows=%d s3://%s/%s", len(df), self.bucket, os_key)
    except Exception as e:
      logger.error("upload failed: %s", e)
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
    logger.info("flush tranche buffer: (%d,%d) rows=%d", row, col, remaining)
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
            "appended to open chunk: (%d,%d) idx=%d added=%d total=%d",
            row,
            col,
            state["next_chunk_idx"],
            take,
            state["open_chunk_rows"],
          )
        if state["open_chunk_rows"] >= self.max_rows:
          state["next_chunk_idx"] += 1
          state["open_chunk_key"] = None
          state["open_chunk_rows"] = 0
          logger.info("open chunk closed, next idx=%d", state["next_chunk_idx"])
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
          "created new open chunk: (%d,%d) idx=%d rows=%d",
          row,
          col,
          state["next_chunk_idx"],
          take,
        )
      else:
        state["next_chunk_idx"] += 1
        state["open_chunk_key"] = None
        state["open_chunk_rows"] = 0
        logger.info(
          "created full chunk: (%d,%d) idx=%d rows=%d",
          row,
          col,
          state["next_chunk_idx"] - 1,
          take,
        )
      remaining -= take
      start += take
    self.buffers[tkey].clear()

  def write(self):
    total = 0
    skipped = 0
    with open(self.input_csv, newline="", encoding="utf-8") as f:
      reader = csv.DictReader(f)
      for row in reader:
        smi = (row.get("input") or "").strip()
        if not smi:
          continue
        if self.checker.seen(smi):
          skipped += 1
          continue
        try:
          r, c, mw, lp = tranche_coordinates(smi)
        except Exception as e:
          logger.warning("invalid smiles skipped: %s", e)
          continue
        self.buffers[(r, c)].append(dict(row))
        self.checker.register(smi)
        total += 1
        if len(self.buffers[(r, c)]) >= self.max_rows:
          self._flush_tranche_buffer(r, c)
    for r, c in list(self.buffers.keys()):
      if self.buffers[(r, c)]:
        self._flush_tranche_buffer(r, c)
    logger.info("write done: new=%d skipped=%d", total, skipped)

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
