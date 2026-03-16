import math
import os
import psutil
import tempfile
import time

from isaura.const import (
  ACCESS_FILE, INDEX_FILE, STORE_DIRECTORY,
  DEFAULT_BUCKET_NAME, DEFAULT_PRIVATE_BUCKET_NAME,
)
from isaura.logging import logger

proc = psutil.Process(os.getpid())


def rss_mb():
  return proc.memory_info().rss / (1024 * 1024)


def log(msg):
  logger.info(f"[{time.strftime('%H:%M:%S')}] {msg} | RSS={rss_mb():.1f} MB")


def avail_mem():
  return int(psutil.virtual_memory().available)


def mem_gb_lim(ratio=0.8, floor_gb=1):
  return max(floor_gb, int(avail_mem() * ratio / 1024**3))


def cpu_cnt(ratio=0.6):
  return max(1, int(math.floor((os.cpu_count() or 1) * ratio)))


def make_temp(pref):
  return tempfile.mkdtemp(prefix=pref, dir=STORE_DIRECTORY)


def get_base(mdi, ver):
  return f"{get_pref(mdi, ver)}/tranches"


def get_pref(mdi, ver):
  return f"{mdi}/{ver}"


def get_coll(mdi, ver):
  return f"{mdi}_{ver}"


def hive_prefix(base):
  return f"{base}/data"


def get_files_glob(bucket, base):
  return f"s3://{bucket}/{base}/*/chunk_*.parquet"


def get_keys(file, base):
  return f"{base}/{file}"


def get_idx_key(base):
  return get_keys(INDEX_FILE, base)


def get_acc_key(base):
  return get_keys(ACCESS_FILE, base)


def get_desc(pref, wanted):
  return f"Fetching hive partitions {pref} ({len(wanted)} inputs)"


def split_csv(df):
  paths = []
  output_dir = make_temp("isaura_push_")
  for bucket in [DEFAULT_BUCKET_NAME, DEFAULT_PRIVATE_BUCKET_NAME]:
    if bucket in df["bucket"].unique():
      path = os.path.join(output_dir, f"{bucket.replace('-', '_')}.csv")
      df[df["bucket"] == bucket].to_csv(path, index=False)
      paths.append(str(path))
  return paths
