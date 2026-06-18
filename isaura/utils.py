import atexit
import math
import os
import psutil
import shutil
import tempfile
import threading
import time
import weakref

from isaura.const import (
  ACCESS_FILE, INDEX_FILE, STORE_DIRECTORY,
  DEFAULT_BUCKET_NAME, DEFAULT_PRIVATE_BUCKET_NAME,
)
from isaura.logging import logger

proc = psutil.Process(os.getpid())
_TEMP_DIRS = set()
_TEMP_DIRS_LOCK = threading.Lock()


def rss_mb():
  """Return the current process RSS memory usage in megabytes."""
  return proc.memory_info().rss / (1024 * 1024)


def log(msg):
  """Log a message at INFO level with a timestamp and current RSS."""
  logger.info(f"[{time.strftime('%H:%M:%S')}] {msg} | RSS={rss_mb():.1f} MB")


def avail_mem():
  """Return the number of bytes of available system memory."""
  return int(psutil.virtual_memory().available)


def mem_gb_lim(ratio=0.8, floor_gb=1, total=False):
  """Return a memory limit in GB as a fraction of RAM, with a minimum floor.

  Args:
      ratio: Fraction of memory to use (default 0.8).
      floor_gb: Minimum value to return in GB (default 1).
      total: Base the fraction on TOTAL RAM instead of currently-available RAM.
          macOS underreports "available" (aggressive compression/caching), which
          can collapse the limit to the floor; use total for a stable ceiling
          when the consumer (e.g. DuckDB) spills to disk on its own.
  """
  base = psutil.virtual_memory().total if total else avail_mem()
  return max(floor_gb, int(base * ratio / 1024**3))


def cpu_cnt(ratio=0.6):
  """Return a CPU thread count as a fraction of available cores, minimum 1.

  Args:
      ratio: Fraction of CPU cores to use (default 0.6).
  """
  return max(1, int(math.floor((os.cpu_count() or 1) * ratio)))


def cleanup_temp_dir(path):
  """Delete a temporary directory and remove it from the global tracking set."""
  if not path:
    return
  try:
    shutil.rmtree(path, ignore_errors=True)
  finally:
    with _TEMP_DIRS_LOCK:
      _TEMP_DIRS.discard(path)


def bind_temp_dirs(owner, *paths):
  """Register temp directories to be deleted when owner is garbage collected.

  Uses weakref finalizers so cleanup happens automatically even if close()
  is never called explicitly.

  Args:
      owner: Object whose lifetime governs the temp directories.
      *paths: Temp directory paths to register.

  Returns:
      List of weakref.finalize objects.
  """
  finalizers = [weakref.finalize(owner, cleanup_temp_dir, path) for path in paths if path]
  owner._temp_dir_finalizers = finalizers
  return finalizers


def release_temp_dirs(owner):
  """Immediately trigger and clear all registered temp dir finalizers for owner."""
  finalizers = list(getattr(owner, "_temp_dir_finalizers", []) or [])
  owner._temp_dir_finalizers = []
  for finalizer in finalizers:
    try:
      finalizer()
    except Exception:
      pass


@atexit.register
def _cleanup_registered_temp_dirs():
  with _TEMP_DIRS_LOCK:
    paths = list(_TEMP_DIRS)
  for path in paths:
    cleanup_temp_dir(path)


def make_temp(pref):
  """Create a uniquely named temporary directory inside STORE_DIRECTORY and track it for cleanup.

  Args:
      pref: Prefix for the directory name (e.g. "isaura_reader_").

  Returns:
      Absolute path to the created directory.
  """
  os.makedirs(STORE_DIRECTORY, exist_ok=True)
  path = tempfile.mkdtemp(prefix=pref, dir=STORE_DIRECTORY)
  with _TEMP_DIRS_LOCK:
    _TEMP_DIRS.add(path)
  return path


def get_base(mdi, ver):
  """Return the tranches key prefix for a model (e.g. "eos1234/v1/tranches")."""
  return f"{get_pref(mdi, ver)}/tranches"


def get_pref(mdi, ver):
  """Return the top-level key prefix for a model (e.g. "eos1234/v1")."""
  return f"{mdi}/{ver}"


def get_coll(mdi, ver):
  """Return the NNS collection name for a model (e.g. "eos1234_v1")."""
  return f"{mdi}_{ver}"


def hive_prefix(base):
  """Return the Hive-style data prefix for a model's Parquet chunks (base/data)."""
  return f"{base}/data"


def get_files_glob(bucket, base):
  """Return an S3 glob URI matching all Parquet chunk files for a model."""
  return f"s3://{bucket}/{base}/*/chunk_*.parquet"


def get_keys(file, base):
  """Return the full MinIO key for a metadata file under a model's base prefix."""
  return f"{base}/{file}"


def get_idx_key(base):
  """Return the MinIO key for a model's JSON index file."""
  return get_keys(INDEX_FILE, base)


def get_acc_key(base):
  """Return the MinIO key for a model's access metadata file."""
  return get_keys(ACCESS_FILE, base)


def get_desc(pref, wanted):
  """Return a human-readable description string for a batch fetch operation."""
  return f"Fetching hive partitions {pref} ({len(wanted)} inputs)"


def split_csv(df):
  """Split a DataFrame by bucket column and write each subset to a separate CSV file.

  Returns:
      List of paths to the written CSV files.
  """
  paths = []
  output_dir = make_temp("isaura_push_")
  for bucket in [DEFAULT_BUCKET_NAME, DEFAULT_PRIVATE_BUCKET_NAME]:
    if bucket in df["bucket"].unique():
      path = os.path.join(output_dir, f"{bucket.replace('-', '_')}.csv")
      df[df["bucket"] == bucket].to_csv(path, index=False)
      paths.append(str(path))
  return paths
