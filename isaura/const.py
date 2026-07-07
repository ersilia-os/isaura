import os
from pathlib import Path

import psutil

try:
  from dotenv import load_dotenv
  _user_env = Path.home() / ".isaura" / ".env"
  if not _user_env.exists():
    _user_env.parent.mkdir(parents=True, exist_ok=True)
    _user_env.write_text(
      "MINIO_ENDPOINT=http://127.0.0.1:9000\n"
      "NNS_ENDPOINT=http://127.0.0.1:8080\n"
      "DEFAULT_BUCKET_NAME=isaura-public\n"
      "DEFAULT_PRIVATE_BUCKET_NAME=isaura-private\n"
      "MINIO_LOCAL_AK=minioadmin123\n"
      "MINIO_LOCAL_SK=minioadmin1234\n"
    )
  load_dotenv(dotenv_path=_user_env, override=True)
except Exception:
  pass

ACCESS_FILE = "access.json"
CATALOG_FILE = "catalog.json"
INDEX_FILE = "index.json"
# Workstream 4: real molecule->(chunk, row_group) index for wide models. Presence of
# this file with meta.format == INDEX_FORMAT marks a "real" index; absence => legacy full scan.
INDEX_SQLITE_FILE = "index.sqlite"
INDEX_FORMAT = "loc-v1"
INDEX_GRANULARITY = "chunk_rowgroup"
MIN_NNS_RESULT_SIZE = 1000
MAX_ROWS = 2000000

BUILD_STATUS_TIMEOUT = 5
BUILD_START_TIMEOUT = 30
BUILD_POLL_INTERVAL = 0.5
BUILD_MAX_WAIT = 1.0

GITHUB_ORG = "ersilia-os"
GITHUB_CONTENT_URL = f"https://raw.githubusercontent.com/{GITHUB_ORG}"
METADATA_JSON = "metadata.json"
METADATA_YML = "metadata.yml"
# Authoritative per-column output types: name,type,direction,description
RUN_COLUMNS_FILE = os.getenv("RUN_COLUMNS_FILE", "model/framework/columns/run_columns.csv")

COLLECTION = os.getenv("COLLECTION", "eos3b5e")
MINIO_ENDPOINT = os.getenv("MINIO_ENDPOINT", "http://127.0.0.1:9000")
TIMEOUT = os.getenv("TIMEOUT", 3600)
MINIO_ENDPOINT_CLOUD = os.getenv("MINIO_ENDPOINT_CLOUD") or "http://83.48.73.209:8080"
NNS_ENDPOINT_BASE = os.getenv("NNS_ENDPOINT") or "http://127.0.0.1:8080"

MINIO_LOCAL_AK = os.getenv("MINIO_LOCAL_AK", "minioadmin123")
MINIO_LOCAL_SK = os.getenv("MINIO_LOCAL_SK", "minioadmin1234")
MINIO_CLOUD_AK = os.getenv("MINIO_CLOUD_AK", None)
MINIO_CLOUD_SK = os.getenv("MINIO_CLOUD_SK", None)
MINIO_PRIV_CLOUD_AK = os.getenv("MINIO_PRIV_CLOUD_AK", None)
MINIO_PRIV_CLOUD_SK = os.getenv("MINIO_PRIV_CLOUD_SK", None)

isaura_temp = os.path.join(Path.home(), "isaura", "isaura-temp")
if not os.path.exists(isaura_temp):
  os.makedirs(isaura_temp)
STORE_DIRECTORY = os.getenv("STORE_DIRECTORY", isaura_temp)

# Wide-model chunk size. 10k rows keeps per-chunk RAM bounded (e.g. 10k x 3200
# float64 cols ~= 250MB vs ~2.5GB at 100k), which the read/pull path relies on.
MAX_ROWS_PER_FILE = int(os.getenv("MAX_ROWS_PER_FILE", "10000"))
IMMUTABLE_CHUNK_COLS_THRESHOLD = int(os.getenv("IMMUTABLE_CHUNK_COLS_THRESHOLD", "128"))
WIDE_OUTPUT_DIM_THRESHOLD = int(os.getenv("WIDE_OUTPUT_DIM_THRESHOLD", "100"))

DEFAULT_WRITE_BATCH_ROWS = int(os.getenv("DEFAULT_WRITE_BATCH_ROWS", "100000"))
WIDE_WRITE_BATCH_ROWS = int(os.getenv("WIDE_WRITE_BATCH_ROWS", "50000"))
WRITE_INPUT_CHUNK_ROWS = int(os.getenv("WRITE_INPUT_CHUNK_ROWS", "50000"))

STREAM_DENSE_FILE_RATIO = float(os.getenv("STREAM_DENSE_FILE_RATIO", "0.25"))
STREAM_DENSE_BATCH_ROWS = int(os.getenv("STREAM_DENSE_BATCH_ROWS", "16384"))
STREAM_PARQUET_THRESHOLD = int(os.getenv("STREAM_PARQUET_THRESHOLD", "200000"))
STREAM_PREFETCH_FILES = int(os.getenv("STREAM_PREFETCH_FILES", "6"))
STREAM_DOWNLOAD_WORKERS = int(os.getenv("STREAM_DOWNLOAD_WORKERS", "4"))

# --- Wide-model read/pull speed knobs (auto-scaled to the machine; env wins) ---
# These tune ONLY the wide-model streaming path (output dim >= WIDE_OUTPUT_DIM_
# THRESHOLD). Narrow models keep the STREAM_* defaults above and the DuckDB path
# untouched. Defaults are computed from CPU/RAM so a server speeds up while a
# low-resource laptop stays conservative; any explicit env var overrides them.
# (const.py is imported by isaura.utils, so we can't import cpu_cnt/avail_mem
# from there without a circular import — the formulas below mirror those.)
def _cpu_frac(ratio):
  return max(1, int((os.cpu_count() or 1) * ratio))


def _avail_gb():
  try:
    return psutil.virtual_memory().available / 1024 ** 3
  except Exception:
    return 4.0


def _wide_dense_batch_default():
  g = _avail_gb()
  if g >= 8:
    return 32768
  if g >= 4:
    return 16384
  return 8192


# Downloads are IO-bound, so cap at 8 to avoid diminishing returns / churn.
WIDE_STREAM_DOWNLOAD_WORKERS = int(os.getenv("WIDE_STREAM_DOWNLOAD_WORKERS", str(min(8, max(2, _cpu_frac(0.75))))))
# Keep the prefetch window ahead of the workers without ballooning tmpdir use.
WIDE_STREAM_PREFETCH_FILES = int(os.getenv("WIDE_STREAM_PREFETCH_FILES", str(min(16, WIDE_STREAM_DOWNLOAD_WORKERS * 2))))
WIDE_STREAM_DENSE_BATCH_ROWS = int(os.getenv("WIDE_STREAM_DENSE_BATCH_ROWS", str(_wide_dense_batch_default())))
# 0/unset => compute adaptively at read time from available RAM and the model's
# output width (see IsauraReader); a positive value pins the wide read batch.
WIDE_READ_BATCH_ROWS = int(os.getenv("WIDE_READ_BATCH_ROWS", "0"))
# Per-batch RAM target (fraction of available memory) used for that adaptive
# calc, plus floor/ceiling clamps on the resulting row count.
WIDE_READ_BATCH_MEM_FRACTION = float(os.getenv("WIDE_READ_BATCH_MEM_FRACTION", "0.10"))
WIDE_READ_BATCH_ROWS_MIN = int(os.getenv("WIDE_READ_BATCH_ROWS_MIN", "5000"))
WIDE_READ_BATCH_ROWS_MAX = int(os.getenv("WIDE_READ_BATCH_ROWS_MAX", "40000"))
# Hard RAM ceiling per decoded batch, derived from each file's ACTUAL bytes/row
# (parquet metadata) rather than an assumed cell size. This is the safety net
# that bounds reads regardless of dtype: a model whose numeric outputs were
# stored as text (e.g. eos4u6p, 3200 string cols) has a huge bytes/row, so the
# stream auto-shrinks its batches; correctly-typed numeric models are unaffected.
WIDE_READ_TARGET_BYTES = int(os.getenv("WIDE_READ_TARGET_BYTES", str(128 * 1024 * 1024)))

WIDE_READ_SLICE = int(os.getenv("WIDE_READ_SLICE", "2000"))

# Memory-bounded ordered wide-read reorder (external bucket-sort). Rows are
# scattered to on-disk "bucket" parquet files keyed by input position, then
# gathered one bucket at a time in input order so peak RAM is ~one bucket.
# Pass B does read_table(bucket).to_pandas(), so the bucket row count is the
# gather-phase RAM unit. Keep it small: for a string-typed wide model (e.g.
# eos4u6p, 3200 cols) a 50k-row bucket is ~8GB in pandas, a 5k-row one ~1GB.
# (Effective span is max(this, ceil(n_inputs / WIDE_REORDER_MAX_OPEN_BUCKETS)).)
WIDE_REORDER_BUCKET_ROWS = int(os.getenv("WIDE_REORDER_BUCKET_ROWS", "5000"))
WIDE_REORDER_MAX_OPEN_BUCKETS = int(os.getenv("WIDE_REORDER_MAX_OPEN_BUCKETS", "256"))
WIDE_REORDER_MIN_FREE_MB = int(os.getenv("WIDE_REORDER_MIN_FREE_MB", "2048"))

CHECKPOINT_EVERY = int(os.getenv("CHECKPOINT_EVERY", "50000"))

PARQUET_ROW_GROUP_SIZE = int(os.getenv("PARQUET_ROW_GROUP_SIZE", "500000"))
# Fallback wide row-group size; used only if byte calibration can't run (empty batch).
WIDE_PARQUET_ROW_GROUP_SIZE = int(os.getenv("WIDE_PARQUET_ROW_GROUP_SIZE", "10000"))

# Wide chunk file / row-group byte targets (Workstream 3). Row counts are derived
# at write time from the first typed table's measured bytes/row (see ChunkState),
# so these replace the fixed MAX_ROWS_PER_FILE / WIDE_PARQUET_ROW_GROUP_SIZE for
# wide models. 256 MB files: the eos4u6p (3200-col) benchmark showed ~16 files
# saturate the download pool/prefetch window, keeping the heavy single-threaded
# decode fed — ~17% faster than 1 GB (4 files) on the widest models, and finer
# for W4 pruning. 64 MB row groups = lowest read RAM (safest at high width).
WIDE_TARGET_FILE_BYTES = int(os.getenv("WIDE_TARGET_FILE_BYTES", str(256 * 1024 * 1024)))
WIDE_TARGET_ROWGROUP_BYTES = int(os.getenv("WIDE_TARGET_ROWGROUP_BYTES", str(64 * 1024 * 1024)))
PARQUET_COMPRESSION = os.getenv("PARQUET_COMPRESSION", "zstd")
PARQUET_COMPRESSION_LEVEL = int(os.getenv("PARQUET_COMPRESSION_LEVEL", "1"))
PARQUET_DATA_PAGE_SIZE = int(os.getenv("PARQUET_DATA_PAGE_SIZE", "1048576"))
WIDE_PARQUET_DATA_PAGE_SIZE = int(os.getenv("WIDE_PARQUET_DATA_PAGE_SIZE", "2097152"))
WIDE_PARQUET_USE_DICTIONARY = os.getenv("WIDE_PARQUET_USE_DICTIONARY", "false").lower() in ("1", "true", "yes")
WIDE_PARQUET_WRITE_STATISTICS = os.getenv("WIDE_PARQUET_WRITE_STATISTICS", "false").lower() in ("1", "true", "yes")

MINIO_MULTIPART_THRESHOLD = int(os.getenv("MINIO_MULTIPART_THRESHOLD", str(8 * 1024 * 1024)))
MINIO_MULTIPART_CHUNKSIZE = int(os.getenv("MINIO_MULTIPART_CHUNKSIZE", str(8 * 1024 * 1024)))
MINIO_MAX_CONCURRENCY = int(os.getenv("MINIO_MAX_CONCURRENCY", "10"))
# urllib3 connection-pool size for the boto3 S3 client. Default is botocore's 10,
# which throttles concurrent transfers (each wide download worker can open up to
# MINIO_MAX_CONCURRENCY part connections). Size it off the busiest path — the
# wide download workers — plus headroom for the prefetch boundary.
MINIO_MAX_POOL_CONNECTIONS = int(os.getenv(
  "MINIO_MAX_POOL_CONNECTIONS",
  str(WIDE_STREAM_DOWNLOAD_WORKERS * MINIO_MAX_CONCURRENCY + 8),
))

BLOOM_FILENAME = os.getenv("BLOOM_FILENAME", "bloom.pkl")
INPUT_C = ["input", "smiles"]
DEFAULT_BUCKET_NAME = os.getenv("DEFAULT_BUCKET_NAME", "isaura-public")
DEFAULT_PRIVATE_BUCKET_NAME = os.getenv("DEFAULT_PRIVATE_BUCKET_NAME", "isaura-private")
BATCH = int(os.getenv("BATCH", 10000))
FLUSH_EVERY = os.getenv("FLUSH_EVERY", 10000)

MKEYS = [
  "Status", "Deployment", "Source", "Source Type", "Task", "Subtask",
  "Output", "Output Dimension", "Tag", "Biomedical Area",
  "Target Organism", "Publication Type", "Publication Year",
]
_INT_KEYS = {"Output Dimension", "Publication Year"}
_KEYMAP = {
  "Source Type": "SourceType",
  "Output Dimension": "OutputDimension",
  "Biomedical Area": "BiomedicalArea",
  "Target Organism": "TargetOrganism",
  "Publication Type": "PublicationType",
  "Publication Year": "PublicationYear",
}

inspect_table = [
  {"key": "model", "name": "model/version", "justify": "left", "style": "bold"},
  {"key": "entries", "name": "entry count", "justify": "right"},
  {"key": "size", "name": "size", "justify": "right"},
  {"key": "chunks", "name": "chunks", "justify": "right"},
]
inspect_table_cloud = [
  {"key": "model", "name": "model/version", "justify": "left", "style": "bold"},
  {"key": "entries", "name": "rows", "justify": "right"},
  {"key": "size", "name": "size", "justify": "right"},
  {"key": "chunks", "name": "chunks", "justify": "right"},
]
