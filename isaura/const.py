import os
from pathlib import Path

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
MIN_NNS_RESULT_SIZE = 1000
MAX_ROWS = 2000000

BUILD_STATUS_TIMEOUT = 5
BUILD_START_TIMEOUT = 30
BUILD_POLL_INTERVAL = 0.5
BUILD_MAX_WAIT = 1.0

MW_BINS = [200, 500]
LOGP_BINS = [-1, 5]

GITHUB_ORG = "ersilia-os"
GITHUB_CONTENT_URL = f"https://raw.githubusercontent.com/{GITHUB_ORG}"
METADATA_JSON = "metadata.json"
METADATA_YML = "metadata.yml"

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

MAX_ROWS_PER_FILE = int(os.getenv("MAX_ROWS_PER_FILE", "100000"))
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

WIDE_READ_SLICE = int(os.getenv("WIDE_READ_SLICE", "2000"))

CHECKPOINT_EVERY = int(os.getenv("CHECKPOINT_EVERY", "50000"))

PARQUET_ROW_GROUP_SIZE = int(os.getenv("PARQUET_ROW_GROUP_SIZE", "500000"))
WIDE_PARQUET_ROW_GROUP_SIZE = int(os.getenv("WIDE_PARQUET_ROW_GROUP_SIZE", "50000"))
PARQUET_COMPRESSION = os.getenv("PARQUET_COMPRESSION", "zstd")
PARQUET_COMPRESSION_LEVEL = int(os.getenv("PARQUET_COMPRESSION_LEVEL", "1"))
PARQUET_DATA_PAGE_SIZE = int(os.getenv("PARQUET_DATA_PAGE_SIZE", "1048576"))
WIDE_PARQUET_DATA_PAGE_SIZE = int(os.getenv("WIDE_PARQUET_DATA_PAGE_SIZE", "2097152"))
WIDE_PARQUET_USE_DICTIONARY = os.getenv("WIDE_PARQUET_USE_DICTIONARY", "false").lower() in ("1", "true", "yes")
WIDE_PARQUET_WRITE_STATISTICS = os.getenv("WIDE_PARQUET_WRITE_STATISTICS", "false").lower() in ("1", "true", "yes")

MINIO_MULTIPART_THRESHOLD = int(os.getenv("MINIO_MULTIPART_THRESHOLD", str(8 * 1024 * 1024)))
MINIO_MULTIPART_CHUNKSIZE = int(os.getenv("MINIO_MULTIPART_CHUNKSIZE", str(8 * 1024 * 1024)))
MINIO_MAX_CONCURRENCY = int(os.getenv("MINIO_MAX_CONCURRENCY", "10"))

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
  {"key": "chunks", "name": "chunks", "justify": "right"},
]
inspect_table_cloud = [
  {"key": "model", "name": "model/version", "justify": "left", "style": "bold"},
  {"key": "entries", "name": "rows", "justify": "right"},
  {"key": "size", "name": "size", "justify": "right"},
  {"key": "chunks", "name": "chunks", "justify": "right"},
]
