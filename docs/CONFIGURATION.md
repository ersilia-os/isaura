# Configuration

Isaura is configured primarily via environment variables. You can set them either:

- 🧾 in a `.env` file (recommended for local development), or
- 💻 by exporting variables in your shell (`export ...`) for CI / production.

> ℹ️ If `python-dotenv` is installed, Isaura can load `.env` automatically. Exported shell variables typically take precedence (depending on `load_dotenv(override=...)`).

---

## 🚀 Minimal `.env` (local)

Create a `.env` file in the repository root:

```bash
# Local MinIO (S3)
MINIO_ENDPOINT=http://127.0.0.1:9000

# NN / ANN service
NNS_ENDPOINT=http://127.0.0.1:8080

# Default buckets
DEFAULT_BUCKET_NAME=isaura-public
DEFAULT_PRIVATE_BUCKET_NAME=isaura-private

# Local MinIO credentials
MINIO_LOCAL_AK=minioadmin123
MINIO_LOCAL_SK=minioadmin1234
````

---

## ☁️ Cloud configuration (optional)

If you want read/write access to Ersilia’s cloud buckets, set:

```bash
# Cloud endpoints
MINIO_ENDPOINT_CLOUD=http://83.48.73.209:8080

# Public bucket access
MINIO_CLOUD_AK="<access_key>"
MINIO_CLOUD_SK="<secret_key>"

# Private bucket access
MINIO_PRIV_CLOUD_AK="<access_key>"
MINIO_PRIV_CLOUD_SK="<secret_key>"
```

> 🔒 Keep credentials out of git. Use `.env` locally and secret managers in production.

---

## 🧾 Environment variables

### 🧭 Core

| Variable               |                    Default | What it does                                                        |
| ---------------------- | -------------------------: | ------------------------------------------------------------------- |
| `MINIO_ENDPOINT`       |    `http://127.0.0.1:9000` | S3 endpoint for the local MinIO store.                              |
| `MINIO_ENDPOINT_CLOUD` | `http://83.48.73.209:8080` | S3 endpoint for the cloud store (used by pull/push and cloud mode). |
| `NNS_ENDPOINT`         |    `http://127.0.0.1:8080` | Base URL for the NN/ANN service (index + approximate retrieval).    |
| `TIMEOUT`              |                     `3600` | Request timeout (seconds) for long-running network calls.           |

### 🪣 Buckets / namespaces

| Variable                      |          Default | What it does                                |
| ----------------------------- | ---------------: | ------------------------------------------- |
| `DEFAULT_BUCKET_NAME`         |  `isaura-public` | Default bucket for public precalculations.  |
| `DEFAULT_PRIVATE_BUCKET_NAME` | `isaura-private` | Default bucket for private precalculations. |

### 🔑 Credentials

| Variable              |          Default | What it does                                    |
| --------------------- | ---------------: | ----------------------------------------------- |
| `MINIO_LOCAL_AK`      |  `minioadmin123` | Local MinIO access key.                         |
| `MINIO_LOCAL_SK`      | `minioadmin1234` | Local MinIO secret key.                         |
| `MINIO_CLOUD_AK`      |           `None` | Cloud access key for public bucket operations.  |
| `MINIO_CLOUD_SK`      |           `None` | Cloud secret key for public bucket operations.  |
| `MINIO_PRIV_CLOUD_AK` |           `None` | Cloud access key for private bucket operations. |
| `MINIO_PRIV_CLOUD_SK` |           `None` | Cloud secret key for private bucket operations. |

### 🗂️ Storage & temp

| Variable          |                Default | What it does                                                               |
| ----------------- | ---------------------: | -------------------------------------------------------------------------- |
| `STORE_DIRECTORY` | `~/isaura/isaura-temp` | Local temp directory used for staging writes/reads and intermediate files. |

### ⚡ Performance (read/write)

| Variable            |   Default | What it does                                                   |
| ------------------- | --------: | -------------------------------------------------------------- |
| `MAX_ROWS_PER_FILE` |  `100000` | Max rows per parquet “chunk” written to the store.             |
| `CHECKPOINT_EVERY`  |   `50000` | How frequently the writer persists bloom/index checkpoints.    |
| `BATCH`             |   `10000` | Batch size for streaming inputs to the NN service.             |
| `FLUSH_EVERY`       |   `10000` | Flush frequency for NN streaming inserts.                      |
| `MAX_ROWS`          | `2000000` | Hard cap used by the writer for buffer/flush logic (internal). |

### 🔎 Indexing & deduplication

| Variable              |     Default | What it does                                                              |
| --------------------- | ----------: | ------------------------------------------------------------------------- |
| `BLOOM_FILENAME`      | `bloom.pkl` | Bloom filter filename used for deduplication and “seen” tracking.         |
| `MIN_NNS_RESULT_SIZE` |      `1000` | Minimum molecules required before ANN retrieval is enabled (approx mode). |

### 🧪 Input headers

| Variable  |              Default | What it does                                                         |
| --------- | -------------------: | -------------------------------------------------------------------- |
| `INPUT_C` | `["input","smiles"]` | Accepted input columns. Uses `input` if present, otherwise `smiles`. |

---

## 🧰 Common setups

### 🏠 Local-only development

* ▶️ Start services: `isaura engine --start`
* ✅ Use defaults for `MINIO_ENDPOINT`, buckets, and local credentials
* 🔍 Prefer exact search first (don’t use `-nn`) until enough data exists

### 🌍 Cloud read/write

* ☁️ Set `MINIO_ENDPOINT_CLOUD`
* 🔐 Export the correct cloud credentials (public and/or private)
* 🧭 Use `--cloud` where supported (or cloud-enabled Python classes)

---

## 🧷 Notes

* 🔎 ANN mode depends on the NN service and an index being built; if ANN fails, fall back to exact search.

---

