<div align="center">

<img src="./isaura/assets/isaura_v2.png" height="160" alt="Isaura logo" />

### Ersilia’s Precalculation Store

Fast, reproducible access to **precalculated model outputs** from the **Ersilia Model Hub** — with a CLI and Python API built for batch workflows.

<br/>

[![Python](https://img.shields.io/badge/Python-%3E%3D3.10-3776AB?style=flat-square&logo=python&logoColor=white)](#)
[![uv](https://img.shields.io/badge/uv-supported-111111?style=flat-square&logo=astral&logoColor=white)](https://docs.astral.sh/uv/)
[![Docker](https://img.shields.io/badge/Docker-required-2496ED?style=flat-square&logo=docker&logoColor=white)](https://www.docker.com/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000?style=flat-square&logo=python&logoColor=white)](https://github.com/psf/black)
[![License](https://img.shields.io/badge/license-MIT-green?style=flat-square)](#license)

<br/>

[Installation](#installation) ·
[CLI](#cli) ·
[Python API](#python-api) ·
[Configuration](#configuration) ·
[Docs](#docs) ·
[Contributing](#contributing)

</div>

---

## Why Isaura?
Isaura is Ersilia’s precalculation store: it **persistently stores model outputs** so researchers can retrieve results instantly instead of repeatedly runningtime-consuming inference. This delivers a major research speed-up—especially in low-resource settings where compute, bandwidth, or infrastructure are limited—by turning repeated calculations into reusable shared artifacts. To support equitable access, Ersilia also provides **free access to public precalculations**, making high-value model outputs available even when local compute isn’t.

Isaura provides a structured store for model results so you can:

- ⚡ **Skip recomputation** by reusing precalculated outputs
- 🧱 Keep artifacts **versioned and organized** (model → version → bucket/project)
- 📦 Store and retrieve results via **S3-compatible object storage (MinIO)**  
- 🔎 Enable **fast retrieval** using its fast engine developed on top of duckdb and for ANN uses vector search / indexing components (Milvus + NN service)

---
## Architecture (high level)

* 📝 **Write:** `CLI / Python API → MinIO`
  Precomputed outputs are stored as chunked artifacts (e.g., Parquet) under `model_id/version`, and Isaura updates lightweight registries (index/metadata/bloom) for deduplication and fast lookup.

* 📥 **Read(exact):** `CLI / Python API → DuckDB query on MinIO → results`
  Inputs are matched against the index, then the corresponding rows are fetched directly from the stored chunks.

* ⚡ **Read (approx / ANN, optional):** `CLI / Python API → NN service (+ Milvus) → nearest match → exact fetch from MinIO`
  For unseen inputs, the NN service finds the closest indexed compound(s); Isaura then retrieves the corresponding stored result from MinIO.


See the deep dive: **[How it works →](docs/HOW_IT_WORKS.md)**

---

## Installation

### Prerequisites

Before installing Isaura, make sure you have the following:

- **Python 3.10+** — [download here](https://www.python.org/downloads/)
- **Git** — used to download the project ([download here](https://git-scm.com/downloads))
- **Docker** — required to run local services like MinIO. [Download Docker Desktop](https://www.docker.com/products/docker-desktop/)
- **Docker Compose** — use docker-compose v2
- **Write permissions**  — make sure you have write permissions in your filesystem

---

### Option A: Standard install *(recommended for most users)*

This is the simplest path. Open a terminal and run the following commands one by one.

**1. Clone the repository** — this downloads the project to your computer:

```bash
git clone https://github.com/ersilia-os/isaura.git
```

**2. Navigate into the project folder:**

```bash
cd isaura
```

**3a. Install Isaura through pip:**

```bash
conda activate <your_env>
pip install -e .
```

The `-e` flag installs it in "editable" mode, meaning any changes you make to the source code are reflected immediately without reinstalling.

---

### Option B: Install with uv *(recommended for developers)*
[uv](https://docs.astral.sh/uv/) is a faster alternative to pip. If you don't have it yet, [install it first](https://docs.astral.sh/uv/getting-started/installation/).

```bash
git clone https://github.com/ersilia-os/isaura.git
cd isaura
uv sync #creates an isolated virtual environment and installs all dependencies automatically.
source .venv/bin/activate  # on Windows: .venv\Scripts\activate
```

---

### Verify installation

Once installed, confirm everything is working by running:

```bash
isaura --help
```

You should see the list of available commands printed to your terminal.

---

### Start local services

Isaura relies on local infrastructure (MinIO for storage, and optionally Milvus + NNS for approximate search). Make sure Docker is running, then start everything with:

```bash
isaura engine --start
```

Local dashboards once running:
- MinIO Console: `http://localhost:9001`

Default local credentials:
```
Username: minioadmin123
Password: minioadmin1234
```
If you plan to upload/download large volumes of data, we recommend disabling Milvus if you will not use the NNS search, as the indexing can use a lot of memory:

```
docker stop milvus-standalone
```

---

## CLI

### Common commands

#### Write (store outputs)

```bash
isaura write -i data/ersilia_output.csv -m eos8a4x -v v2 -pn myproject --access public
```

#### Read (retrieve outputs)

```bash
isaura read -i data/inputs.csv -m eos8a4x -v v2 -pn myproject -o data/outputs.csv
```

#### Copy artifacts to local directory

```bash
isaura copy -m eos8a4x -v v1 -pn myproject -o ~/Documents/isaura-backup/
```

#### Inspect available entries

```bash
isaura inspect -m eos8a4x -v v1 -o reports/available.csv
```

#### Upload to cloud store

The cloud only hosts two canonical buckets: `isaura-public` and `isaura-private`. If your local work uses a custom project name, you need to copy (or move) it into the appropriate canonical bucket first, then push to cloud.

**Step 1 — write outputs to your local project:**

```bash
isaura write -i data/ersilia_output.csv -m eos8a4x -v v1 -pn myproject --access public
```

**Step 2 — copy (or move) into the canonical bucket:**

Isaura routes each entry automatically based on the `--access` tag set during write: `public` → `isaura-public`, `private` → `isaura-private`.

```bash
# copy (keeps data in myproject as well)
isaura copy -m eos8a4x -v v1 -pn myproject

# or move (removes data from myproject after copying)
isaura move -m eos8a4x -v v1 -pn myproject
```

**Step 3 — push the canonical bucket to cloud:**

```bash
isaura push -m eos8a4x -v v1 -pn isaura-public
# or for private data:
isaura push -m eos8a4x -v v1 -pn isaura-private
```

Cloud credentials must be set beforehand (in `.env` or exported in the terminal in each session):

```bash
export MINIO_ENDPOINT_CLOUD="<cloud-endpoint>"
export MINIO_CLOUD_AK="<access-key>"        # public bucket
export MINIO_CLOUD_SK="<secret-key>"
export MINIO_PRIV_CLOUD_AK="<access-key>"   # private bucket
export MINIO_PRIV_CLOUD_SK="<secret-key>"
```

> See [CONFIGURATION](docs/CONFIGURATION.md) for the full list of env vars.


---

## Python API

```python
from isaura.manage import IsauraWriter, IsauraReader
```
Write the precalculation
```python
writer = IsauraWriter(
    input_csv="data/input.csv",
    model_id="eos8a4x",
    model_version="v1",
    bucket="my-project",
    access="public",
)
writer.write()
```
Read the stored calculation
```python
reader = IsauraReader(
    model_id="eos8a4x",
    model_version="v1",
    bucket="my-project",
    input_csv="data/query.csv",
    approximate=False,
)
reader.read(output_csv="results.csv")
```

More examples for CLI and API usage: **[API and CLI usage](docs/API_AND_CLI_USAGE.md)**

---

## Configuration

Isaura reads configuration from environment variables.

### Recommended: `.env`

Create a `.env` file in the repo root:

```bash
MINIO_ENDPOINT=http://127.0.0.1:9000
NNS_ENDPOINT=http://127.0.0.1:8080
DEFAULT_BUCKET_NAME=isaura-public
DEFAULT_PRIVATE_BUCKET_NAME=isaura-private
```

### Cloud credentials (optional)

```bash
export MINIO_CLOUD_AK="<access_key>"
export MINIO_CLOUD_SK="<secret_key>"

export MINIO_PRIV_CLOUD_AK="<access_key>"
export MINIO_PRIV_CLOUD_SK="<secret_key>"
```
> You can define those credentials in the .env as well

See the full list: **[CONFIGURATION](docs/CONFIGURATION.md)**

---

## MinIO Client (optional but recommended)

Install `mc` to manage buckets:

```bash
brew install minio/stable/mc   # macOS
# or Linux:
curl -O https://dl.min.io/client/mc/release/linux-amd64/mc && chmod +x mc && sudo mv mc /usr/local/bin/
```

Configure alias:

```bash
mc alias set local http://localhost:9000 minioadmin123 minioadmin1234
mc ls local
```

---

## Docs

* 📘 **How it works**: [here](docs/HOW_IT_WORKS.md)
* ⚙️ **Configuration**: [here](docs/CONFIGURATION.md)
* 🧰 **CLI and API reference**: [here](docs/API_AND_CLI_USAGE.md)
* 🧪 **Benchmark**: [here](docs/BENCHMARK.md)
* 🩹 **Troubleshooting / recovery**: [here](docs/TROUBLESHOOTING.md)

---

## Contributing

PRs are welcome. Please run format + lint before pushing:

```bash
uv run ruff format .
```

If you’re changing CLI behavior, please update **[here](docs/API_AND_CLI_USAGE.md)**.

---

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](isaura/assets/Ersilia_Brand.png)
