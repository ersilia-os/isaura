<div align="center">

<img src="./isaura/assets/isaura_v2.png" height="160" alt="Isaura logo" />

### Ersilia's Precalculation Store

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
[Configuration](#configuration) ·
[Docs](#docs)

</div>

---

## Why Isaura?
Isaura is Ersilia's precalculation store: it **persistently stores model outputs** so researchers can retrieve results instantly instead of repeatedly running time-consuming inference. This delivers a major research speed-up — especially in low-resource settings where compute, bandwidth, or infrastructure are limited — by turning repeated calculations into reusable shared artifacts. To support equitable access, Ersilia also provides **free access to public precalculations**, making high-value model outputs available even when local compute isn't.

Isaura provides a structured store for model results so you can:

- ⚡ **Skip recomputation** by reusing precalculated outputs
- 🧱 Keep artifacts **versioned and organized** (model → version → bucket/project)
- 📦 Store and retrieve results via **S3-compatible object storage (MinIO)**
- 🔎 Enable **fast retrieval** using its engine built on top of DuckDB

---

## Installation

### Prerequisites

- **Python 3.10+** — [download here](https://www.python.org/downloads/)
- **Git** — [download here](https://git-scm.com/downloads)
- **Docker Desktop** — [download here](https://www.docker.com/products/docker-desktop/) — must be open before starting local services

---

### Step 1 — Clone and install

```bash
git clone https://github.com/ersilia-os/isaura.git
cd isaura
pip install -e .
```

> **Using uv?** Run `uv sync` instead and activate the environment with `source .venv/bin/activate`.

A local configuration file is created automatically at `~/.isaura/.env` with sensible defaults the first time you run any `isaura` command.

---

### Step 2 — Start local services

Make sure Docker Desktop is open, then run:

```bash
isaura engine --start
```

This starts a local MinIO instance and automatically creates the `isaura-public` and `isaura-private` buckets. You can explore the MinIO console at `http://localhost:9001` (user: `minioadmin123`, password: `minioadmin1234`).

---

### Step 3 — (Optional) Set up remote credentials

If you have access to Ersilia's remote store, add your cloud credentials interactively:

```bash
isaura configure --remote
```

You will be prompted for the cloud endpoint and access keys. Credentials are saved locally to `~/.isaura/.env` — nothing is sent anywhere.

---

### Step 4 — Verify your setup

```bash
isaura configure --test-credentials
```

This checks connectivity for local and cloud (if configured) and prints a result table. All rows you care about should show `✓ connected`.

```
┏━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┓
┃ Target        ┃ Bucket         ┃ Result      ┃
┡━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━┩
│ Local         │ isaura-public  │ ✓ connected │
│ Cloud public  │ isaura-public  │ ✓ connected │
│ Cloud private │ isaura-private │ ✓ connected │
└───────────────┴────────────────┴─────────────┘
```

---

## CLI

### Configuration

```bash
isaura configure                     # show current configuration
isaura configure --remote            # add or update remote/cloud credentials
isaura configure --update            # update a single credential interactively
isaura configure --show-secrets      # show all credential values unmasked
isaura configure --test-credentials  # test local and cloud connectivity
```

### Local services

```bash
isaura engine           # show status of Docker and MinIO
isaura engine --start   # start local MinIO
isaura engine --stop    # stop local MinIO
```

### Projects

A **project** is a named MinIO bucket used as a staging area before data reaches the canonical `isaura-public` or `isaura-private` buckets.

```bash
isaura info                  # list local projects with access type and creation date
isaura info --remote         # list remote projects

isaura create -pn myproject --access public   # create a new local project
isaura destroy -pn myproject                  # destroy an entire project
isaura destroy -pn myproject -m eos8a4x -v v1 # destroy a specific model version
```

> `isaura-public` and `isaura-private` are reserved — they cannot be fully destroyed, but individual model versions inside them can be removed with `--model-id` and `--version`.

### Writing outputs

Store model outputs in a project bucket:

```bash
isaura write -i data/ersilia_output.csv -m eos8a4x -v v1 -pn myproject
```

The model ID in `--model-id` is validated against the filename — if you accidentally pass a file from a different model, isaura will warn you before writing anything.

### Reading outputs

Retrieve stored outputs for a set of inputs:

```bash
# explicit version
isaura read -i data/inputs.csv -m eos8a4x -v v1 -pn myproject -o data/outputs.csv

# omit --version to automatically use the latest stored version
isaura read -i data/inputs.csv -m eos8a4x -pn myproject -o data/outputs.csv
```

`--project-name` is required. `--output` is optional — without it results are printed but not saved.

### Inspecting what is stored

Check which molecules from an input file are already cached:

```bash
isaura inspect --model_id eos8a4x -v v1 --access public -i data/inputs.csv -o reports/available.csv
```

Browse all models in a project:

```bash
isaura catalog -pn myproject           # local
isaura catalog -pn isaura-public -r    # remote
```

### Publishing to the cloud

The cloud hosts two canonical buckets: `isaura-public` and `isaura-private`. Write to a local project first, persist it into the canonical bucket, then push to cloud.

**Step 1 — write to your local project:**

```bash
isaura write -i data/ersilia_output.csv -m eos8a4x -v v1 -pn myproject
```

**Step 2 — persist into the canonical bucket:**

```bash
isaura persist -m eos8a4x -v v1 -pn myproject
```

This routes molecules tagged `public` → `isaura-public` and `private` → `isaura-private` based on the project's access setting.

**Step 3 — push to cloud:**

```bash
isaura push -m eos8a4x -v v1 -pn isaura-public
```

> Cloud credentials must be configured first with `isaura configure --remote`.

### Pulling from the cloud

Download precalculations from the remote store into your local MinIO:

```bash
# explicit version
isaura pull -i data/inputs.csv -m eos8a4x -v v1 -pn isaura-public

# omit --version to automatically pull the latest stored version
isaura pull -i data/inputs.csv -m eos8a4x -pn isaura-public
```

### Storage statistics

Generate a JSON inventory of all models in a bucket:

```bash
isaura stats -pn isaura-public -o ./reports
isaura stats -pn myproject --remote -o ./reports
```

---

## Configuration

Configuration is stored in `~/.isaura/.env` and created automatically on first run. You can view or update it at any time with:

```bash
isaura configure                 # view current config
isaura configure --update        # update any value interactively
isaura configure --remote        # add cloud credentials
```

See the full list of available variables: **[CONFIGURATION →](docs/CONFIGURATION.md)**

---

## Docs

* 📘 **How it works**: [here](docs/HOW_IT_WORKS.md)
* ⚙️ **Configuration**: [here](docs/CONFIGURATION.md)
* 🧰 **CLI and API reference**: [here](docs/API_AND_CLI_USAGE.md)
* 🧪 **Benchmark**: [here](docs/BENCHMARK.md)
* 🩹 **Troubleshooting / recovery**: [here](docs/TROUBLESHOOTING.md)

---

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](isaura/assets/Ersilia_Brand.png)
