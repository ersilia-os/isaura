<div id="top"></div>

<p align="center">
  <img src="./isaura/assets/isaura.png" height="110" alt="Isaura logo"><br><br>
  <img src="https://img.shields.io/badge/python-%3E%3D3.8-blue.svg?style=flat-square&logo=python&logoColor=white" alt="Python >=3.8">
  <img src="https://img.shields.io/badge/code%20style-black-000000.svg?style=flat-square&logo=python&logoColor=white" alt="Code style: Black">
</p>


<h3 align="center">The Ersilia's Precalculation Store</h3>




This repository provides an interface and a CLI to the precalculated data available from the Ersilia Model Hub. Initial benchmark is made and can be found [here](BENCHMARK.md). And the detailed mechanism on to how it works [here](HOW_IT_WORKS.md).

---

## Quick start guide
First install a high-performance python project manager called `uv` based on this [manual](https://docs.astral.sh/uv/getting-started/installation/)
### 1. Clone the repository

```bash
git clone https://github.com/ersilia-os/isaura.git
cd isaura
uv sync
source .venv/bin/activate
```
### 2. Install all isaura services
#### Prerequisites

- [Docker](https://www.docker.com/get-started) installed and running
- `docker-compose` ubuntu we need to have them installed. Use this instruction for more detail [here](https://docs.docker.com/engine/install/ubuntu/).
- `docker-compose` macOS as `brew install docker-compose`

---
#### Fastest way to start all the services
```bash
isaura engine --start
```

#### 3. Install MinIO Client (mc) for fine grained control and management over the object store

The MinIO Client (`mc`) is a command-line tool to interact with MinIO or any S3-compatible storage.

#### Install (Linux/macOS)

```bash
curl -O https://dl.min.io/client/mc/release/linux-amd64/mc
chmod +x mc
sudo mv mc /usr/local/bin/
```

#### Or with Homebrew (macOS)

```bash
brew install minio/stable/mc
```

---

#### Configure the MinIO Client

```bash
mc alias set local http://localhost:9000 minioadmin123 minioadmin1234
```
- Example command to list the projects for `local`:

```bash
mc ls local
```

You can find more detailed docs [here](https://github.com/minio/mc?tab=readme-ov-file) on how to use `mc`.

#### Access the Web Console

Open your browser and go to:
👉 [http://localhost:9001](http://localhost:9001)

Login using:

```
Username: minioadmin123
Password: minioadmin1234
```

---
### Cloud functionality usage
You can export the following env varibles:
- For read/write the could public data
```bash
export MINIO_CLOUD_AK = <Key here> # access key
export MINIO_CLOUD_SK = <Key here> # secrete key
```
- For read/write the could private data
```bash
export MINIO_PRIV_CLOUD_AK = <Key here> # access key
export MINIO_PRIV_CLOUD_SK = <Key here> # secrete key
```



### Command at a Glance

| Command   | Alias | Required Options                                                      | Optional Options                                                                                                                                     | What it does                                                                                                                                       |                        |                                                                           |
| --------- | ----- | --------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------- | ------------------------------------------------------------------------- |
| `write`   | —     | `-i/--input-file`, `-m/--model`                                       | `-pn/--project-name`, `--access [public                                                                                                              | private                                                                                                                                            | both]`, `-v/--version` | Upload/write outputs for the given model & version using rows from a CSV. The output csv file should have an input header named `input` not `smiles`. This prevent collision when uploading the results for the second time and better for standardize.  |
| `read`    | —     | `-i/--input-file`, `-m/--model`                                       | `-pn/--project-name`, `--access`, `-v/--version`, `-o/--output-file`, `-nn`                                                                                 | Read/download results for inputs in a CSV and optionally save as CSV/HDF5.                                                                         |                        |                                                                           |
| `copy`    | `cp`  | `-m/--model`, `-v/--version`, `-pn/--project-name`, `-o/--output-dir` | —                                                                                                                                                    | Copy all artifacts for a model/version from a project to a local directory. If `-o` is omitted in code, it logs counts; with `-o` it writes files. |                        |                                                                           |
| `move`    | `mv`  | `-m/--model`, `-v/--version`, `-pn/--project-name`                    | —                                                                                                                                                    | Move/relocate server-side artifacts for a model/version within the project space.                                                                  |                        |                                                                           |
| `remove`  | `rm`  | `-m/--model`, `-v/--version`, `-pn/--project-name`, `-y/--yes`        | —                                                                                                                                                    | Permanently delete artifacts for a model/version from a project. Safety-guarded by `--yes`.                                                        |                        |                                                                           |
| `inspect` | —     | `-m/--model`, `-v/--version`, `-o/--output-file`                       | `-pn/--project-name`, `--access`, `-i/--input-file`, `--cloud`| Inspect available items or validate inputs. With `-i`, validates inputs and writes a report; without `-i`, lists available entries.                |                        |                                                                           |
| `catalog` | —     | `-pn/--project-name`                                                                    |  `--cloud`                                                                                                            | List models present in a project (bucket), optionally filtered by an id prefix.                                                                    |                        |                                                                           |

### Brief CLI usage examples
> Buckets: are just a storage directory for model calcultaion (just a term used by the minio). 

| **Simple Example**                      | **Command**                                                                                                   | **Description**                                                                                           |
| -------------------------------- | ------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| 🧾 **Write results calculation**           | `isaura write -i data/ersilia_output.csv -m eos8a4x -v v2 -pn myproject --access public`                        | Upload/write outputs(needs to have the input column named as `input`) for the given model and version using a CSV as input.                                |
| 📥 **Read results**              | `isaura read -i data/inputs.csv -m eos8a4x -v v2 -pn myproject -o data/outputs.csv`                     | Read results corresponding to inputs and save them to an output CSV file.                                 |
| 📂 **Copy buckets**    | `isaura copy -m eos8a4x -v v1 -pn myproject-private -o ~/Documents/files/`                                     | Copy all model artifacts from a project to a local directory.                                             |
| 🚚 **Move buckets**            | `isaura move -m eos9876 -v v1 -pn myproject-private`                                                           | Move or relocate artifacts for a specific model/version within the project.                               |
| 🗑️ **Remove buckets**         | `isaura remove -m eos8a4x -v v1 -pn myproject-private --yes`                                                   | Permanently delete all artifacts for a model/version from a project (requires confirmation with `--yes`). |
| 🔍 **Inspect inputs (validate)** | `isaura inspect inputs -m eos8a4x -v v1 -pn myproject -i data/inputs.csv -o reports/inspect_report.csv` | Validate input data for a model and output a report.                                                      |
| 📋 **List available model results**      | `isaura inspect -m eos8a4x -v v1  -o reports/available.csv`                                      | List all available inputs or files related to a model/version.                                            |
| 📚 **Catalog project models**    | `isaura catalog -pn myproject`                                                                 | Display all models within a project                                    |

### Troubleshooting Initial Setup Failures (Ersilia ↔ Isaura Integration)

This guide lists practical cleanup + restart steps to resolve failures during the **initial setup stage** of the **Ersilia + Isaura integration**, especially when retrieval stays slow after setup.

---

#### When to use this

Use these steps if you see any of the following:

* Setup fails during first-time initialization
* Models don’t load correctly
* Retrieval is slow even after a successful run
* You suspect stale buckets/volumes/containers are causing inconsistent behavior

---

### Prerequisites

* You have access to:

  * Local MinIO console at `http://localhost:9000/`
  * Docker (for containers/images)
  * Isaura CLI (`isaura`)
  * Ersilia CLI (`ersilia`)
* You know the model ID you’re testing (e.g. `eosxxxx`)

---

### Step-by-step troubleshooting

#### 1) Clean stale model artifacts from MinIO

1. Open: `http://localhost:9000/`
2. In each of these buckets:

   * `isaura-public`
   * `isaura-private`
   * `ersilia`
3. Locate the **model of interest** and **remove its stored artifacts** (select → remove).

This helps eliminate stale cached model files and index artifacts that can break or slow setup.

---

#### 2) Remove old Milvus volumes (local data reset)

If you’re using local Milvus storage, remove old volumes:

```bash
sudo rm -rf ~/isaura
```

⚠️ **Warning:** This deletes local Isaura/Milvus persisted state.

---

#### 3) Remove Isaura-related Docker images and containers

Remove all containers and images related to:

* Milvus
* `ersiliaos/nns`

Use your usual Docker cleanup approach (stop/remove containers, then remove images).
This is meant to guarantee you’re not running with stale layers or a broken container state.

---

#### 4) Restart the Milvus container

After cleanup, restart Milvus so it comes up fresh and healthy.

(If you manage Milvus via `docker compose`, restart through that; otherwise restart the container directly.)

---

#### 5) Reinstall Isaura and reinitialize the engine

If there were any changes (or you suspect version/config mismatch), reinstall Isaura, then run:

```bash
isaura engine -s
```

This reboots the Isaura engine setup path and recreates needed runtime state.

---

#### 6) Re-serve the model and run a retrieval test

Serve the model:

```bash
ersilia serve eosxxxx -rs -ws -a public
```

Run a batch test:

```bash
ersilia run -i input.csv -o output.csv -b 10000
```

---

#### 7) Run it twice and compare performance

Run the same `ersilia run ...` command **two times**.

Expected behavior:

* **First run:** may be slower (initial indexing/warmup)
* **Second run:** should show **fast retrieval**

If the **second run is still slow**, the issue is likely not resolved.

---

### If it’s still slow after the second run

At that point:

* Contact the admin
* Open an issue with:

  * the model ID (`eosxxxx`)
  * logs from Milvus + Isaura + Ersilia (startup + run)
  * confirmation you cleaned MinIO buckets and removed `~/isaura`

#### API usage examples
```python
from isaura.manage import (
    IsauraWriter,
    IsauraReader,
    IsauraMover,
    IsauraCopy,
    IsauraRemover,
    IsauraInspect,
    IsauraPull,
    IsauraPush,
)


writer = IsauraWriter(
    input_csv="data/input.csv",
    model_id="eos8a4x",
    model_version="v1",
    bucket="my-project",
    access="public",  # can be 'public', 'private', or 'both'
)
writer.write()


reader = IsauraReader(
    model_id="eos8a4x",
    model_version="v1",
    bucket="my-project",
    input_csv="data/query.csv",
    approximate=False,  # use ANN if True
)
reader.read(output_csv="results.csv")


puller = IsauraPull(
    model_id="eos8a4x",
    model_version="v1",
    bucket="my-project",
    input_csv="data/ids.csv",
)
puller.pull()


pusher = IsauraPush(
    model_id="eos8a4x",
    model_version="v1",
    bucket="my-project",
)
pusher.push()


copier = IsauraCopy(
    model_id="eos8a4x",
    model_version="v1",
    bucket="my-project",
    output_dir="backups/",
)
copier.copy()


mover = IsauraMover(
    model_id="eos8a4x",
    model_version="v1",
    bucket="my-project",
)
mover.move()


remover = IsauraRemover(
    model_id="eos8a4x",
    model_version="v1",
    bucket="my-project",
)
remover.remove()


inspector = IsauraInspect(
    model_id="eos8a4x",
    model_version="v1",
    project_name="my-project",
    access="public",
    cloud=False,
)

# List available inputs
df_inputs = inspector.list_available(output_file="inputs.csv")

# Inspect specific input CSV
df_inspected = inspector.inspect_inputs("data/input.csv", "inspected_results.csv")

# Inspect all models in a project
df_models = inspector.inspect_models("my-project")
```
