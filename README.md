<div id="top"></div>
<p align="center">
  <img src="/assets/isaura.png" height="220" alt="Isuara logo">
</p>
<h2 align="center"> The Isaura data store</h2>

This repository provides an interface and a CLI to the precalculated data available from the Ersilia Model Hub. Initial benchmark is made and can be found [here](BENCHMARK.md).

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
| `inspect` | —     | (none strictly; behavior changes with flags)                          | `-m/--model`, `-v/--version`, `-pn/--project-name`, `--access`, `-i/--input-file`, `-o/--output-file`, `--cloud` (optional) | Inspect available items or validate inputs. With `-i`, validates inputs and writes a report; without `-i`, lists available entries.                |                        |                                                                           |
| `catalog` | —     | —                                                                     | `-pn/--project-name`, , `--cloud` (optional)                                                                                                                | List models present in a project (bucket), optionally filtered by an id prefix.                                                                    |                        |                                                                           |

### Brief CLI usage examples
> Buckets: are just a storage directory for model calcultaion (just a term used by the minio). 

| **Example**                      | **Command**                                                                                                   | **Description**                                                                                           |
| -------------------------------- | ------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| 🧾 **Write results calculation**           | `isaura write -i data/inputs.csv -m eos9876 -v v2 -pn myproject --access public`                        | Upload/write outputs for the given model and version using a CSV as input.                                |
| 📥 **Read results**              | `isaura read -i data/inputs.csv -m eos9876 -v v2 -pn myproject -o data/outputs.csv`                     | Read results corresponding to inputs and save them to an output CSV file.                                 |
| 📂 **Copy buckets**    | `isaura copy -m eos9876 -v v1 -pn myproject-private -o ~/Documents/files/`                                     | Copy all model artifacts from a project to a local directory.                                             |
| 🚚 **Move buckets**            | `isaura move -m eos9876 -v v1 -pn myproject-private`                                                           | Move or relocate artifacts for a specific model/version within the project.                               |
| 🗑️ **Remove buckets**         | `isaura remove -m eos9876 -v v1 -pn myproject-private --yes`                                                   | Permanently delete all artifacts for a model/version from a project (requires confirmation with `--yes`). |
| 🔍 **Inspect inputs (validate)** | `isaura inspect inputs -m eos9876 -v v1 -pn myproject -i data/inputs.csv -o reports/inspect_report.csv` | Validate input data for a model and output a report.                                                      |
| 📋 **List available model results**      | `isaura inspect -m eos9876 -v v1  -o reports/available.csv`                                      | List all available inputs or files related to a model/version.                                            |
| 📚 **Catalog project models**    | `isaura catalog -pn myproject`                                                                 | Display all models within a project                                    |

