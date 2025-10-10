# The Isaura data store

This repository provides an interface and a CLI to the precalculated data available from the Ersilia Model Hub.

## Quick start guide
First install a high-performance python project manager called `uv` based on this [manual](https://docs.astral.sh/uv/getting-started/installation/)
### 1. Clone the repository

```bash
git clone https://github.com/ersilia-os/isaura.git
cd isaura
uv sync
source .venv/bin/activate
```
### 2. Install Minio (high performance object store) as below
### Prerequisites

- [Docker](https://www.docker.com/get-started) installed and running
- A terminal (Linux, macOS, or WSL on Windows)

---

### Step 1: Pull the MinIO Docker Image

```bash
docker pull minio/minio
````

---

### Step 2: Create a Data Directory

Create a local folder to persist data (optional but recommended):

```bash
mkdir -p ~/minio-data
```

---

### Step 3: Run MinIO Server with Environment Variables

Replace the example keys with your own secure values.

```bash
docker run -d \
  -p 9000:9000 \
  -p 9001:9001 \
  --name minio \
  -v ~/minio-data:/data \
  -e "MINIO_ROOT_USER=minioadmin" \
  -e "MINIO_ROOT_PASSWORD=minioadmin" \
  minio/minio server /data --console-address ":9001"
```

**Explanation:**

* `-v ~/minio-data:/data`: mounts your local folder
* `--console-address ":9001"`: enables the web console at port 9001
* Access UI: [http://localhost:9001](http://localhost:9001)
* S3 endpoint: [http://localhost:9000](http://localhost:9000)

---

### Step 4: Access the Web Console

Open your browser and go to:
👉 [http://localhost:9001](http://localhost:9001)

Login using:

```
Username: minioadmin
Password: minioadmin
```

---

### Step 5: Install MinIO Client (mc) for fine grained control and management over the object store

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
mc alias set local http://localhost:9000 minioadmin minioadmin
```
- Example command to list the projects for `local`:

```bash
mc ls local
```

You can find more detailed docs [here](https://github.com/minio/mc?tab=readme-ov-file) on how to use `mc`.

### Command at a Glance

| Command   | Alias | Required Options                                                      | Optional Options                                                                                                                                     | What it does                                                                                                                                       |                        |                                                                           |
| --------- | ----- | --------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------- | ------------------------------------------------------------------------- |
| `write`   | —     | `-i/--input-file`, `-m/--model`                                       | `-pn/--project-name`, `--access [public                                                                                                              | private                                                                                                                                            | both]`, `-v/--version` | Upload/write outputs for the given model & version using rows from a CSV. |
| `read`    | —     | `-i/--input-file`, `-m/--model`                                       | `-pn/--project-name`, `--access`, `-v/--version`, `-o/--output-file`                                                                                 | Read/download results for inputs in a CSV and optionally save as CSV/HDF5.                                                                         |                        |                                                                           |
| `copy`    | `cp`  | `-m/--model`, `-v/--version`, `-pn/--project-name`, `-o/--output-dir` | —                                                                                                                                                    | Copy all artifacts for a model/version from a project to a local directory. If `-o` is omitted in code, it logs counts; with `-o` it writes files. |                        |                                                                           |
| `move`    | `mv`  | `-m/--model`, `-v/--version`, `-pn/--project-name`                    | —                                                                                                                                                    | Move/relocate server-side artifacts for a model/version within the project space.                                                                  |                        |                                                                           |
| `remove`  | `rm`  | `-m/--model`, `-v/--version`, `-pn/--project-name`, `-y/--yes`        | —                                                                                                                                                    | Permanently delete artifacts for a model/version from a project. Safety-guarded by `--yes`.                                                        |                        |                                                                           |
| `inspect` | —     | (none strictly; behavior changes with flags)                          | `what` argument (`inputs`, default: `inputs`), `-m/--model`, `-v/--version`, `-pn/--project-name`, `--access`, `-i/--input-file`, `-o/--output-file` | Inspect available items or validate inputs. With `-i`, validates inputs and writes a report; without `-i`, lists available entries.                |                        |                                                                           |
| `catalog` | —     | —                                                                     | `-pn/--project-name`, `-f/--filter`                                                                                                                  | List models present in a project (bucket), optionally filtered by an id prefix.                                                                    |                        |                                                                           |

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

