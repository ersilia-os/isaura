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
Username: admin
Password: supersecretkey123
```

---

### Step 5: Install MinIO Client (mc)

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

## Configure the MinIO Client

Add your local MinIO server as a host:

```bash
mc alias set local http://localhost:9000 minioadmin minioadmin
```