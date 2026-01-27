## How it wotks

<p align="center">
  <img src="../isaura/assets/isaura_mechanism.png" alt="Isaura logo"><br><br>
</p>



### Write / ingest
- Accepts **Ersilia outputs format** for any model
- Uses **DuckDB** as the **writing engine** to produce Parquet chunks
- Stores everything in **MinIO** (object storage)

### Read / fetch
You provide a **CSV file** containing a list of inputs.

- **Exact search**
  - Isaura uses **DuckDB** to query Parquet directly from MinIO
  - Returns results in **Ersilia outputs format**

- **Approximate search**
  - Isaura uses **Milvus** (via **NNS**) to find the **top-1 most similar stored input**
  - Then Isaura fetches the cached output for that matched input using **DuckDB + MinIO**
  - Returns results in **Ersilia outputs format**

---

## Architecture

- User interacts through **Isaura CLI**
- **DuckDB** is the query + write engine over Parquet
- **MinIO** stores Parquet data + indexes
- **Milvus** stores input fingerprints for approximate search
- **NNS** is a Go-based REST API sitting in front of Milvus for high-performance ingest/query from the Milvus:
  - https://github.com/ersilia-os/nn-search/tree/main/api

---

## Components (quick summary)

### DuckDB
Isaura’s compute engine for:
- writing Parquet chunks during ingestion
- querying Parquet during exact retrieval

DuckDB reads/writes data in MinIO via its S3/HTTPFS support.

**DuckDB “URL” pattern used by Isaura (logical path):**
- Parquet data:
  - `s3://<bucket>/<eos_id>/<version>/data/chunk_*.parquet`
- Bloom filter index:
  - `s3://<bucket>/<eos_id>/<version>/bloom.pkl`
- Access control metadata:
  - `s3://<bucket>/<eos_id>/<version>/access.json`

> Note: With MinIO, DuckDB is typically configured with an S3 endpoint (e.g. `http://<minio-host>:9000`) plus credentials; the object paths remain `s3://...`.

### MinIO
Object storage for:
- cached outputs as Parquet chunk files
- `bloom.pkl` membership index
- `access.json` access metadata per input

### Milvus (approximate-only)
Milvus is used **only** for approximate search.

- You query with a list of inputs
- Milvus returns, for each input, the **top-1 nearest** stored input
- Similarity: **Jaccard similarity**
- Current representation: **1024-bit Morgan fingerprint** (CPD inputs)
- Collection naming convention:
  - `{ersilia_eos_id}_{version}`

Milvus stores *input fingerprints and identifiers*, not full Ersilia outputs.

### NNS (Milvus REST gateway)
NNS is a Go-based REST API server to ingest/query Milvus with high performance:
- https://github.com/ersilia-os/nn-search/tree/main/api

---

## Projects (buckets) and defaults

Isaura organizes data in **projects**, backed by MinIO buckets.

### Default projects
- `isaura-public`
- `isaura-private`

You can create custom projects/buckets as needed (see CLI commands below / `--help`).

---

## Storage layout

Example structure for a single model (`eosid`) and version:

```
isaura-public/
eosid/
version/
bloom.pkl
access.json
data/
chunk_{idx}.parquet

```

### Parquet chunking rule
Each `chunk_{idx}.parquet` contains **up to 2,000,000 rows** (2M max) for performance reasons (DuckDB scan efficiency / row grouping).

---

## Fast “does this input exist?” checks (Bloom filter)

Isaura keeps a Bloom filter per `(project, model, version)` at:

- `<bucket>/<eos_id>/<version>/bloom.pkl`

This allows Isaura to quickly rule out missing inputs before doing heavier Parquet scans.

> Bloom filters may produce false positives, but not false negatives.

---

## Access control with `access.json`

Each `(project, model, version)` also includes:

- `<bucket>/<eos_id>/<version>/access.json`

This file stores **inputs and their access classification** (e.g. `public` or `private`).

Isaura uses this to decide where results should live in the default buckets.

---

## Copying calculations into default buckets (and Milvus registration)

When you **copy calculations** for a given:
- `project + model (eos_id) + version`

Isaura performs the following:

1. Reads `<custom_project>/<eos_id>/<version>/access.json`
2. For each input:
   - writes/stores it into `isaura-public` **or** `isaura-private` based on access metadata
3. Updates the Bloom filter(s) accordingly (`bloom.pkl`)
4. **Registers inputs into Milvus** (this happens during copy from a custom project to the default buckets)

This means:
- Inputs are **indexed in Milvus** at the moment they become part of the default store
- Approximate search operates over what has been copied/registered into the default store

