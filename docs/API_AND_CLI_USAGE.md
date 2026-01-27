
# Python API

Isaura exposes a small set of high-level managers for writing, reading, and moving precalculated artifacts.

> Tip: Most workflows are “write once, read many times”. After the first write, reads become fast.

---

## Install & import

```python
from isaura.manage import (
  IsauraWriter,
  IsauraReader,
  IsauraCopy,
  IsauraMover,
  IsauraRemover,
  IsauraInspect,
  IsauraPull,
  IsauraPush,
  IsauraStat,
)
````

---

## Data format

### Inputs

Isaura expects a column named:

* `input` (preferred), or
* `smiles` (fallback)

### Writer: output schema

`IsauraWriter` uploads full rows, preserving all CSV columns.
It uses `input`/`smiles` as the lookup key and stores the rest as payload columns.

---

## Write (upload/store results)

Writes a CSV (or DataFrame) into the Isaura store for a given:

* `model_id` (e.g. `eos8a4x`)
* `model_version` (e.g. `v1`)
* `bucket` (project namespace; default is `isaura-public`)

```python
from isaura.manage import IsauraWriter

with IsauraWriter(
  input_csv="data/ersilia_output.csv",
  model_id="eos8a4x",
  model_version="v1",
  bucket="isaura-public",     # optional (defaults to public bucket)
  access="public",            # metadata tag; typically public/private
  endpoint=None,              # optional override MinIO endpoint
  access_key=None,            # optional MinIO access key
  secrete=None,               # optional MinIO secret key
) as w:
  w.write()
```

### Write from a DataFrame

```python
import pandas as pd
from isaura.manage import IsauraWriter

df = pd.read_csv("data/ersilia_output.csv")

with IsauraWriter(
  input_csv="data/ersilia_output.csv",  
  model_id="eos8a4x",
  model_version="v1",
  bucket="isaura-public",
) as w:
  w.write(df=df)
```

### Deduplication behavior

* Isaura uses a bloom filter + index to skip duplicates.
* If an `input` is already present, it will be skipped.

---

## Read (retrieve results)

Reads rows matching inputs from a CSV (or DataFrame).
You can run either:

* **Exact** search: `approximate=False`
* **Approximate/ANN** search: `approximate=True` (requires enough indexed data)

```python
from isaura.manage import IsauraReader

r = IsauraReader(
  model_id="eos8a4x",
  model_version="v1",
  input_csv="data/inputs.csv",
  approximate=False,
  bucket="isaura-public",
  endpoint=None,
  access_key=None,
  secrete=None,
)

df = r.read(output_csv="data/outputs.csv")
```

### Approximate retrieval (ANN)

```python
r = IsauraReader(
  model_id="eos8a4x",
  model_version="v1",
  input_csv="data/inputs.csv",
  approximate=True,
)

df = r.read()
```

**Notes**

* ANN is only enabled when the index size is above a minimum (`MIN_NNS_RESULT_SIZE`).
* ANN contacts an NN service using `NNS_ENDPOINT` and may trigger/ensure index building.

---

## Inspect (discover what exists)

Use `IsauraInspect` to list available inputs and models across buckets.

```python
from isaura.manage import IsauraInspect

insp = IsauraInspect(
  cloud=False,
  project_name=None,
  access="both",     # "public" | "private" | "both"
  endpoint=None,
)

# List unique inputs that exist and which bucket owns them
df_available = insp.list_available(output_csv="reports/available_inputs.csv")

# Check if a given input CSV is available
df_check = insp.inspect_inputs(
  input_csv="data/inputs.csv",
  output_csv="reports/inspect_report.csv",
)

# Inspect models in a bucket (model_id/version + entry counts)
models = insp.inspect_models(bucket="isaura-public")
```

---

## Copy / Move / Remove

These are “server-side” operations on stored artifacts for a model/version/bucket.

```python
from isaura.manage import IsauraCopy, IsauraMover, IsauraRemover

IsauraCopy(model_id="eos8a4x", model_version="v1", bucket="isaura-public", output_dir="backup/").copy()

IsauraMover(model_id="eos8a4x", model_version="v1", bucket="isaura-public", output_dir="newplace/").move()

IsauraRemover(model_id="eos8a4x", model_version="v1", bucket="isaura-public", output_dir=None).remove()
```

---

## Pull / Push (cloud sync workflows)

### Pull

Pulls from cloud into a local workspace by reading from the cloud endpoint and materializing objects.

```python
from isaura.manage import IsauraPull

p = IsauraPull(
  model_id="eos8a4x",
  model_version="v1",
  bucket="isaura-public",
  input_csv="data/inputs.csv",
)
p.pull()
```

### Push

Pushes local public/private data to cloud by:

1. listing available inputs,
2. splitting by bucket,
3. reading locally,
4. writing to cloud.

```python
from isaura.manage import IsauraPush

IsauraPush(model_id="eos8a4x", model_version="v1", bucket="my-project").push()
```

> Cloud credentials must be provided via env vars (see `docs/CONFIGURATION.md`).

---

## Stats (inventory report)

Generate a JSON report of what’s stored, including optional column counts.

```python
from isaura.manage import IsauraStat

st = IsauraStat(
  project_name=None,
  access="both",
  cloud=False,
  include_columns=True,
  include_column_names=False,
)
out_path = st.write_json("reports/isaura_stats.json")
print(out_path)
```

The stats report includes:

* buckets scanned
* model counts by bucket
* molecules per model
* histogram of “how many models each molecule appears in”
* optional parquet column info (best effort)

<br>

---


# CLI usage

Isaura ships with a CLI for writing/reading/store operations and local engine management.

> The CLI mirrors the Python API: write/read/inspect/copy/move/remove + engine start.

---

## Quick start

```bash
isaura engine --start
````

Local dashboards:

* MinIO Console: `http://localhost:9001`
* Milvus API: `http://localhost:8080`

---

## Common commands

### Write (upload/store results)

Uploads outputs for a model/version from a CSV.
Your CSV should contain an input identifier column (`input` preferred; `smiles` supported).

```bash
isaura write \
  -i data/ersilia_output.csv \
  -m eos8a4x \
  -v v1 \
  -pn myproject \
  --access public
```

**Behavior**

* Skips duplicates using the bloom/index registry.
* Persists metadata to `access.json`.

---

### Read (retrieve results)

Fetches stored rows corresponding to inputs in a CSV.

```bash
isaura read \
  -i data/inputs.csv \
  -m eos8a4x \
  -v v1 \
  -pn myproject \
  -o data/outputs.csv
```

### Read with approximate matching (ANN)

```bash
isaura read \
  -i data/inputs.csv \
  -m eos8a4x \
  -v v1 \
  -pn myproject \
  -o data/outputs.csv \
  -nn
```

**Notes**

* ANN requires a sufficiently large index.
* Uses `NNS_ENDPOINT` and may trigger/build ANN index server-side.

---

## Inspect & catalog

### Inspect inputs (validate availability)

```bash
isaura inspect inputs \
  -m eos8a4x \
  -v v1 \
  -pn myproject \
  -i data/inputs.csv \
  -o reports/inspect_report.csv
```

### List available entries for a model/version

```bash
isaura inspect \
  -m eos8a4x \
  -v v1 \
  -o reports/available.csv
```

### Catalog models in a project (bucket)

```bash
isaura catalog -pn myproject
```

---

## Transfer operations

### Copy artifacts (server → local)

```bash
isaura copy \
  -m eos8a4x \
  -v v1 \
  -pn myproject \
  -o ~/Documents/isaura-backup/
```

### Move artifacts (server-side relocate)

```bash
isaura move \
  -m eos8a4x \
  -v v1 \
  -pn myproject
```

### Remove artifacts (permanent)

```bash
isaura remove \
  -m eos8a4x \
  -v v1 \
  -pn myproject \
  --yes
```

---

## Cloud sync

### Pull from cloud → local

```bash
isaura pull \
  -m eos8a4x \
  -v v1 \
  -pn isaura-public \
  -i data/inputs.csv
```

### Push local → cloud

```bash
isaura push \
  -m eos8a4x \
  -v v1 \
  -pn myproject
```

Cloud auth is via env vars:

* `MINIO_CLOUD_AK`, `MINIO_CLOUD_SK`
* `MINIO_PRIV_CLOUD_AK`, `MINIO_PRIV_CLOUD_SK`
* `MINIO_ENDPOINT_CLOUD`

---

## Configuration (env vars)

Common:

* `MINIO_ENDPOINT` (default `http://127.0.0.1:9000`)
* `NNS_ENDPOINT` (default `http://127.0.0.1:8080`)
* `DEFAULT_BUCKET_NAME` (default `isaura-public`)
* `DEFAULT_PRIVATE_BUCKET_NAME` (default `isaura-private`)
* `BATCH`, `FLUSH_EVERY`, `TIMEOUT`

Recommended: define them in a `.env` file (loaded automatically if `python-dotenv` is installed).

---

## Reference: flags (most used)

| Flag                  | Meaning                            |
| --------------------- | ---------------------------------- |
| `-i, --input-file`    | Input CSV                          |
| `-o, --output-file`   | Output CSV                         |
| `-m, --model`         | Model ID (e.g. `eos8a4x`)          |
| `-v, --version`       | Model version (e.g. `v1`)          |
| `-pn, --project-name` | Bucket/project name                |
| `--access`            | `public` / `private` / `both`      |
| `-nn`                 | Enable approximate (ANN) retrieval |
| `--cloud`             | Use cloud endpoint/credentials     |

---

