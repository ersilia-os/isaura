
# Python API

Isaura exposes a small set of high-level managers for writing, reading, and moving precalculated artifacts.

> Tip: Most workflows are "write once, read many times". After the first write, reads become fast.

---

## Install & import

```python
from isaura.manage import (
  IsauraWriter,
  IsauraReader,
  IsauraCopy,
  IsauraMolRemover,
  IsauraInspect,
  IsauraPull,
  IsauraPush,
  IsauraStat,
)
```

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

Writes a CSV (or DataFrame) into the Isaura store for a given model ID, version, and project bucket.

```python
from isaura.manage import IsauraWriter

with IsauraWriter(
  input_csv="data/ersilia_output.csv",
  model_id="eos8a4x",
  model_version="v1",
  bucket="myproject",
  endpoint=None,        # optional MinIO endpoint override
  access_key=None,      # optional MinIO access key
  secrete=None,         # optional MinIO secret key
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
  bucket="myproject",
) as w:
  w.write(df=df)
```

### Deduplication behavior

* Isaura uses a bloom filter + index to skip duplicates.
* If an `input` is already present, it will be skipped.

---

## Read (retrieve results)

Reads rows matching inputs from a CSV (or DataFrame).

```python
from isaura.manage import IsauraReader

with IsauraReader(
  model_id="eos8a4x",
  model_version="v1",
  input_csv="data/inputs.csv",
  approximate=False,
  bucket="myproject",
  endpoint=None,
  access_key=None,
  secrete=None,
) as r:
  df = r.read(output_csv="data/outputs.csv")
```

---

## Inspect (discover what exists)

Use `IsauraInspect` to list available inputs and models across buckets.

```python
from isaura.manage import IsauraInspect

insp = IsauraInspect(
  cloud=False,
  project_name="myproject",
  access="both",     # "public" | "private" | "both"
  endpoint=None,
)

# Check if a given input CSV is available
df_check = insp.inspect_inputs(
  input_csv="data/inputs.csv",
  output_csv="reports/inspect_report.csv",
)

# Inspect models in a bucket (model_id/version + entry counts)
models = insp.inspect_models(bucket="myproject")
```

---

## Persist / Remove

These are server-side operations on stored artifacts for a model/version/bucket.

```python
from isaura.manage import IsauraCopy, IsauraMolRemover

# Persist: route rows from a project bucket into isaura-public / isaura-private
IsauraCopy(model_id="eos8a4x", model_version="v1", bucket="myproject").copy()

# Remove specific molecules (listed in a CSV) from a model/version in a bucket
with IsauraMolRemover(
  model_id="eos8a4x", model_version="v1", bucket="myproject", input_csv="data/to_remove.csv"
) as r:
  n_removed, n_not_found = r.remove()
```

> To delete an entire model/version or a whole project bucket, use the CLI:
> `isaura destroy -pn myproject -m eos8a4x -v v1` (or `isaura destroy -pn myproject`).

---

## Pull / Push (cloud sync)

### Pull

Downloads outputs from the cloud store into your local MinIO.

```python
from isaura.manage import IsauraPull

with IsauraPull(
  model_id="eos8a4x",
  model_version="v1",
  bucket="isaura-public",
  input_csv="data/inputs.csv",
) as p:
  p.pull()
```

### Push

Uploads local public/private data to the cloud.

```python
from isaura.manage import IsauraPush

IsauraPush(model_id="eos8a4x", model_version="v1", bucket="myproject").push()
```

> Cloud credentials must be provided via environment variables (see `docs/CONFIGURATION.md`).

---

## Stats (inventory report)

Generate a JSON report of what's stored, including optional column counts.

```python
from isaura.manage import IsauraStat

st = IsauraStat(
  project_name="isaura-public",
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
* histogram of "how many models each molecule appears in"
* optional parquet column info (best effort)

<br>

---


# CLI reference

---

## Configuration

```bash
isaura configure                     # show current configuration
isaura configure --remote            # add or update remote/cloud credentials
isaura configure --update            # update a single credential interactively
isaura configure --show-secrets      # show all credential values unmasked
isaura configure --test-credentials  # test local and cloud connectivity
```

---

## Local services

```bash
isaura engine           # show Docker and MinIO status
isaura engine --start   # start local MinIO
isaura engine --stop    # stop local MinIO
```

---

## Projects

A **project** is a named MinIO bucket used as a staging area.

```bash
isaura info                  # list local projects (name, access, created, path)
isaura info --remote         # list remote projects

isaura create -pn myproject --access public    # create a new local project
isaura destroy -pn myproject                   # destroy entire project
isaura destroy -pn myproject -m eos8a4x -v v1  # destroy a specific model version
```

> `isaura-public` and `isaura-private` are reserved — they cannot be fully destroyed, but individual model versions inside them can be removed.

---

## Write

Store model outputs in a project bucket:

```bash
isaura write -i data/ersilia_output.csv -m eos8a4x -v v1 -pn myproject
```

| Flag | Description |
|---|---|
| `-i, --input` | Input CSV (required) |
| `-m, --model-id` | Ersilia model ID (required) |
| `-v, --version` | Model version, default `v1` |
| `-pn, --project-name` | Target bucket (required) |

> If the filename contains a model ID that doesn't match `--model-id`, isaura will error before writing anything.

---

## Read

Retrieve stored outputs for a set of inputs:

```bash
# explicit version
isaura read -i data/inputs.csv -m eos8a4x -v v1 -pn myproject -o data/outputs.csv

# auto-resolve to latest stored version
isaura read -i data/inputs.csv -m eos8a4x -pn myproject -o data/outputs.csv
```

| Flag | Description |
|---|---|
| `-i, --input` | Input CSV (required) |
| `-m, --model-id` | Ersilia model ID (required) |
| `-v, --version` | Model version — defaults to latest stored version if omitted |
| `-pn, --project-name` | Source bucket (required) |
| `-o, --output` | Output CSV path (optional — prints row count if omitted) |
| `-V, --verbose` | Show detailed internal logs |

---

## Inspect

Check which molecules from an input CSV are already cached:

```bash
isaura inspect --model_id eos8a4x -v v1 --access public -i data/inputs.csv -o reports/available.csv
```

| Flag | Description |
|---|---|
| `--model_id, -m` | Ersilia model ID (required) |
| `-v, --version` | Model version, default `v1` |
| `-pn, --project-name` | Restrict search to this bucket |
| `--access` | `public`, `private`, or `both` — ignored when `--project-name` is set |
| `-i, --input` | Input CSV to check (required) |
| `-o, --output` | Output CSV path (required) |
| `-r, --remote` | Query the remote store (default: local) |

---

## Catalog

List all models stored in a project:

```bash
isaura catalog -pn myproject           # local
isaura catalog -pn isaura-public -r    # remote (requires cloud credentials)
```

---

## Publishing to the cloud

**Step 1 — write to a local project:**

```bash
isaura write -i data/ersilia_output.csv -m eos8a4x -v v1 -pn myproject
```

**Step 2 — persist into the canonical bucket:**

```bash
isaura persist -m eos8a4x -v v1 -pn myproject
```

Routes molecules to `isaura-public` or `isaura-private` based on the project's access setting.

**Step 3 — push to cloud:**

```bash
isaura push -m eos8a4x -v v1 -pn isaura-public
```

> Cloud credentials must be configured first with `isaura configure --remote`.

---

## Pull

Download precalculations from the remote store into local MinIO:

```bash
# explicit version
isaura pull -i data/inputs.csv -m eos8a4x -v v1 -pn isaura-public

# auto-resolve to latest stored version
isaura pull -i data/inputs.csv -m eos8a4x -pn isaura-public
```

| Flag | Description |
|---|---|
| `-i, --input` | Input CSV (required) |
| `-m, --model-id` | Ersilia model ID (required) |
| `-v, --version` | Model version — defaults to latest stored version if omitted |
| `-pn, --project-name` | Remote bucket to pull from |
| `-V, --verbose` | Show detailed internal logs |

---

## Stats

Generate a JSON inventory of all models in a bucket:

```bash
isaura stats -pn isaura-public -o ./reports
isaura stats -pn myproject --remote -o ./reports
```

| Flag | Description |
|---|---|
| `-pn, --project-name` | Bucket to scan (required) |
| `-o, --output-dir` | Folder to write the JSON report (required) |
| `-r, --remote` | Scan the remote store (default: local) |
