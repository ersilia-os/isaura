# Troubleshooting

This guide focuses on **local setup failures**, **stale artifacts**, and **slow retrieval** when using **Ersilia + Isaura**.

> ⚠️ Warning: Some steps delete local data (buckets, volumes, files). Proceed only if you’re okay resetting local state.

---

## Quick symptoms → likely cause

### Setup fails early / won’t complete cleanly
- Docker services not running / wrong compose setup
- MinIO or Milvus ports already in use
- Stale volumes from a previous run

### NN / ANN retrieval fails
- NNS container not started
- NNS endpoint misconfigured (`NNS_ENDPOINT`)
- ANN index not built yet (or build failed)

### Old artifacts keep getting reused unexpectedly
- Objects still exist in MinIO buckets
- Local state directory not cleared (e.g. `~/isaura`)
- Index/metadata left behind (`index.json`, `access.json`, bloom)

### Second run is still slow
- ANN index not ready
- Query is scanning too much (missing `input` column match)
- DuckDB is memory constrained / temp directory issues

---

## Health checks (do these first)

### 1) Confirm services are running

```bash
docker ps
````

You should see containers for (names may vary):

* MinIO
* Milvus
* NNS service

### 2) Confirm ports are reachable

MinIO (API): `http://localhost:9000`
MinIO (Console): `http://localhost:9001`
Milvus API: `http://localhost:8080`

Quick checks:

```bash
curl -s http://localhost:9000/minio/health/live
curl -s http://localhost:9000/minio/health/ready
curl -s http://localhost:8080/health || true
```

### 3) Verify environment configuration

```bash
env | egrep "MINIO_|NNS_|DEFAULT_BUCKET|CHECKPOINT|BLOOM|STORE_DIRECTORY" || true
```

If you rely on a `.env` file, confirm it is being loaded (or exported vars are set).

---

## Reset workflow (recommended if things are weird)

### Step 0) Optional: snapshot logs first

```bash
docker ps
docker logs <container_name_or_id> --tail 200
```

---

### 1) Delete stale objects in MinIO buckets

Buckets commonly involved:

* `isaura-public`
* `isaura-private`
* `ersilia` (if used)

Use either the UI or `mc`.

#### Option A: MinIO Console (UI)

1. Open `http://localhost:9001`
2. Open each bucket
3. Remove the model folder you’re troubleshooting (e.g. `eosxxxx/v1/...`)

#### Option B: MinIO client (`mc`)

```bash
mc alias set local http://localhost:9000 minioadmin123 minioadmin1234
mc ls local
mc rm -r --force local/isaura-public/eosxxxx/
mc rm -r --force local/isaura-private/eosxxxx/
```

> Tip: delete only the affected model/version prefix to avoid wiping everything.

---

### 2) Reset local persisted state

Some deployments persist local state under:

```bash
rm -rf ~/isaura
```

If your setup uses a different path, check `STORE_DIRECTORY` or your compose volumes.

---

### 3) Restart the engine services

Fast path:

```bash
isaura engine --start
```

If you need a full container refresh:

```bash
docker compose -f isaura/configs/docker-compose.yml down -v
docker compose -f isaura/configs/docker-compose.yml up -d
```

> The `-v` wipes volumes (use only if you really want a clean slate).

---

## Milvus / ANN troubleshooting

### Inspect Milvus collection metadata

```bash
curl -X POST "http://localhost:8080/info?collection=eosxxxx_v1" | jq
```

Replace `eosxxxx_v1` with your actual collection name.

If the collection exists but looks stale, rebuild your artifacts and try again.

---

### ANN index not ready / stuck building

ANN search depends on the NN service index being built.

#### 1) Check NN service status (if exposed)

If your NN service has a build endpoint:

```bash
curl -s "http://127.0.0.1:8080/build_index?collection=eosxxxx_v1" | jq
```

You may see fields like:

* `exists`
* `is_building`
* `is_finished`
* `is_failed`
* `progress_pct`

#### 2) Rebuild / trigger index build

If supported:

```bash
curl -X POST "http://127.0.0.1:8080/build_index?collection=eosxxxx_v1&rebuild=1" | jq
```

#### 3) If ANN keeps failing

* Confirm `NNS_ENDPOINT` is correct
* Ensure the NNS container is running
* Use exact search until the index is ready

---

## Read/write issues

### “inputs not indexed” error on read

This happens when you request inputs that are missing from `index.json`.

Fix:

* Ensure you wrote results for those inputs first (`isaura write`)
* Or validate inputs using inspect:

  ```bash
  isaura inspect --model_id eosxxxx -v v1 --access public -i data/inputs.csv -o reports/inspect_report.csv
  ```

### Wrong input header

Reader logic prefers `input`, otherwise it falls back to `smiles`.

Make sure your input CSV has one of:

* `input`
* `smiles`

---

## Performance debugging

### Second run still slow

Expected behavior:

* 1st run: slower (initialization, index warm-up)
* 2nd run: fast (cached / indexed retrieval)

If the second run is still slow:

1. Ensure you’re using **exact matches** (header correct).
2. Ensure the dataset is not tiny if using ANN (ANN may be disabled for small sets).
3. Check machine resources:

   ```bash
   free -h || true
   df -h
   ```
4. Verify DuckDB temp dir is writable and has space (`STORE_DIRECTORY` / temp path).
5. If ANN is enabled, ensure index is ready (see ANN section above).

---

## Ersilia integration validation

### Minimal validation loop

Serve:

```bash
ersilia serve eosxxxx -rs -ws -a public
```

Run inference (twice):

```bash
ersilia run -i input.csv -o output.csv -b 10000
ersilia run -i input.csv -o output.csv -b 10000
```

Expected:

* First run: slower
* Second run: significantly faster

---

## When opening an issue (copy/paste checklist)

Include:

* Model + version:

  * `eosxxxx`
  * `v1`
  * collection name (if relevant): `eosxxxx_v1`
* Environment:

  ```bash
  env | egrep "MINIO_|NNS_|DEFAULT_BUCKET|CHECKPOINT|BLOOM|STORE_DIRECTORY" || true
  ```
* Running containers:

  ```bash
  docker ps
  ```
* Logs (last 200 lines):

  ```bash
  docker logs <minio_container> --tail 200
  docker logs <milvus_container> --tail 200
  docker logs <nns_container> --tail 200
  ```
* Milvus info:

  ```bash
  curl -X POST "http://localhost:8080/info?collection=eosxxxx_v1" | jq
  ```
* CLI outputs:

  * `isaura engine --start`
  * `isaura write ...`
  * `isaura read ...`

---

