import os, pandas as pd, numpy as np, json, math, os, psutil, pyarrow as pa, pyarrow.csv as pa_csv, requests, subprocess, sys, tempfile, time, yaml
import pyarrow.compute as pc
from io import StringIO
from collections import defaultdict
from loguru import logger
from typing import TypeVar
from rich.progress import (
  Progress,
  SpinnerColumn,
  TextColumn,
  BarColumn,
  TimeRemainingColumn,
  TimeElapsedColumn,
)
from rich.text import Text
from rich.table import Table
from rich.logging import RichHandler
from rich.progress import Progress
from rich.console import Console
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
from pathlib import Path

try:
  from dotenv import load_dotenv

  load_dotenv(override=False)
except Exception:
  pass
logger.remove()
console = Console()
logger.level("DEBUG", color="<cyan><bold>")
logger.level("INFO", color="<blue><bold>")
logger.level("WARNING", color="<white><bold><bg yellow>")
logger.level("ERROR", color="<white><bold><bg red>")
logger.level("CRITICAL", color="<white><bold><bg red>")
logger.level("SUCCESS", color="<black><bold><bg green>")
ACCESS_FILE = "access.json"
INDEX_FILE = "index.json"
MIN_NNS_RESULT_SIZE = 1000
MAX_ROWS = 2000000
BUILD_STATUS_TIMEOUT = 5
BUILD_START_TIMEOUT = 30
BUILD_POLL_INTERVAL = 0.5
BUILD_MAX_WAIT = 1.0
MW_BINS = [200, 500]
LOGP_BINS = [-1, 5]
GITHUB_ORG = "ersilia-os"
GITHUB_CONTENT_URL = f"https://raw.githubusercontent.com/{GITHUB_ORG}"
METADATA_JSON = "metadata.json"
METADATA_YML = "metadata.yml"
COLLECTION = os.getenv("COLLECTION", "eos3b5e")
MINIO_ENDPOINT = os.getenv("MINIO_ENDPOINT", "http://127.0.0.1:9000")
TIMEOUT = os.getenv("TIMEOUT", 3600)
MINIO_ENDPOINT_CLOUD = os.getenv("MINIO_ENDPOINT_CLOUD") or "http://83.48.73.209:8080"
NNS_ENDPOINT_BASE = os.getenv("NNS_ENDPOINT") or "http://127.0.0.1:8080"
MINIO_LOCAL_AK = os.getenv("MINIO_LOCAL_AK", "minioadmin123")
MINIO_LOCAL_SK = os.getenv("MINIO_LOCAL_SK", "minioadmin1234")
MINIO_CLOUD_AK = os.getenv("MINIO_CLOUD_AK", None)
MINIO_CLOUD_SK = os.getenv("MINIO_CLOUD_SK", None)
MINIO_PRIV_CLOUD_AK = os.getenv("MINIO_PRIV_CLOUD_AK", None)
MINIO_PRIV_CLOUD_SK = os.getenv("MINIO_PRIV_CLOUD_SK", None)
isaura_temp = os.path.join(Path.home(), "isaura", "isaura-temp")
if not os.path.exists(isaura_temp):
  os.makedirs(isaura_temp)
STORE_DIRECTORY = os.getenv("STORE_DIRECTORY", isaura_temp)
MAX_ROWS_PER_FILE = int(os.getenv("MAX_ROWS_PER_FILE", "100000"))
IMMUTABLE_CHUNK_COLS_THRESHOLD = int(os.getenv("IMMUTABLE_CHUNK_COLS_THRESHOLD", "128"))
WIDE_OUTPUT_DIM_THRESHOLD = int(os.getenv("WIDE_OUTPUT_DIM_THRESHOLD", "100"))
CHECKPOINT_EVERY = int(os.getenv("CHECKPOINT_EVERY", "50000"))
BLOOM_FILENAME = os.getenv("BLOOM_FILENAME", "bloom.pkl")
INPUT_C = ["input", "smiles"]
DEFAULT_BUCKET_NAME = os.getenv("DEFAULT_BUCKET_NAME", "isaura-public")
DEFAULT_PRIVATE_BUCKET_NAME = os.getenv("DEFAULT_PRIVATE_BUCKET_NAME", "isaura-private")
BATCH = int(os.getenv("BATCH", 10000))
FLUSH_EVERY = os.getenv("FLUSH_EVERY", 10000)
proc = psutil.Process(os.getpid())
MKEYS = [
  "Status",
  "Deployment",
  "Source",
  "Source Type",
  "Task",
  "Subtask",
  "Output",
  "Output Dimension",
  "Tag",
  "Biomedical Area",
  "Target Organism",
  "Publication Type",
  "Publication Year",
]
_INT_KEYS = {"Output Dimension", "Publication Year"}
_KEYMAP = {
  "Source Type": "SourceType",
  "Output Dimension": "OutputDimension",
  "Biomedical Area": "BiomedicalArea",
  "Target Organism": "TargetOrganism",
  "Publication Type": "PublicationType",
  "Publication Year": "PublicationYear",
}


def get(mid, file):
  return requests.get(f"{GITHUB_CONTENT_URL}/{mid}/main/{file}")


def get_base(mdi, ver):
  return f"{get_pref(mdi, ver)}/tranches"


def get_desc(pref, wanted):
  return f"Fetching hive partitions {pref} ({len(wanted)} inputs)"


def get_files_glob(bucket, base):
  return f"s3://{bucket}/{base}/*/chunk_*.parquet"


def get_keys(file, base):
  return f"{base}/{file}"


def get_idx_key(base):
  return get_keys(INDEX_FILE, base)


def get_acc_key(base):
  return get_keys(ACCESS_FILE, base)


def get_pref(mdi, ver):
  return f"{mdi}/{ver}"


def get_coll(mdi, ver):
  return f"{mdi}_{ver}"


def get_params(collection):
  return {"collection": collection, "batch": str(BATCH), "flush_every": str(FLUSH_EVERY)}


def get_header():
  return {"Content-Type": "text/plain"}


def hive_prefix(base):
  return f"{base}/data"


def make_temp(pref):
  return tempfile.mkdtemp(prefix=pref, dir=STORE_DIRECTORY)


def rss_mb():
  return proc.memory_info().rss / (1024 * 1024)


def log(msg):
  logger.info(f"[{time.strftime('%H:%M:%S')}] {msg} | RSS={rss_mb():.1f} MB")


def avail_mem():
  return int(psutil.virtual_memory().available)


def mem_gb_lim(ratio=0.8, floor_gb=1):
  return max(floor_gb, int(avail_mem() * ratio / 1024**3))


def cpu_cnt(ratio=0.6):
  return max(1, int(math.floor((os.cpu_count() or 1) * ratio)))


def list_parquet_keys(store, bucket, base):
  prefix = hive_prefix(base) + "/"
  t0 = time.perf_counter()
  try:
    keys = []
    for obj in store.list_keys(bucket, prefix):
      k = obj["Key"]
      if k.endswith(".parquet") and "/chunk_" in k:
        keys.append(k)
    keys.sort()
    uris = [f"s3://{bucket}/{k}" for k in keys]
    dt = time.perf_counter() - t0
    logger.info(f"[list_parquet_keys] listed {len(uris)} files in {dt:.3f}s bucket={bucket} prefix={prefix}")
    return uris if uris else get_files_glob(bucket, base)
  except Exception as e:
    logger.warning(f"[list_parquet_keys] boto3 listing failed ({e}), falling back to glob")
    return get_files_glob(bucket, base)


def pick_meta(d: dict) -> dict:
  def first(v):
    if v is None:
      return None
    if isinstance(v, list):
      return None if not v else v[0]
    return v

  out = {}
  for k in MKEYS:
    v = first(d.get(k))
    kk = _KEYMAP.get(k, k)
    if v is None:
      out[kk] = None
      continue
    if k in _INT_KEYS:
      try:
        out[kk] = int(v)
        continue
      except:
        pass
    out[kk] = str(v)
  return out


def output_dimension_from_metadata(meta):
  if not isinstance(meta, dict):
    return None
  value = meta.get("OutputDimension")
  if value is None:
    value = meta.get("Output Dimension")
  try:
    return int(value)
  except Exception:
    return None


def chunk_row_limit(output_dimension):
  if output_dimension is not None and int(output_dimension) >= WIDE_OUTPUT_DIM_THRESHOLD:
    return MAX_ROWS_PER_FILE
  return MAX_ROWS


def fetch_schema_from_github(model_id):
  try:
    r = get(model_id, METADATA_JSON)
    data = json.load(StringIO(r.text))
  except:
    r = get(model_id, METADATA_YML)
    data = yaml.safe_load(StringIO(r.text))
  return pick_meta(data)


def run_docker_compose(up=True):
  try:
    path = Path(__file__).parent / "configs" / "docker-compose.yml"
    cmd = ["docker", "compose", "-f", path, "up", "-d"] if up else ["docker", "compose", "-f", path, "down"]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    logger.info(f"Docker Compose {('started' if up else 'stopped')} successfully.")
    logger.debug(result.stdout.strip())
    return True
  except subprocess.CalledProcessError as e:
    logger.error(f"Docker Compose failed: {e.stderr.strip()}")
  except FileNotFoundError:
    logger.error("Docker is not installed or not in PATH.")
  except Exception as e:
    logger.error(f"Unexpected error: {e}")
  return False


def show_figlet():
  path = Path(__file__).parent / "assets" / "figlet.txt"
  text = Path(path).read_text(encoding="utf-8")
  start_color = (0, 255, 255)
  end_color = (255, 0, 255)
  content = "".join(text)
  gradient = Text()
  for i, ch in enumerate(content):
    ratio = i / max(1, len(content) - 1)
    r = int(start_color[0] + (end_color[0] - start_color[0]) * ratio)
    g = int(start_color[1] + (end_color[1] - start_color[1]) * ratio)
    b = int(start_color[2] + (end_color[2] - start_color[2]) * ratio)
    gradient.append(ch, style=f"rgb({r},{g},{b})")
  print()
  console.print(gradient, justify="center")
  console.print(Text(f"[🚀] [New] Version 2.1.16", style="bold bright_black"), justify="center")
  print()


def split_csv(df):
  paths = []
  output_dir = make_temp("isaura_push_")
  for bucket in [DEFAULT_BUCKET_NAME, DEFAULT_PRIVATE_BUCKET_NAME]:
    if bucket in df["bucket"].unique():
      path = os.path.join(output_dir, f"{bucket.replace('-', '_')}.csv")
      df[df["bucket"] == bucket].to_csv(path, index=False)
      paths.append(str(path))
  return paths


def query(conn, header, wanted, file_glob, columns="*", tmpdir="/tmp", preserve_order=False):
  chunks = list(
    query_batched(
      conn, header, wanted, file_glob, columns=columns, tmpdir=tmpdir, preserve_order=preserve_order
    )
  )
  if not chunks:
    return pd.DataFrame()
  return pd.concat(chunks, ignore_index=True)


def query_batched(
  conn, header, wanted, file_glob, batch_size=10000, columns="*", tmpdir="/tmp", preserve_order=False
):
  if not wanted:
    logger.debug("[query_batched] empty wanted list — nothing to query")
    return
  mem_lim = mem_gb_lim()
  threads = cpu_cnt()
  try:
    conn.execute(f"SET memory_limit='{mem_lim}GB'")
    conn.execute(f"SET temp_directory='{tmpdir}'")
    conn.execute("PRAGMA enable_object_cache")
    conn.execute(f"SET threads TO {threads}")
  except Exception:
    pass
  wanted_list = list(wanted)
  if isinstance(file_glob, list):
    escaped = ", ".join((f"'{u}'" for u in file_glob))
    src_expr = f"read_parquet([{escaped}])"
    src_desc = f"{len(file_glob)} explicit files"
  else:
    src_expr = f"read_parquet('{file_glob}')"
    src_desc = f"glob={file_glob}"
  logger.info(
    f"[query_batched] inputs={len(wanted_list)} arrow_batch_size={batch_size} mem_limit={mem_lim}GB threads={threads} source={src_desc}"
  )
  if preserve_order:
    sql = f"\n          WITH p AS (\n            SELECT {columns}\n            FROM {src_expr}\n            WHERE {header} IN (SELECT {header} FROM __wanted_inputs)\n          )\n          SELECT p.*\n          FROM p\n          JOIN __wanted_inputs w\n            ON p.{header} = w.{header}\n          ORDER BY w.__o\n      "
    wdf = pd.DataFrame({header: wanted_list, "__o": np.arange(len(wanted_list), dtype=np.int64)})
  else:
    sql = f"\n          SELECT {columns}\n          FROM {src_expr}\n          WHERE {header} IN (SELECT {header} FROM __wanted_inputs)\n      "
    wdf = pd.DataFrame({header: wanted_list})
  conn.register("__wanted_inputs", wdf)
  total_rows = 0
  n_arrow_batches = 0
  sql_t0 = time.perf_counter()
  try:
    rb_reader = conn.execute(sql).fetch_record_batch(batch_size)
    for batch in rb_reader:
      n_arrow_batches += 1
      df = batch.to_pandas(split_blocks=True, self_destruct=True)
      total_rows += len(df)
      logger.debug(
        f"[query_batched] arrow_batch #{n_arrow_batches} rows={len(df)} total_rows={total_rows} rss={rss_mb():.0f}MB"
      )
      yield df
  finally:
    conn.unregister("__wanted_inputs")
  sql_dt = time.perf_counter() - sql_t0
  logger.info(
    f"[query_batched] finished total_rows={total_rows} arrow_batches={n_arrow_batches} elapsed={sql_dt:.2f}s rss={rss_mb():.0f}MB"
  )


STREAM_PARQUET_THRESHOLD = int(os.getenv("STREAM_PARQUET_THRESHOLD", "200000"))


def stream_parquet_filtered(store, bucket, prefix, wanted, header="input", batch_size=10000):
  import pyarrow.parquet as pq

  wanted_set = set(wanted) if not isinstance(wanted, set) else wanted
  remaining = len(wanted_set)
  wanted_arr = pa.array(list(wanted_set))
  tmpdir = make_temp("isaura_stream_")
  keys = []
  for obj in store.list_keys(bucket, prefix):
    k = obj["Key"]
    if k.endswith(".parquet") and "/chunk_" in k:
      keys.append(k)
  keys.sort()
  logger.info(
    f"[stream] starting: {len(keys)} parquet files, {remaining} wanted inputs, batch_size={batch_size} bucket={bucket} prefix={prefix}"
  )
  total_rg_scanned = 0
  total_rg_hit = 0
  total_rows_yielded = 0
  n_chunks_yielded = 0
  for ki, key in enumerate(keys):
    if remaining <= 0:
      break
    local = os.path.join(tmpdir, f"s_{ki}.parquet")
    dl_t0 = time.perf_counter()
    try:
      store.download_file(bucket, key, local)
    except Exception as e:
      logger.warning(f"[stream] skip {key}: {e}")
      continue
    dl_dt = time.perf_counter() - dl_t0
    file_size = os.path.getsize(local)
    try:
      pf = pq.ParquetFile(local)
      n_rg = pf.metadata.num_row_groups
      logger.debug(
        f"[stream] file {ki + 1}/{len(keys)} key={key} size={file_size / (1024 * 1024):.1f}MB row_groups={n_rg} download={dl_dt:.2f}s"
      )
      file_hits = 0
      for rg_idx in range(n_rg):
        if remaining <= 0:
          break
        total_rg_scanned += 1
        key_col = pf.read_row_group(rg_idx, columns=[header]).column(header)
        mask = pc.is_in(key_col, wanted_arr)
        if pc.any(mask).as_py() is not True:
          del key_col, mask
          logger.debug(f"[stream] file {ki + 1} rg {rg_idx}/{n_rg} — no matches, skipped")
          continue
        total_rg_hit += 1
        filtered = pf.read_row_group(rg_idx).filter(mask)
        del key_col, mask
        if filtered.num_rows == 0:
          del filtered
          continue
        file_hits += filtered.num_rows
        logger.debug(
          f"[stream] file {ki + 1} rg {rg_idx}/{n_rg} — matched {filtered.num_rows} rows, cols={filtered.num_columns}"
        )
        for start in range(0, filtered.num_rows, batch_size):
          chunk = filtered.slice(start, batch_size).to_pandas(split_blocks=True, self_destruct=True)
          matched = set(chunk[header].astype(str).str.strip())
          remaining -= len(matched & wanted_set)
          wanted_set -= matched
          n_chunks_yielded += 1
          total_rows_yielded += len(chunk)
          yield chunk
          del chunk
        del filtered
      if file_hits > 0:
        logger.debug(f"[stream] file {ki + 1}/{len(keys)} yielded {file_hits} rows, remaining={remaining}")
      del pf
    except Exception as e:
      logger.warning(f"[stream] error reading {key}: {e}")
    finally:
      try:
        os.remove(local)
      except Exception:
        pass
    if (ki + 1) % 10 == 0:
      logger.info(
        f"[stream] progress files={ki + 1}/{len(keys)} remaining={remaining} rows_yielded={total_rows_yielded} rg_scanned={total_rg_scanned} rg_hit={total_rg_hit} rss={rss_mb():.0f}MB"
      )
  try:
    os.rmdir(tmpdir)
  except Exception:
    pass
  logger.info(
    f"[stream] done files={len(keys)} rg_scanned={total_rg_scanned} rg_hit={total_rg_hit} chunks_yielded={n_chunks_yielded} rows_yielded={total_rows_yielded} unmatched={remaining} rss={rss_mb():.0f}MB"
  )


def stream_parquet_filtered_ordered(store, bucket, prefix, wanted, header="input", batch_size=10000):
  wanted_list = [str(v).strip() for v in wanted if str(v).strip()]
  if not wanted_list:
    logger.info("[stream-ordered] empty wanted list")
    return
  unresolved = set(wanted_list)
  resolved = {}
  schema_cols = None
  next_idx = 0
  emitted = 0

  def missing_row(key):
    base = {header: key}
    if schema_cols:
      for col in schema_cols:
        if col != header:
          base[col] = None
    return base

  def make_frame(rows):
    if not rows:
      return pd.DataFrame(columns=schema_cols or [header])
    if schema_cols:
      return pd.DataFrame.from_records(rows, columns=schema_cols)
    return pd.DataFrame.from_records(rows)

  def flush_ready():
    nonlocal next_idx, emitted
    rows = []
    while next_idx < len(wanted_list):
      key = wanted_list[next_idx]
      row = resolved.get(key)
      if row is None:
        if key in unresolved:
          break
        rows.append(missing_row(key))
      else:
        rows.append(dict(row))
      next_idx += 1
      emitted += 1
      if len(rows) >= batch_size:
        yield make_frame(rows)
        rows = []
    if rows:
      yield make_frame(rows)

  for chunk in stream_parquet_filtered(
    store, bucket, prefix, set(unresolved), header=header, batch_size=batch_size
  ):
    if chunk is None or chunk.empty:
      continue
    if schema_cols is None:
      schema_cols = list(chunk.columns)
    for row in chunk.to_dict("records"):
      key = str(row.get(header) or "").strip()
      if key and key not in resolved:
        resolved[key] = row
        unresolved.discard(key)
    for ready in flush_ready():
      logger.debug(
        f"[stream-ordered] yield rows={len(ready)} emitted={emitted}/{len(wanted_list)} unresolved={len(unresolved)}"
      )
      yield ready
    if not unresolved:
      break
  unresolved_before_fill = len(unresolved)
  unresolved.clear()
  if next_idx < len(wanted_list):
    tail = []
    for key in wanted_list[next_idx:]:
      row = resolved.get(key)
      if row is None:
        tail.append(missing_row(key))
      else:
        tail.append(dict(row))
      if len(tail) >= batch_size:
        yield make_frame(tail)
        tail = []
    if tail:
      yield make_frame(tail)
    if unresolved_before_fill:
      logger.warning(
        f"[stream-ordered] missing={unresolved_before_fill} rows after full scan; emitted blank placeholders"
      )


def group_inputs(wanted, index, bloom=None, force=False):
  logger.debug(f"Checking {len(wanted)} inputs in the {len(index)} index!")
  try:
    if bloom:
      miss = [s for s in wanted if not bloom.seen(s)]
    g = defaultdict(set)
    if miss and (not force):
      logger.error(
        f"inputs not indexed: {miss[:5]}{('...' if len(miss) > 5 else '')} total_missing={len(miss)}"
      )
      sys.exit(1)
    if not force:
      for s in wanted:
        r, c = index[s]
        g[int(r), int(c)].add(s)
      return g
    if miss and force:
      for s in miss:
        r, c, _, _ = tranche_coordinates(s)
        g[int(r), int(c)].add(s)
      return g
  except Exception as e:
    logger.error(e)
    return None


def line_gen(df, chunksize=100000):
  for start in range(0, len(df), chunksize):
    chunk = df.iloc[start : start + chunksize]
    for s in chunk["input"].astype(str):
      yield (s + "\n").encode("utf-8")
    log(f"sent chunk rows={len(chunk):,}")
  log("finished streaming")


def build_index_status(collection):
  r = requests.get(
    f"{NNS_ENDPOINT_BASE}/build_index", params={"collection": collection}, timeout=BUILD_STATUS_TIMEOUT
  )
  r.raise_for_status()
  return r.json()


def start_build_index(collection, nlist=None, rebuild=False, wait=False):
  params = {"collection": collection}
  if nlist is not None:
    params["nlist"] = str(nlist)
  if rebuild:
    params["rebuild"] = "1"
  if wait:
    params["wait"] = "1"
  r = requests.post(f"{NNS_ENDPOINT_BASE}/build_index", params=params, timeout=BUILD_START_TIMEOUT)
  r.raise_for_status()
  return r.json()


def ensure_index_ready(collection, max_wait_s=BUILD_MAX_WAIT):
  st = build_index_status(collection)
  if not st.get("exists", False):
    try:
      start_build_index(collection, wait=False)
    except requests.RequestException as e:
      return (False, {"error": f"start_build_index_failed: {e}", "status": st})
    st = build_index_status(collection)
  if st.get("is_failed", False):
    return (False, {"error": "index_failed", "status": st})
  if st.get("is_finished", False) and (not st.get("is_building", False)):
    return (True, {"status": st})
  t_end = time.time() + max_wait_s
  while time.time() < t_end:
    time.sleep(BUILD_POLL_INTERVAL)
    st = build_index_status(collection)
    if st.get("is_failed", False):
      return (False, {"error": "index_failed", "status": st})
    if st.get("is_finished", False) and (not st.get("is_building", False)):
      return (True, {"status": st})
  return (False, {"error": "index_building", "status": st})


def post_apprx(df, collection, nlist=None):
  t0 = time.time()
  try:
    with requests.Session() as s:
      log("start streaming to nns api. This process sometimes appears to be slow. Please have some patience!")
      resp = s.post(
        f"{NNS_ENDPOINT_BASE}/insert",
        params=get_params(collection),
        data=line_gen(df),
        headers=get_header(),
        timeout=None,
      )
    dt = (time.time() - t0) * 1000
    log(f"done total_ms={dt:.0f}. Body: {resp.text[:1000]}")
    try:
      start_build_index(collection, nlist=nlist, rebuild=False, wait=False)
      st = build_index_status(collection)
      logger.info(
        f"index build triggered collection={collection} exists={st.get('exists')} building={st.get('is_building')} finished={st.get('is_finished')} progress={st.get('progress_pct')}"
      )
    except requests.RequestException as e:
      logger.error(f"index build trigger failed collection={collection}: {e}")
    return resp
  except requests.RequestException as e:
    logger.error(f"approx NN search failed. The NNS server container may not be sarted!: {e}")
    return None


def get_apprx(inputs, collection, fallback_search=None):
  try:
    logger.info(f"Sending {len(inputs)} inputs for ANN search server to get top 1 similar compounds")
    ok, meta = ensure_index_ready(collection, max_wait_s=BUILD_MAX_WAIT)
    if not ok:
      st = meta.get("status", {})
      err = meta.get("error", "index_not_ready")
      msg = f"ANN index not ready collection={collection} error={err} exists={st.get('exists')} building={st.get('is_building')} finished={st.get('is_finished')} failed={st.get('is_failed')} progress={st.get('progress_pct')}"
      logger.error(msg)
      if callable(fallback_search):
        logger.info("Please use conventional/exact search")
        return fallback_search(inputs, collection)
      return []
    r = requests.post(
      f"{NNS_ENDPOINT_BASE}/search", json={"collection": collection, "smiles": inputs}, timeout=TIMEOUT
    )
    r.raise_for_status()
    payload = r.json() or {}
    results = payload.get("results", [])
    return [x.get("match") for x in results if "match" in x]
  except requests.RequestException as e:
    logger.error(f"approx NN search failed. The NNS server container may not be sarted!: {e}")
    if callable(fallback_search):
      logger.info("Falling back to conventional/exact search")
      return fallback_search(inputs, collection)
    return []


def write_access_file(existed, data, access, dir):
  try:
    if data:
      m = [{"input": i, "access": access} for i in data]
      if existed:
        m = m + existed
      with open(dir, "w") as f:
        json.dump(m, f, indent=2)
  except Exception as e:
    logger.error(e)


def tranche_coordinates(smiles):
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    raise ValueError("Invalid SMILES")
  mw, logp = (Descriptors.MolWt(mol), Crippen.MolLogP(mol))
  for i, edge in enumerate(MW_BINS):
    if mw <= edge:
      col = i + 1
      break
  else:
    col = len(MW_BINS) + 1
  for j, edge in enumerate(LOGP_BINS):
    if logp <= edge:
      row = j + 1
      break
  else:
    row = len(LOGP_BINS) + 1
  return (row, col, mw, logp)


def make_table(title, cols, rows):
  t = Table(title=title)
  for c in cols:
    t.add_column(c["name"], justify=c.get("justify", "left"), style=c.get("style", ""))
  for r in rows:
    t.add_row(*[str(r.get(c["key"], "")) for c in cols])
  return t


inspect_table = [
  {"key": "model", "name": "model/version", "justify": "left", "style": "bold"},
  {"key": "entries", "name": "entry count", "justify": "right"},
  {"key": "chunks", "name": "chunks", "justify": "right"},
]
T = TypeVar("T")


def track_write_progress(rows, total=None, description="Writing rows", console=None):
  if total is None:
    progress = Progress(
      SpinnerColumn(),
      TextColumn("[progress.description]{task.description}"),
      TextColumn("{task.completed} rows"),
      TimeElapsedColumn(),
      console=console,
      transient=True,
    )
  else:
    progress = Progress(
      SpinnerColumn(),
      TextColumn("[progress.description]{task.description}"),
      BarColumn(),
      TextColumn("{task.completed}/{task.total}"),
      TimeElapsedColumn(),
      TimeRemainingColumn(),
      console=console,
      transient=True,
    )
  with progress:
    task_id = progress.add_task(description, total=total)
    for row in rows:
      yield row
      progress.advance(task_id, 1)


def spinner(message, fn, *args, **kwargs):
  c = Console()
  with c.status(message, spinner="dots"):
    result = fn(*args, **kwargs)
  return result


class StreamingCsvSink:
  def __init__(self, path):
    self._path = path
    self._fp = None
    self._header_written = False
    self.rows_written = 0

  def __enter__(self):
    self._fp = open(self._path, "wb")
    return self

  def __exit__(self, et, ev, tb):
    self.close()

  def close(self):
    if self._fp is not None and (not self._fp.closed):
      self._fp.close()

  def write_table(self, table):
    if table is None or table.num_rows == 0:
      return
    opts = pa_csv.WriteOptions(include_header=not self._header_written)
    pa_csv.write_csv(table, self._fp, write_options=opts)
    self.rows_written += table.num_rows
    self._header_written = True

  def write_batch(self, batch):
    self.write_table(pa.Table.from_batches([batch]))

  def write_df(self, df):
    if df is None or df.empty:
      return
    self.write_table(pa.Table.from_pandas(df, preserve_index=False))


class Logger:
  def __init__(self):
    self.logger = logger
    self._console = None
    self._file = None
    self._log_to_console()

  def _log_to_console(self):
    if self._console is None:
      rich_handler = RichHandler(
        rich_tracebacks=True, markup=True, log_time_format="%H:%M:%S", show_path=False
      )
      self._rich_console = rich_handler.console
      self._console = self.logger.add(rich_handler, format="{message}", colorize=True)

  @property
  def console(self):
    return self._rich_console

  def _unlog_from_console(self):
    if self._console is not None:
      try:
        self.logger.remove(self._console)
      except Exception:
        pass
      self._console = None

  def set_verbosity(self, verbose):
    if verbose:
      self._log_to_console()
    else:
      self._unlog_from_console()

  def debug(self, text):
    self.logger.debug(text)

  def info(self, text):
    self.logger.info(text)

  def warning(self, text):
    self.logger.warning(text)

  def error(self, text):
    self.logger.error(text)

  def critical(self, text):
    self.logger.critical(text)

  def success(self, text):
    self.logger.success(text)


logger = Logger()
