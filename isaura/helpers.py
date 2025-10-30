import os, pandas as pd, numpy as np, json, math, os, psutil, requests, subprocess, sys, tempfile, time
from contextlib import contextmanager
from collections import defaultdict
from loguru import logger
from typing import TypeVar, Optional
from rich.progress import (
  Progress,
  SpinnerColumn,
  TextColumn,
  BarColumn,
  TimeRemainingColumn,
  TimeElapsedColumn,
  ProgressColumn,
  Task,
)
from rich.text import Text
from rich.table import Table
from rich.logging import RichHandler
from rich.progress import Progress
from rich.console import Console
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
from pathlib import Path

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
MAX_ROWS = 2_000_000

MW_BINS = [200, 500]
LOGP_BINS = [-1, 5]

COLLECTION = os.getenv("COLLECTION", "eos3b5e")
MINIO_ENDPOINT = os.getenv("MINIO_ENDPOINT", "http://127.0.0.1:9000")
TIMEOUT = os.getenv("TIMEOUT", 3600)
MINIO_ENDPOINT_CLOUD = os.getenv("MINIO_ENDPOINT_CLOUD") or "http://83.48.73.209:8080"
NNS_ENDPOINT_BASE = os.getenv("NNS_ENDPOINT") or "http://127.0.0.1:8080/"
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
CHECKPOINT_EVERY = int(os.getenv("CHECKPOINT_EVERY", "50000"))
BLOOM_FILENAME = os.getenv("BLOOM_FILENAME", "bloom.pkl")
INPUT_C = ["input", "smiles"]
DEFAULT_BUCKET_NAME = os.getenv("DEFAULT_BUCKET_NAME", "isaura-public")
DEFAULT_PRIVATE_BUCKET_NAME = os.getenv("DEFAULT_PRIVATE_BUCKET_NAME", "isaura-private")
BATCH = int(os.getenv("BATCH", 10_000))
FLUSH_EVERY = os.getenv("FLUSH_EVERY", 10_000)
proc = psutil.Process(os.getpid())


# fmt: off
def get_base(mdi, ver): return f"{get_pref(mdi, ver)}/tranches"
def get_desc(pref, wanted): return f"Fetching hive partitions {pref} ({len(wanted)} inputs)"
def get_files_glob(bucket, base): return f"s3://{bucket}/{base}/*/chunk_*.parquet"
def get_keys(file, base): return f"{base}/{file}"
def get_idx_key(base): return get_keys(INDEX_FILE, base)
def get_acc_key(base): return get_keys(ACCESS_FILE,base)
def get_pref(mdi, ver): return f"{mdi}/{ver}"
def get_coll(mdi, ver): return f"{mdi}_{ver}"
def get_params(collection): return {"collection": collection, "batch": str(BATCH), "flush_every": str(FLUSH_EVERY)}
def get_header(): return {"Content-Type": "text/plain"}
def hive_prefix(base): return f"{base}/data"
def make_temp(pref): return tempfile.mkdtemp(prefix=pref, dir=STORE_DIRECTORY)
def rss_mb(): return proc.memory_info().rss / (1024*1024)
def log(msg): logger.info(f"[{time.strftime('%H:%M:%S')}] {msg} | RSS={rss_mb():.1f} MB")
def avail_mem(): return int(psutil.virtual_memory().available)
def mem_gb_lim(ratio=0.8, floor_gb=1): return max(floor_gb, int(avail_mem() * ratio / (1024**3)))
def cpu_cnt(ratio=0.6): return max(1, int(math.floor((os.cpu_count() or 1) * ratio)))
# fmt: on


def run_docker_compose(up=True):
  try:
    path = Path(__file__).parent / "configs" / "docker-compose.yml"
    cmd = ["docker", "compose", "-f", path, "up", "-d"] if up else ["docker", "compose", "-f", path, "down"]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    logger.info(f"Docker Compose {'started' if up else 'stopped'} successfully.")
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
  path = Path(__file__).parent.parent / "assets" / "figlet.txt"
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
  console.print(Text(f"Version 2.0.1", style="bold bright_black"), justify="center")
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


def query(conn, header, wanted, file_glob, columns="*", tmpdir="/tmp"):
  if not wanted:
    return pd.DataFrame()
  try:
    conn.execute(f"SET memory_limit='{mem_gb_lim()}GB'")
    conn.execute(f"SET temp_directory='{tmpdir}'")
    conn.execute("PRAGMA enable_object_cache")
    conn.execute(f"SET threads TO {cpu_cnt()}")
  except Exception:
    pass
  wanted_list = list(wanted)
  order = np.arange(len(wanted_list), dtype=np.int64)
  wdf = pd.DataFrame({header: wanted_list, "__o": order})
  conn.register("wanted_inputs", wdf)
  sql = f"""
        WITH p AS (
          SELECT {columns}
          FROM read_parquet('{file_glob}')
          WHERE {header} IN (SELECT {header} FROM wanted_inputs)
        )
        SELECT p.*
        FROM p
        JOIN wanted_inputs w
          ON p.{header} = w.{header}
        ORDER BY w.__o
    """
  try:
    out = conn.execute(sql).fetchdf()
  finally:
    conn.unregister("wanted_inputs")
  return out


def group_inputs(wanted, index, force=False):
  logger.debug(f"Checking {len(wanted)} inputs in the index!")
  try:
    miss = [s for s in wanted if s not in index]
    g = defaultdict(set)
    if miss and not force:
      logger.error(
        f"inputs not indexed: {miss[:5]}{'...' if len(miss) > 5 else ''} total_missing={len(miss)}"
      )
      sys.exit(1)
    if not force:
      for s in wanted:
        r, c = index[s]
        g[(int(r), int(c))].add(s)
      return g
    if miss and force:
      for s in miss:
        r, c, _, _ = tranche_coordinates(s)
        g[(int(r), int(c))].add(s)
      return g
  except Exception as e:
    logger.error(e)
    return None


def line_gen(df, chunksize=100_000):
  for start in range(0, len(df), chunksize):
    chunk = df.iloc[start : start + chunksize]
    for s in chunk["input"].astype(str):
      yield (s + "\n").encode("utf-8")
    log(f"sent chunk rows={len(chunk):,}")
  log("finished streaming")


def post_apprx(df, colelction):
  t0 = time.time()
  try:
    with requests.Session() as s:
      log("start streaming to nns api. This process sometimes appears to be slow. Please have some patience!")
      resp = s.post(
        f"{NNS_ENDPOINT_BASE}/insert",
        params=get_params(colelction),
        data=line_gen(df),
        headers=get_header(),
        timeout=None,
      )
    dt = (time.time() - t0) * 1000
    log(f"done total_ms={dt:.0f}. Body: {resp.text[:1000]}")
  except requests.RequestException as e:
    logger.error(f"approx NN search failed. The NNS server container may not be sarted!: {e}")
    return []


def get_apprx(inputs, collection):
  try:
    logger.info(f"Sending {len(inputs)} inputs for ANN search server to get top 1 similar compounds")
    r = requests.post(
      f"{NNS_ENDPOINT_BASE}/search", json={"collection": collection, "smiles": inputs}, timeout=TIMEOUT
    )
    r.raise_for_status()
    return [x["input"] for x in r.json().get("results", [])]
  except requests.RequestException as e:
    logger.error(f"approx NN search failed. The NNS server container may not be sarted!: {e}")
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
  mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)

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

  return row, col, mw, logp


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


class RowCountColumn(ProgressColumn):
  def render(self, task: Task) -> Text:
    return Text(f"{int(task.completed):,} rows")


class RowSpeedColumn(ProgressColumn):
  def render(self, task: Task) -> Text:
    return Text(f"{task.speed:,.0f} rows/s" if task.speed else "– rows/s", style="dim")


def make_row_progress(transient: bool = True) -> Progress:
  return Progress(
    SpinnerColumn(),
    TextColumn("[bold cyan]{task.fields[desc]}[/]"),
    BarColumn(),
    RowCountColumn(),
    RowSpeedColumn(),
    TimeElapsedColumn(),
    TimeRemainingColumn(),
    transient=transient,
  )


@contextmanager
def progress(desc: str, total_rows: Optional[int] = None, transient: bool = True):
  with make_row_progress(transient=transient) as prog:
    task_id = prog.add_task("rows", total=total_rows, desc=desc)
    yield prog, task_id


def spinner(message, fn, *args, **kwargs):
  c = Console()
  with c.status(message, spinner="dots"):
    result = fn(*args, **kwargs)
  return result


# fmt: off
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
      self._console = self.logger.add(
        rich_handler,
        format="{message}",
        colorize=True,
      )

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

  def debug(self, text): self.logger.debug(text)
  def info(self, text): self.logger.info(text)
  def warning(self, text): self.logger.warning(text)
  def error(self, text): self.logger.error(text)
  def critical(self, text): self.logger.critical(text)
  def success(self, text): self.logger.success(text)


logger = Logger()
