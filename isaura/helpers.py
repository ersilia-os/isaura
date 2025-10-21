import csv, json, os, psutil, requests, tempfile, time
import pandas as pd
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
  DownloadColumn,
  TransferSpeedColumn,
)
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
MAX_ROWS = 100_000

MW_BINS = [200, 250, 300, 325, 350, 375, 400, 425, 450, 500]
LOGP_BINS = [-1, 0, 1, 2, 2.5, 3, 3.5, 4, 4.5, 5]

COLLECTION = os.getenv("COLLECTION", "eos3b5e")
MINIO_ENDPOINT = os.getenv("MINIO_ENDPOINT") or "http://127.0.0.1:9000"
TIMEOUT = os.getenv("TIMEOUT", 3600)
MINIO_ENDPOINT_CLOUD = os.getenv("MINIO_ENDPOINT_CLOUD") or "http://83.48.73.209:8080"
NNS_ENDPOINT_BASE = os.getenv("NNS_ENDPOINT") or "http://127.0.0.1:8080/"
MINIO_ACCESS_KEY = os.getenv("MINIO_ACCESS_KEY", "minioadmin")
MINIO_SECRET_KEY = os.getenv("MINIO_SECRET_KEY", "minioadmin")
isaura_temp = os.path.join(Path.home(), "eos", "isaura-temp")
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
def get_files(bucket, r, c, base): return f"s3://{bucket}/{hive_prefix(r, c, base)}/chunk_*.parquet"
def get_keys(file, base): return f"{base}/{file}"
def get_idx_key(base): return get_keys(INDEX_FILE, base)
def get_acc_key(base): return get_keys(ACCESS_FILE,base)
def get_pref(mdi, ver): return f"{mdi}/{ver}"
def get_coll(mdi, ver): return f"{mdi}_{ver}"
def get_params(collection): return {"collection": collection, "batch": str(BATCH), "flush_every": str(FLUSH_EVERY)}
def get_header(): return {"Content-Type": "text/plain"}
def hive_prefix(r, c, base): return f"{base}/row={r}/col={c}"
def make_temp(pref): return tempfile.mkdtemp(prefix=pref, dir=STORE_DIRECTORY)
def query(files): return f"SELECT * FROM read_parquet('{files}', hive_partitioning=1)"
def rss_mb(): return proc.memory_info().rss / (1024*1024)
def log(msg): logger.info(f"[{time.strftime('%H:%M:%S')}] {msg} | RSS={rss_mb():.1f} MB")
# fmt: on


def filter_out(out, objc, wanted, header):
  return (
    out[out[header].isin(wanted)]
    .assign(__o=lambda d: d[header].map({s: i for i, s in enumerate(wanted)}))
    .sort_values(objc)
    .drop(columns=[objc, "row", "col"], errors="ignore")
    .reset_index(drop=True)
  )


def group_inputs(wanted, index, force=False):
  try:
    miss = [s for s in wanted if s not in index]
    print(len(miss), len(index))
    g = defaultdict(set)
    if miss and not force:
      raise RuntimeError(
        f"inputs not indexed: {miss[:5]}{'...' if len(miss) > 5 else ''} total_missing={len(miss)}"
      )
    if not force:
      print("Or This happeded")
      for s in wanted:
        r, c = index[s]
        g[(int(r), int(c))].add(s)
      return g
    if miss and force:
      print("This happeded")
      for s in miss:
        r, c, _, _ = tranche_coordinates(s)
        g[(int(r), int(c))].add(s)
      return g
  except Exception as e:
    print(e)


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
    logger.error(f"approx NN search failed: {e}")
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
    logger.error(f"approx NN search failed: {e}")
    return []


def write_access_file(data, access, dir):
  with open(data, "r") as f:
    data = csv.DictReader(f)
    input = [(d.get("input") or d.get("smiles")).strip() for d in data]
    m = [{"input": i, "access": access} for i in input]
  with open(dir, "w") as f:
    json.dump(m, f, indent=2)


def tranche_coordinates(smiles):
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    raise ValueError("Invalid SMILES")
  mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
  for i, edge in enumerate(MW_BINS):
    if mw <= edge:
      col = i
      break
  else:
    col = len(MW_BINS)
  for j, edge in enumerate(LOGP_BINS):
    if logp <= edge:
      row = j
      break
  else:
    row = len(LOGP_BINS)
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
  {"key": "tranches", "name": "tranches", "justify": "right"},
  {"key": "chunks", "name": "chunks", "justify": "right"},
]

T = TypeVar("T")


def make_download_progress(transient: bool = True) -> Progress:
  return Progress(
    SpinnerColumn(),
    TextColumn("[bold cyan]{task.fields[desc]}[/]"),
    BarColumn(),
    DownloadColumn(binary_units=True),
    TransferSpeedColumn(),
    TimeRemainingColumn(),
    transient=transient,
  )


@contextmanager
def progress(desc: str, total_bytes: Optional[int] = None, transient: bool = True):
  with make_download_progress(transient=transient) as progress:
    task_id = progress.add_task("download", total=total_bytes or 0, desc=desc)
    yield progress, task_id


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
