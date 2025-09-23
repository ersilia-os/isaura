import csv, gzip, io, os, requests, time
from contextlib import contextmanager
from io import StringIO
from loguru import logger
from typing import TypeVar, Optional, Iterable, Dict, List, Any, Iterator, Tuple
from rich.progress import (
  Progress,
  SpinnerColumn,
  TextColumn,
  BarColumn,
  TimeRemainingColumn,
  DownloadColumn,
  TransferSpeedColumn,
)
from rich.logging import RichHandler
from rich.progress import Progress
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

ROTATION = "10 MB"

logger.remove()

logger.level("DEBUG", color="<cyan><bold>")
logger.level("INFO", color="<blue><bold>")
logger.level("WARNING", color="<white><bold><bg yellow>")
logger.level("ERROR", color="<white><bold><bg red>")
logger.level("CRITICAL", color="<white><bold><bg red>")
logger.level("SUCCESS", color="<black><bold><bg green>")


def get_schema(model_id):
  st = time.perf_counter()
  try:
    response = requests.get(
      f"{GITHUB_CONTENT_URL}/{model_id}/main/{PREDEFINED_COLUMN_FILE}"
    )
  except requests.RequestException:
    logger.warning("Couldn't fetch column name from github!")
    return None

  csv_data = StringIO(response.text)
  reader = csv.DictReader(csv_data)

  if "name" not in reader.fieldnames:
    logger.warning("Couldn't fetch column name from github. Column name not found.")
    return None

  if "type" not in reader.fieldnames:
    logger.warning("Couldn't fetch data type from github. Column name not found.")
    return None

  rows = list(reader)
  col_name = [row["name"] for row in rows if row["name"]]
  col_dtype = [row["type"] for row in rows if row["type"]]
  shape = len(col_dtype)
  if len(col_name) == 0 and len(col_dtype) == 0:
    return None
  et = time.perf_counter()
  logger.info(f"Column metadata fetched in {et - st:.2} seconds!")
  return col_name, col_dtype, shape


def resolve_dtype(dtype):
  unique_dtypes = list(set(dtype))
  dtype = "float" if "float" in unique_dtypes else unique_dtypes[0]
  if dtype == "integer":
    return int
  if dtype == "float":
    return float
  return str


def is_gz_shard(shard: Dict[str, Any]) -> bool:
  key = shard.get("key", "")
  url = shard.get("url", "")
  return key.endswith(".gz") or url.endswith(".gz")


def get_size(client, url: str, default: int = 0) -> int:
  try:
    print(url)
    headers = client.head(url)
    return int(headers.get("Content-Length", default) or default)
  except Exception:
    return int(default or 0)


def collect_gz_shards(
  shards: Iterable[Dict[str, Any]], client
) -> Tuple[List[Dict[str, Any]], int]:
  gz = [s for s in shards if is_gz_shard(s)]
  total = 0
  for s in gz:
    size = int(s.get("size", 0) or 0)
    url = s.get("url", "")
    total += get_size(client, url, size)
  return gz, total


def iter_csv(text: str) -> Iterator[List[str]]:
  reader = csv.reader(io.StringIO(text))
  for row in reader:
    if row:
      yield row


def parse_row(row: List[str]) -> Optional[Tuple[str, List[float]]]:
  if len(row) < 3:
    return None
  smi = row[1].strip()
  try:
    values = [float(r) for r in row[2:]]
  except Exception:
    return None
  return smi, values


def flush(smiles: List[str], values: List[List[float]]) -> None:
  if smiles:
    smiles.clear()
    values.clear()


def process_shard(
  shard: Dict[str, Any], client, batch_size: int, progress=None, task=None
) -> None:
  buf = io.BytesIO()
  for chunk in client.download_stream(shard["url"]):
    if not chunk:
      continue
    buf.write(chunk)
    if progress and task:
      progress.update(task, advance=len(chunk))
  text = gzip.decompress(buf.getvalue()).decode("utf-8", errors="replace")

  smiles, values, fps = [], [], []
  for row in iter_csv(text):
    rec = parse_row(row)
    if rec is None:
      continue
    smi, vals = rec
    smiles.append(smi)
    values.append(vals)
    if len(smiles) >= batch_size:
      pass
    # we do something here


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
def download_progress(
  desc: str, total_bytes: Optional[int] = None, transient: bool = True
):
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


# Enums
class JobStatus:
  PENDING = "PENDING"
  RUNNING = "RUNNING"
  FAILED = "FAILED"
  SUCCEEDED = "SUCCEEDED"


REDIS_EXPIRATION = 3600 * 24 * 7
REDIS_PORT = 6379
REDIS_CONTAINER_NAME = "redis"
REDIS_IMAGE = "redis:latest"
REDIS_HOST = "127.0.0.1"
DEFAULT_API_NAME = "run"
S3_BUCKET_URL = "https://ersilia-models.s3.eu-central-1.amazonaws.com"
S3_BUCKET_URL_ZIP = "https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com"
INFERENCE_STORE_API_URL = (
  "https://5x2fkcjtei.execute-api.eu-central-1.amazonaws.com/dev/precalculations"
)
API_BASE = "https://hov95ejni7.execute-api.eu-central-1.amazonaws.com/dev/predict"
GITHUB_ORG = "ersilia-os"
GITHUB_CONTENT_URL = f"https://raw.githubusercontent.com/{GITHUB_ORG}"
GITHUB_ERSILIA_REPO = "ersilia"
PREDEFINED_COLUMN_FILE = "model/framework/columns/run_columns.csv"
TIMEOUT = 3600


MW_BINS = [200, 250, 300, 325, 350, 375, 400, 425, 450, 500]
LOGP_BINS = [-1, 0, 1, 2, 2.5, 3, 3.5, 4, 4.5, 5]


def tranche_coordinates(smiles):
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    raise ValueError("Invalid SMILES")
  mw = Descriptors.MolWt(mol)
  logp = Crippen.MolLogP(mol)
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


MINIO_ENDPOINT = os.getenv("MINIO_ENDPOINT", "http://127.0.0.1:9000")
MINIO_ACCESS_KEY = os.getenv("MINIO_ACCESS_KEY", "minioadmin")
MINIO_SECRET_KEY = os.getenv("MINIO_SECRET_KEY", "minioadmin")
STORE_DIRECTORY = os.getenv("STORE_DIRECTORY", ".")
MAX_ROWS_PER_FILE = int(os.getenv("MAX_ROWS_PER_FILE", "100000"))
CHECKPOINT_EVERY = int(os.getenv("CHECKPOINT_EVERY", "50000"))
BLOOM_FILENAME = os.getenv("BLOOM_FILENAME", "bloom.pkl")
DEFAULT_BUCKET_NAME = os.getenv("DEFAULT_BUCKET_NAME", "isaura-public")
