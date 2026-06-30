import csv
import json
import subprocess

import requests
import yaml
from io import StringIO
from pathlib import Path

from isaura.const import (
  DEFAULT_BUCKET_NAME, DEFAULT_PRIVATE_BUCKET_NAME,
  GITHUB_CONTENT_URL, METADATA_JSON, METADATA_YML,
  MINIO_ENDPOINT, MINIO_LOCAL_AK, MINIO_LOCAL_SK,
  MKEYS, RUN_COLUMNS_FILE, _INT_KEYS, _KEYMAP,
)
from isaura.logging import logger


def _github_get(mid, file):
  """Fetch a raw file from the Ersilia model GitHub repository."""
  return requests.get(f"{GITHUB_CONTENT_URL}/{mid}/main/{file}")


def fetch_run_columns(model_id):
  """Fetch and parse a model's run_columns.csv from GitHub.

  run_columns.csv is the authoritative declaration of a model's output columns
  (header: name,type,direction,description). Returns an ordered dict mapping each
  declared output column name to its lower-cased declared type (e.g. "float",
  "integer", "string"), or None if the file is missing / unreadable so the caller
  can log loudly and fall back to inference.

  Args:
      model_id: Ersilia model identifier (e.g. "eos4u6p").

  Returns:
      dict[str, str] of {column_name: declared_type}, or None if unavailable.
  """
  try:
    r = _github_get(model_id, RUN_COLUMNS_FILE)
  except Exception as e:
    logger.warning(f"[run_columns] fetch failed for {model_id}: {e}")
    return None
  if r.status_code != 200:
    logger.warning(f"[run_columns] {model_id}: HTTP {r.status_code} for {RUN_COLUMNS_FILE}")
    return None
  out = {}
  try:
    reader = csv.DictReader(StringIO(r.text))
    if not reader.fieldnames or "name" not in reader.fieldnames or "type" not in reader.fieldnames:
      logger.warning(f"[run_columns] {model_id}: unexpected header {reader.fieldnames}")
      return None
    for row in reader:
      name = (row.get("name") or "").strip()
      typ = (row.get("type") or "").strip().lower()
      if name:
        out[name] = typ
  except Exception as e:
    logger.warning(f"[run_columns] {model_id}: parse error: {e}")
    return None
  if not out:
    logger.warning(f"[run_columns] {model_id}: no columns parsed from {RUN_COLUMNS_FILE}")
    return None
  logger.debug(f"[run_columns] {model_id}: {len(out)} declared columns")
  return out


def pick_meta(d):
  """Extract and normalise a standard set of metadata fields from a raw model metadata dict.

  Handles both JSON and YAML source formats. List values are collapsed to their
  first element. Integer fields are cast to int where possible.

  Args:
      d: Raw metadata dict from GitHub (JSON or YAML parsed).

  Returns:
      Dict with normalised keys defined in MKEYS / _KEYMAP.
  """
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
      except Exception:
        pass
    out[kk] = str(v)
  return out


def output_dimension_from_metadata(meta):
  """Return the OutputDimension integer from a metadata dict, or None if absent or unparseable."""
  if not isinstance(meta, dict):
    return None
  value = meta.get("OutputDimension")
  if value is None:
    value = meta.get("Output Dimension")
  try:
    return int(value)
  except Exception:
    return None


def fetch_schema_from_github(model_id):
  """Fetch and normalise model metadata from GitHub, trying JSON then YAML.

  Args:
      model_id: Ersilia model identifier (e.g. "eos1234").

  Returns:
      Normalised metadata dict from pick_meta().
  """
  try:
    r = _github_get(model_id, METADATA_JSON)
    data = json.load(StringIO(r.text))
  except Exception:
    r = _github_get(model_id, METADATA_YML)
    data = yaml.safe_load(StringIO(r.text))
  return pick_meta(data)


def write_access_file(existed, data, access, dir):
  """Write or merge an access metadata JSON file mapping inputs to their access level.

  Args:
      existed: Existing list of access records (or None if file is new).
      data: List of new input SMILES strings to add.
      access: Access level string ("public" or "private").
      dir: Local file path to write the JSON to.
  """
  try:
    if data:
      m = [{"input": i, "access": access} for i in data]
      if existed:
        m = m + existed
      with open(dir, "w") as f:
        json.dump(m, f, indent=2)
  except Exception as e:
    logger.error(e)


def docker_is_installed() -> bool:
  """Return True if the docker CLI is available on PATH."""
  try:
    subprocess.run(["docker", "--version"], capture_output=True, check=True, timeout=5)
    return True
  except (FileNotFoundError, subprocess.CalledProcessError):
    return False


def docker_is_running() -> bool:
  """Return True if the Docker daemon is active and reachable."""
  try:
    result = subprocess.run(["docker", "info"], capture_output=True, timeout=5)
    return result.returncode == 0
  except Exception:
    return False


def get_engine_status() -> dict:
  """Return the status of Docker and the isaura containers.

  Returns:
      Dict with keys "docker" (bool) and "containers" (list of {name, status}).
  """
  running = docker_is_installed() and docker_is_running()
  containers = []
  if running:
    for name in ["minio"]:
      try:
        result = subprocess.run(
          ["docker", "ps", "-a", "--filter", f"name=^{name}$", "--format", "{{.Status}}"],
          capture_output=True, text=True, timeout=5,
        )
        status = result.stdout.strip() or "not created"
      except Exception:
        status = "unknown"
      containers.append({"name": name, "status": status})
  return {"docker": running, "containers": containers}


def ensure_default_buckets(timeout: int = 30) -> None:
  """Wait for local MinIO to be ready, then create the default buckets if they don't exist."""
  import time
  import boto3
  from botocore.config import Config

  url = f"{MINIO_ENDPOINT.rstrip('/')}/minio/health/ready"
  for _ in range(timeout):
    try:
      if requests.get(url, timeout=2).status_code == 200:
        break
    except Exception:
      pass
    time.sleep(1)
  else:
    logger.warning("MinIO did not become ready in time — skipping default bucket creation.")
    return

  client = boto3.client(
    "s3",
    endpoint_url=MINIO_ENDPOINT,
    aws_access_key_id=MINIO_LOCAL_AK,
    aws_secret_access_key=MINIO_LOCAL_SK,
    region_name="us-east-1",
    config=Config(signature_version="s3v4", s3={"addressing_style": "path"}),
  )
  for bucket in [DEFAULT_BUCKET_NAME, DEFAULT_PRIVATE_BUCKET_NAME]:
    try:
      client.head_bucket(Bucket=bucket)
    except Exception:
      client.create_bucket(Bucket=bucket)
      logger.info(f"Created default bucket: {bucket}")


def run_docker_compose(up=True) -> bool:
  """Start or stop the MinIO service via Docker Compose.

  Args:
      up: If True, runs `docker compose up -d`. If False, runs `docker compose down`.

  Returns:
      True on success, False on failure.
  """
  if not docker_is_installed():
    logger.error("Docker is not installed. Please install Docker Desktop: https://docs.docker.com/get-docker/")
    return False
  if not docker_is_running():
    logger.error("Docker is not running. Please open Docker Desktop and try again.")
    return False

  path = Path(__file__).parent / "configs" / "docker-compose.yml"

  if up:
    try:
      check = subprocess.run(
        ["docker", "image", "inspect", "minio/minio:latest"],
        capture_output=True, timeout=5,
      )
      if check.returncode != 0:
        logger.info("MinIO image not found locally — it will be downloaded now. This may take a few minutes on first run.")
    except Exception:
      pass
    logger.info("Spinning up Docker containers...")
    cmd = ["docker", "compose", "-f", str(path), "up", "-d"]
  else:
    logger.info("Stopping Docker containers...")
    cmd = ["docker", "compose", "-f", str(path), "down"]

  result = subprocess.run(cmd)
  if result.returncode == 0:
    logger.info(f"Containers {'started' if up else 'stopped'} successfully.")
    if up:
      ensure_default_buckets()
    return True
  logger.error(f"Docker Compose {'up' if up else 'down'} failed (exit code {result.returncode}).")
  return False


def show_figlet():
  """Print the isaura ASCII art banner with a cyan-to-magenta gradient to the console."""
  from rich.text import Text
  from isaura.logging import console

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
  try:
    from importlib.metadata import version
    _version = version("isaura")
  except Exception:
    _version = "unknown"
  console.print(Text(f"Version {_version}", style="bold bright_black"), justify="center")
  print()
