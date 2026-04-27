import json
import subprocess

import requests
import yaml
from io import StringIO
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

from isaura.const import (
  GITHUB_CONTENT_URL, LOGP_BINS, METADATA_JSON, METADATA_YML,
  MKEYS, MW_BINS, _INT_KEYS, _KEYMAP,
)
from isaura.logging import logger


def _github_get(mid, file):
  """Fetch a raw file from the Ersilia model GitHub repository."""
  return requests.get(f"{GITHUB_CONTENT_URL}/{mid}/main/{file}")


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


def tranche_coordinates(smiles):
  """Compute the (row, col, mw, logp) tranche coordinates for a SMILES string.

  Bins molecular weight and LogP into grid coordinates used to assign a
  molecule to a specific Parquet tranche. Used as a fallback when a molecule
  is missing from the JSON index.

  Args:
      smiles: A valid SMILES string.

  Returns:
      Tuple of (row, col, mw, logp).

  Raises:
      ValueError: If the SMILES string cannot be parsed by RDKit.
  """
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    raise ValueError("Invalid SMILES")
  mw, logp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
  col = next((i + 1 for i, edge in enumerate(MW_BINS) if mw <= edge), len(MW_BINS) + 1)
  row = next((j + 1 for j, edge in enumerate(LOGP_BINS) if logp <= edge), len(LOGP_BINS) + 1)
  return (row, col, mw, logp)


def run_docker_compose(up=True):
  """Start or stop the MinIO/Milvus/NNS services via Docker Compose.

  Args:
      up: If True, runs `docker compose up -d`. If False, runs `docker compose down`.

  Returns:
      True on success, False on failure.
  """
  try:
    path = Path(__file__).parent / "configs" / "docker-compose.yml"
    cmd = ["docker", "compose", "-f", path, "up", "-d"] if up else ["docker", "compose", "-f", path, "down"]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    logger.info(f"Docker Compose {'started' if up else 'stopped'} successfully.")
    return True
  except subprocess.CalledProcessError as e:
    logger.error(f"Docker Compose failed: {e.stderr.strip()}")
  except FileNotFoundError:
    logger.error("Docker is not installed or not in PATH.")
  except Exception as e:
    logger.error(f"Unexpected error: {e}")
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
  console.print(Text("[New] Version 2.1.16", style="bold bright_black"), justify="center")
  print()
