import csv
import subprocess
import tempfile

import pandas as pd
import pytest

from isaura.helpers import DEFAULT_BUCKET_NAME as PUB



def cli(*args):
  r = subprocess.run(["isaura", *args], capture_output=True, text=True)
  assert r.returncode == 0, f"isaura {' '.join(args)}\n--- stderr ---\n{r.stderr.strip()}"
  return r


def load_inputs(path):
  with open(path, newline="") as f:
    return {(row.get("input") or row.get("smiles", "")).strip() for row in csv.DictReader(f)} - {""}


def input_col(df):
  return "input" if "input" in df.columns else "smiles"




def test_1_write(cfg, state):
  cli(
    "write",
    "-m",
    cfg["model"],
    "-v",
    cfg["version"],
    "-pn",
    cfg["src"],
    "--access",
    "public",
    "-i",
    cfg["input"],
  )
  state["wrote"] = True


def test_2_read(cfg, state):
  if not state.get("wrote"):
    pytest.skip("write did not pass")

  out = tempfile.mktemp(suffix=".csv")
  state["read_out"] = out

  cli(
    "read",
    "-m",
    cfg["model"],
    "-v",
    cfg["version"],
    "-pn",
    cfg["src"],
    "-i",
    cfg["input"],
    "-o",
    out,
  )

  df = pd.read_csv(out)
  assert len(df) > 0, "read returned 0 rows"

  wanted = load_inputs(cfg["input"])
  got = set(df[input_col(df)].astype(str))
  assert wanted.issubset(got), f"inputs missing from output: {wanted - got}"

  state["read_rows"] = len(df)


def test_3_read_batched(cfg, state):
  if not state.get("wrote"):
    pytest.skip("write did not pass")

  from isaura.manage import IsauraReader

  reader = IsauraReader(
    model_id=cfg["model"],
    model_version=cfg["version"],
    bucket=cfg["src"],
    input_csv=cfg["input"],
    approximate=False,
  )

  chunks = list(reader.read_batched(batch_size=3))
  assert chunks, "read_batched yielded no chunks"

  total = sum(len(c) for c in chunks)
  assert total == state.get("read_rows", total), f"batched total {total} != read total {state['read_rows']}"




def test_4_inspect(cfg, state):
  if not state.get("wrote"):
    pytest.skip("write did not pass")

  out = tempfile.mktemp(suffix=".csv")
  cli(
    "inspect",
    "-m",
    cfg["model"],
    "-v",
    cfg["version"],
    "-pn",
    cfg["src"],
    "--access",
    "public",
    "-o",
    out,
  )

  df = pd.read_csv(out)
  assert len(df) > 0, "inspect returned empty result"



def test_5_copy(cfg, state):
  if not state.get("wrote"):
    pytest.skip("write did not pass")

  cli(
    "copy",
    "-m",
    cfg["model"],
    "-v",
    cfg["version"],
    "-pn",
    cfg["src"],
  )
  state["copied"] = True


def test_6_move(cfg, state):
  if not state.get("wrote"):
    pytest.skip("write did not pass")

  cli(
    "move",
    "-m",
    cfg["model"],
    "-v",
    cfg["version"],
    "-pn",
    cfg["src"],
  )
  state["moved"] = True



def test_7_remove_model(cfg, state):
  if not (state.get("copied") or state.get("moved")):
    pytest.skip("copy/move did not pass")

  cli(
    "remove",
    "-m",
    cfg["model"],
    "-v",
    cfg["version"],
    "-pn",
    PUB,
    "--yes",
  )



def test_8_remove_project(cfg, state):
  cli("remove", "-pn", cfg["src"], "--yes")
