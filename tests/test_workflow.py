import csv
import json
import subprocess
import tempfile
import pandas as pd
import pytest
from click.testing import CliRunner
from isaura.helpers import DEFAULT_BUCKET_NAME as PUB
from pandas.testing import assert_frame_equal


def cli(*args):
  r = subprocess.run(["isaura", *args], capture_output=True, text=True)
  assert r.returncode == 0, f"isaura {' '.join(args)}\n--- stderr ---\n{r.stderr.strip()}"
  return r


def load_inputs(path):
  with open(path, newline="") as f:
    return {(row.get("input") or row.get("smiles", "")).strip() for row in csv.DictReader(f)} - {""}


def input_col(df):
  return "input" if "input" in df.columns else "smiles"


def normalized(df):
  ic = input_col(df)
  return df.sort_values(ic).reset_index(drop=True)


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
  cli("read", "-m", cfg["model"], "-v", cfg["version"], "-pn", cfg["src"], "-i", cfg["input"], "-o", out)
  df = pd.read_csv(out)
  original = pd.read_csv(cfg["input"])
  assert len(df) > 0, "read returned 0 rows"
  wanted = load_inputs(cfg["input"])
  got = set(df[input_col(df)].astype(str))
  assert wanted.issubset(got), f"inputs missing from output: {wanted - got}"
  cols = [c for c in original.columns if c in df.columns]
  assert cols, "read output does not share columns with input"
  assert_frame_equal(normalized(df[cols]), normalized(original[cols]), check_dtype=False, check_like=False)
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
  cli("inspect", "-m", cfg["model"], "-v", cfg["version"], "-pn", cfg["src"], "--access", "public", "-o", out)
  df = pd.read_csv(out)
  assert len(df) > 0, "inspect returned empty result"


def test_5_copy(cfg, state):
  if not state.get("wrote"):
    pytest.skip("write did not pass")
  cli("copy", "-m", cfg["model"], "-v", cfg["version"], "-pn", cfg["src"])
  state["copied"] = True


def test_6_move(cfg, state):
  if not state.get("wrote"):
    pytest.skip("write did not pass")
  cli("move", "-m", cfg["model"], "-v", cfg["version"], "-pn", cfg["src"])
  state["moved"] = True


def test_7_remove_model(cfg, state):
  if not (state.get("copied") or state.get("moved")):
    pytest.skip("copy/move did not pass")
  cli("remove", "-m", cfg["model"], "-v", cfg["version"], "-pn", PUB, "--yes")


def test_8_remove_project(cfg, state):
  cli("remove", "-pn", cfg["src"], "--yes")


def test_9_stats_writes_model_sizes(monkeypatch, tmp_path):
  from isaura.cli import cli as isaura_cli

  class FakeInspect:
    def __init__(self, *args, **kwargs):
      self.cloud = kwargs.get("cloud", False)
      self.heavy_index = kwargs.get("heavy_index", False)

    def buckets(self):
      return ["isaura-test"]

    def iter_models(self, bucket):
      yield "eos3b5e", "v1"

    def load_index(self, bucket, model_id, model_version):
      return {"A": [1, 1], "B": [1, 1]}

    def iter_object_meta(self, bucket, prefix=""):
      return iter([
        {"Key": f"{prefix}/tranches/data/chunk_1.parquet", "Size": 1073741824},
        {"Key": f"{prefix}/tranches/index.json", "Size": 512},
      ])

    def parquet_columns(self, bucket, parquet_key):
      return ["input", "mol_weight"]

  monkeypatch.setattr("isaura.manage.IsauraInspect", FakeInspect)
  monkeypatch.setattr("isaura.manage.fetch_schema_from_github", lambda mid: {"id": mid})

  out_dir = tmp_path / "stats"
  runner = CliRunner()
  result = runner.invoke(
    isaura_cli,
    ["stats", "-pn", "isaura-test", "--access", "public", "-o", str(out_dir), "-d", "."],
  )

  assert result.exit_code == 0, result.output
  files = list(out_dir.glob("isaura_stats_*.json"))
  assert len(files) == 1
  data = json.loads(files[0].read_text(encoding="utf-8"))
  assert data["models_total"] == 1
  assert data["models"][0]["total_bytes"] == 1073742336
  assert data["models"][0]["total_gb"] == 1.000000
