import os

import pandas as pd
import pytest

duckdb = pytest.importorskip("duckdb")
pytest.importorskip("pyarrow")

from isaura.query import chunked_query_batched, query, query_batched


@pytest.fixture
def store(tmp_path):
  # Two local parquet chunks, 300 rows; no MinIO needed.
  paths = []
  for c in range(2):
    p = os.path.join(str(tmp_path), f"chunk_{c}.parquet")
    ids = range(c * 150, (c + 1) * 150)
    pd.DataFrame({"input": [f"s{i}" for i in ids], "x": [float(i) for i in ids]}).to_parquet(p)
    paths.append(p)
  return paths, str(tmp_path)


def test_query_batched_returns_matching_rows(store):
  # Issue #31: duckdb 2 removed fetch_record_batch, so every narrow read failed.
  paths, tmp = store
  wanted = ["s3", "s149", "s150", "s299", "missing"]
  out = pd.concat(query_batched(duckdb.connect(), "input", wanted, paths, batch_size=2, tmpdir=tmp))
  assert sorted(out["input"]) == sorted(["s3", "s149", "s150", "s299"])
  assert dict(zip(out["input"], out["x"]))["s150"] == 150.0


def test_query_batched_preserve_order(store):
  paths, tmp = store
  wanted = ["s299", "s0", "s150", "s7"]
  out = query(duckdb.connect(), "input", wanted, paths, tmpdir=tmp, preserve_order=True)
  assert list(out["input"]) == wanted


@pytest.mark.parametrize("preserve_order", [False, True])
def test_chunked_query_batched(store, preserve_order):
  paths, tmp = store
  wanted = [f"s{i}" for i in range(0, 300, 7)]
  out = pd.concat(
    chunked_query_batched(
      duckdb.connect(), "input", wanted, paths, slice_size=5, batch_size=3, tmpdir=tmp, preserve_order=preserve_order
    )
  )
  got = list(out["input"])
  assert got == wanted if preserve_order else sorted(got) == sorted(wanted)


def test_connection_reusable_after_early_stop(store):
  # duckdb 2 keeps a connection busy while a result is open; stopping early must release it.
  paths, tmp = store
  conn = duckdb.connect()
  gen = query_batched(conn, "input", [f"s{i}" for i in range(300)], paths, batch_size=10, tmpdir=tmp)
  next(gen)
  gen.close()
  out = query(conn, "input", ["s1"], paths, tmpdir=tmp)
  assert list(out["input"]) == ["s1"]


def test_cleanup_failure_does_not_mask_query_error(store, monkeypatch):
  # Issue #31: a failing unregister in `finally` replaced the real error with a misleading one.
  paths, tmp = store
  conn = duckdb.connect()

  def boom(*args, **kwargs):
    raise RuntimeError("real query error")

  class Conn:
    def __getattr__(self, name):
      return getattr(conn, name)

    def unregister(self, name):
      raise duckdb.Error("cleanup failure")

  monkeypatch.setattr("isaura.query._arrow_reader", boom)
  with pytest.raises(RuntimeError, match="real query error"):
    list(query_batched(Conn(), "input", ["s1"], paths, tmpdir=tmp))
  with pytest.raises(RuntimeError, match="real query error"):
    list(chunked_query_batched(Conn(), "input", ["s1"], paths, tmpdir=tmp))
