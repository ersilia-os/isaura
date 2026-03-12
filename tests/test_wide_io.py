import os
import pandas as pd
import pytest

pq = pytest.importorskip("pyarrow.parquet")

from isaura.base import TrancheState, _BaseTransfer
from isaura.manage import IsauraReader


class LocalStore:
  def __init__(self, root):
    self.root = root

  def _path(self, bucket, key):
    return os.path.join(self.root, bucket, key)

  def upload_file(self, local, bucket, key, extra_args=None):
    target = self._path(bucket, key)
    os.makedirs(os.path.dirname(target), exist_ok=True)
    with open(local, "rb") as src, open(target, "wb") as dst:
      dst.write(src.read())

  def download_file(self, bucket, key, local):
    source = self._path(bucket, key)
    os.makedirs(os.path.dirname(local), exist_ok=True)
    with open(source, "rb") as src, open(local, "wb") as dst:
      dst.write(src.read())

  def list_keys(self, bucket, prefix):
    root = self._path(bucket, prefix)
    if not os.path.exists(root):
      return
    for current, _, files in os.walk(root):
      for name in files:
        rel = os.path.relpath(os.path.join(current, name), os.path.join(self.root, bucket))
        yield {"Key": rel}


def test_tranche_flush_handles_wide_rows_without_full_frame(tmp_path):
  store = LocalStore(str(tmp_path))
  tranche = TrancheState(store, "bucket", "model/v1/tranches", str(tmp_path), 4)
  rows = []
  for i in range(6):
    row: dict[str, object] = {"input": f"s{i}"}
    for j in range(128):
      row[f"f{j}"] = float(i + j)
    rows.append(row)
  cols = list(rows[0].keys())

  tranche.flush(rows[:3], cols)
  tranche.flush(rows[3:], cols)

  key = os.path.join(str(tmp_path), "bucket", "model/v1/tranches/data/chunk_1.parquet")
  out = pq.read_table(key).to_pandas()
  assert len(out) == 3
  assert list(out.columns) == cols

  key2 = os.path.join(str(tmp_path), "bucket", "model/v1/tranches/data/chunk_2.parquet")
  out2 = pq.read_table(key2).to_pandas()
  assert len(out2) == 3
  assert list(out2.columns) == cols


def test_tranche_flush_df_handles_wide_rows_with_immutable_chunks(tmp_path):
  store = LocalStore(str(tmp_path))
  tranche = TrancheState(store, "bucket", "model/v1/tranches", str(tmp_path), 4)
  rows = []
  for i in range(6):
    row: dict[str, object] = {"input": f"s{i}"}
    for j in range(128):
      row[f"f{j}"] = float(i + j)
    rows.append(row)
  df = pd.DataFrame(rows)
  cols = list(df.columns)

  tranche.flush_df(df.iloc[:3], cols)
  tranche.flush_df(df.iloc[3:], cols)

  key1 = os.path.join(str(tmp_path), "bucket", "model/v1/tranches/data/chunk_1.parquet")
  key2 = os.path.join(str(tmp_path), "bucket", "model/v1/tranches/data/chunk_2.parquet")
  out1 = pq.read_table(key1).to_pandas()
  out2 = pq.read_table(key2).to_pandas()

  assert len(out1) == 3
  assert len(out2) == 3
  assert list(out1.columns) == cols
  assert list(out2.columns) == cols


def test_reader_read_batched_to_csv_skips_index_load(monkeypatch, tmp_path):
  class FakeStore:
    def __init__(self, *args, **kwargs):
      pass

  class FakeDuck:
    def __init__(self, *args, **kwargs):
      self.con = object()

  class FakeBloom:
    def __init__(self, *args, **kwargs):
      pass

    def seen(self, v):
      return True

  monkeypatch.setattr("isaura.manage.MinioStore", FakeStore)
  monkeypatch.setattr("isaura.manage.DuckDBMinio", FakeDuck)
  monkeypatch.setattr("isaura.manage.BloomIndex", FakeBloom)

  called = {"query_batched": 0}

  def fake_query_batched(*args, **kwargs):
    called["query_batched"] += 1
    yield pd.DataFrame([
      {"input": "a", "x": 1.5},
      {"input": "b", "x": 2.5},
    ])

  monkeypatch.setattr("isaura.manage.query_batched", fake_query_batched)

  reader = IsauraReader(model_id="m", model_version="v1", bucket="bucket", input_csv="unused.csv", approximate=False)
  reader._load_index = lambda: (_ for _ in ()).throw(AssertionError("index should not load"))

  out = tmp_path / "out.csv"
  df = pd.DataFrame([{"input": "a"}, {"input": "b"}])
  chunks = list(reader.read_batched(output_csv=str(out), df=df))

  assert called["query_batched"] == 1
  assert list(pd.read_csv(out)["input"]) == ["a", "b"]
  assert list(chunks[0]["input"]) == ["a", "b"]


def test_reader_read_streaming_mode_skips_duckdb(monkeypatch):
  class FakeStore:
    def __init__(self, *args, **kwargs):
      pass

  class FakeDuck:
    def __init__(self, *args, **kwargs):
      self.con = object()

  class FakeBloom:
    def __init__(self, *args, **kwargs):
      pass

    def seen(self, v):
      return True

  monkeypatch.setattr("isaura.manage.MinioStore", FakeStore)
  monkeypatch.setattr("isaura.manage.DuckDBMinio", FakeDuck)
  monkeypatch.setattr("isaura.manage.BloomIndex", FakeBloom)
  monkeypatch.setattr("isaura.manage.STREAM_PARQUET_THRESHOLD", 1)
  monkeypatch.setattr("isaura.manage.query", lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("duckdb path should not run")))

  def fake_stream(*args, **kwargs):
    yield pd.DataFrame([{"input": "a", "x": 1.5}, {"input": "b", "x": 2.5}])

  monkeypatch.setattr("isaura.manage.stream_parquet_filtered", fake_stream)

  reader = IsauraReader(model_id="m", model_version="v1", bucket="bucket", input_csv="unused.csv", approximate=False)
  df = pd.DataFrame([{"input": "a"}, {"input": "b"}])
  out = reader.read(df=df)

  assert list(out["input"]) == ["a", "b"]


def test_reader_read_batched_streaming_mode_skips_duckdb(monkeypatch, tmp_path):
  class FakeStore:
    def __init__(self, *args, **kwargs):
      pass

  class FakeDuck:
    def __init__(self, *args, **kwargs):
      self.con = object()

  class FakeBloom:
    def __init__(self, *args, **kwargs):
      pass

    def seen(self, v):
      return True

  monkeypatch.setattr("isaura.manage.MinioStore", FakeStore)
  monkeypatch.setattr("isaura.manage.DuckDBMinio", FakeDuck)
  monkeypatch.setattr("isaura.manage.BloomIndex", FakeBloom)
  monkeypatch.setattr("isaura.manage.STREAM_PARQUET_THRESHOLD", 1)
  monkeypatch.setattr(
    "isaura.manage.query_batched",
    lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("duckdb path should not run")),
  )

  def fake_stream(*args, **kwargs):
    yield pd.DataFrame([{"input": "a", "x": 1.5}])
    yield pd.DataFrame([{"input": "b", "x": 2.5}])

  monkeypatch.setattr("isaura.manage.stream_parquet_filtered", fake_stream)

  reader = IsauraReader(model_id="m", model_version="v1", bucket="bucket", input_csv="unused.csv", approximate=False)
  out = tmp_path / "out.csv"
  df = pd.DataFrame([{"input": "a"}, {"input": "b"}])
  chunks = list(reader.read_batched(output_csv=str(out), df=df))

  assert list(pd.read_csv(out)["input"]) == ["a", "b"]
  assert [list(chunk["input"]) for chunk in chunks] == [["a"], ["b"]]


def test_base_transfer_select_rows_batched_streaming_mode_skips_duckdb(monkeypatch):
  class FakeStore:
    endpoint = "http://127.0.0.1:9000"
    access = "a"
    secret = "s"

    def __init__(self, *args, **kwargs):
      pass

    def download_file(self, *args, **kwargs):
      raise AssertionError("download_file should not run in this unit test")

  class FakeDuck:
    def __init__(self, *args, **kwargs):
      self.con = object()

  class FakeBloom:
    def __init__(self, *args, **kwargs):
      pass

  monkeypatch.setattr("isaura.base.MinioStore", FakeStore)
  monkeypatch.setattr("isaura.base.DuckDBMinio", FakeDuck)
  monkeypatch.setattr("isaura.base.BloomIndex", FakeBloom)
  monkeypatch.setattr("isaura.base.STREAM_PARQUET_THRESHOLD", 1)
  monkeypatch.setattr(
    "isaura.base.query_batched",
    lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("duckdb path should not run")),
  )

  def fake_stream(*args, **kwargs):
    yield pd.DataFrame([{"input": "a", "x": 1.5}])

  monkeypatch.setattr("isaura.base.stream_parquet_filtered", fake_stream)

  transfer = _BaseTransfer("m", "v1", "bucket")
  chunks = list(transfer.select_rows_batched(["a"]))

  assert len(chunks) == 1
  assert list(chunks[0]["input"]) == ["a"]
