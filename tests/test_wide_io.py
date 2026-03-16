import os
import pandas as pd
import pytest

pq = pytest.importorskip("pyarrow.parquet")

from isaura.base import TrancheState, _BaseTransfer
from isaura.helpers import chunk_row_limit, stream_parquet_filtered_ordered
from isaura.manage import IsauraReader, IsauraWriter


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
  monkeypatch.setattr(
    "isaura.manage.list_parquet_keys", lambda *args, **kwargs: ["s3://bucket/fake/chunk_1.parquet"]
  )

  called = {"query_batched": 0}

  def fake_query_batched(*args, **kwargs):
    called["query_batched"] += 1
    yield pd.DataFrame([
      {"input": "a", "x": 1.5},
      {"input": "b", "x": 2.5},
    ])

  monkeypatch.setattr("isaura.manage.query_batched", fake_query_batched)

  reader = IsauraReader(
    model_id="m", model_version="v1", bucket="bucket", input_csv="unused.csv", approximate=False
  )
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
  monkeypatch.setattr(
    "isaura.manage.query",
    lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("duckdb path should not run")),
  )

  def fake_ordered_stream(*args, **kwargs):
    yield pd.DataFrame([{"input": "a", "x": 1.5}, {"input": "b", "x": 2.5}])

  monkeypatch.setattr("isaura.manage.stream_parquet_filtered_ordered", fake_ordered_stream)

  reader = IsauraReader(
    model_id="m", model_version="v1", bucket="bucket", input_csv="unused.csv", approximate=False
  )
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

  def fake_ordered_stream(*args, **kwargs):
    yield pd.DataFrame([{"input": "a", "x": 1.5}])
    yield pd.DataFrame([{"input": "b", "x": 2.5}])

  monkeypatch.setattr("isaura.manage.stream_parquet_filtered_ordered", fake_ordered_stream)

  reader = IsauraReader(
    model_id="m", model_version="v1", bucket="bucket", input_csv="unused.csv", approximate=False
  )
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


def test_chunk_row_limit_uses_output_dimension_threshold():
  assert chunk_row_limit(150) == 100000
  assert chunk_row_limit(99) == 2000000
  assert chunk_row_limit(None) == 2000000


def test_reader_wide_model_uses_fast_stream_then_reorders(monkeypatch):
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
  monkeypatch.setattr("isaura.manage.fetch_schema_from_github", lambda model_id: {"OutputDimension": 150})

  called = {"stream": 0}

  def fake_stream(*args, **kwargs):
    called["stream"] += 1
    yield pd.DataFrame([{"input": "b", "x": 2.5}, {"input": "a", "x": 1.5}])

  monkeypatch.setattr("isaura.manage.stream_parquet_filtered", fake_stream)

  reader = IsauraReader(
    model_id="m", model_version="v1", bucket="bucket", input_csv="unused.csv", approximate=False
  )
  df = pd.DataFrame([{"input": "a"}, {"input": "b"}])
  chunks = list(reader.read_batched(df=df))

  assert called["stream"] == 1
  all_rows = pd.concat(chunks, ignore_index=True)
  assert list(all_rows["input"]) == ["a", "b"]


def test_writer_uses_small_chunks_for_wide_models(monkeypatch):
  class FakeStore:
    def __init__(self, *args, **kwargs):
      pass

    def ensure_bucket(self, bucket):
      return None

  class FakeBloom:
    def __init__(self, *args, **kwargs):
      self._added = 0

    def seen(self, v):
      return False

    def register(self, v, rc=None):
      self._added += 1

    def persist(self):
      return None

  class FakeChunkState:
    def __init__(self, *args, **kwargs):
      self.max_rows = args[-1]

    def flush_df(self, df, schema_cols):
      return None

  monkeypatch.setattr("isaura.manage.MinioStore", FakeStore)
  monkeypatch.setattr("isaura.manage.BloomIndex", FakeBloom)
  monkeypatch.setattr("isaura.manage.TrancheState", FakeChunkState)
  monkeypatch.setattr("isaura.manage.fetch_schema_from_github", lambda model_id: {"OutputDimension": 150})

  writer = IsauraWriter(input_csv="unused.csv", model_id="m", model_version="v1", bucket="bucket")
  assert writer.max_rows == 100000


def test_writer_buffers_raw_rows_until_flush(monkeypatch):
  flushed = {}

  class FakeStore:
    def __init__(self, *args, **kwargs):
      pass

    def ensure_bucket(self, bucket):
      return None

    def upload_file(self, *args, **kwargs):
      return None

    def download_file(self, *args, **kwargs):
      raise FileNotFoundError

  class FakeBloom:
    def __init__(self, *args, **kwargs):
      self._added = 0

    def seen(self, v):
      return False

    def register(self, v, rc=None):
      self._added += 1

    def persist(self):
      return None

  class FakeChunkState:
    def __init__(self, *args, **kwargs):
      pass

    def flush(self, rows, schema_cols):
      flushed["rows"] = list(rows)
      flushed["schema_cols"] = list(schema_cols)

  monkeypatch.setattr("isaura.manage.MinioStore", FakeStore)
  monkeypatch.setattr("isaura.manage.BloomIndex", FakeBloom)
  monkeypatch.setattr("isaura.manage.TrancheState", FakeChunkState)
  monkeypatch.setattr("isaura.manage.fetch_schema_from_github", lambda model_id: {"OutputDimension": 10})

  writer = IsauraWriter(input_csv="unused.csv", model_id="m", model_version="v1", bucket="bucket")
  writer.max_rows = 2
  writer.write(df=pd.DataFrame([{"input": "a", "x": 1}, {"input": "b", "x": 2}]), show_progress=False)

  assert writer.buffers == []
  assert [row["input"] for row in flushed["rows"]] == ["a", "b"]
  assert flushed["schema_cols"] == ["input", "x"]


def test_stream_parquet_filtered_ordered_handles_single_missing(monkeypatch):
  def fake_stream(*args, **kwargs):
    yield pd.DataFrame([{"input": "b", "x": 2.0}, {"input": "a", "x": 1.0}])
    yield pd.DataFrame([{"input": "d", "x": 4.0}])

  monkeypatch.setattr("isaura.stream.stream_parquet_filtered", fake_stream)

  chunks = list(
    stream_parquet_filtered_ordered(
      store=None,
      bucket="bucket",
      prefix="prefix",
      wanted=["a", "missing", "b", "d"],
      header="input",
      batch_size=10,
    )
  )

  out = pd.concat(chunks, ignore_index=True)
  assert list(out["input"]) == ["a", "missing", "b", "d"]
  assert pd.isna(out.loc[1, "x"])
