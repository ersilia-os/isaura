import gc
import os
import pandas as pd
import pytest

pq = pytest.importorskip("pyarrow.parquet")

from isaura.base import TrancheState, _BaseTransfer
from isaura.helpers import chunk_row_limit, stream_parquet_filtered_ordered
from isaura.manage import IsauraReader, IsauraWriter
from isaura.stream import stream_parquet_filtered


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


def test_tranche_flush_supports_heterogeneous_column_types(tmp_path):
  store = LocalStore(str(tmp_path))
  tranche = TrancheState(store, "bucket", "model/v1/tranches", str(tmp_path), 10)
  rows = [
    {
      "input": "a",
      "mixed": 1.25,
      "payload": {"x": 1},
      "blob": b"abc",
      "flag": True,
    },
    {
      "input": "b",
      "mixed": "C1=CC=CC=C1",
      "payload": ["u", "v"],
      "blob": float("nan"),
      "flag": None,
    },
  ]
  cols = ["input", "mixed", "payload", "blob", "flag"]

  tranche.flush(rows, cols)

  key = os.path.join(str(tmp_path), "bucket", "model/v1/tranches/data/chunk_1.parquet")
  out = pq.read_table(key).to_pandas()
  assert list(out["input"]) == ["a", "b"]
  assert list(out["mixed"]) == ["1.25", "C1=CC=CC=C1"]
  assert list(out["payload"]) == ['{"x": 1}', '["u", "v"]']
  assert list(out["blob"]) == [b"abc", None]
  assert list(out["flag"]) == [True, None]


def test_reader_read_batched_to_csv_skips_index_load(monkeypatch, tmp_path):
  class FakeStore:
    def __init__(self, *args, **kwargs):
      pass

    def require_bucket(self, bucket):
      return None

    def ensure_bucket(self, bucket):
      return None

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

    def require_bucket(self, bucket):
      return None

    def ensure_bucket(self, bucket):
      return None

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

    def require_bucket(self, bucket):
      return None

    def ensure_bucket(self, bucket):
      return None

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

    def require_bucket(self, bucket):
      return None

    def ensure_bucket(self, bucket):
      return None

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


def test_reader_wide_model_ordered_uses_duckdb_keep_missing(monkeypatch):
  class FakeStore:
    def __init__(self, *args, **kwargs):
      pass

    def require_bucket(self, bucket):
      return None

    def ensure_bucket(self, bucket):
      return None

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

  called = {"query": 0, "stream": 0, "kwargs": None}

  def fake_query_batched(con, header, wanted, files, **kwargs):
    called["query"] += 1
    called["kwargs"] = kwargs
    # Simulate the function's contract: rows in wanted order.
    yield pd.DataFrame({header: list(wanted), "x": [float(i) for i in range(len(wanted))]})

  def fake_stream(*args, **kwargs):
    called["stream"] += 1
    yield pd.DataFrame()

  monkeypatch.setattr("isaura.manage.query_batched", fake_query_batched)
  monkeypatch.setattr("isaura.manage.stream_parquet_filtered", fake_stream)

  reader = IsauraReader(
    model_id="m", model_version="v1", bucket="bucket", input_csv="unused.csv", approximate=False
  )
  df = pd.DataFrame([{"input": "a"}, {"input": "b"}])
  chunks = list(reader.read_batched(df=df))  # ordered=True default (read command)

  # Ordered wide reads go through the single-pass DuckDB path with order + miss
  # preservation, not the unordered streaming path.
  assert called["query"] == 1
  assert called["stream"] == 0
  assert called["kwargs"].get("preserve_order") is True
  assert called["kwargs"].get("keep_missing") is True
  all_rows = pd.concat(chunks, ignore_index=True)
  assert list(all_rows["input"]) == ["a", "b"]


def test_reader_wide_model_unordered_streams_without_reordering(monkeypatch):
  class FakeStore:
    def __init__(self, *args, **kwargs):
      pass

    def require_bucket(self, bucket):
      return None

    def ensure_bucket(self, bucket):
      return None

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
    yield pd.DataFrame([{"input": "b", "x": 2.5}])
    yield pd.DataFrame([{"input": "a", "x": 1.5}])

  monkeypatch.setattr("isaura.manage.stream_parquet_filtered", fake_stream)

  reader = IsauraReader(
    model_id="m", model_version="v1", bucket="bucket", input_csv="unused.csv", approximate=False
  )
  # ordered=False must never materialize the full frame or reorder it.
  reader._reorder = lambda *a, **k: (_ for _ in ()).throw(AssertionError("should not reorder"))
  reader._reorder_batched = lambda *a, **k: (_ for _ in ()).throw(AssertionError("should not reorder"))

  df = pd.DataFrame([{"input": "a"}, {"input": "b"}])
  chunks = list(reader.read_batched(df=df, ordered=False))

  assert called["stream"] == 1
  # Chunks are emitted incrementally in store order, not reordered to wanted order.
  assert [list(chunk["input"]) for chunk in chunks] == [["b"], ["a"]]


def test_reader_wide_model_read_preserves_input_descriptor_pairing(monkeypatch, tmp_path):
  """Critical: read output must be in input order with each row's descriptors
  staying matched to its own SMILES, even when the store returns rows shuffled
  and some inputs are missing."""
  class FakeStore:
    def __init__(self, *args, **kwargs):
      pass

    def require_bucket(self, bucket):
      return None

    def ensure_bucket(self, bucket):
      return None

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

  # Stored descriptors (each value encodes its own input so a mismatch is
  # visible). The chunked DuckDB path is responsible for returning rows in
  # wanted order with placeholders for misses; here we mock that contract and
  # assert the reader streams it to CSV faithfully (the reorder logic itself is
  # tested against a real DuckDB in test_chunked_query_batched_*).
  stored = {
    "smiC": {"d0": "C-0", "d1": "C-1"},
    "smiA": {"d0": "A-0", "d1": "A-1"},
    "smiE": {"d0": "E-0", "d1": "E-1"},
  }

  def fake_query_batched(con, header, wanted, files, **kwargs):
    assert kwargs.get("preserve_order") is True
    assert kwargs.get("keep_missing") is True
    rows = []
    for k in wanted:
      desc = stored.get(k, {"d0": None, "d1": None})
      rows.append({header: k, **desc})
    yield pd.DataFrame(rows)

  monkeypatch.setattr("isaura.manage.query_batched", fake_query_batched)

  reader = IsauraReader(
    model_id="m", model_version="v1", bucket="bucket", input_csv="unused.csv", approximate=False
  )
  # Input order includes a missing molecule ("smiB") in the middle.
  df = pd.DataFrame([{"input": s} for s in ["smiA", "smiB", "smiC", "smiE"]])
  out = tmp_path / "out.csv"
  list(reader.read_batched(output_csv=str(out), df=df))  # default ordered=True

  result = pd.read_csv(out)
  # Rows in exact input order, one row per input, missing input kept as a row.
  assert list(result["input"]) == ["smiA", "smiB", "smiC", "smiE"]
  # Every present input keeps its OWN descriptors.
  for s, prefix in [("smiA", "A"), ("smiC", "C"), ("smiE", "E")]:
    row = result[result["input"] == s].iloc[0]
    assert row["d0"] == f"{prefix}-0"
    assert row["d1"] == f"{prefix}-1"
  # Missing input has blank descriptors, not another molecule's values.
  missing = result[result["input"] == "smiB"].iloc[0]
  assert pd.isna(missing["d0"]) and pd.isna(missing["d1"])


def test_query_batched_ordered_keep_missing_preserves_pairing(tmp_path):
  """Core guard: the single-pass ordered query (preserve_order + keep_missing)
  must return rows in exact wanted order (including duplicates), each input
  paired with its own descriptors, and missing inputs as placeholder (NaN) rows.
  Uses a real DuckDB over a local parquet file (no S3)."""
  duckdb = pytest.importorskip("duckdb")
  from isaura.query import query_batched

  # Parquet rows are deliberately shuffled; each descriptor encodes its input.
  stored = pd.DataFrame([
    {"input": "smiC", "d0": "C-0", "d1": "C-1"},
    {"input": "smiA", "d0": "A-0", "d1": "A-1"},
    {"input": "smiE", "d0": "E-0", "d1": "E-1"},
    {"input": "smiD", "d0": "D-0", "d1": "D-1"},
  ])
  pq_path = tmp_path / "chunk_1.parquet"
  stored.to_parquet(pq_path, index=False)

  # Wanted order: a miss ("smiMISS") in the middle, and a duplicate ("smiA").
  wanted = ["smiA", "smiMISS", "smiC", "smiD", "smiE", "smiA"]

  con = duckdb.connect(":memory:")
  chunks = list(
    query_batched(
      con, "input", wanted, [str(pq_path)],
      batch_size=10, tmpdir=str(tmp_path), preserve_order=True, keep_missing=True,
    )
  )
  out = pd.concat(chunks, ignore_index=True)

  # Exact input order, exactly one row per wanted entry (incl. the duplicate).
  assert list(out["input"]) == wanted
  assert len(out) == len(wanted)
  # Every present input carries its OWN descriptors.
  for i, k in enumerate(wanted):
    if k == "smiMISS":
      continue
    prefix = k[-1]  # smiA -> A
    assert out.loc[i, "d0"] == f"{prefix}-0"
    assert out.loc[i, "d1"] == f"{prefix}-1"
  # Missing input is a placeholder row: input preserved, descriptors NaN.
  miss_idx = wanted.index("smiMISS")
  assert pd.isna(out.loc[miss_idx, "d0"])
  assert pd.isna(out.loc[miss_idx, "d1"])


def test_writer_uses_small_chunks_for_wide_models(monkeypatch):
  class FakeStore:
    def __init__(self, *args, **kwargs):
      pass

    def require_bucket(self, bucket):
      return None

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

    def require_bucket(self, bucket):
      return None

    def ensure_bucket(self, bucket):
      return None

    def upload_file(self, *args, **kwargs):
      return None

    def download_file(self, *args, **kwargs):
      raise FileNotFoundError

  class FakeBloom:
    def __init__(self, *args, **kwargs):
      self._added = 0
      self.index = {}

    def seen(self, v):
      return False

    def register(self, v, rc=None):
      self._added += 1

    def persist(self):
      return None

  class FakeChunkState:
    def __init__(self, *args, **kwargs):
      pass

    def _list_chunks(self):
      return []

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


def test_stream_parquet_filtered_ordered_closes_inner_stream_on_early_completion(monkeypatch):
  closed = {"value": False}

  def fake_stream(*args, **kwargs):
    try:
      yield pd.DataFrame([{"input": "a", "x": 1.0}])
      yield pd.DataFrame([{"input": "b", "x": 2.0}])
    finally:
      closed["value"] = True

  monkeypatch.setattr("isaura.stream.stream_parquet_filtered", fake_stream)

  chunks = list(
    stream_parquet_filtered_ordered(
      store=None,
      bucket="bucket",
      prefix="prefix",
      wanted=["a"],
      header="input",
      batch_size=10,
    )
  )

  assert list(pd.concat(chunks, ignore_index=True)["input"]) == ["a"]
  assert closed["value"] is True


def test_reader_tempdir_is_cleaned_when_instance_is_collected(monkeypatch, tmp_path):
  created = []

  class FakeStore:
    def __init__(self, *args, **kwargs):
      pass

    def require_bucket(self, bucket):
      return None

    def ensure_bucket(self, bucket):
      return None

  class FakeDuck:
    def __init__(self, *args, **kwargs):
      self.con = object()

  class FakeBloom:
    def __init__(self, *args, **kwargs):
      pass

    def seen(self, v):
      return True

  def fake_make_temp(prefix):
    path = tmp_path / f"{prefix}{len(created)}"
    path.mkdir()
    created.append(path)
    return str(path)

  monkeypatch.setattr("isaura.manage.MinioStore", FakeStore)
  monkeypatch.setattr("isaura.manage.DuckDBMinio", FakeDuck)
  monkeypatch.setattr("isaura.manage.BloomIndex", FakeBloom)
  monkeypatch.setattr("isaura.manage.fetch_schema_from_github", lambda model_id: {})
  monkeypatch.setattr("isaura.manage.make_temp", fake_make_temp)

  def build_reader():
    reader = IsauraReader(
      model_id="m", model_version="v1", bucket="bucket", input_csv="unused.csv", approximate=False
    )
    assert os.path.isdir(reader.tmpdir)
    return reader.tmpdir

  tempdir = build_reader()
  gc.collect()

  assert not os.path.exists(tempdir)


def test_base_transfer_close_removes_owned_tempdirs(monkeypatch, tmp_path):
  created = []

  class FakeStore:
    endpoint = "http://127.0.0.1:9000"
    access = "a"
    secret = "s"

    def __init__(self, *args, **kwargs):
      pass

    def require_bucket(self, bucket):
      return None

    def ensure_bucket(self, bucket):
      return None

  class FakeDuck:
    def __init__(self, *args, **kwargs):
      self.con = object()

  class FakeBloom:
    def __init__(self, *args, **kwargs):
      pass

  def fake_make_temp(prefix):
    path = tmp_path / f"{prefix}{len(created)}"
    path.mkdir()
    created.append(path)
    return str(path)

  monkeypatch.setattr("isaura.base.MinioStore", FakeStore)
  monkeypatch.setattr("isaura.base.DuckDBMinio", FakeDuck)
  monkeypatch.setattr("isaura.base.BloomIndex", FakeBloom)
  monkeypatch.setattr("isaura.base.make_temp", fake_make_temp)

  transfer = _BaseTransfer("m", "v1", "bucket")
  owned = [transfer.tmpdir, transfer.tmpdir_sinkw]

  assert all(os.path.isdir(path) for path in owned)

  transfer.close()

  assert all(not os.path.exists(path) for path in owned)


def test_stream_parquet_filtered_cleans_tmpdir_when_closed_early(monkeypatch, tmp_path):
  source = tmp_path / "chunk_1.parquet"
  pd.DataFrame([{"input": "a", "x": 1.0}]).to_parquet(source, index=False)

  created = []

  class FakeStore:
    def list_keys(self, bucket, prefix):
      yield {"Key": f"{prefix}/chunk_1.parquet"}

    def download_file(self, bucket, key, local):
      with open(source, "rb") as src, open(local, "wb") as dst:
        dst.write(src.read())

  def fake_make_temp(prefix):
    path = tmp_path / f"{prefix}{len(created)}"
    path.mkdir()
    created.append(path)
    return str(path)

  monkeypatch.setattr("isaura.stream.make_temp", fake_make_temp)

  gen = stream_parquet_filtered(FakeStore(), "bucket", "prefix", ["a"], header="input", batch_size=1)
  chunk = next(gen)
  tempdir = created[0]

  assert list(chunk["input"]) == ["a"]
  assert tempdir.exists()

  gen.close()

  assert not tempdir.exists()
