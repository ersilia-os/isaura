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


def test_reader_wide_model_uses_fast_stream_then_reorders(monkeypatch):
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
    yield pd.DataFrame([{"input": "b", "x": 2.5}, {"input": "a", "x": 1.5}])

  monkeypatch.setattr("isaura.manage.stream_parquet_filtered", fake_stream)

  reader = IsauraReader(
    model_id="m", model_version="v1", bucket="bucket", input_csv="unused.csv", approximate=False
  )
  df = pd.DataFrame([{"input": "a"}, {"input": "b"}])
  chunks = list(reader.read_batched(df=df))

  # Ordered wide reads stream chunk files, then reorder to wanted order.
  assert called["stream"] == 1
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
  # visible). The wide read streams chunk files (here shuffled, with one input
  # missing) and is responsible for reordering to wanted order with blank
  # placeholders for misses; assert the reader streams that to CSV faithfully.
  stored_rows = [
    {"input": "smiC", "d0": "C-0", "d1": "C-1"},
    {"input": "smiA", "d0": "A-0", "d1": "A-1"},
    {"input": "smiE", "d0": "E-0", "d1": "E-1"},
  ]

  def fake_stream(*args, **kwargs):
    # Rows come back shuffled and missing "smiB" — reorder must fix both.
    yield pd.DataFrame(stored_rows)

  monkeypatch.setattr("isaura.manage.stream_parquet_filtered", fake_stream)

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


def _make_wide_reader(monkeypatch):
  """Build an IsauraReader with mocked store/duck/bloom for unit-testing the
  in-process _reorder_external (which only needs self.tmpdir)."""
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
  return IsauraReader(
    model_id="m", model_version="v1", bucket="bucket", input_csv="unused.csv", approximate=False
  )


def test_reorder_external_worst_case_order_is_bounded(monkeypatch):
  # wanted[0]'s row arrives in the LAST chunk — the case that OOMs an in-memory
  # ordered buffer. Force a tiny bucket span so we exercise multiple buckets and
  # assert the gather never loads more than one span of rows at once.
  monkeypatch.setattr("isaura.manage.WIDE_REORDER_BUCKET_ROWS", 2)
  monkeypatch.setattr("isaura.manage.WIDE_REORDER_MAX_OPEN_BUCKETS", 256)
  reader = _make_wide_reader(monkeypatch)

  wanted = ["a", "b", "c", "d", "e"]
  # Store order is reversed vs wanted, split across chunks (a is last).
  chunks = [
    pd.DataFrame([{"input": "e", "x": 4.0}, {"input": "d", "x": 3.0}]),
    pd.DataFrame([{"input": "c", "x": 2.0}, {"input": "b", "x": 1.0}]),
    pd.DataFrame([{"input": "a", "x": 0.0}]),
  ]

  import isaura.manage as mgr
  real_read_table = mgr.pq.read_table
  max_bucket_rows = {"n": 0}

  def spy_read_table(path, *a, **k):
    tbl = real_read_table(path, *a, **k)
    max_bucket_rows["n"] = max(max_bucket_rows["n"], tbl.num_rows)
    return tbl

  monkeypatch.setattr("isaura.manage.pq.read_table", spy_read_table)

  out = pd.concat(list(reader._reorder_external(wanted, "input", iter(chunks), batch_size=10)), ignore_index=True)

  assert list(out["input"]) == wanted
  assert list(out["x"]) == [0.0, 1.0, 2.0, 3.0, 4.0]
  # No single bucket load exceeded the span (2).
  assert max_bucket_rows["n"] <= 2


def test_reorder_external_pairing_guard_raises_on_mismatch(monkeypatch):
  reader = _make_wide_reader(monkeypatch)
  wanted = ["a", "b"]

  # Corrupt the bucket on read-back so a stored row's header no longer matches
  # the input position it was filed under — the fatal mispairing the guard
  # must catch instead of silently emitting wrong descriptors.
  import pyarrow as pa
  import isaura.manage as mgr
  real_read_table = mgr.pq.read_table

  def corrupt_read_table(path, *a, **k):
    tbl = real_read_table(path, *a, **k)
    df = tbl.to_pandas()
    df["input"] = "WRONG"
    return pa.Table.from_pandas(df, preserve_index=False)

  monkeypatch.setattr("isaura.manage.pq.read_table", corrupt_read_table)

  chunks = [pd.DataFrame([{"input": "a", "x": 1.0}, {"input": "b", "x": 2.0}])]
  with pytest.raises(ValueError, match="FATAL"):
    list(reader._reorder_external(wanted, "input", iter(chunks), batch_size=10))


def test_reorder_external_preserves_duplicates(monkeypatch):
  reader = _make_wide_reader(monkeypatch)
  wanted = ["smiA", "smiB", "smiA", "smiC", "smiA"]
  chunks = [
    pd.DataFrame([
      {"input": "smiC", "d0": "C0"},
      {"input": "smiA", "d0": "A0"},
      {"input": "smiB", "d0": "B0"},
    ]),
  ]
  out = pd.concat(list(reader._reorder_external(wanted, "input", iter(chunks), batch_size=10)), ignore_index=True)
  assert list(out["input"]) == wanted
  # All three smiA positions carry smiA's own descriptor.
  assert list(out["d0"]) == ["A0", "B0", "A0", "C0", "A0"]


def test_reorder_external_all_misses_returns_blank_rows(monkeypatch):
  reader = _make_wide_reader(monkeypatch)
  wanted = ["x", "y", "z"]
  # Source yields nothing found.
  out = pd.concat(
    list(reader._reorder_external(wanted, "input", iter([pd.DataFrame()]), batch_size=10)),
    ignore_index=True,
  )
  assert list(out["input"]) == wanted
  assert len(out) == 3


def test_reorder_external_cleans_spill_dir(monkeypatch):
  reader = _make_wide_reader(monkeypatch)
  before = set(os.listdir(reader.tmpdir)) if os.path.isdir(reader.tmpdir) else set()
  chunks = [pd.DataFrame([{"input": "a", "x": 1.0}])]
  list(reader._reorder_external(["a"], "input", iter(chunks), batch_size=10))
  after = set(os.listdir(reader.tmpdir)) if os.path.isdir(reader.tmpdir) else set()
  # No leftover reorder_* spill directories.
  assert not [d for d in after - before if d.startswith("reorder_")]


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


class _RecordingParquetStore:
  """Serves real parquet chunk files from disk and records downloaded keys."""

  def __init__(self, files):
    # files: list of (key, path)
    self.files = list(files)
    self.downloaded = []

  def list_keys(self, bucket, prefix):
    for key, _ in self.files:
      yield {"Key": key}

  def download_file(self, bucket, key, local):
    self.downloaded.append(key)
    source = dict(self.files)[key]
    with open(source, "rb") as src, open(local, "wb") as dst:
      dst.write(src.read())


def _build_chunk_files(tmp_path, n_files, rows_per_file):
  files = []
  for ci in range(n_files):
    df = pd.DataFrame(
      [{"input": f"m{ci}_{ri}", "x": float(ci * 1000 + ri)} for ri in range(rows_per_file)]
    )
    path = tmp_path / f"chunk_{ci}.parquet"
    df.to_parquet(path, index=False)
    files.append((f"prefix/chunk_{ci}.parquet", str(path)))
  return files


@pytest.mark.parametrize("prefetch", [1, 3, 10])
def test_stream_parquet_filtered_results_invariant_to_prefetch_depth(monkeypatch, tmp_path, prefetch):
  files = _build_chunk_files(tmp_path, n_files=6, rows_per_file=4)
  wanted = [f"m{ci}_{ri}" for ci in range(6) for ri in range(4)]

  monkeypatch.setattr("isaura.stream.STREAM_PREFETCH_FILES", prefetch)
  monkeypatch.setattr("isaura.stream.STREAM_DOWNLOAD_WORKERS", 4)

  store = _RecordingParquetStore(files)
  got = []
  for chunk in stream_parquet_filtered(store, "bucket", "prefix", wanted, header="input", batch_size=2):
    got.extend(chunk["input"].tolist())

  # Unordered stream: every wanted row returned exactly once, regardless of prefetch depth.
  assert sorted(got) == sorted(wanted)


def test_stream_parquet_filtered_early_exit_bounds_downloads(monkeypatch, tmp_path):
  # All wanted inputs live in the first chunk; later files must not all be downloaded.
  files = _build_chunk_files(tmp_path, n_files=20, rows_per_file=4)
  wanted = ["m0_0", "m0_1", "m0_2", "m0_3"]

  monkeypatch.setattr("isaura.stream.STREAM_PREFETCH_FILES", 3)
  monkeypatch.setattr("isaura.stream.STREAM_DOWNLOAD_WORKERS", 2)

  store = _RecordingParquetStore(files)
  got = []
  for chunk in stream_parquet_filtered(store, "bucket", "prefix", wanted, header="input", batch_size=2):
    got.extend(chunk["input"].tolist())

  assert sorted(got) == sorted(wanted)
  # Early-exit caps wasted downloads at the prefetch window past the resolving file.
  assert len(store.downloaded) < len(files)
  assert len(store.downloaded) <= 1 + 3
