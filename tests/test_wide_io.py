import gc
import os
import pandas as pd
import pytest

pq = pytest.importorskip("pyarrow.parquet")
pa = pytest.importorskip("pyarrow")

from isaura.base import TrancheState, _BaseTransfer, LocationSink, build_location_index
from isaura.const import MAX_ROWS, MAX_ROWS_PER_FILE
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


def _wide_chunk_rows(n, ncols=150, start=0):
  rows = []
  for i in range(start, start + n):
    row = {"input": f"s{i}"}
    for j in range(ncols):
      row[f"d{j}"] = float(i * 1000 + j)
    rows.append(row)
  return rows


def _read_all_chunks(tmp_path, bucket="bucket", base="model/v1/tranches"):
  import glob
  d = os.path.join(str(tmp_path), bucket, base, "data")
  paths = sorted(
    glob.glob(os.path.join(d, "chunk_*.parquet")),
    key=lambda p: int(os.path.splitext(os.path.basename(p))[0].split("_")[1]),
  )
  return [pq.read_table(p).to_pandas() for p in paths]


def test_wide_chunk_writes_full_chunks_then_one_partial(tmp_path):
  # output_dimension >= 100 => wide byte-sized accumulation path. Rows from many
  # separate flushes must ACCUMULATE into full chunks, not each be abandoned as
  # its own tiny partial (the pre-W3 ensure()-resets-open bug -> ~2x files).
  store = LocalStore(str(tmp_path))
  tranche = TrancheState(
    store, "bucket", "model/v1/tranches", str(tmp_path), 999999, output_dimension=150
  )
  tranche._calibrated = True  # pin file_rows deterministically; skip byte calibration
  tranche.max_rows = 4
  cols = ["input"] + [f"d{j}" for j in range(150)]
  for i in range(13):
    tranche.flush(_wide_chunk_rows(1, start=i), cols)
  tranche.finalize_chunks()

  chunks = _read_all_chunks(tmp_path)
  # Full chunks of 4, then a single trailing partial — never alternating full/tiny.
  assert [len(c) for c in chunks] == [4, 4, 4, 1]
  out = pd.concat(chunks, ignore_index=True)
  assert list(out["input"]) == [f"s{i}" for i in range(13)]


def test_wide_chunk_calibrates_chunk_size_from_bytes(tmp_path, monkeypatch):
  # Tiny byte targets so file_rows is derived small and files roll over; assert
  # every non-final chunk is exactly the calibrated file_rows (bytes, not a fixed
  # row count) with a single trailing partial, and all rows survive in order.
  monkeypatch.setattr("isaura.base.WIDE_TARGET_FILE_BYTES", 200_000)
  monkeypatch.setattr("isaura.base.WIDE_TARGET_ROWGROUP_BYTES", 50_000)
  store = LocalStore(str(tmp_path))
  tranche = TrancheState(
    store, "bucket", "model/v1/tranches", str(tmp_path), 999999, output_dimension=150
  )
  cols = ["input"] + [f"d{j}" for j in range(150)]
  n = 300
  for start in range(n):  # one row per flush -> exact chunk boundaries
    tranche.flush(_wide_chunk_rows(1, start=start), cols)
  tranche.finalize_chunks()

  assert tranche._calibrated
  file_rows = tranche.max_rows
  chunks = _read_all_chunks(tmp_path)
  sizes = [len(c) for c in chunks]
  assert len(chunks) >= 2  # calibration made file_rows small enough to roll over
  assert all(s == file_rows for s in sizes[:-1])  # full chunks, none abandoned early
  assert 0 < sizes[-1] <= file_rows  # single trailing partial
  assert sum(sizes) == n
  out = pd.concat(chunks, ignore_index=True)
  assert list(out["input"]) == [f"s{i}" for i in range(n)]


def _all_row_group_sizes(tmp_path, bucket="bucket", base="model/v1/tranches"):
  import glob
  d = os.path.join(str(tmp_path), bucket, base, "data")
  paths = sorted(
    glob.glob(os.path.join(d, "chunk_*.parquet")),
    key=lambda p: int(os.path.splitext(os.path.basename(p))[0].split("_")[1]),
  )
  sizes = []
  for p in paths:
    md = pq.ParquetFile(p).metadata
    sizes += [md.row_group(i).num_rows for i in range(md.num_row_groups)]
  return sizes


def test_wide_chunk_row_groups_coalesced_to_target(tmp_path, monkeypatch):
  # Row groups must be coalesced to the calibrated write_batch_rows regardless of how
  # small each flush is; pre-fix each tiny flush became its own row group (the 64 MB
  # WIDE_TARGET_ROWGROUP_BYTES knob was dead). Assert uniform rg-sized groups + one partial.
  monkeypatch.setattr("isaura.base.WIDE_TARGET_FILE_BYTES", 400_000)
  monkeypatch.setattr("isaura.base.WIDE_TARGET_ROWGROUP_BYTES", 100_000)
  store = LocalStore(str(tmp_path))
  tranche = TrancheState(
    store, "bucket", "model/v1/tranches", str(tmp_path), 999999, output_dimension=150
  )
  cols = ["input"] + [f"d{j}" for j in range(150)]
  n = 500
  for start in range(n):  # one row per flush -> must NOT yield 500 one-row row groups
    tranche.flush(_wide_chunk_rows(1, start=start), cols)
  tranche.finalize_chunks()

  rg = tranche.write_batch_rows
  assert rg > 1  # calibration derived a real row-group size
  assert tranche.max_rows % rg == 0  # file is a whole number of row groups
  sizes = _all_row_group_sizes(tmp_path)
  assert sum(sizes) == n  # no rows lost
  assert all(s == rg for s in sizes[:-1])  # coalesced to the target, not per-flush
  assert 0 < sizes[-1] <= rg  # single trailing partial group


def test_wide_loc_capture_matches_actual_layout(tmp_path, monkeypatch):
  # W4 step 1: as wide chunks are written, the LocationSink must record each key's
  # ACTUAL (chunk, row_group). Cross-check against the real layout by reading the input
  # column back per row group from every chunk. Many one-row flushes exercise coalescing.
  monkeypatch.setattr("isaura.base.WIDE_TARGET_FILE_BYTES", 400_000)
  monkeypatch.setattr("isaura.base.WIDE_TARGET_ROWGROUP_BYTES", 100_000)
  store = LocalStore(str(tmp_path))
  tranche = TrancheState(
    store, "bucket", "model/v1/tranches", str(tmp_path), 999999, output_dimension=150
  )
  loc_path = os.path.join(str(tmp_path), "loc.sqlite")
  tranche.loc_sink = LocationSink(loc_path)
  cols = ["input"] + [f"d{j}" for j in range(150)]
  n = 500
  for start in range(n):
    tranche.flush(_wide_chunk_rows(1, start=start), cols)
  tranche.finalize_chunks()
  tranche.loc_sink.close()

  import glob, sqlite3
  truth = {}  # key -> (chunk_idx, rg_ordinal), read back from the actual parquet layout
  d = os.path.join(str(tmp_path), "bucket", "model/v1/tranches", "data")
  for p in glob.glob(os.path.join(d, "chunk_*.parquet")):
    cidx = int(os.path.splitext(os.path.basename(p))[0].split("_")[1])
    pf = pq.ParquetFile(p)
    for rg in range(pf.metadata.num_row_groups):
      for k in pf.read_row_group(rg, columns=["input"]).column("input").to_pylist():
        truth[k] = (cidx, rg)

  conn = sqlite3.connect(loc_path)
  captured = {k: (c, g) for k, c, g in conn.execute("SELECT key, chunk, rg FROM loc")}
  conn.close()

  assert len(captured) == n  # every written key captured exactly once
  assert captured == truth  # captured (chunk, row_group) == the real layout


def test_wide_loc_capture_guard_rejects_unaligned_write(tmp_path):
  # A mislocated key returns a silent blank row on read (the full-scan fallback only fires
  # on keys fully ABSENT), so a non-final write off a row-group boundary must fail loud.
  store = LocalStore(str(tmp_path))
  tranche = TrancheState(
    store, "bucket", "model/v1/tranches", str(tmp_path), 999999, output_dimension=150
  )
  tranche.loc_sink = LocationSink(os.path.join(str(tmp_path), "loc.sqlite"))
  tranche.write_batch_rows = 4  # row-group size
  tbl = pa.table({"input": ["a", "b"], "d0": [1.0, 2.0]})
  with pytest.raises(RuntimeError, match="mislocate"):
    tranche._capture_locs(tbl, chunk_idx=0, start_row=1, final=False)


def test_wide_location_index_built_and_uploaded(tmp_path, monkeypatch):
  # W4 step 2: finalize the captured locations into index.sqlite (unique key index + meta
  # marker) and upload it. Its loc rows must match the actual parquet layout.
  import glob, sqlite3
  monkeypatch.setattr("isaura.base.WIDE_TARGET_FILE_BYTES", 400_000)
  monkeypatch.setattr("isaura.base.WIDE_TARGET_ROWGROUP_BYTES", 100_000)
  store = LocalStore(str(tmp_path))
  base = "model/v1/tranches"
  tranche = TrancheState(store, "bucket", base, str(tmp_path), 999999, output_dimension=150)
  tranche.loc_sink = LocationSink(os.path.join(str(tmp_path), "loc.sqlite"))
  cols = ["input"] + [f"d{j}" for j in range(150)]
  n = 400
  for start in range(n):
    tranche.flush(_wide_chunk_rows(1, start=start), cols)
  tranche.finalize_chunks()
  entries = build_location_index(tranche.loc_sink, store, "bucket", base)
  assert entries == n

  idx_path = os.path.join(str(tmp_path), "bucket", base, "index.sqlite")
  assert os.path.exists(idx_path)  # uploaded to {base}/index.sqlite
  conn = sqlite3.connect(idx_path)
  meta = dict(conn.execute("SELECT k, v FROM meta"))
  assert meta["format"] == "loc-v1"
  assert meta["granularity"] == "chunk_rowgroup"
  assert int(meta["entries"]) == n
  captured = {k: (c, g) for k, c, g in conn.execute("SELECT key, chunk, rg FROM loc")}
  conn.close()

  truth = {}
  for p in glob.glob(os.path.join(str(tmp_path), "bucket", base, "data", "chunk_*.parquet")):
    cidx = int(os.path.splitext(os.path.basename(p))[0].split("_")[1])
    pf = pq.ParquetFile(p)
    for rg in range(pf.metadata.num_row_groups):
      for k in pf.read_row_group(rg, columns=["input"]).column("input").to_pylist():
        truth[k] = (cidx, rg)
  assert len(captured) == n
  assert captured == truth  # index locations == real parquet layout


def test_wide_location_index_skips_upload_when_empty(tmp_path):
  # No captured rows must NOT overwrite an existing index with an empty one.
  store = LocalStore(str(tmp_path))
  base = "model/v1/tranches"
  sink = LocationSink(os.path.join(str(tmp_path), "loc.sqlite"))
  entries = build_location_index(sink, store, "bucket", base)
  assert entries == 0
  assert not os.path.exists(os.path.join(str(tmp_path), "bucket", base, "index.sqlite"))


def _build_wide_store_with_index(tmp_path, monkeypatch, n=400, ncols=150):
  # Small byte targets so the data rolls into several chunks (pruning is observable).
  monkeypatch.setattr("isaura.base.WIDE_TARGET_FILE_BYTES", 200_000)
  monkeypatch.setattr("isaura.base.WIDE_TARGET_ROWGROUP_BYTES", 50_000)
  store = LocalStore(str(tmp_path))
  base = "model/v1/tranches"
  tr = TrancheState(store, "bucket", base, str(tmp_path), 999999, output_dimension=ncols)
  tr.loc_sink = LocationSink(os.path.join(str(tmp_path), "loc.sqlite"))
  cols = ["input"] + [f"d{j}" for j in range(ncols)]
  for start in range(n):
    tr.flush(_wide_chunk_rows(1, start=start), cols)
  tr.finalize_chunks()
  build_location_index(tr.loc_sink, store, "bucket", base)
  return store, base


def _reader_stub(store, base, tmp_path):
  import types
  return types.SimpleNamespace(store=store, bucket="bucket", base=base, tmpdir=str(tmp_path))


def _concat_stream(store, base, wanted, locations):
  parts = list(
    stream_parquet_filtered(store, "bucket", base + "/data/", wanted, header="input", locations=locations)
  )
  df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
  return df.sort_values("input").reset_index(drop=True) if not df.empty else df


def test_wide_index_prune_matches_full_scan(tmp_path, monkeypatch):
  # A clustered read must (a) select fewer chunks than the total and (b) return output
  # byte-identical to a full scan (locations=None). Correctness first, speed second.
  from isaura.manage import IsauraReader
  store, base = _build_wide_store_with_index(tmp_path, monkeypatch, n=400)
  import glob
  n_chunks = len(glob.glob(os.path.join(str(tmp_path), "bucket", base, "data", "chunk_*.parquet")))
  assert n_chunks >= 2  # multiple chunks so pruning is meaningful

  wanted = [f"s{i}" for i in range(8)]  # contiguous -> clustered into few chunks
  locs = IsauraReader._compute_wide_locations(_reader_stub(store, base, tmp_path), wanted)
  assert locs is not None
  assert len(locs) < n_chunks  # pruned away at least one chunk

  full = _concat_stream(store, base, wanted, None)
  pruned = _concat_stream(store, base, wanted, locs)
  assert list(pruned["input"]) == sorted(wanted)  # found them all
  assert full.equals(pruned)  # identical output despite skipping chunks


def test_wide_index_coverage_fallback_returns_none(tmp_path, monkeypatch):
  # A wanted key absent from the index (stale/partial) must force a full scan, never a
  # silent miss. _compute_wide_locations returns None so the streamer scans everything.
  from isaura.manage import IsauraReader
  store, base = _build_wide_store_with_index(tmp_path, monkeypatch, n=200)
  stub = _reader_stub(store, base, tmp_path)
  assert IsauraReader._compute_wide_locations(stub, ["s1", "s2"]) is not None
  assert IsauraReader._compute_wide_locations(stub, ["s1", "NOT_IN_INDEX"]) is None


def test_rewrite_chunk_preserves_zstd_compression(tmp_path):
  # A delete rewrites chunks; it must re-compress with the write path's codec (zstd), not
  # ParquetWriter's snappy default — otherwise removing rows grows the model on disk.
  import types
  from isaura.manage import IsauraMolRemover
  store = LocalStore(str(tmp_path))
  base = "model/v1/tranches"
  chunk_key = f"{base}/data/chunk_1.parquet"
  df = pd.DataFrame(_wide_chunk_rows(50, ncols=150))  # wide (150 output cols)
  tbl = pa.Table.from_pandas(df, preserve_index=False)
  src = str(tmp_path / "src.parquet")
  pq.write_table(tbl, src, compression="snappy")  # simulate the wrong-codec starting point
  store.upload_file(src, "bucket", chunk_key)
  assert pq.ParquetFile(src).metadata.row_group(0).column(0).compression == "SNAPPY"

  stub = types.SimpleNamespace(store=store, bucket="bucket", tmpdir=str(tmp_path))
  before, after = IsauraMolRemover._rewrite_chunk(stub, chunk_key, {"s5", "s6", "s7"})
  assert (before, after) == (50, 47)  # 3 removed

  out = str(tmp_path / "back.parquet")
  store.download_file("bucket", chunk_key, out)
  md = pq.ParquetFile(out).metadata
  assert md.row_group(0).column(0).compression == "ZSTD"  # re-compressed with zstd, not snappy
  assert md.num_rows == 47


def test_wide_index_absent_returns_none(tmp_path):
  # No index.sqlite at all -> full scan (legacy models).
  from isaura.manage import IsauraReader
  store = LocalStore(str(tmp_path))
  base = "model/v1/tranches"
  os.makedirs(os.path.join(str(tmp_path), "bucket", base, "data"), exist_ok=True)
  stub = _reader_stub(store, base, tmp_path)
  assert IsauraReader._compute_wide_locations(stub, ["s1"]) is None


def test_wide_push_merge_unions_bloom_from_sqlite(tmp_path, monkeypatch):
  # New-wide push (4d): source has index.sqlite + no index.json. The merge must union the
  # source's keys into the cloud bloom (so pushed molecules are found) and NOT create a
  # cloud index.json.
  import types
  from isaura.base import BloomIndex
  from isaura.manage import IsauraPush
  store, base = _build_wide_store_with_index(tmp_path, monkeypatch, n=200)  # source bucket="bucket"
  assert not os.path.exists(os.path.join(str(tmp_path), "bucket", base, "index.json"))
  push = types.SimpleNamespace(model_id="model", model_version="v1")
  tmpd = str(tmp_path / "merge_tmp")
  os.makedirs(tmpd, exist_ok=True)
  entries, _ = IsauraPush._merge_bloom_index(push, store, "bucket", store, "cloud", tmpd)
  assert entries == 200
  assert not os.path.exists(os.path.join(str(tmp_path), "cloud", base, "index.json"))  # stays JSON-free
  cloud_bi = BloomIndex(store, "cloud", base, str(tmp_path / "cloud_reload"), load_index=False)
  assert all(cloud_bi.seen(f"s{i}") for i in range(200))  # every pushed key now in cloud bloom


def test_inspect_load_index_falls_back_to_sqlite(tmp_path, monkeypatch):
  # Inspect (4e): a wide model with no index.json must enumerate its keys via index.sqlite.
  from isaura.manage import IsauraInspect
  store, base = _build_wide_store_with_index(tmp_path, monkeypatch, n=150)
  insp = IsauraInspect(cloud=False)
  insp._clients = lambda b: (store,)
  insp.get_json = lambda *a, **k: None  # simulate absent JSON index
  out = insp.load_index("bucket", "model", "v1")
  assert set(out.keys()) == {f"s{i}" for i in range(150)}
  insp.get_json = lambda *a, **k: {"X": [1, 1]}  # JSON present -> use it, don't touch sqlite
  assert insp.load_index("bucket", "model", "v1") == {"X": [1, 1]}
  insp_cloud = IsauraInspect(cloud=True)  # cloud + non-heavy + non-force -> lightweight empty
  assert insp_cloud.load_index("bucket", "model", "v1") == {}


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
  assert chunk_row_limit(150) == MAX_ROWS_PER_FILE
  assert chunk_row_limit(99) == MAX_ROWS
  assert chunk_row_limit(None) == MAX_ROWS


def test_resolve_column_types_maps_declared_to_arrow():
  import pyarrow as pa
  from isaura.parquet import resolve_column_types
  rc = {"d0": "float", "d1": "integer", "d2": "string"}
  t = resolve_column_types(rc, ["key", "input", "d0", "d1", "d2"])
  assert t["key"] == pa.string() and t["input"] == pa.string()
  assert t["d0"] == pa.float32() and t["d1"] == pa.int32() and t["d2"] == pa.string()
  # no run_columns -> empty (caller falls back to inference)
  assert resolve_column_types(None, ["key", "d0"]) == {}


def test_build_typed_array_casts_and_hard_fails():
  import pyarrow as pa
  from isaura.parquet import build_typed_array
  # numeric-as-text -> float32, blanks/None -> null
  a = build_typed_array(["0.5", "", None, "1.25"], pa.float32())
  assert a.type == pa.float32() and a.to_pylist()[1] is None and a.to_pylist()[0] == 0.5
  # integers
  assert build_typed_array(["3", "-7"], pa.int32()).to_pylist() == [3, -7]
  # hard-fail: non-numeric in float, float in int, int32 overflow
  for vals, t in [(["x"], pa.float32()), (["1.5"], pa.int32()), (["3000000000"], pa.int32())]:
    with pytest.raises(ValueError):
      build_typed_array(vals, t)


def test_validate_columns_contract():
  from isaura.parquet import validate_columns
  rc = {"d0": "float", "d1": "float"}
  ok, _ = validate_columns(rc, ["key", "input", "d0", "d1"])
  assert ok
  bad, msg = validate_columns(rc, ["key", "input", "d0"])  # missing d1
  assert not bad and "mismatch" in msg
  assert validate_columns(None, ["key", "d0"])[0]  # no contract -> ok


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


def test_writer_requires_explicit_bucket():
  # API must not silently default to the canonical isaura-public bucket (the CLI requires -pn).
  # The check fires before any store/metadata work, so no mocks are needed.
  from isaura.manage import IsauraWriter
  with pytest.raises(ValueError, match="bucket"):
    IsauraWriter(input_csv="x.csv", model_id="m", model_version="v1")
  with pytest.raises(ValueError, match="bucket"):
    IsauraWriter(input_csv="x.csv", model_id="m", model_version="v1", bucket=None)


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
  assert writer.max_rows == MAX_ROWS_PER_FILE


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

    def finalize_chunks(self):
      return None

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
