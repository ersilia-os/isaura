import json
import os
import time

import pyarrow as pa

from isaura.const import (
  DEFAULT_WRITE_BATCH_ROWS,
  MAX_ROWS, MAX_ROWS_PER_FILE,
  PARQUET_COMPRESSION, PARQUET_COMPRESSION_LEVEL,
  PARQUET_DATA_PAGE_SIZE, PARQUET_ROW_GROUP_SIZE,
  WIDE_OUTPUT_DIM_THRESHOLD,
  WIDE_PARQUET_DATA_PAGE_SIZE, WIDE_PARQUET_ROW_GROUP_SIZE,
  WIDE_PARQUET_USE_DICTIONARY, WIDE_PARQUET_WRITE_STATISTICS,
  WIDE_WRITE_BATCH_ROWS,
)
from isaura.logging import logger
from isaura.utils import get_files_glob, hive_prefix


def is_wide(output_dimension):
  """Return True if the output dimension meets the threshold for wide-table mode."""
  return output_dimension is not None and int(output_dimension) >= WIDE_OUTPUT_DIM_THRESHOLD


def chunk_row_limit(output_dimension):
  """Return the max rows per Parquet chunk file, using a lower limit for wide models."""
  return MAX_ROWS_PER_FILE if is_wide(output_dimension) else MAX_ROWS


def chunk_write_batch_rows(output_dimension):
  """Return the number of rows to write per batch, using a smaller batch for wide models."""
  return WIDE_WRITE_BATCH_ROWS if is_wide(output_dimension) else DEFAULT_WRITE_BATCH_ROWS


def parquet_writer_kwargs(output_dimension):
  """Return (constructor_kwargs, write_table_kwargs) for PyArrow ParquetWriter.

  Wide models (100+ output columns) use smaller row groups, disabled
  dictionary encoding, and no statistics to keep file sizes manageable.
  Narrow models use standard compression and row group settings.

  Args:
      output_dimension: Number of output columns, or None if unknown.

  Returns:
      Tuple of (ctor_kw dict, wt_kw dict) to pass to ParquetWriter.
  """
  wide = is_wide(output_dimension)
  ctor_kw = {
    "compression": PARQUET_COMPRESSION,
    "compression_level": PARQUET_COMPRESSION_LEVEL,
  }
  wt_kw = {}
  if wide:
    wt_kw["row_group_size"] = WIDE_PARQUET_ROW_GROUP_SIZE
    ctor_kw["use_dictionary"] = WIDE_PARQUET_USE_DICTIONARY
    ctor_kw["write_statistics"] = WIDE_PARQUET_WRITE_STATISTICS
    ctor_kw["write_batch_size"] = WIDE_PARQUET_DATA_PAGE_SIZE
  else:
    wt_kw["row_group_size"] = PARQUET_ROW_GROUP_SIZE
    ctor_kw["write_batch_size"] = PARQUET_DATA_PAGE_SIZE
  return ctor_kw, wt_kw


def list_parquet_keys(store, bucket, base):
  """List all Parquet chunk keys for a model and return them as s3:// URIs.

  Falls back to a glob URI if the boto3 listing fails.

  Args:
      store: MinioStore used for listing.
      bucket: Source bucket.
      base: Model tranches prefix (e.g. "eos1234/v1/tranches").

  Returns:
      List of s3:// URIs, or a glob URI string if listing fails.
  """
  prefix = hive_prefix(base) + "/"
  t0 = time.perf_counter()
  try:
    keys = []
    for obj in store.list_keys(bucket, prefix):
      k = obj["Key"]
      if k.endswith(".parquet") and "/chunk_" in k:
        keys.append(k)
    keys.sort()
    uris = [f"s3://{bucket}/{k}" for k in keys]
    dt = time.perf_counter() - t0
    logger.debug(f"[list_parquet_keys] listed {len(uris)} files in {dt:.3f}s bucket={bucket} prefix={prefix}")
    return uris if uris else get_files_glob(bucket, base)
  except Exception as e:
    logger.warning(f"[list_parquet_keys] boto3 listing failed ({e}), falling back to glob")
    return get_files_glob(bucket, base)


# --- Column dtype handling (run_columns.csv -> parquet types) -----------------
# Declared run_columns type -> arrow type. float->float32, integer->int32 keep
# files compact; numeric outputs stored as text are the source of the bloat/OOM.
_DECLARED_ARROW_TYPES = {
  "float": pa.float32(),
  "integer": pa.int32(),
  "int": pa.int32(),
  "string": pa.string(),
  "text": pa.string(),
}
# Columns isaura adds itself (not in run_columns); always stored as string.
PREFIX_STRING_COLS = ("key", "input", "smiles")


def arrow_type_for(declared_type):
  """Map a run_columns declared type (e.g. 'float') to an arrow type, or None if unknown."""
  return _DECLARED_ARROW_TYPES.get((declared_type or "").strip().lower())


def resolve_column_types(run_columns, schema_cols):
  """Build {column_name: arrow type} for a model's stored schema.

  Prefix columns (key/input/smiles) are always string; output columns take their
  declared type from run_columns (float->float32, integer->int32, string->string).
  A column absent from run_columns (and not a prefix) defaults to string. Returns
  {} when run_columns is falsy so callers fall back to inference.
  """
  if not run_columns:
    return {}
  types = {}
  for col in schema_cols:
    if col in PREFIX_STRING_COLS:
      types[col] = pa.string()
      continue
    at = arrow_type_for(run_columns.get(col))
    types[col] = at if at is not None else pa.string()
  return types


def validate_columns(run_columns, schema_cols):
  """Check the stored schema's output columns exactly match run_columns (by name).

  schema_cols includes isaura's prefix (key + input/smiles); the rest must equal the
  run_columns declaration, order-independent. Returns (ok, message); the caller
  decides whether to hard-fail. (ok, "") when run_columns is unavailable.
  """
  if not run_columns:
    return True, ""
  output_cols = [c for c in schema_cols if c not in PREFIX_STRING_COLS]
  declared, actual = set(run_columns), set(output_cols)
  if declared == actual:
    return True, ""
  missing, extra = declared - actual, actual - declared
  msg = (
    f"column contract mismatch vs run_columns.csv: {len(actual)} output cols present "
    f"vs {len(declared)} declared; missing={sorted(missing)[:5]} extra={sorted(extra)[:5]}"
  )
  return False, msg


def _stringify(v):
  """Coerce a scalar to str for a string column (None/NaN -> None, collections -> JSON)."""
  if v is None or (isinstance(v, float) and v != v):
    return None
  if isinstance(v, bytes):
    return v.decode("utf-8", errors="replace")
  if isinstance(v, (list, dict)):
    return json.dumps(v, ensure_ascii=False)
  if isinstance(v, (tuple, set)):
    return json.dumps(list(v), ensure_ascii=False)
  return v if isinstance(v, str) else str(v)


def build_typed_array(values, target_type):
  """Build an arrow array of target_type from python values, hard-failing on bad data.

  String targets stringify (collections->JSON). Numeric targets treat blanks/NaN as
  null and then do a SAFE cast that RAISES on a value that cannot be represented (a
  non-numeric in a float/int column, or an integer outside the target's range) — we
  never silently coerce bad values to null.
  """
  if target_type is None or pa.types.is_string(target_type):
    return pa.array([_stringify(v) for v in values], type=pa.string())
  norm = []
  for v in values:
    if v is None or (isinstance(v, float) and v != v) or (isinstance(v, str) and v.strip() == ""):
      norm.append(None)
    else:
      norm.append(v)
  try:
    return pa.array(norm).cast(target_type)  # safe cast: raises on overflow/parse failure
  except (pa.ArrowInvalid, pa.ArrowNotImplementedError, pa.ArrowTypeError) as e:
    raise ValueError(f"failed to cast column to {target_type}: {e}") from e


def resolve_write_types(model_id, schema_cols):
  """Resolve {column_name: arrow type} for a write, enforcing the run_columns contract.

  Fetches the model's run_columns.csv, validates the schema against it, and returns the
  type map. On a column mismatch it raises RuntimeError (hard-fail) unless
  ISAURA_SKIP_COLUMN_CHECK is set. If run_columns.csv is unavailable it logs loudly and
  returns {} so the caller falls back to inferred types.
  """
  from isaura.metadata import fetch_run_columns  # local import avoids any import-cycle risk

  run_columns = fetch_run_columns(model_id)
  if not run_columns:
    logger.error(
      f"[dtype] run_columns.csv unavailable for {model_id} — falling back to inferred types"
    )
    return {}
  ok, msg = validate_columns(run_columns, schema_cols)
  if not ok:
    if os.getenv("ISAURA_SKIP_COLUMN_CHECK", "").strip().lower() in ("1", "true", "yes"):
      logger.warning(f"[dtype] {model_id}: {msg} — proceeding (ISAURA_SKIP_COLUMN_CHECK set)")
    else:
      raise RuntimeError(f"[dtype] {model_id}: {msg}")
  types = resolve_column_types(run_columns, schema_cols)
  n_typed = sum(1 for t in types.values() if not pa.types.is_string(t))
  logger.debug(f"[dtype] {model_id}: applied run_columns types ({n_typed} numeric columns)")
  return types
