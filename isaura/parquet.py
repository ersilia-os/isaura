import time

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
    logger.info(f"[list_parquet_keys] listed {len(uris)} files in {dt:.3f}s bucket={bucket} prefix={prefix}")
    return uris if uris else get_files_glob(bucket, base)
  except Exception as e:
    logger.warning(f"[list_parquet_keys] boto3 listing failed ({e}), falling back to glob")
    return get_files_glob(bucket, base)
