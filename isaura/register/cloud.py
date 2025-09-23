import logging, os, subprocess, sys
import boto3, awswrangler as wr, pandas as pd
from pathlib import Path
from config.app import DataLakeConfig, WorkerConfig
from precalculator.models import Prediction, validate_dataframe_schema

INPUT_NAME = "reference_library"
INPUT_FILE_NAME = f"{INPUT_NAME}.csv"
PROCESSED_FILE_NAME = "input.csv"
OUTPUT_FILE_NAME = "output.csv"


class IsauraWriter:
  pass


class IsauraReader:
  pass


class LocalWriter(IsauraWriter):
  pass


class CloudWriter(IsauraWriter):
  pass


class PredictionWriter:
  def __init__(
    self,
    data_config: DataLakeConfig,
    worker_config: WorkerConfig,
    model_id: str,
    dev: bool,
  ):
    self.data_config = data_config
    self.worker_config = worker_config
    self.model_id = model_id
    self.s3 = boto3.client("s3")
    self.dev = dev
    self.logger = logging.getLogger("PredictionWriter")
    self.logger.setLevel(logging.INFO)
    if self.dev:
      logging.basicConfig(stream=sys.stdout, level=logging.INFO)
      logging.getLogger("botocore").setLevel(logging.WARNING)

  def fetch(self) -> str:
    logger = self.logger
    input_filename = (
      f"{INPUT_NAME}_{self.worker_config.sample}.csv"
      if self.worker_config.sample
      else INPUT_FILE_NAME
    )
    self.s3.download_file(
      self.data_config.s3_bucket_name, input_filename, INPUT_FILE_NAME
    )
    logger.info(f"Downloaded {input_filename} from S3")
    start_row, end_row = self._split_csv()
    logger.info(f"Partitioned rows {start_row} to {end_row}")
    return PROCESSED_FILE_NAME

  def predict(self, input_file_path: str) -> str:
    logger = self.logger
    exe = os.environ.get("ERSILIA_EXE", "ersilia")
    logger.info(f"Calling Ersilia CLI for model {self.model_id}")
    subprocess.run([exe, "-v", "fetch", self.model_id, "--from_github"], check=True)
    subprocess.run(
      [exe, "-v", "serve", self.model_id, "--disable-local-cache"], check=True
    )
    subprocess.run(
      [
        exe,
        "-v",
        "run",
        "-i",
        input_file_path,
        "-o",
        OUTPUT_FILE_NAME,
        "--batch_size",
        "1000",
      ],
      check=True,
    )
    return OUTPUT_FILE_NAME

  def postprocess(self, ersilia_output_path: str) -> pd.DataFrame:
    df = pd.read_csv(ersilia_output_path)
    output_cols = [c for c in df.columns if c not in ("key", "input")]
    df["output"] = df[output_cols].apply(
      lambda row: ",".join(str(v) for v in row.values), axis=1
    )
    df["model_id"] = self.model_id
    return df[["key", "input", "output", "model_id"]]

  def write_to_lake(self, outputs: pd.DataFrame) -> None:
    validate_dataframe_schema(outputs, Prediction)  # type: ignore[arg-type]
    wr.s3.to_parquet(
      df=outputs,
      path=os.path.join(
        "s3://",
        self.data_config.s3_bucket_name,
        self.data_config.athena_prediction_table,
      ),
      dataset=True,
      database=self.data_config.athena_database,
      table=self.data_config.athena_prediction_table,
      partition_cols=["model_id"],
    )

  def _split_csv(self) -> tuple[int, int]:
    df = pd.read_csv(INPUT_FILE_NAME)
    n, d = len(df), max(1, int(self.worker_config.denominator))
    k = max(1, int(self.worker_config.numerator))
    base = n // d
    rem = n % d
    sizes = [base + (1 if i < rem else 0) for i in range(d)]
    starts = [0]
    for s in sizes[:-1]:
      starts.append(starts[-1] + s)
    idx = min(max(1, k), d) - 1
    start_row = starts[idx]
    end_row = start_row + sizes[idx]
    df.iloc[start_row:end_row].to_csv(PROCESSED_FILE_NAME, index=False)
    return start_row, end_row
