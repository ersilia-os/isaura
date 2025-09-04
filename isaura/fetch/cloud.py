import requests, time, uuid
from typing import Iterable, Dict, Any
from isaura.helpers import (
  download_progress,
  get_schema,
  collect_gz_shards,
  process_shard,
  logger,
  API_BASE,
  TIMEOUT,
  JobStatus,
)


class GetPrecalCloud:
  def __init__(self, model_id: str, file_name: str = None, nsmaples: int = 1000):
    super().__init__()
    self.nsamples = nsmaples
    self.model_id = model_id
    self.file_name = file_name
    self.request_id = None
    self.fetch_type = "all"
    self.schema = get_schema(self.model_id)
    assert self.schema is not None, "Model schema can not be fetched from github."
    self.header = ["key", "input"] + self.schema[0]
    self.api = ApiClient()

  def fetch(self) -> None:
    shards = self._submit_and_get_shards()
    self.merge_shards(shards, self.model_id)

  def _submit_and_get_shards(self) -> list:
    self.request_id = str(uuid.uuid4())

    logger.info(f"Submitting job for model_id={self.model_id} with all samples!")

    payload = {
      "requestid": self.request_id,
      "modelid": self.model_id,
      "fetchtype": self.fetch_type,
      "nsamples": self.nsamples,
      "dim": len(self.schema[0]),
    }
    job_id = self.api.post_json(f"{API_BASE}/submit", json=payload)["jobId"]

    logger.info(f"Job submitted successfully. job_id={job_id}")

    start = time.time()
    while True:
      status = self.api.get_json(f"{API_BASE}/status", params={"jobId": job_id})["status"]

      logger.info(f"Job status: {status}")

      if status == JobStatus.SUCCEEDED:
        logger.info("Job succeeded.")
        break

      if status == JobStatus.FAILED:
        error = self.api.get_json(f"{API_BASE}/status", params={"jobId": job_id}).get(
          "errorMessage", "Unknown error"
        )
        logger.error(f"Job failed: {error}")
        raise RuntimeError(f"Job failed: {error}")

      if time.time() - start > 3600:
        logger.error("Job polling timed out.")
        raise RuntimeError("Job polling timed out")

      time.sleep(5)

    shards = self.api.get_json(f"{API_BASE}/result", params={"jobId": job_id})["files"]

    if not shards:
      logger.error(f"No shards returned for job_id={job_id}")
      raise RuntimeError("No shards returned")

    return shards

  def merge_shards(
    self,
    shards: Iterable[Dict[str, Any]],
    model_id: str,
    batch_size: int = 5000,
  ):
    gz_shards, total_size = collect_gz_shards(shards, self.api)
    logger.info(f"⬇ Total download size: {total_size / 1024:.1f} KB")

    logger.info("⬇ No reorder list provided; writing to Milvus")

    with download_progress(desc="⬇ downloading shards", total_bytes=total_size, transient=True) as (
      progress,
      task,
    ):
      for shard in gz_shards:
        process_shard(shard, self.api, batch_size=batch_size, progress=progress, task=task)


class ApiClient:
  @staticmethod
  def get_json(endpoint: str, params: dict = None):
    r = requests.get(endpoint, params=params, timeout=TIMEOUT)
    r.raise_for_status()
    return r.json()

  @staticmethod
  def post_json(endpoint: str, json: dict = None):
    r = requests.post(endpoint, json=json, timeout=TIMEOUT)
    r.raise_for_status()
    return r.json()

  @staticmethod
  def head(url: str):
    r = requests.head(url, timeout=TIMEOUT)
    r.raise_for_status()
    return r.headers

  @staticmethod
  def download_stream(url: str, chunk_size: int = 64 * 1024):
    r = requests.get(url, stream=True, timeout=TIMEOUT)
    r.raise_for_status()
    for chunk in r.iter_content(chunk_size):
      yield chunk
