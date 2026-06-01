import sys
import time
from collections import defaultdict

import requests

from isaura.const import (
  BATCH, BUILD_MAX_WAIT, BUILD_POLL_INTERVAL, BUILD_START_TIMEOUT,
  BUILD_STATUS_TIMEOUT, FLUSH_EVERY, NNS_ENDPOINT_BASE, TIMEOUT,
)
from isaura.logging import logger
from isaura.utils import log


def _get_params(collection):
  """Return query params for NNS API requests (collection, batch size, flush interval)."""
  return {"collection": collection, "batch": str(BATCH), "flush_every": str(FLUSH_EVERY)}


def _get_header():
  """Return HTTP headers for plain-text NNS streaming requests."""
  return {"Content-Type": "text/plain"}


def _line_gen(df, chunksize=100000):
  """Yield newline-terminated SMILES strings from a DataFrame as encoded byte chunks.

  Used to stream molecule inputs to the NNS /insert endpoint without loading
  the full DataFrame into a single request body.

  Args:
      df: DataFrame with an "input" column of SMILES strings.
      chunksize: Number of rows to process per iteration.

  Yields:
      Encoded byte strings, one SMILES per line.
  """
  for start in range(0, len(df), chunksize):
    chunk = df.iloc[start:start + chunksize]
    for s in chunk["input"].astype(str):
      yield (s + "\n").encode("utf-8")
    log(f"sent chunk rows={len(chunk):,}")
  log("finished streaming")


def build_index_status(collection):
  """Return the current index build status dict for a collection from the NNS service."""
  r = requests.get(f"{NNS_ENDPOINT_BASE}/build_index", params={"collection": collection}, timeout=BUILD_STATUS_TIMEOUT)
  r.raise_for_status()
  return r.json()


def start_build_index(collection, nlist=None, rebuild=False, wait=False):
  """Trigger an index build for a collection on the NNS service.

  Args:
      collection: Collection name in the NNS service.
      nlist: Optional number of IVF clusters for the index.
      rebuild: If True, force a full rebuild even if an index already exists.
      wait: If True, block until the build is complete.

  Returns:
      Response JSON dict from the NNS service.
  """
  params = {"collection": collection}
  if nlist is not None:
    params["nlist"] = str(nlist)
  if rebuild:
    params["rebuild"] = "1"
  if wait:
    params["wait"] = "1"
  r = requests.post(f"{NNS_ENDPOINT_BASE}/build_index", params=params, timeout=BUILD_START_TIMEOUT)
  r.raise_for_status()
  return r.json()


def ensure_index_ready(collection, max_wait_s=BUILD_MAX_WAIT):
  """Wait until the NNS index for a collection is built and ready for search.

  Triggers a build if no index exists, then polls until it finishes or times out.

  Args:
      collection: Collection name in the NNS service.
      max_wait_s: Maximum seconds to wait before returning a timeout failure.

  Returns:
      Tuple of (success: bool, metadata: dict). On failure, metadata contains
      an "error" key describing what went wrong.
  """
  st = build_index_status(collection)
  if not st.get("exists", False):
    try:
      start_build_index(collection, wait=False)
    except requests.RequestException as e:
      return (False, {"error": f"start_build_index_failed: {e}", "status": st})
    st = build_index_status(collection)
  if st.get("is_failed", False):
    return (False, {"error": "index_failed", "status": st})
  if st.get("is_finished", False) and not st.get("is_building", False):
    return (True, {"status": st})
  t_end = time.time() + max_wait_s
  while time.time() < t_end:
    time.sleep(BUILD_POLL_INTERVAL)
    st = build_index_status(collection)
    if st.get("is_failed", False):
      return (False, {"error": "index_failed", "status": st})
    if st.get("is_finished", False) and not st.get("is_building", False):
      return (True, {"status": st})
  return (False, {"error": "index_building", "status": st})


def post_apprx(df, collection, nlist=None):
  """Stream molecule SMILES to the NNS /insert endpoint and trigger an index rebuild.

  Inserts are streamed line-by-line to avoid large request bodies. After
  insertion, a non-blocking index build is triggered so future searches
  reflect the new data.

  Args:
      df: DataFrame with an "input" column of SMILES strings.
      collection: Target collection in the NNS service.
      nlist: Optional number of IVF clusters for the triggered index build.

  Returns:
      Response object from the /insert request, or None on failure.
  """
  t0 = time.time()
  try:
    with requests.Session() as s:
      logger.debug("start streaming to nns api")
      resp = s.post(
        f"{NNS_ENDPOINT_BASE}/insert",
        params=_get_params(collection),
        data=_line_gen(df),
        headers=_get_header(),
        timeout=None,
      )
    dt = (time.time() - t0) * 1000
    logger.debug(f"done total_ms={dt:.0f}. Body: {resp.text[:1000]}")
    try:
      start_build_index(collection, nlist=nlist, rebuild=False, wait=False)
      st = build_index_status(collection)
      logger.debug(
        f"index build triggered collection={collection} exists={st.get('exists')} "
        f"building={st.get('is_building')} finished={st.get('is_finished')}"
      )
    except requests.RequestException as e:
      logger.debug(f"index build trigger failed collection={collection}: {e}")
    return resp
  except requests.RequestException as e:
    logger.debug(f"approx NN search failed: {e}")
    return None


def get_apprx(inputs, collection, fallback_search=None):
  """Find the nearest stored molecule for each input using the NNS service.

  Ensures the index is ready before searching. For each input molecule,
  returns the SMILES of the closest molecule already in the store (by
  chemical similarity). If the index is not ready or the search fails,
  calls fallback_search if provided.

  Args:
      inputs: List of molecule SMILES strings to search for.
      collection: Collection name in the NNS service.
      fallback_search: Optional callable(inputs, collection) used if NNS fails.

  Returns:
      List of nearest-neighbor SMILES strings (one per input), or [] on failure.
  """
  try:
    logger.info(f"Sending {len(inputs)} inputs for ANN search")
    ok, meta = ensure_index_ready(collection, max_wait_s=BUILD_MAX_WAIT)
    if not ok:
      st = meta.get("status", {})
      err = meta.get("error", "index_not_ready")
      logger.error(f"ANN index not ready collection={collection} error={err}")
      if callable(fallback_search):
        return fallback_search(inputs, collection)
      return []
    r = requests.post(
      f"{NNS_ENDPOINT_BASE}/search",
      json={"collection": collection, "smiles": inputs},
      timeout=TIMEOUT,
    )
    r.raise_for_status()
    payload = r.json() or {}
    results = payload.get("results", [])
    return [x.get("match") for x in results if "match" in x]
  except requests.RequestException as e:
    logger.debug(f"approx NN search failed: {e}")
    if callable(fallback_search):
      return fallback_search(inputs, collection)
    return []


def group_inputs(wanted, index, bloom=None, force=False):
  """Group wanted molecule inputs by their (row, chunk) coordinates from the JSON index.

  Used in approximate read mode to map NNS-resolved SMILES to the right
  Parquet chunk file. If force=True, molecules missing from the index have
  their coordinates computed from SMILES via tranche_coordinates().

  Args:
      wanted: List of molecule SMILES strings to group.
      index: JSON index dict mapping SMILES to (row, chunk) tuples.
      bloom: Optional BloomIndex used to detect missing inputs.
      force: If True, compute coordinates for missing inputs instead of exiting.

  Returns:
      Dict mapping (row, chunk) tuples to sets of SMILES strings, or None on error.
  """
  try:
    miss = []
    if bloom:
      miss = [s for s in wanted if not bloom.seen(s)]
    g = defaultdict(set)
    if miss and not force:
      logger.error(f"inputs not indexed: {miss[:5]} total_missing={len(miss)}")
      sys.exit(1)
    if not force:
      for s in wanted:
        r, c = index[s]
        g[int(r), int(c)].add(s)
      return g
    if miss and force:
      from isaura.metadata import tranche_coordinates
      for s in miss:
        r, c, _, _ = tranche_coordinates(s)
        g[int(r), int(c)].add(s)
      return g
  except Exception as e:
    logger.error(e)
    return None
