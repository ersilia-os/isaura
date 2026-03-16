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
  return {"collection": collection, "batch": str(BATCH), "flush_every": str(FLUSH_EVERY)}


def _get_header():
  return {"Content-Type": "text/plain"}


def _line_gen(df, chunksize=100000):
  for start in range(0, len(df), chunksize):
    chunk = df.iloc[start:start + chunksize]
    for s in chunk["input"].astype(str):
      yield (s + "\n").encode("utf-8")
    log(f"sent chunk rows={len(chunk):,}")
  log("finished streaming")


def build_index_status(collection):
  r = requests.get(f"{NNS_ENDPOINT_BASE}/build_index", params={"collection": collection}, timeout=BUILD_STATUS_TIMEOUT)
  r.raise_for_status()
  return r.json()


def start_build_index(collection, nlist=None, rebuild=False, wait=False):
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
  t0 = time.time()
  try:
    with requests.Session() as s:
      log("start streaming to nns api")
      resp = s.post(
        f"{NNS_ENDPOINT_BASE}/insert",
        params=_get_params(collection),
        data=_line_gen(df),
        headers=_get_header(),
        timeout=None,
      )
    dt = (time.time() - t0) * 1000
    log(f"done total_ms={dt:.0f}. Body: {resp.text[:1000]}")
    try:
      start_build_index(collection, nlist=nlist, rebuild=False, wait=False)
      st = build_index_status(collection)
      logger.info(
        f"index build triggered collection={collection} exists={st.get('exists')} "
        f"building={st.get('is_building')} finished={st.get('is_finished')}"
      )
    except requests.RequestException as e:
      logger.error(f"index build trigger failed collection={collection}: {e}")
    return resp
  except requests.RequestException as e:
    logger.error(f"approx NN search failed: {e}")
    return None


def get_apprx(inputs, collection, fallback_search=None):
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
    logger.error(f"approx NN search failed: {e}")
    if callable(fallback_search):
      return fallback_search(inputs, collection)
    return []


def group_inputs(wanted, index, bloom=None, force=False):
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
