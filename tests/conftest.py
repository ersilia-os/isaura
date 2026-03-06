import os
import pytest

MODEL = os.getenv("ISAURA_TEST_MODEL", "eos3b5e")
VERSION = os.getenv("ISAURA_TEST_VERSION", "v1")
SRC = os.getenv("ISAURA_TEST_BUCKET", "isaura-test")
INPUT = os.getenv("ISAURA_TEST_INPUT", "tests/data.csv")


@pytest.fixture(scope="session")
def cfg():
  return {"model": MODEL, "version": VERSION, "src": SRC, "input": INPUT}


@pytest.fixture(scope="session")
def state():
  return {}
