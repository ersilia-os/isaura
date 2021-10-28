# imports
from pathlib import Path
import os

# system variables
ISAURA = os.path.join(str(Path.home()), "eos", "isaura")
if not os.path.exists(ISAURA):
    os.makedirs(ISAURA)

# internal variables

## logging
LOGGING_FILE = "console.log"
SILENCE_FILE = ".silence.json"
VERBOSE_FILE = ".verbose.json"

## filesystem

REPOSITORY_PATH = os.path.join(ISAURA, "lake")
if not os.path.exists(REPOSITORY_PATH):
    os.makedirs(REPOSITORY_PATH)

HDF5_EXTENSION = "h5"

VALUES = "Values"
KEYS = "Keys"
FEATURES = "Features"
