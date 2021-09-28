# imports
from pathlib import Path
import os

# system variables
ISAURA = os.path.join(str(Path.home()), "isaura")
if not os.path.exists(ISAURA):
    os.makedirs(ISAURA)

# internal variables

## logging
LOGGING_FILE = "console.log"
SILENCE_FILE = ".silence.json"
VERBOSE_FILE = ".verbose.json"

HDF5_EXTENSION = "h5"
