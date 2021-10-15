import os
import h5py
from ..default import REPOSITORY_PATH
from ..default import HDF5_EXTENSION
from .. import logger
import subprocess


class IsauraBase(object):
    def __init__(self, model_id, verbose=True):
        self.logger = logger
        if verbose:
            self.logger.set_verbosity(1)
        else:
            self.logger.set_verbosity(0)
        self.model_id = model_id
        self.repository_path = REPOSITORY_PATH

    @property
    def local_file_name(self):
        return self.model_id + "_local." + HDF5_EXTENSION

    @property
    def public_file_name(self):
        return self.model_id + "_public." + HDF5_EXTENSION

    @property
    def local_data_path(self):
        return os.path.join(self.repository_path, self.local_file_name)

    @property
    def public_data_path(self):
        return os.path.join(self.repository_path, self.public_file_name)

    def avail_data_files(self):
        data_files = []
        if self._check_h5_exists(self.local_data_path):
            data_files.append(self.local_data_path)
        if self._check_h5_exists(self.public_data_path):
            data_files.append(self.public_data_path)
        return data_files

    def _check_h5_exists(self, data_path):
        return os.path.isfile(data_path)

    def _check_api_exists(self, data_path, api_name):
        with h5py.File(data_path, "r") as f:
            if api_name in f.keys():
                return True
            return False

    def run_cmd(self, cmd):
        subprocess.Popen(cmd, shell=True, env=os.environ).wait()
