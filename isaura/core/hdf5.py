from isaura.core.base import IsauraBase
from isaura.handle.io import Hdf5
import h5py
import os
import numpy as np


class Hdf5Explorer(Hdf5):
    def __init__(self, model_id):
        Hdf5.__init__(self, model_id=model_id)

    def local_h5_exists(self):
        return self._check_h5_exists(self.local_data_path)

    def public_h5_exists(self):
        return self._check_h5_exists(self.public_data_path)

    def apis(self):
        return self.list_apis()


class Hdf5ApiExplorer(Hdf5):
    def __init__(self, model_id, api_name):
        Hdf5.__init__(self, model_id=model_id)
        self.set_curr_api(api_name)

    def api_exists(self):
        return self._check_api_exists(self.api_name)

    @property
    def dtype(self):
        for data_path in self.avail_data_files():
            with h5py.File(data_path, "r") as f:
                g = f.get(self.api_name)
                return g["Values"].dtype
