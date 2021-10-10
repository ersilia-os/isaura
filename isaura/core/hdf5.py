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
            
    @property
    def size(self):
        return os.path.getsize(self.data_path)


class Hdf5ApiExplorer(Hdf5):

    def __init__(self, model_id, api_name):
        Hdf5.__init__(self, model_id=model_id)
        self.set_curr_api(api_name)

    def api_exists(self):
        return self._check_api_exists(self.api_name)

    @property
    def shape(self):
        with h5py.File(self.data_path, "rb") as f:
            return f[self.api_name][data].shape

    @property
    def dtype(self):
        with h5py.File(self.data_path, "rb") as f:
            return f[self.api_name][data].dtype

if __name__ == "__main__":  #TESTING
    h = Hdf5Explorer("eos4e40")
    h.set_curr_api("Predict")
    h.write_api(["e","d"], [[90, np.nan, 70], [61, 51, 41]])
    for i in h.read_api():
        print(i)