from isaura.core.base import IsauraBase
from isaura.handle.io import Hdf5
import h5py
import os
import numpy as np


class Hdf5Explorer(Hdf5):
    def __init__(self, model_id):
        Hdf5.__init__(self, model_id=model_id)
        self.model_id = model_id

    def local_h5_exists(self):
        return self._check_h5_exists(self.local_data_path)

    def public_h5_exists(self):
        return self._check_h5_exists(self.public_data_path)

    def apis(self):
        return self.list_apis()

    def info(self):
        i = {}
        for api_name in self.apis():
            api = Hdf5ApiExplorer(model_id=self.model_id, api_name=api_name)
            i[api_name] = api.info()
        return i


class Hdf5ApiExplorer(Hdf5):
    def __init__(self, model_id, api_name):
        Hdf5.__init__(self, model_id=model_id)
        self.set_curr_api(api_name)

    def api_exists(self):
        return self._check_api_exists(self.api_name)

    def info(self, head=10):
        i = {}
        for data_path in self.avail_data_files():
            fn = filename = os.path.basename(data_path)
            with h5py.File(data_path, "r") as f:
                g = f.get(self.api_name)
                if g is None:
                    res = None
                else:
                    res = {
                        "shape": g["Values"].shape,
                        "keys": g["Keys"].shape,
                        "features": g["Features"].shape,
                        "keys_head": [x.decode("utf-8") for x in g["Keys"][:head]],
                        "features_head": [
                            x.decode("utf-8") for x in g["Features"][:head]
                        ],
                        "dtype": g["Values"].dtype.name,
                    }
            i[fn] = res
        return i
