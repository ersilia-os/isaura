from ..core.base import IsauraBase
from..dtypes.numeric import Ranges
import h5py
import numpy as np

class Reader(IsauraBase):

    def __init__(self, data_path, model_id):
        IsauraBase.__init__(self, model_id)
        self.ranges = Ranges()
        self.data_path = data_path

    def read_by_idx(self, api_name, idxs):
        with h5py.File(self.data_path, "r") as f:
            for i in idxs:      #Optimise this with batching?
                yield self._decode(f.get(api_name)["Values"][i])

    def read_by_key(self, api_name, keys):
        with h5py.File(self.data_path, "r") as f:
            key_index_dict = self._index_keys(api_name)
            for k in keys:      #Optimise this with batching?
                if k in key_index_dict.keys():
                    pos = key_index_dict[k]
                    yield self._decode(f.get(api_name)["Values"][pos])

    def yield_api(self, api_name):
        if self._check_api_exists(self.data_path, api_name):
            with h5py.File(self.data_path, "r") as f:
                for key, data in zip(f.get(api_name)["Keys"].asstr(), f.get(api_name)["Values"]):
                    yield key, self._decode(data)
        return False

    def get_apis(self):
        with h5py.File(self.data_path, "r") as f:
            return list(f.keys())

    def _get_keys(self, api_name):
        keys = []
        with h5py.File(self.data_path, "r") as f:
            if self._check_api_exists(self.data_path, api_name):
                keys = list(f.get(api_name)["Keys"].asstr())
        return keys

    def _index_keys(self, api_name):
        indices = {}
        with h5py.File(self.data_path, "r") as f:
            if self._check_api_exists(self.data_path, api_name):
                keys = f.get(api_name)["Keys"].asstr()
                indices = {k:i for i,k in enumerate(keys)}
        return indices

    def _decode(self, data):
        d = data
        for index, element in enumerate(d):
            if element == self.ranges.max_of_type(element):
                d[index] = np.nan
        return d