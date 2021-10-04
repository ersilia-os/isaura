from ..core.base import IsauraBase
from..dtypes.numeric import Ranges
import h5py
import numpy as np

class Reader(IsauraBase):

    def __init__(self, model_id):
        IsauraBase.__init__(self, model_id)
        self.ranges = Ranges()

    def read_by_idx(self, api_name, idxs):
        pass

    def read_by_key(self, api_name):
        pass

    def yield_api(self, api_name):
        if self._check_api_exists(api_name):
            with h5py.File(self.data_path, "r") as f:
                for key, data in zip(f.get(api_name)["Keys"].asstr(), f.get(api_name)["Values"]):
                    yield key, self._decode(data)
            return True
        return False

    def _get_keys(self, api_name):
        with h5py.File(self.data_path, "r") as f:
            keys = list(f.get(api_name)["Keys"].asstr())
        return keys

    def _decode(self, data):
        d = data
        for index, element in enumerate(d):
            if element == self.ranges.max_of_type(element):
                d[index] = np.nan
        return d