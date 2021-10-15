from ..core.base import IsauraBase
from ..dtypes.numeric import Ranges
import h5py
import numpy as np


class Reader(IsauraBase):
    def __init__(self, data_path, model_id):
        IsauraBase.__init__(self, model_id)
        self.ranges = Ranges()
        self.data_path = data_path

    def read_by_idx(self, api_name, idxs):
        idxs = np.array(idxs)
        mapping = np.argsort(idxs)
        mapping_inv = np.argsort(mapping)
        idxs_sorted = idxs[mapping]
        with h5py.File(self.data_path, "r") as f:
            grp = f.get(api_name)
            values = self._decode(grp["Values"][idxs_sorted])
            values = values[mapping_inv]
            return values

    def read_by_key(self, api_name, keys):
        key_index_dict = self._index_keys(api_name)
        found_keys = []
        idxs = []
        for k in keys:
            if k in key_index_dict:
                idxs += [key_index_dict[k]]
                found_keys += [k]
        values = self.read_by_idx(api_name, idxs)
        result = {"keys": found_keys, "values": values}
        return result

    def _get_features(self, api_name):
        with h5py.File(self.data_path, "r") as f:
            features = list(f.get(api_name)["Features"].asstr())
        return features

    def yield_api(self, api_name):
        if self._check_api_exists(self.data_path, api_name):
            with h5py.File(self.data_path, "r") as f:
                for key, data in zip(
                    f.get(api_name)["Keys"].asstr(), f.get(api_name)["Values"]
                ):
                    yield key, self._decode(data)
        return False

    def read_api(self, api_name):
        if self._check_api_exists(self.data_path, api_name):
            with h5py.File(self.data_path, "r") as f:
                keys, values = (
                    np.array(f.get(api_name)["Keys"].asstr()),
                    np.array(f.get(api_name)["Values"]),
                )
                return keys, self._decode(values)
        return None, None

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
                indices = {k: i for i, k in enumerate(keys)}
        return indices

    def _decode(self, data):
        max_of_type = self.ranges.max_of_type(data)
        data = data.astype(np.float32)
        data[data == max_of_type] = np.nan
        return data

    def get_dtype(self, api_name):
        if self._check_api_exists(self.data_path, api_name):
            with h5py.File(self.data_path, "r") as f:
                value = f.get(api_name)["Values"][0]
                return value.dtype.type
        return None
