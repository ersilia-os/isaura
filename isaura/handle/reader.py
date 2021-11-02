from ..core.base import IsauraBase
from ..dtypes.numeric import Ranges, TypeMask
import h5py
import numpy as np


class Reader(IsauraBase):
    def __init__(self, data_path, model_id):
        IsauraBase.__init__(self, model_id)
        self.ranges = Ranges()
        self.data_path = data_path

    def _asstr(self, grp, name):
        v = grp[name][:]
        return [x.decode("utf-8") for x in v]

    def read_by_idx(self, api_name, idxs):
        idxs = np.array(idxs)
        mapping = np.argsort(idxs)
        mapping_inv = np.argsort(mapping)
        idxs_sorted = idxs[mapping]
        with h5py.File(self.data_path, "r") as f:
            grp = f.get(api_name)
            if grp is None:
                return None
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
        if values is None:
            return None
        result = {"keys": found_keys, "values": values}
        return result

    def _get_features(self, api_name):
        with h5py.File(self.data_path, "r") as f:
            grp = f.get(api_name)
            if grp is None:
                return None
            features = self._asstr(grp, "Features")
        return features

    def read_api(self, api_name):
        if self._check_api_exists(self.data_path, api_name):
            with h5py.File(self.data_path, "r") as f:
                grp = f.get(api_name)
                keys, values = (
                    np.array(self._asstr(grp, "Keys")),
                    np.array(grp["Values"]),
                )
                return keys, self._decode(values)
        return None, None

    def get_apis(self):
        with h5py.File(self.data_path, "r") as f:
            return list(f.keys())

    def _get_keys(self, api_name):
        keys = []
        if self._check_api_exists(self.data_path, api_name):
            with h5py.File(self.data_path, "r") as f:
                grp = f.get(api_name)
                keys = self._asstr(grp, "Keys")
        return keys

    def _index_keys(self, api_name):
        indices = {}
        keys = self._get_keys(api_name)
        indices = {k: i for i, k in enumerate(keys)}
        return indices

    def _decode(self, data):
        mask_value = TypeMask(data.dtype).mask
        data = data.astype(np.float32)
        data[data == mask_value] = np.nan
        return data

    def get_dtype(self, api_name):
        if self._check_api_exists(self.data_path, api_name):
            with h5py.File(self.data_path, "r") as f:
                value = f.get(api_name)["Values"][0]
                return value.dtype.type
        return None
