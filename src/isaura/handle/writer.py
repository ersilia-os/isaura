from ..core.base import IsauraBase
from ..handle.mapper import Mapper

# from isaura.dtypes.numeric import NumericDataTyper
from ..dtypes.numeric import TypeMask
import h5py
import numpy as np


class Writer(IsauraBase):
    def __init__(self, model_id, path=""):
        IsauraBase.__init__(self, model_id)
        if path == "":
            self.path = self.local_data_path
        else:
            self.path = path

        self._create_h5(self.path)
        self.set_dtype(np.float32)

    def _unique_keys(self, keys):
        _, idxs = np.unique(keys, return_index=True)
        return idxs

    def write(self, api_name, keys, input):
        arr_keys, arr_values = self._format_encode(keys, input)
        new_keys, new_values = self._filter_keys(api_name, arr_keys, arr_values)
        if self._check_api_exists(self.path, api_name):
            if new_keys.shape != (0,):
                self._append_api(api_name, new_keys, new_values)
        elif new_keys.shape[0] > 0:
            self._write_new_api(api_name, new_keys, new_values)

    def write_append(self, api_name, keys, input):
        arr_keys, arr_values = self._format_encode(keys, input)
        new_keys, new_values = self._filter_single_file(
            self.path, api_name, arr_keys, arr_values
        )
        if self._check_api_exists(self.path, api_name):
            if new_keys.shape != (0,):
                self._append_api(api_name, new_keys, new_values)
        elif new_keys.shape[0] > 0:
            self._write_new_api(api_name, new_keys, new_values)

    def _append_api(self, api_name, keys, values):
        with h5py.File(self.path, "a") as f:
            grp = f.get(api_name)
            grp["Keys"].resize((grp["Keys"].shape[0] + keys.shape[0]), axis=0)
            grp["Keys"][-keys.shape[0] :] = keys

            grp["Values"].resize((grp["Values"].shape[0] + values.shape[0]), axis=0)
            grp["Values"][-values.shape[0] :] = values

    def _write_features_api(self, api_name, features):
        with h5py.File(self.path, "a") as f:
            grp = f.get(api_name)
            keys = grp.keys()
            if "Features" not in keys:
                arr_features = np.array(list(features), h5py.string_dtype())
                grp.create_dataset("Features", data=arr_features)

    def _write_new_api(self, api_name, keys, values):
        with h5py.File(self.path, "a") as f:
            grp = f.create_group(api_name)
            grp.create_dataset(
                "Keys", shape=keys.shape, data=keys, maxshape=(None,), chunks=True
            )
            grp.create_dataset(
                "Values",
                shape=values.shape,
                data=values,
                maxshape=(None, values.shape[1]),
                chunks=True,
                dtype=self.dtype,
            )

    def _create_h5(self, path):
        if not self._check_h5_exists(path):
            f = h5py.File(path, "w")
            f.close()

    def _format_encode(self, keys, input):
        arr_keys = np.array(list(keys), h5py.string_dtype())
        arr_values = np.array(list(input))

        idxs = self._unique_keys(arr_keys)
        arr_keys = arr_keys[idxs]
        arr_values = arr_values[idxs]

        if len(arr_values.shape) > 1:
            for record in arr_values:
                self._encode(record, self.dtype)
        else:
            self._encode(arr_values, self.dtype)
        return arr_keys, arr_values

    def _filter_keys(self, api_name, arr_keys, arr_values):
        m = Mapper(self.model_id)
        new_keys, new_values = [], []
        filter = m.check_keys(api_name, arr_keys)["unavailable_keys"]
        for k, v in zip(filter.keys(), filter.values()):
            if k is not None:
                new_keys.append(k)
                new_values.append(arr_values[v])

        return np.array(new_keys, h5py.string_dtype()), np.array(new_values)

    def _filter_single_file(self, file_path, api_name, arr_keys, arr_values):
        m = Mapper(self.model_id)
        new_keys, new_values = [], []
        filter = m.filter_file(file_path, api_name, arr_keys)["unavailable_keys"]
        for k, v in zip(filter.keys(), filter.values()):
            if k is not None:
                new_keys.append(k)
                new_values.append(arr_values[v])

        return np.array(new_keys, h5py.string_dtype()), np.array(new_values)

    def _encode(self, data, dtype):
        for index, element in enumerate(data):
            if np.isnan(element):
                data[index] = self.mask_value

    def _get_dtype_info(self, dtype):
        dtype_char = np.sctype2char(dtype)
        if dtype_char == "f":
            return np.finfo(dtype)
        elif dtype_char == "i":
            return np.iinfo(dtype)

    def change_api_name(self, curr_api, new_api):
        with h5py.File(self.path, "a") as f:
            f.move(curr_api, new_api)

    def remove_api(self, api_name):
        with h5py.File(self.path, "a") as f:
            del f[api_name]

    def set_dtype(self, new_dtype):
        self.dtype = new_dtype
        self.mask_value = TypeMask(new_dtype).mask
