from ..core.base import IsauraBase
from .reader import Reader
from .writer import Writer
from .retyper import Retyper
from .repacker import Repacker
from ..dtypes.numeric import NumericDataTyper
from ..default import HDF5_EXTENSION
import numpy as np
import os


class Appender(IsauraBase):
    def __init__(
        self, curr_path, curr_model_id, curr_api, split=True
    ):  # File intended to add to
        IsauraBase.__init__(self, model_id=curr_model_id)
        self.curr_api = curr_api
        self.curr_path = curr_path
        self.curr_reader = Reader(self.curr_path, curr_model_id)
        self.curr_dtype = self.curr_reader.get_dtype(self.curr_api)
        self.split = split
        self.is_empty = False

        if type(self.curr_dtype) == type(None):  # If file empty, choose smallest dtype
            self.curr_dtype = np.uint8
            self.is_empty = True

    def append_from(
        self, path, model_id, api_name, exclude_list
    ):  # File intended to add from
        new_data_reader = Reader(path, model_id)
        keys, values = new_data_reader.read_api(api_name)
        new_dtype = self._check_new_dtype(values)

        if type(new_dtype) != type(None):  # If there is data to append
            (
                filtered_keys,
                filtered_values,
                secret_keys,
                secret_values,
            ) = self.filter_secret_keys(keys, values, exclude_list)

            if self.split == False:

                if not self.is_empty and self._dtype_clash(
                    self.curr_dtype, new_dtype
                ):  # If append_to file has data and incoming data is incompatible with that file
                    self.retype_data(new_dtype)
                w = Writer(self.model_id, path=self.curr_path)
                w.set_dtype(new_dtype)
                w.write_append(self.curr_api, filtered_keys, filtered_values)

                if self.is_empty:
                    self.append_features_list(new_data_reader)

            else:
                w = Writer(
                    self.model_id, path=self.curr_path[:-3] + "_add." + HDF5_EXTENSION
                )
                w.set_dtype(new_dtype)
                filtered_pub_keys, filtered_pub_values = w._filter_single_file(
                    self.curr_path,
                    self.curr_api,
                    list(filtered_keys),
                    list(filtered_values),
                )
                w.write_append(self.curr_api, filtered_pub_keys, filtered_pub_values)

    def append_features_list(self, reader):
        features = reader._get_features(self.curr_api)
        w = Writer(self.model_id, path=self.curr_path)
        w._write_features_api(self.curr_api, features)

    def _check_new_dtype(self, data):
        if type(data) != type(None):
            dt = NumericDataTyper(data)
            return dt.best()
        return None

    def _dtype_clash(self, curr_dtype, new_dtype):
        final_dtype = self._find_new_dtype(curr_dtype, new_dtype)
        if final_dtype == curr_dtype:
            return False
        return True

    def _find_new_dtype(self, dtype1, dtype2):
        return np.promote_types(dtype1, dtype2)

    def retype_data(self, new_dtype):
        keys, values = self.curr_reader.read_api(self.curr_api)
        rt = Retyper(self.curr_path, self.model_id, self.curr_api)
        rt.recast(keys, values, new_dtype)

    def filter_secret_keys(self, new_keys, new_values, secret_keys):
        secret_set = set(secret_keys)
        to_append, keep_secret = {}, {}
        for k, v in zip(new_keys, new_values):
            if k not in secret_set:
                to_append[k] = v
            else:
                keep_secret[k] = v
        return (
            to_append.keys(),
            to_append.values(),
            keep_secret.keys(),
            keep_secret.values(),
        )
