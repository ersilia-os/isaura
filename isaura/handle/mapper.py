from ..core.base import IsauraBase
from ..handle.reader import Reader
import h5py
import numpy as np


class Mapper(IsauraBase):
    def __init__(self, model_id):
        IsauraBase.__init__(self, model_id)

    def check_keys(
        self, api_name, new_keys
    ):  # Change to display which file keys are in?
        curr_api_keys = set()
        for file in self.avail_data_files():
            r = Reader(file, self.model_id)
            keys = r._get_keys(api_name)
            curr_api_keys.update(keys)
        return self._filter(curr_api_keys, new_keys)

    def filter_file(self, file_path, api_name, new_keys):
        r = Reader(file_path, self.model_id)
        curr_api_keys = r._get_keys(api_name)
        return self._filter(curr_api_keys, new_keys)

    def _filter(self, curr_api_keys, new_keys):
        avail_keys, unavail_keys = {}, {}
        for i, k in enumerate(new_keys):
            if k in curr_api_keys:
                avail_keys[k] = i
            else:
                unavail_keys[k] = i
        result = {"available_keys": avail_keys, "unavailable_keys": unavail_keys}
        return result
