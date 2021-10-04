from ..core.base import IsauraBase
from ..handle.reader import Reader
import h5py
import numpy as np

class Mapper(IsauraBase):

    def __init__(self, model_id, api_name):
        IsauraBase.__init__(self, model_id)
        self.reader = Reader(model_id)
        self.api_name = api_name

    def check_keys(self, new_keys):
        keys_set, keys_list = set(new_keys), list(new_keys)
        curr_api_keys = set(self.reader._get_keys(self.api_name))
        avail_keys, unavail_keys = {}, {}

        for k in keys_set:
            if k not in curr_api_keys:
                avail_keys[k] = keys_list.index(k)
        unavail_keys = keys_set.difference(avail_keys)
        result = {
            "available_keys": avail_keys,
            "unavailable_keys": unavail_keys
        }

        return result

    
