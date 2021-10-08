from ..core.base import IsauraBase
from ..handle.reader import Reader
import h5py
import numpy as np

class Mapper(IsauraBase):

    def __init__(self, model_id):
        IsauraBase.__init__(self, model_id)
        self.reader = Reader(model_id)

    def check_keys(self, api_name, new_keys):
        curr_api_keys = set(self.reader._get_keys(api_name))
        avail_keys, unavail_keys = {}, {}
        for i,k in enumerate(new_keys):
            if k in curr_api_keys:
                avail_keys[k] = i
            else:
                unavail_keys[k] = i
        result = {
            "available_keys": avail_keys,
            "unavailable_keys": unavail_keys
        }
        return result
