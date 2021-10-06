from isaura.core.base import IsauraBase
from isaura.handle.reader import Reader
from isaura.handle.writer import Writer
from isaura.handle.mapper import Mapper
import numpy as np

class Hdf5(IsauraBase):

    def __init__(self, model_id):
        IsauraBase.__init__(self, model_id)
        self.r = Reader(self.model_id)
        self.w = Writer(self.model_id)
        self.m = Mapper(self.model_id)
        self.api_name = ""

    def set_curr_api(self, api):
        self.api_name = api

    def get_curr_api(self):
        return self.api_name

    def read_api(self):
        return self.r.yield_api(self.api_name)

    def list_apis(self):
        return self.r.get_apis()

    def list_keys(self):
        return self.r._get_keys(self.api_name)

    def read_by_key(self, key_list):
        return self.r.read_by_key(self.api_name, key_list)

    def read_by_index(self, index_list):
        return self.r.read_by_idx(api_name, self.index_list)

    def check_keys_exist(self, api_name, key_list):
        return self.m.check_keys(api_name, key_list)

    def filter_keys(self, keys, values):
        arr_keys = np.array(keys)
        arr_values = np.array(list(values))
        return self.w._filter_keys(self.api_name, keys, values)

    def write_api(self, keys, values):
        self.w.write(self.api_name, keys, values)

    def merge(self, old_model, old_api_name):
        reader = Reader(old_model)
        for key, value in reader.yield_api(old_api_name):   #Inefficient
            self.w.write(self.api_name, key, value)

    def _resolve_dtype_clash(self, curr_h5, new_data):
        pass