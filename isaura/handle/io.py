import numpy as np
from isaura.core.base import IsauraBase
from isaura.handle.reader import Reader
from isaura.handle.writer import Writer
from isaura.handle.mapper import Mapper


class Hdf5(IsauraBase):

    def __init__(self, model_id):
        IsauraBase.__init__(self, model_id)
        self.r = Reader(self.model_id)
        self.w = Writer(self.model_id)
        self.m = Mapper(self.model_id)

    def read_api(self, api_name):
        return self.r.yield_api(api_name)

    def list_apis(self):
        return self.r.get_apis()

    def list_keys(self, api_name):
        return self.r._get_keys(api_name)

    def read_by_key(self, api_name, key_list):
        return self.r.read_by_key(api_name, key_list)

    def read_by_index(self, api_name, index_list):
        return self.r.read_by_idx(api_name, index_list)

    def check_exists(self, api_name, key_list):
        return self.m.check_keys(api_name, key_list)

    def filter_keys(self, api_name, keys, values):
        arr_keys = np.array(keys)
        arr_values = np.array(list(values))
        return self.w._filter_keys(api_name, keys, values)

    def write_api(self, api_name, keys, values):
        self.w.write(api_name, keys, values)

    def merge(self, old_model, old_api_name, curr_api):
        reader = Reader(old_model)
        for key, value in reader.yield_api(old_api_name):   #Inefficient
            self.w.write(curr_api, key, value)

    def _resolve_dtype_clash(self, curr_h5, new_data):
        pass

if __name__ == "__main__":  #TESTING
    h = Hdf5("eos4e40")
    h.write_api("Predict", ["y","z"], [[90, np.nan, 70], [61, 51, 41]])
    h.merge("eos4e40_1", "Predict", "Predict")
    for line in h.read_api("Predict"):
        print(line)
