from ..core.base import IsauraBase
from ..handle.reader import Reader
from ..handle.writer import Writer
from ..handle.mapper import Mapper
from ..handle.retyper import Retyper
from ..handle.appender import Appender
import numpy as np


class Hdf5(IsauraBase):
    def __init__(self, model_id):
        IsauraBase.__init__(self, model_id)
        self.w = Writer(self.model_id)
        self.m = Mapper(self.model_id)
        self.api_name = ""

    # TO DO optimisiation: if keys found in one file then remove them from search in next file

    def set_curr_api(self, api):
        self.api_name = api

    def get_curr_api(self):
        return self.api_name

    def read_api(self):
        for file in self._get_readers():
            for entry in file.yield_api(self.api_name):
                yield entry

    def list_apis(self):
        api_set = set()
        for reader in self._get_readers():
            api_set.update(reader.get_apis())
        return api_set

    def list_keys(self):
        keys = []
        for reader in self._get_readers():
            if reader._check_api_exists(reader.data_path, self.api_name):
                keys.append(reader._get_keys(self.api_name))
        return keys

    def read_by_key(self, key_list):
        results_dict = {}
        for reader in self._get_readers():
            result = reader.read_by_key(self.api_name, key_list)
            keys = result["keys"]
            values = result["values"]
            for i in range(len(keys)):
                results_dict[keys[i]] = values[i]
        for k in key_list:
            yield results_dict[k]

    # TO DO Fix read_by_index with file indexing

    # def read_by_index(self, index_list):
    # for reader in self._get_readers():
    # return reader.read_by_idx(api_name, self.index_list)
    # return False

    def check_keys_exist(self, key_list):
        return self.m.check_keys(self.api_name, key_list)

    def filter_keys(self, keys, values):
        arr_keys = np.array(keys)
        arr_values = np.array(list(values))
        return self.w._filter_keys(self.api_name, keys, values)

    def write_api(self, keys, values):
        # Batch here with generator input
        self.w.write(self.api_name, keys, values)

    def retype_local(self):
        rt = Retyper(self.local_data_path, self.model_id, self.api_name)
        rt.retype()

    def retype_public(self):
        rt = Retyper(self.public_data_path, self.model_id, self.api_name)
        rt.retype()

    def retype(self, path):
        rt = Retyper(path, self.model_id, self.api_name)
        rt.retype()

    def append(self, append_from_path, append_to_path, secret_keys):
        a = Appender(append_to_path, self.model_id, self.api_name, split=False)
        a.append_from(append_from_path, self.model_id, self.api_name, secret_keys)

    def append_split(self, append_from_path, append_to_path, secret_keys):
        a = Appender(append_to_path, self.model_id, self.api_name, split=True)
        a.append_from(append_from_path, self.model_id, self.api_name, secret_keys)

    def append_excl_file(self, append_from_path, append_to_path, secrets_path):
        r = Reader(secrets_path, self.model_id)
        a = Appender(append_to_path, self.model_id, self.api_name, split=False)
        a.append_from(append_from_path, self.model_id, self.api_name, r._get_keys(self.api_name))

    def get_features(self):
        for reader in self._get_readers():
            return reader._get_features(self.api_name)

    def write_features(self, features):
        self.w._write_features_api(self.api_name, features)

    def _get_readers(self):
        readers = []
        for path in self.avail_data_files():
            readers.append(Reader(path, self.model_id))
        return readers
