from ..core.base import IsauraBase
from .reader import Reader
from .writer import Writer
from ..dtypes.numeric import NumericDataTyper
from ..default import HDF5_EXTENSION
import numpy as np
import h5py
import os


class Retyper(IsauraBase):
    def __init__(self, data_path, model_id, api_name):
        IsauraBase.__init__(self, model_id=model_id)
        self.reader = Reader(data_path, model_id)
        self.data_path = data_path
        self.api_name = api_name
        self.suffix = "_old"

    def retype(self):
        curr_dtype = self.reader.get_dtype(self.api_name)
        keys, data = self.reader.read_api(self.api_name)

        if type(data) != type(None):
            new_dtype = self._best_dtype(data)
            if new_dtype != curr_dtype:
                self.recast(keys, data, new_dtype)
                self._clean(self.data_path[:-3] + "_backup.h5")

    def retype_with_dtype(self, new_dtype):
        curr_dtype = self.reader.get_dtype(self.api_name)
        keys, data = self.reader.read_api(self.api_name)

        if type(data) != type(None):
            if new_dtype != curr_dtype:
                self.recast(keys, data, new_dtype)
                self._clean(self.data_path[:-3] + "_backup.h5")

    def recast(self, keys, data, new_dtype):
        self.run_cmd("cp " + self.data_path + " " + self._backup_path())
        writer = Writer(self.model_id, path=self.data_path)
        writer.change_api_name(self.api_name, self.api_name + self.suffix)
        writer.set_dtype(new_dtype)
        writer.write(self.api_name, keys, data)
        writer.remove_api(self.api_name + self.suffix)

        os.rename(self.data_path, self._temp_path())
        self.run_cmd("h5repack " + self._temp_path() + " " + self.data_path)

    def _clean(self, old_path):
        try:
            old_r = Reader(old_path, self.model_id)
            new_r = Reader(self.data_path, self.model_id)
            old_keys, old_values = old_r.read_api(self.api_name)
            new_keys, new_values = new_r.read_api(self.api_name)
            if (
                old_keys.shape == new_keys.shape
                and old_values.shape == new_values.shape
            ):
                os.remove(self._temp_path())
                os.remove(self._backup_path())

        except:
            os.remove(self.data_path)
            os.remove(self._temp_path())
            self.run_cmd("cp " + self._backup_path() + " " + self.data_path)

    def restore_backup(self):
        pass

    def _best_dtype(self, data):
        dtyper = NumericDataTyper(data)
        return dtyper.best()

    def _resolve_dtype_clash(self, curr_h5, new_data):
        pass

    def _backup_path(self):
        return self.data_path[:-3] + "_backup." + HDF5_EXTENSION

    def _temp_path(self):
        return self.data_path[:-3] + self.suffix + "." + HDF5_EXTENSION
