from .base import IsauraBase
import h5py
import os


class Hdf5Explorer(IsauraBase):

    def __init__(self, model_id):
        IsauraBase.__init__(self, model_id=model_id)

    def apis(self):
        pass

    def size(self):
        return os.path.getsize(self.data_path)


class Hdf5ApiExplorer(IsauraBase):

    def __init__(self, model_id, api_name):
        IsauraBase.__init__(self, model_id=model_id)
        self.api_name = api_name

    def exists(self):
        pass

    def shape(self):
        with h5py.File(self.data_path, "rb") as f:
            return f[self.api_name][data].shape

    def dtype(self):
        with h5py.File(self.data_path, "rb") as f:
            return f[self.api_name][data].dtype
