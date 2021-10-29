from ..core.base import IsauraBase
from .writer import Writer
from ..default import HDF5_EXTENSION
import os


class Repacker(IsauraBase):
    def __init__(self, data_path, model_id):
        IsauraBase.__init__(self, model_id)
        self.data_path = data_path
        self.suffix = "_old"

    def repack_data(self, api_name, keys, data):
        self.run_cmd("cp " + self.data_path + " " + self._backup_path(self.data_path))
        writer = Writer(self.model_id, path=self.data_path)
        writer.change_api_name(api_name, api_name + self.suffix)
        if len(keys) > 0:
            writer.write_append(api_name, keys, data)
        writer.remove_api(api_name + self.suffix)
        self.repack()

    def repack(self):
        os.rename(self.data_path, self._temp_path(self.data_path))
        self.run_cmd(
            "h5repack " + self._temp_path(self.data_path) + " " + self.data_path
        )
        os.remove(self._temp_path(self.data_path))
        os.remove(self._backup_path(self.data_path))

    def _backup_path(self, data_path):
        return data_path[:-3] + "_backup." + HDF5_EXTENSION

    def _temp_path(self, data_path):
        return data_path[:-3] + self.suffix + "." + HDF5_EXTENSION
