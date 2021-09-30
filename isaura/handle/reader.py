from ..core.base import IsauraBase
import h5py


class Reader(IsauraBase):

    def __init__(self, model_id):
        IsauraBase.__init__(self, model_id)

    def read_by_idx(self, api_name, idxs):
        pass

    def read_by_key(self, api_name):
        pass

    def yield_api(self, api_name):
        if self._check_api_exists(api_name):
            with h5py.File(self.data_path, "r") as f:
                for data in f.get(api_name)["Values"]:
                    yield data
            return True
        return False

    def _get_indices_by_key(self):
        with h5py.File(self.data_path, "r") as f:
            keys = f.get(self.api)["Keys"]
            indices = {k:i for i,k in enumerate(keys)}
        return indices

