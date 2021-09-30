from ..core.base import IsauraBase
from isaura.dtypes.numeric import NumericDataTyper
import h5py
import numpy as np


class Writer(IsauraBase):

    def __init__(self, model_id):
        IsauraBase.__init__(self, model_id)
        self._create_h5(model_id)
        #self.model = ErsiliaModel(model_id)
        #self.model.serve()

    def write(self, api_name, keys, input):
        new_keys = np.array(keys)
        new_values = list(input)

        #TO DO manage empty values here
        #TO DO filter out values that already exist

        dtypes = ''
        for i in np.transpose(new_values):  #TO DO check all Nan cols
            dtyper = NumericDataTyper(i)
            dtypes += str(np.dtype(dtyper.best())) + ","
        dtypes = dtypes[:-1]

        #Loop through new datatypes against current h5 file datatypes
        #Upcast data if necessary

        new_values = [tuple(x) for x in new_values]
        np_arr = np.array(new_values, dtype=np.dtype(dtypes))

        if self._check_api_exists(api_name):
            self._append_api(api_name, new_keys, np_arr)
        else:
            self._write_new_api(api_name, new_keys, np_arr)


    def _append_api(self, api_name, keys, values):
        with h5py.File(self.data_path, "a") as f:
            grp = f.get(api_name)
            grp["Keys"].resize((grp["Keys"].shape[0] + keys.shape[0]), axis=0)
            grp["Keys"][-keys.shape[0]:] = keys

            grp["Values"].resize((grp["Values"].shape[0] + values.shape[0]), axis=0)
            grp["Values"][-values.shape[0]:] = values

    def _write_new_api(self, api_name, keys, values):
        with h5py.File(self.data_path, "a") as f:
            grp = f.create_group(api_name)
            grp.create_dataset("Keys", shape=keys.shape, data=keys, maxshape=(None,), chunks=True)
            grp.create_dataset("Values", shape=values.shape, data=values, maxshape=(None,), chunks=True)

    def _create_h5(self, model_id):
        if not self._check_h5_exists():
            f = h5py.File(self.data_path, "w")
            f.close()