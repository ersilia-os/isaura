from ..core.base import IsauraBase
from ..handle.mapper import Mapper
from isaura.dtypes.numeric import NumericDataTyper
import h5py
import numpy as np


class Writer(IsauraBase):

    def __init__(self, model_id):
        IsauraBase.__init__(self, model_id)
        self._create_h5(model_id)

    def write(self, api_name, keys, input):
        arr_keys = np.array(list(keys), h5py.string_dtype())     #dtype='<S27')   #InChi Keys??
        arr_values = np.array(list(input))

        #TO DO manage empty values here
        #TO DO manage batching/large datasets before here


        ### Storage Optimisation Code - Implement Later ###
        #new_values = list(input)
        #dtypes = ''
        #for i in np.transpose(new_values):  #TO DO check all Nan cols
            #dtyper = NumericDataTyper(i)
            #dtypes += str(np.dtype(dtyper.best())) + ","
        #dtypes = dtypes[:-1]

        #ew_values = [tuple(x) for x in new_values]
        #np_arr = np.array(new_values, dtype=np.dtype(dtypes))


        if self._check_api_exists(api_name):
            new_keys, new_values = self._filter_keys(api_name, arr_keys, arr_values)
            if new_keys.shape != (0,):
                self._append_api(api_name, new_keys, new_values)
        else:
            self._write_new_api(api_name, arr_keys, arr_values)


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
            grp.create_dataset("Values", shape=values.shape, data=values, maxshape=(None, values.shape[1]), chunks=True, dtype='f')

    def _create_h5(self, model_id):
        if not self._check_h5_exists():
            f = h5py.File(self.data_path, "w")
            f.close()

    def _filter_keys(self, api_name, arr_keys, arr_values):
        m = Mapper(self.model_id, api_name)
        new_keys, new_values = [], []
        filter = m.check_keys(arr_keys)["available_keys"]
        for k, i in zip(filter.keys(), filter.values()):
            if k is not None:
                new_keys.append(k)
                new_values.append(arr_values[i])

        return np.array(new_keys), np.array(new_values)

    def _encode(self):
        pass