import h5py
import os.path
import numpy as np
from isaura.dtypes.numeric import NumericDataTyper

class Data(object):

    def __init__(self, keys, inps, txts, vals):
        self.keys = keys
        self.inps = inps
        self.txts = txts
        self.vals = vals


class Hdf5(object):


    def __init__(self, model_id, file_path, api): #Pass through model object here instead?
        self.model_id = model_id
        self.path = os.path.join(file_path, model_id + ".h5")
        self._check_h5_exists()
        self.api = api
        self.api_exists = self._check_api_exists()
        if self.api_exists:
            self.key_indices = self._get_indices_by_key()


    def write_api(self, keys, iter):
        with h5py.File(self.path, "a") as f:
            new_keys = np.array(keys)
            new_values = list(iter)

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
            print(dtypes)


            if self.api_exists:
                grp = f.get(self.api)
                grp["Keys"].resize((grp["Keys"].shape[0] + new_keys.shape[0]), axis = 0)
                grp["Keys"][-new_keys.shape[0]:] = new_keys

                grp["Values"].resize((grp["Values"].shape[0] + np_arr.shape[0]), axis = 0)
                grp["Values"][-np_arr.shape[0]:] = np_arr
            else:
                grp = f.create_group(self.api)
                grp.create_dataset("Keys", shape=new_keys.shape, data=new_keys, maxshape=(None, ), chunks=True)
                grp.create_dataset("Values", shape=np_arr.shape, data=np_arr, maxshape=(None, ), chunks=True)
                self.api_exists = True

    def read_api(self):
        if self.api_exists:
            with h5py.File(self.path, "r") as f:
                for data in f.get(self.api)["Values"]:
                    yield data
        return False

    def _filter_existing(self):
        pass

    def _resolve_dtype_clash(self, curr_h5, new_data):
        pass

    def _get_indices_by_key(self):
        with h5py.File(self.path, "r") as f:
            keys = f.get(self.api)["Keys"]
            indices = {k:i for i,k in enumerate(keys)}
        return indices


    def _get_indices_by_inp(self, inps):
        #with open(self.h5_file, "rb") as f:
        pass


    def _get_indices_by_txt(self, txts):
        #with open(self.h5_file, "rb") as f:
        pass

    def insert(self, data):
        pass

    def query(self):
        pass

    def _check_h5_exists(self):
        if not os.path.isfile(self.path):
            f = h5py.File(os.path.join(self.path), "w")
            f.close()

    def _check_api_exists(self):
        with h5py.File(self.path, "r") as f:
            if self.api in f.keys():
                return True
            return False

if __name__ == "__main__":  #TESTING
    h = Hdf5("eos4e40", "/home/jason/", "Predict")
    h.write_api([0,1], iter([[10, 20, 30], [1.98455484, 120, -130.2]]))
    for i in h.read_api():
        print(i, np.dtype(i))