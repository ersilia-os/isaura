import h5py


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
        self._check_h5_exists(module_id, self.path)
        self.api_exists = self._check_api_exists(self.path, api)
        if self.api_exists:
            self.key_indices = self._get_indices_by_key(self.path, api)


    #def write_api(self, api, iter):
        #with h5py.File() as f:

    #def read_api(self, api):
        #with h5py.File(self.path, "r") as f:


    # Change to h5py file management
    def _get_indices_by_key(self, path, api):
        with h5py.File(path, "r") as f:
            keys = f.get(api).keys()
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

    def _check_h5_exists(self, module_id, file_path):
        if not os.path.isfile(file_path):
            f = h5py.File(module_id + ".h5", "w")
            f.close() #with

    def _check_api_exists(self, path, api):
        with h5py.File(path, "w") as f:
            if api in f.keys():
                return True
            return False
