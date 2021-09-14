import h5py


class Data(object):

    def __init__(self, keys, inps, txts, vals):
        self.keys = keys
        self.inps = inps
        self.txts = txts
        self.vals = vals


class Hdf5(object):

    def __init__(self):
        pass

    def _get_indices_by_key(self, keys):
        with open(self.h5_file, "rb") as f:


    def _get_indices_by_inp(self, inps):
        with open(self.h5_file, "rb") as f:


    def _get_indices_by_txt(self, txts):
        with open(self.h5_file, "rb") as f:


    def insert(self, data):
        pass

    def query(self):
        pass

    def 
