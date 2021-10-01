import numpy as np
from isaura.handle.reader import Reader
from isaura.handle.writer import Writer

class Data(object):

    def __init__(self, keys, inps, txts, vals):
        self.keys = keys
        self.inps = inps
        self.txts = txts
        self.vals = vals


class Hdf5(object):


    def __init__(self, model_id): #Pass through model object here instead?
        self.model_id = model_id

    def write_api(self, keys, iter):
        pass

    def _filter_existing(self):
        pass

    def _resolve_dtype_clash(self, curr_h5, new_data):
        pass

    def insert(self, data):
        pass

    def query(self):
        pass

if __name__ == "__main__":  #TESTING
    w = Writer("eos4e40")
    w.write("Predict", iter(["c","b"]), iter([[10, 20, 30], [1.98455484, 120, -130.2]]))

    r = Reader("eos4e40")
    for key, record in r.yield_api("Predict"):
        print(str(key), record)