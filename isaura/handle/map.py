

class Mapper(IsauraBase):

    def __init__(self, model_id):
        IsauraBase.__init__(self, model_id)
        self.reader = Hdf5Reader(self.data_path)

    def map_keys(self, keys):
        keys_set = set(keys)
        avail_keys = {}
        for i, key in enumerate(self.reader):
            if key in keys_set:
                avail_keys[key] = i
        unavail_keys = keys_set.difference(avail_keys.keys())
        result = {
            "available_keys": avail_keys,
            "unavailable_keys": unavail_keys
        }
        return result

    
