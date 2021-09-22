import numpy as np


class Ranges(object):
    def __init__(self):
        self.uints, self.ints, self.floats = self.data_lists()
     
    def data_lists(self):
        uints = ( self.uint8(), self.uint16(), self.uint32(), self.uint64() )
        ints = ( self.int8(), self.int16(), self.int32(), self.int64() )
        floats = ( self.float16(), self.float32(), self.float64() )
        return uints, ints, floats

    def uint8(self):
        dtype = np.uint8
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max-3)

    def uint16(self):
        dtype = np.uint16
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max-3)

    def uint32(self):
        dtype = np.uint32
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max-3)

    def uint64(self):
        dtype = np.uint64
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max-3)

    def int8(self):
        dtype = np.int8
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max-3)

    def int16(self):
        dtype = np.int16
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max-3)

    def int32(self):
        dtype = np.int32
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max-3)

    def int64(self):
        dtype = np.int64
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max-3)

    def float16(self):
        dtype = np.float16
        info = np.finfo(dtype)
        return (dtype, info.min, info.max-3)

    def float32(self):
        dtype = np.float32
        info = np.finfo(dtype)
        return (dtype, info.min, info.max-3)

    def float64(self):
        dtype = np.float64
        info = np.finfo(dtype)
        return (dtype, info.min, info.max-3)

    def is_in_range(self, type, min_value, max_value):
        if min_value < type[1]:
            return False
        if max_value > type[2]:
            return False
        return True

    def smallest_dtype(self, data_type, min_value, max_value):  #Determine smallest applicable container
        dtype = data_type[-1][0]    #Default to largest size
        for d in data_type:
            if self.is_in_range(d, min_value, max_value):
                dtype = d[0]
                break
        return dtype

    def best(self, is_integer, min_value, max_value):   #Split into data types then make call for sizing
        if is_integer:
            if min_value > 0:   #Integer and only positive
                dtype = self.smallest_dtype(self.uints, min_value, max_value)
            else:               #Integer and negative values
                dtype = self.smallest_dtype(self.ints, min_value, max_value)
        else:                   #Floats
            dtype = self.smallest_dtype(self.floats, min_value, max_value)
        return dtype





class NumericDataTyper(object):

    def __init__(self, ary):
        self.ary_flat = ary.ravel()
        self.min = np.min(ary)
        self.max = np.max(ary)
        self.is_integer = self._is_integer()
        self._ranges = Ranges()
        self._dtype_orig = np.dtype(ary[0])

    def _is_integer(self):
        for x in self.ary_flat:
            if not np.equal(np.mod(x, 1), 0):
            #if not isinstance(x, int):
                return False
        return True

    def actual(self):
        return self._dtype_orig

    def best(self):
        return self._ranges.best(self.is_integer, self.min, self.max)

    def data_size(self, dt):    #Byte size of data type
        return np.dtype(dt).itemsize
