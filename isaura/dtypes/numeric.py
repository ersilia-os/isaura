import numpy as np


class Ranges(object):

    def __init__(self):
        pass

    def int8(self):
        dtype = np.int8
        info = np.iinfo(dtype)
        return (0, dtype, info.min, info.max)

    def int16(self):
        dtype = np.int16
        info = np.iinfo(dtype)
        return (1, dtype, info.min, info.max)

    def int32(self):
        dtype = np.int32
        info = np.iinfo(dtype)
        return (2, dtype, info.min, info.max)

    def int64(self):
        dtype = np.int64
        info = np.iinfo(dtype)
        return (3, dtype, info.min, info.max)

    def float8(self):
        dtype = np.float8
        info = np.iinfo(dtype)
        return (4, dtype, info.min, info.max)

    def float16(self):
        dtype = np.float16
        info = np.iinfo(dtype)
        return (5, dtype, info.min, info.max)

    def float32(self):
        dtype = np.float32
        info = np.iinfo(dtype)
        return (6, dtype, info.min, info.max)

    def float64(self):
        dtype = np.float64
        info = np.info(dtype)
        return (7, dtype, info.min, info.max)

    def is_in_range(self, func, min_value, max_value):
        if min_value < x:
            return False
        if max_value > x:
            return False
        return True

    def best(self, is_integer, min_value, max_value):
        if is_integer:

        else:




class NumericDataTyper(object):

    def __init__(self, ary):
        self.ary_flat = ary.ravel()
        self.min = np.min
        self.max = np.max
        self.is_integer = self._is_integer()
        self._ranges = Ranges()
        self._dtype_orig = np.dtype(ary)

    def _is_integer(self):
        for x in self.ary_flat:
            if not x.is_integer():
                return False
        return True

    def actual(self):
        return self._dtype_orig

    def best(self):
        return self._ranges.best(self.is_integer, self.min, self.max)
