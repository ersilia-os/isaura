import numpy as np


class Ranges(object):
    def __init__(self):
        self.uints, self.ints, self.floats = self.data_lists()

    def data_lists(self):
        uints = (self.uint8(), self.uint16(), self.uint32(), self.uint64())
        ints = (self.int8(), self.int16(), self.int32(), self.int64())
        floats = (self.float16(), self.float32(), self.float64())
        return uints, ints, floats

    def uint8(self):
        dtype = np.uint8
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max - 1)

    def uint16(self):
        dtype = np.uint16
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max - 1)

    def uint32(self):
        dtype = np.uint32
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max - 1)

    def uint64(self):
        dtype = np.uint64
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max - 1)

    def int8(self):
        dtype = np.int8
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max - 1)

    def int16(self):
        dtype = np.int16
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max - 1)

    def int32(self):
        dtype = np.int32
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max - 1)

    def int64(self):
        dtype = np.int64
        info = np.iinfo(dtype)
        return (dtype, info.min, info.max - 1)

    def float16(self):
        dtype = np.float16
        info = np.finfo(dtype)
        precision = info.precision
        return (dtype, info.min, info.max - 1, precision)

    def float32(self):
        dtype = np.float32
        info = np.finfo(dtype)
        precision = info.precision
        return (dtype, info.min, info.max - 1, precision)

    def float64(self):
        dtype = np.float64
        info = np.finfo(dtype)
        precision = info.precision
        return (dtype, info.min, info.max - 1, precision)

    def is_in_range(self, type, min_value, max_value):
        if min_value < type[1]:
            return False
        if max_value > type[2]:
            return False
        return True

    def smallest_dtype(
        self, type_list, min_value, max_value
    ):  # Determine smallest applicable container
        dtype = type_list[-1][0]  # Default to largest size
        for d in type_list:
            if self.is_in_range(d, min_value, max_value):
                dtype = d[0]
                break
        return dtype

    def best_int(
        self, min_value, max_value
    ):  # Split into data types then make call for sizing
        if min_value > 0:
            dtype = self.smallest_dtype(self.uints, min_value, max_value)
        else:
            dtype = self.smallest_dtype(self.ints, min_value, max_value)
        return dtype

    def best_float(self, min_value, max_value, max_decimals):
        float_list = self.floats
        dtype = float_list[-1][0]  # Default to largest size
        for f in float_list:
            if (
                self.is_in_range(f, min_value, max_value) and max_decimals < f[3]
            ):  # Check range and precision
                dtype = f[0]
                break
        return dtype

    def max_of_type(self, data):
        for dt in self.floats:
            if data.dtype == dt[0]:
                return dt[2]


class NumericDataTyper(object):
    def __init__(self, ary):
        self.ary_flat = ary.ravel()
        self.min = np.nanmin(ary)
        self.max = np.nanmax(ary)
        self.is_integer = self._is_integer()
        self._ranges = Ranges()
        self._dtype_orig = ary[0].dtype

    def _is_integer(self):
        for x in self.ary_flat:
            if not np.isnan(x) and not np.equal(
                np.mod(x, 1), 0
            ):  # Not NaN and 0 modulo
                return False
        return True

    def _count_decimals(self):
        most_decimals = 0
        for x in self.ary_flat:
            x_str = str(x)
            decimals = len(x_str) - x_str.index(".") - 1
            if decimals > most_decimals:
                most_decimals = decimals
        return most_decimals

    def actual(self):
        return self._dtype_orig

    def best(self):
        if self.is_integer:
            return self._ranges.best_int(self.min, self.max)
        else:
            return self._ranges.best_float(self.min, self.max, self._count_decimals())

    def data_size(self, dt):  # Byte size of data type
        return np.dtype(dt).itemsize


class TypeMask(object):
    def __init__(self, dtype):
        self.dtype = dtype
        try:
            info = np.iinfo(self.dtype)
        except:
            info = np.finfo(self.dtype)
        self.mask = info.max
