'''Binary IO designed to be compatible with C code
'''

import numpy as np
import struct


STRUCT_INT_CODE = 'i'
assert struct.calcsize(STRUCT_INT_CODE) == 4


class BinaryWriter(object):

    def __init__(self, fp):
        self.fp = fp

    def write(self, frmt, *args):
        self.fp.write(struct.pack(frmt, *args))

    def write_int(self, i):
        self.write(STRUCT_INT_CODE, i)

    def write_double(self, d):
        self.write('d', d)

    def write_array(self, a):
        self.fp.write(a.tostring())


class BinaryReader(object):

    def __init__(self, fp):
        self.fp = fp

    def pull(self, nbytes):
        bytes = self.fp.read(nbytes)
        if len(bytes) != nbytes:
            raise EOFError
        return bytes

    def read(self, frmt):
        op, = struct.unpack(frmt, self.pull(struct.calcsize(frmt)))
        return op

    def read_int(self):
        return self.read(STRUCT_INT_CODE)

    def read_double(self):
        return self.read('d')

    def read_array(self, dtype, shape):
        bytes = self.pull(dtype(1).itemsize * np.multiply.reduce(shape))
        return np.fromstring(bytes, dtype=dtype).reshape(shape)

    def check_magic(self, magic):
        i = self.read_int()
        if (i & 0xfffffff) != (magic & 0xfffffff):
            raise ValueError("bad magic %X; expected %X" % (i, magic))

