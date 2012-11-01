'''Raw interface to c extension
'''

from __future__ import division
from __future__ import absolute_import

import ctypes

import numpy as np

from .cutil import *
from .cinf import maybe_wire_dylib, select_lib
from .tempfile import temporary_file
from .fingerprint import PairOrientationFingerPrint


LIB_NAME = 'libpointpot'

maybe_wire_dylib(LIB_NAME)

clib = deflib(select_lib(LIB_NAME),
              symbol_translator = lambda s: 'POINTPOT_' + s,
              defs=[
    [c_void_p, 'COUL_load', [c_void_p]],
    [c_void, 'COUL_free', [c_void_p]],
    [c_double, 'BASE_calculate_potential', [c_void_p,
                                            c_double_p, c_double_p,
                                            c_double_p, c_double_p]],
    [c_void, 'BASE_calculate_forces_torques', [c_void_p,
                                               c_double_p, c_double_p,
                                               c_double_p, c_double_p,
                                               c_double_p, c_double_p,
                                               c_double_p, c_double_p]],
    [c_void, 'BASE_calculate_point_charge_set_sizes', [c_void_p,
                                                       c_int_p, c_int_p,
                                                       c_double_p, c_double_p,
                                                       c_double_p, c_double_p]],

    ])


class CCoulombPointPotential(object):

    c_ptr = None

    def __init__(self, c_ptr):
        self.c_ptr = c_ptr

    def __del__(self):
        c_ptr = self.c_ptr
        self.c_ptr = None
        if c_ptr is not None:
            clib.COUL_free(c_ptr)

    @classmethod
    def from_binary_file(cls, fp):
        return cls(clib.COUL_load(ctypes.pythonapi.PyFile_AsFile(fp)))

    @classmethod
    def from_binary_filepath(cls, path):
        with open(path) as fp:
            return cls.from_binary_file(fp)

    @classmethod
    def from_coulomb_point_potential(cls, cpp, truncate=False):
        with temporary_file(suffix='bin') as path:
            cpp.write_binary_filepath(path, truncate=truncate)
            return cls.from_binary_filepath(path)

    def evaluate_potential(self, r, omega, theta_i, theta_j):
        fp = PairOrientationFingerPrint(r=r, omega=omega, theta_i=theta_i, theta_j=theta_j)
        return self.evaluate_pair_potential(fp.create_pair_with_orientation())

    def evaluate_pair_potential(self, pair):
        holder, pos_i, rot_i, pos_j, rot_j = self.unpack_pair(pair)
        return clib.BASE_calculate_potential(self.c_ptr, pos_i, rot_i, pos_j, rot_j)

    def ex_evaluate_pair_potential(self, pos_i, rot_i, pos_j, rot_j):
        holder, pos_i, rot_i, pos_j, rot_j = self.setup_pair(pos_i, rot_i, pos_j, rot_j)
        return clib.BASE_calculate_potential(self.c_ptr, pos_i, rot_i, pos_j, rot_j)

    def evaluate_pair_forces_and_torques(self, pair):
        holder, pos_i, rot_i, pos_j, rot_j = self.unpack_pair(pair)
        [[force_i, torque_i], [force_j, torque_j]] = res = np.empty((2,2,3), dtype=c_double)
        clib.BASE_calculate_forces_torques(self.c_ptr,
                                           force_i.ctypes.data_as(c_double_p),
                                           torque_i.ctypes.data_as(c_double_p),
                                           force_j.ctypes.data_as(c_double_p),
                                           torque_j.ctypes.data_as(c_double_p),
                                           pos_i, rot_i, pos_j, rot_j)
        return res

    def ex_evaluate_pair_forces_and_torques(self, pos_i, rot_i, pos_j, rot_j):
        holder, pos_i, rot_i, pos_j, rot_j = self.setup_pair(pos_i, rot_i, pos_j, rot_j)
        [[force_i, torque_i], [force_j, torque_j]] = res = np.empty((2,2,3), dtype=c_double)
        clib.BASE_calculate_forces_torques(self.c_ptr,
                                           force_i.ctypes.data_as(c_double_p),
                                           torque_i.ctypes.data_as(c_double_p),
                                           force_j.ctypes.data_as(c_double_p),
                                           torque_j.ctypes.data_as(c_double_p),
                                           pos_i, rot_i, pos_j, rot_j)
        return res

    def calculate_charge_set_sizes(self, r, omega, theta_i, theta_j):
        fp = PairOrientationFingerPrint(r=r, omega=omega, theta_i=theta_i, theta_j=theta_j)
        return self.calculate_pair_point_charge_set_sizes(fp.create_pair_with_orientation())

    def calculate_pair_point_charge_set_sizes(self, pair):
        i = ctypes.c_int(0)
        j = ctypes.c_int(0)
        holder, pos_i, rot_i, pos_j, rot_j = self.unpack_pair(pair)
        clib.BASE_calculate_point_charge_set_sizes(self.c_ptr, ctypes.byref(i), ctypes.byref(j),
                                                   pos_i, rot_i, pos_j, rot_j)
        return i.value, j.value

    @classmethod
    def unpack_pair(cls, pair):
        pos_i = np.asarray(pair.p_i.position, dtype=c_double, order='C')
        pos_j = np.asarray(pair.p_j.position, dtype=c_double, order='C')
        assert pos_i.shape == (3,)
        assert pos_j.shape == (3,)

        rot_i = np.asarray(pair.p_i.rotation, dtype=c_double, order='C')
        rot_j = np.asarray(pair.p_j.rotation, dtype=c_double, order='C')
        assert rot_i.shape == (3,3)
        assert rot_j.shape == (3,3)

        return cls.setup_pair(pos_i, rot_i, pos_j, rot_j)

    @staticmethod
    def setup_pair(pos_i, rot_i, pos_j, rot_j):
        return [[pos_i, rot_i, pos_j, rot_j],
                pos_i.ctypes.data_as(c_double_p),
                rot_i.ctypes.data_as(c_double_p),
                pos_j.ctypes.data_as(c_double_p),
                rot_j.ctypes.data_as(c_double_p)]


