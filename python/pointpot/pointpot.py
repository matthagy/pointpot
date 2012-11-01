'''
'''

import numpy as np

from .linutils import vlen
from .fingerprint import PairOrientationFingerPrint
from ..binio import BinaryWriter


class CoulombPointPotential(object):

    def __init__(self, r_colloid, r_min, r_max, base_colloid, potential_scale=1.0):
        self.r_colloid = r_colloid
        self.r_min = r_min
        self.r_max = r_max
        self.base_colloid = base_colloid
        self.potential_scale = potential_scale
        self.r_cutoff = self.r_max - 2*self.r_colloid
        assert self.r_cutoff >= 0.0

    def scale(self, factor):
        return self.__class__(self.r_colloid, self.r_min, self.r_max,
                              self.base_colloid,
                              self.potential_scale * factor)

    def evaluate_potential(self, r, omega, theta_i, theta_j):
        fp = PairOrientationFingerPrint(r=r, omega=omega, theta_i=theta_i, theta_j=theta_j)
        return self.evaluate_pair_potential(fp.create_pair_with_orientation())

    def evaluate_pair_potential(self, pair):
        return self.potential_scale * calculate_energy_coulomb(
            pair.p_i.to_particle_collection(self.base_colloid),
            pair.p_j.to_particle_collection(self.base_colloid),
            self.r_cutoff)

    def evaluate_pair_forces_and_torques(self, pair):
        if not self.check_r_range(vlen(pair.p_j.position - pair.p_i.position)):
            return np.zeros((2,2,3))

        return (self.potential_scale *
                np.array([self.evaluate_1_force_torque(pair.p_i, pair.p_j),
                          self.evaluate_1_force_torque(pair.p_j, pair.p_i)]))

    def evaluate_1_force_torque(self, p_center, p_other):
        return calculate_forces_coulomb(p_center.to_particle_collection(self.base_colloid),
                                        p_other.to_particle_collection(self.base_colloid),
                                        self.r_cutoff)

    def check_r_range(self, r):
        if r < self.r_min:
            raise ValueError("bad r=%.4g < r_min=%.4g" % (r, self.r_min))
        return r <= self.r_max

    binary_magic = 0x734234

    def write_binary_filepath(self, filepath, **kwds):
        with open(filepath, 'w') as fp:
            self.write_binary_file(fp, **kwds)

    def write_binary_file(self, fp, truncate=False):
        bw = BinaryWriter(fp)
        bw.write_int(self.binary_magic)
        bw.write_int(self.base_colloid.N)
        bw.write_array(self.base_colloid.locations.astype(np.float64))

        assert ((self.base_colloid.types == 'A') |
                (self.base_colloid.types == 'B')).all()

        types = np.where(self.base_colloid.types == 'A', 1, -1).astype(np.int32)
        bw.write_array(types)

        bw.write_int(int(truncate))
        bw.write_double(self.r_colloid)
        bw.write_double(self.r_cutoff)
        bw.write_double(self.potential_scale)

        bw.write_int(self.binary_magic)


def calculate_forces_coulomb(p1, p2, r_cutoff):
    raise RuntimeError('use cpointpot for force evaluation')

def calculate_energy_coulomb(p1, p2, r_cutoff):
    raise RuntimeError('use cpointpot for energy evaluation')


