"""Tools for calculating the internal energy of a pairs of SubSphereCollection's.
"""

from __future__ import division
from __future__ import with_statement

import numpy as np
import transformations

from .bases import AutoRepr
from .util import float_fix
from .linutils import (vx,vy,vz,cartesian_basis,
                      vlen, vnorm, vperp,
                      to_rotation_matrix, rotate_positions, times_column3,
                      signed_angle_between_unit_vecs)

from .particle import Particle, Pair

pi = np.pi


class PairOrientationFingerPrint(AutoRepr):
    '''A fingerprint unqiue to the mutal orientation of a pair of particles that
       is invarient their spatial location as classified by six variables.
       See the associated picture pair_orientation_fingerprint.png to see
       the meaning of each variable as well as associated external documentation.
    '''

    r_epsilon = 1e-12
    angle_epsilon = 0.5 * 2*pi / 360.0

    def __init__(self, r=1, omega=0, theta_i=0, phi_i=0, theta_j=0, phi_j=0):
        #assert r>2*radius
        #assert -pi <= omega <= pi, 'bad omega %r' % (omega,)
        #assert 0 <= theta_i <= pi
        #assert -pi <= phi_i <= pi
        #assert 0 <= theta_j <= pi
        #assert -pi <= phi_j <= pi
        self.r = r
        self.omega = omega
        self.theta_i = theta_i
        self.phi_i = phi_i
        self.theta_j = theta_j
        self.phi_j = phi_j

    def repr_args(self):
        return tuple(float_fix(op,4) for op in (self.r, self.omega, self.theta_i, self.phi_i, self.theta_j, self.phi_j))

    def to_integer_tuple(self):
        return tuple(int(round(f/epsilon)) for f,epsilon in [
            (self.r, self.r_epsilon),
            (self.omega, self.angle_epsilon),
            (self.theta_i, self.angle_epsilon),
            (self.phi_i, self.angle_epsilon) ,
            (self.theta_j, self.angle_epsilon),
            (self.phi_j, self.angle_epsilon)])

    def __hash__(self):
        return reduce(lambda a,b: a^b, map(hash, self.to_integer_tuple()))

    def __eq__(self, other):
        if not isinstance(other, PairOrientationFingerPrint):
            return NotImplemented
        return self.to_integer_tuple() == other.to_integer_tuple()

    def to_string_key(self):
        return ','.join(map(str, self.to_integer_tuple()))

    @classmethod
    def from_pair(cls, pair):
        r_ij = pair.p_j.position - pair.p_i.position
        r_ij_unit = vnorm(r_ij)

        orient_vec_i, phi_i = pair.p_i.calculate_orientation(x_axis=r_ij_unit)
        orient_vec_j, phi_j = pair.p_j.calculate_orientation(x_axis=-r_ij_unit)

        theta_i = np.arccos(np.dot(r_ij_unit, orient_vec_i))
        theta_j = np.arccos(np.dot(-r_ij_unit, orient_vec_j))

        n_i = vperp(r_ij_unit, orient_vec_i)
        n_j = vperp(-r_ij_unit, orient_vec_j)
        omega = signed_angle_between_unit_vecs(n_i, n_j, r_ij_unit)

        return cls(vlen(r_ij), omega, theta_i, phi_i, theta_j, phi_j)

    def create_pair_with_orientation(self):
        # pair is constructed with p_i at the origin and p_j translated along
        # the x-axis by distance of r.  z serves as the orientation reference
        # axis.

        def setup_particle_orientation(invert, theta, phi):
            # spin about z-axis vector to create phi spin
            phi_mat = transformations.rotation_matrix(phi + (pi if invert else 0), vz)
            # rotate z-x axis to create orientation vector
            axis_rotation = transformations.rotation_matrix((-1 if invert else 1) * -(theta - pi/2), vy)

            return transformations.concatenate_matrices(axis_rotation, phi_mat)

        mat_p_i = setup_particle_orientation(0, self.theta_i, self.phi_i)
        mat_p_j = setup_particle_orientation(1, self.theta_j, self.phi_j)

        # rotate p_j about separation axis (x-axis) create target omega
        omega_mat = transformations.rotation_matrix(self.omega, (1, 0, 0))
        mat_p_j = transformations.concatenate_matrices(omega_mat, mat_p_j)

        return Pair(Particle((0,0,0),      mat_p_i),
                    Particle((self.r,0,0), mat_p_j))


# # # # # # # # # # # # #
# Potential Calculators #
# # # # # # # # # # # # #

class PotentialCalculator(object):

    def calculate_energy(self, pair, col_i, col_j=None):
        return self.calculate_energy_ex(pair.p_i.to_particle_collection(col_i),
                                        pair.p_j.to_particle_collection(col_j or col_i))


class CoulombPotentialCalculator(PotentialCalculator):

    def __init__(self, scale=1.0, cutoff=1e30):
        self.scale = scale
        self.cutoff = cutoff

    def calculate_energy_ex(self, col_i, col_j):
        return self.scale * cpairenergy.calculate_energy_coulomb(
            col_i, col_j, self.cutoff)


class InversePotentialCalculator(PotentialCalculator):

    def __init__(self, exponent=1, cutoff=None):
        self.exponent = exponent
        self.cutoff = cutoff

    def get_key(self):
        return 'inverse%d%s' % (self.exponent, self.cutoff)

    def calculate_energy_ex(self, col_i, col_j):
        return cpairenergy.inverse_potential_calculate_energy(
                                   col_i, col_j,
                                   exponent=self.exponent, cutoff=self.cutoff)


class YakawaPotentialCalculator(PotentialCalculator):

    def __init__(self, kappa=0.01**-1):
        self.kappa = kappa

    def get_key(self):
        return 'yakawa%.4g' % (self.kappa,)

    def calculate_energy_ex(self, col_i, col_j):
        return cpairenergy.calculate_energy_yakawa(col_i, col_j, self.kappa)


class YukawaCutoffPotentialCalculator(PotentialCalculator):

    def __init__(self, A=1, kappa=5e-9**-1, cutoff=30e-9):
        self.A = A
        self.kappa = kappa
        self.cutoff = cutoff

    def get_key(self):
        return 'yA%.4gkappa%.4gc%.4g' % (self.A, self.kappa, self.cutoff)

    def calculate_energy_ex(self, col_i, col_j):
        return self.A * cpairenergy.calculate_energy_yukawa_cutoff(col_i, col_j, self.kappa, self.cutoff)


class SpunYukawaCutoffPotentialCalculator(YukawaCutoffPotentialCalculator):

    def __init__(self, n_spins=6, **kwds):
        super(SpunYukawaCutoffPotentialCalculator, self).__init__(**kwds)
        self.n_spins = n_spins

    def get_key(self):
        return '%ss%d' % (super(SpunYukawaCutoffPotentialCalculator, self).get_key(),
                          self.n_spins)

    def calculate_energy_ex(self, col_i, col_j):
        raise RuntimeError("shouldn't call this method")

    def calculate_energy(self, pair, col_i, col_j=None):
        fp = PairOrientationFingerPrint.from_pair(pair)
        return self.calculate_spun_energy(fp.r, fp.omega, fp.theta_i, fp.theta_j, col_i, col_j)

    def calculate_spun_energy(self, r, omega, theta_i, theta_j, col_i, col_j=None):

        col_j = col_j or col_i

        spins = np.pi * np.linspace(-1.0, 1.0, num=self.n_spins, endpoint=False)

        theta_i = self.normalize_angle(theta_i)
        theta_j = self.normalize_angle(theta_j)

        acc_spun_energies = []
        for spin_i in spins:
            for spin_j in spins:
                spun_fp = PairOrientationFingerPrint(r=r, omega=omega,
                                                     theta_i=theta_i, theta_j=theta_j,
                                                     phi_i=spin_i, phi_j=spin_j)
                spun_pair = spun_fp.create_pair_with_orientation()
                spun_energy = cpairenergy.calculate_energy_yukawa_cutoff(
                                  spun_pair.p_i.to_particle_collection(col_i),
                                  spun_pair.p_j.to_particle_collection(col_j),
                                  self.kappa, self.cutoff)
                acc_spun_energies.append(spun_energy)

        return self.A * np.array(acc_spun_energies).mean()

    @staticmethod
    def normalize_angle(theta):
        theta = abs(theta) % (2*pi)
        if theta >= pi:
            theta = pi - (theta-pi)
        if theta > 0.99 * pi:
            theta = 0.98 * pi
        assert 0 <= theta < pi, 'bad angle %s' % (theta,)
        return theta


class YukawaLJPotentialCalculator(PotentialCalculator):

    def __init__(self, A_like=1, A_unlike=None, kappa=0.01**-1, yukawa_cutoff=None,
                       sigma=1, epsilon=1, lj_cutoff=None):
        A_unlike = A_like if A_unlike is None else A_unlike
        vars(self).update(locals())
        del self.self

    def get_key(self):
        return 'yukawalj' + ','.join('None' if el is None else '%.4g' % el for el in
                                     (self.A_like, self.A_unlike, self.kappa, self.yukawa_cutoff,
                                      self.sigma, self.epsilon, self.lj_cutoff))

    def calculate_energy_ex(self, col_i, col_j):
        return cpairenergy.calculate_energy_yakawa_lj(
                              col_i, col_j,
                              self.A_like, self.A_unlike, self.kappa,
                              self.yukawa_cutoff,
                              self.sigma, self.epsilon, self.lj_cutoff)


