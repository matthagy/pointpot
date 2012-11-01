'''
'''

from __future__ import division
from __future__ import absolute_import

import numpy as np

from .bases import AutoRepr
from .util import float_fix
from .linutils import (vx,vy,vz,cartesian_basis, vperp,
                       to_rotation_matrix, rotate_positions, times_column3,
                       signed_angle_between_unit_vecs)


class Particle(AutoRepr):
    '''Represents the position and orientation of a particle
    '''

    def __init__(self, position=None, rotation=None):
        self.position = np.array(position) if position is not None else np.zeros(3)
        self.rotation = to_rotation_matrix(rotation)

    def repr_args(self):
        return [tuple(float_fix(el,4) for el in self.position),
                list([float_fix(el,4) for el in row]
                     for row in np.array(self.rotation))]

    def copy(self):
        return self.__class__(self.position, self.rotation)

    def calculate_orientation(self, x_axis=vx):
        '''Rotate basis into coordinate frame of particle.
           Use vz in this frame as a vector characteristic of the particles orientation.
           Additionally calculate a spin about this axis by comparing vx in this frame
           to an x_axis.
        '''
        cin_vx, cin_vy, cin_vz = cartesian_in_rotation = rotate_positions(self.rotation, cartesian_basis)
        vec = cin_vz
        n_x = vperp(vec, x_axis)
        phi = signed_angle_between_unit_vecs(n_x, cin_vx, vec)
        return vec, phi

    def calculate_orientation_vector(self):
        return times_column3(self.rotation, vz)

    def to_particle_collection(self, col):
        return col.rotated_by_matrix(self.rotation).translated(self.position)


class Pair(AutoRepr):
    '''Represents the state of a pair of particles
    '''

    def __init__(self, p_i, p_j):
        self.p_i = p_i
        self.p_j = p_j

    def repr_args(self):
        return self.p_i, self.p_j

