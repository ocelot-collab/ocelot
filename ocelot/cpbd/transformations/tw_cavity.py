from copy import copy

import numpy as np

from ocelot.common.globals import speed_of_light
from ocelot.cpbd.transformations.transfer_map import TransferMap


class TWCavityTM(TransferMap):
    def __init__(self, l=0, v=0, phi=0, freq=0):
        TransferMap.__init__(self)
        self.length = l
        self.dx = 0
        self.dy = 0
        self.tilt = 0
        self.v = v
        self.phi = phi
        self.freq = freq
        self.delta_e_z = lambda z: self.v * np.cos(self.phi * np.pi / 180.) * z / self.length
        self.delta_e = self.v * np.cos(self.phi * np.pi / 180.)
        self.R_z = lambda z, energy: np.dot(
            self.tw_cavity_R_z(z, self.v * z / self.length, energy, self.freq, self.phi),
            self.f_entrance(z, self.v * z / self.length, energy, self.phi))
        self.R = lambda energy: np.dot(self.f_exit(self.length, self.v, energy, self.phi),
                                       self.R_z(self.length, energy))

    def tw_cavity_R_z(self, z, V, E, freq, phi=0.):
        """
        :param z: length
        :param de: delta E
        :param f: frequency
        :param E: initial energy
        :return: matrix
        """
        phi = phi * np.pi / 180.
        de = V * np.cos(phi)
        r12 = z * E / de * np.log(1. + de / E) if de != 0 else z
        r22 = E / (E + de)
        r65 = V * np.sin(phi) / (E + de) * (2 * np.pi / (speed_of_light / freq)) if freq != 0 else 0
        r66 = r22
        cav_matrix = np.array([[1, r12, 0., 0., 0., 0.],
                               [0, r22, 0., 0., 0., 0.],
                               [0., 0., 1, r12, 0., 0.],
                               [0., 0., 0, r22, 0., 0.],
                               [0., 0., 0., 0., 1., 0],
                               [0., 0., 0., 0., r65, r66]]).real
        return cav_matrix

    def f_entrance(self, z, V, E, phi=0.):
        phi = phi * np.pi / 180.
        de = V * np.cos(phi)
        r = np.eye(6)
        r[1, 0] = -de / z / 2. / E
        r[3, 2] = r[1, 0]
        return r

    def f_exit(self, z, V, E, phi=0.):
        phi = phi * np.pi / 180.
        de = V * np.cos(phi)
        r = np.eye(6)
        r[1, 0] = +de / z / 2. / (E + de)
        r[3, 2] = r[1, 0]
        return r

    @classmethod
    def create_from_element(cls, element, params=None):
        return cls(v=element.v, freq=element.freq, phi=element.phi)
