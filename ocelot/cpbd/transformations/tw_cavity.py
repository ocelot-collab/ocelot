from copy import copy

import numpy as np

from ocelot.cpbd.transformations.transfer_map import TransferMap, TMTypes
from ocelot.cpbd.elements.element import Element
from ocelot.common.globals import speed_of_light, m_e_GeV

class TWCavityTM(TransferMap):
    """[summary]
    Implementation of the Cavity Transformation.
    The concrete element atom have to implement:
    create_cavity_tm_main_params(self, energy: float, delta_length: float) -> CavityParams
    create_cavity_tm_entrance_params(self, energy: float, delta_length: float) -> CavityParams
    create_cavity_tm_exit_params(self, energy: float, delta_length: float) -> CavityParams
    """

    def __init__(self, create_tm_param_func, delta_e_func, tm_type: TMTypes, length: float,
                 delta_length: float) -> None:
        super().__init__(create_tm_param_func, delta_e_func, tm_type, length, delta_length=delta_length)

    @classmethod
    def from_element(cls, element: Element, tm_type: TMTypes = TMTypes.MAIN, delta_l=None, **params):
        return cls.create(entrance_tm_params_func=element.create_cavity_tm_entrance_params,
                          delta_e_func=element.create_delta_e,
                          main_tm_params_func=element.create_cavity_tm_main_params,
                          exit_tm_params_func=element.create_cavity_tm_exit_params,
                          tm_type=tm_type, length=element.l, delta_length=delta_l)

    def map4cav(self, X, E, delta_length, length):
        params = self.get_params(E)
        if delta_length is not None:
            V = params.v * delta_length / length if length != 0 else params.v
            z = delta_length
        else:
            V = params.v
            z = length

        beta0 = 1
        igamma2 = 0
        g0 = 1e10
        if E != 0:
            g0 = E / m_e_GeV
            igamma2 = 1. / (g0 * g0)
            beta0 = np.sqrt(1. - igamma2)

        phi = params.phi * np.pi / 180.

        X4 = np.copy(X[4])
        X5 = np.copy(X[5])
        X = self.mul_p_array(X, energy=E)  # t_apply(R, T, X, dx, dy, tilt)
        delta_e = V * np.cos(phi)

        T566 = 1.5 * z * igamma2 / (beta0 ** 3)
        T556 = 0.
        T555 = 0.
        if E + delta_e > 0:
            k = 2. * np.pi * params.freq / speed_of_light
            E1 = E + delta_e
            g1 = E1 / m_e_GeV
            beta1 = np.sqrt(1. - 1. / (g1 * g1))

            X[5] = X5 * E * beta0 / (E1 * beta1) + V * beta0 / (E1 * beta1) * (
                    np.cos(-X4 * beta0 * k + phi) - np.cos(phi))

            dgamma = V / m_e_GeV
            if delta_e > 0:
                T566 = z * (beta0 ** 3 * g0 ** 3 - beta1 ** 3 * g1 ** 3) / (
                        2 * beta0 * beta1 ** 3 * g0 * (g0 - g1) * g1 ** 3)
                T556 = beta0 * k * z * dgamma * g0 * (beta1 ** 3 * g1 ** 3 + beta0 * (g0 - g1 ** 3)) * np.sin(phi) / (
                        beta1 ** 3 * g1 ** 3 * (g0 - g1) ** 2)
                T555 = beta0 ** 2 * k ** 2 * z * dgamma / 2. * (
                        dgamma * (2 * g0 * g1 ** 3 * (beta0 * beta1 ** 3 - 1) + g0 ** 2 + 3 * g1 ** 2 - 2) / (
                        beta1 ** 3 * g1 ** 3 * (g0 - g1) ** 3) * np.sin(phi) ** 2 -
                        (g1 * g0 * (beta1 * beta0 - 1) + 1) / (beta1 * g1 * (g0 - g1) ** 2) * np.cos(phi))
        X[4] += T566 * X5 * X5 + T556 * X4 * X5 + T555 * X4 * X4

        return X

    def map_function(self, X, energy: float):
        if self.tm_type == TMTypes.MAIN:
            return self.map4cav(X, energy, self.delta_length, self.length)
        else:
            return self.mul_p_array(X, energy=energy)


class TWCavityTM_old(TransferMap):
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
