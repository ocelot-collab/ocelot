from copy import copy

import numpy as np

from ocelot.cpbd.high_order import m_e_GeV
from ocelot.common.globals import speed_of_light
from ocelot.cpbd.transformations.first_order import TransferMap


class CavityTM(TransferMap):
    def __init__(self, v=0, freq=0., phi=0.):
        TransferMap.__init__(self)
        self.v = v
        self.freq = freq
        self.phi = phi
        self.vx_up = 0.
        self.vy_up = 0.
        self.vx_down = 0.
        self.vy_down = 0.
        phase_term = np.cos(self.phi * np.pi / 180.)
        self.delta_e_z = lambda z: phase_term * self.v * z / self.length if self.length != 0 else phase_term * self.v
        self.delta_e = self.v * np.cos(self.phi * np.pi / 180.)
        self.map = lambda X, energy: self.map4cav(X, energy, self.v, self.freq, self.phi, self.length)

    def map4cav(self, X, E, V, freq, phi, z=0):
        beta0 = 1
        igamma2 = 0
        g0 = 1e10
        if E != 0:
            g0 = E / m_e_GeV
            igamma2 = 1. / (g0 * g0)
            beta0 = np.sqrt(1. - igamma2)

        phi = phi * np.pi / 180.

        X4 = np.copy(X[4])
        X5 = np.copy(X[5])
        X = self.mul_p_array(X, energy=E)  # t_apply(R, T, X, dx, dy, tilt)
        delta_e = V * np.cos(phi)

        T566 = 1.5 * z * igamma2 / (beta0 ** 3)
        T556 = 0.
        T555 = 0.
        if E + delta_e > 0:
            k = 2. * np.pi * freq / speed_of_light
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

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.delta_e = m.delta_e_z(s)
        v = m.v * s / self.length if self.length != 0 else m.v
        m.map = lambda X, energy: m.map4cav(X, energy, v, m.freq, m.phi, s)
        return m

    @classmethod
    def create_from_element(cls, element, params=None):
        return cls(v=element.v, freq=element.freq, phi=element.phi)
