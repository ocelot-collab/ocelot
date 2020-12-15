from copy import copy

import numpy as np

from ocelot.cpbd.high_order import m_e_GeV
from ocelot.cpbd.transformations.tm_utils import transform_vec_ent, transform_vec_ext
from ocelot.cpbd.transformations.first_order import TransferMap


class KickTM(TransferMap):
    def __init__(self, angle=0., k1=0., k2=0., k3=0., nkick=1):
        TransferMap.__init__(self)
        self.angle = angle
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.nkick = nkick

    @classmethod
    def create_from_element(cls, element, params=None):
        return cls(angle=element.angle, k1=element.k1, k2=element.k2,
                   k3=element.k3 if hasattr(element, 'k3') else 0.,
                   nkick=params['nkick'] if 'nkick' in params else 1)

    def kick(self, X, l, angle, k1, k2, k3, energy, nkick=1):
        """
        does not work for dipole
        """
        gamma = energy / m_e_GeV
        coef = 0.
        beta = 1.
        if gamma != 0:
            gamma2 = gamma * gamma
            beta = 1. - 0.5 / gamma2
            coef = 1. / (beta * beta * gamma2)
        l = l / nkick
        angle = angle / nkick

        dl = l / 2.
        k1 = k1 * l
        k2 = k2 * l / 2.
        k3 = k3 * l / 6.

        for i in range(nkick):
            x = X[0] + X[1] * dl
            y = X[2] + X[3] * dl

            p = -angle * X[5] + 0j
            xy1 = x + 1j * y
            xy2 = xy1 * xy1
            xy3 = xy2 * xy1
            p += k1 * xy1 + k2 * xy2 + k3 * xy3
            X[1] = X[1] - np.real(p)
            X[3] = X[3] + np.imag(p)
            X[4] = X[4] + np.real(angle * xy1) / beta - X[5] * dl * coef

            X[0] = x + X[1] * dl
            X[2] = y + X[3] * dl
            # X[4] -= X[5] * dl * coef
        return X

    def kick_apply(self, X, l, angle, k1, k2, k3, energy, nkick, dx, dy, tilt):
        if dx != 0 or dy != 0 or tilt != 0:
            X = transform_vec_ent(X, dx, dy, tilt)
        self.kick(X, l, angle, k1, k2, k3, energy, nkick=nkick)
        if dx != 0 or dy != 0 or tilt != 0:
            X = transform_vec_ext(X, dx, dy, tilt)

        return X

    def map_function(self, delta_length=None, length=None):
        return lambda X, energy: self.kick_apply(X, delta_length if delta_length else self.length, self.angle, self.k1, self.k2, self.k3, energy, self.nkick, self.dx, self.dy, self.tilt)

