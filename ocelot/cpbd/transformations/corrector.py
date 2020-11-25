from copy import copy
import logging

import numpy as np

from ocelot.cpbd.transformations.second_order import SecondTM

_logger = logging.getLogger(__name__)


class CorrectorTM(SecondTM):
    def __init__(self, angle_x=0., angle_y=0., r_z_no_tilt=None, t_mat_z_e=None):
        SecondTM.__init__(self, r_z_no_tilt, t_mat_z_e)
        self.angle_x = angle_x
        self.angle_y = angle_y
        self.map = lambda X, energy: self.kick(X, self.length, self.length, self.angle_x, self.angle_y, energy)
        self.B_z = lambda z, energy: self.kick_b(z, self.length, angle_x, angle_y)

    def kick_b(self, z, l, angle_x, angle_y):
        if l == 0:
            hx = 0.
            hy = 0.
        else:
            hx = angle_x / l
            hy = angle_y / l

        dx = hx * z * z / 2.
        dy = hy * z * z / 2.
        dx1 = hx * z if l != 0 else angle_x
        dy1 = hy * z if l != 0 else angle_y
        b = np.array([[dx], [dx1], [dy], [dy1], [0.], [0.]])
        return b

    def kick(self, X, z, l, angle_x, angle_y, energy):
        _logger.debug('invoking kick_b')
        b = self.kick_b(z, l, angle_x, angle_y)
        if self.multiplication is not None and self.t_mat_z_e is not None:
            X1 = self.t_apply(R=self.R(energy), T=self.t_mat_z_e(z, energy), X=X, dx=0, dy=0, tilt=0, U5666=0.)
        else:
            X1 = np.dot(self.R(energy), X)
        X1 = np.add(X1, b)
        X[:] = X1[:]
        return X

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.delta_e = m.delta_e_z(s)
        m.map = lambda X, energy: m.kick(X, s, self.length, m.angle_x, m.angle_y, energy)
        return m
