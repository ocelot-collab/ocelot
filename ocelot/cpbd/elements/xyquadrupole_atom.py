import numpy as np

from ocelot.cpbd.high_order import m_e_GeV
from ocelot.cpbd.r_matrix import rot_mtx
from ocelot.cpbd.elements.magnet import Magnet
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.tm_params.second_order_params import SecondOrderParams
from ocelot.cpbd.high_order import t_nnn


class XYQuadrupoleAtom(Magnet):
    """
    quadrupole
    l - length of lens in [m],
    k1 - strength of quadrupole lens in [1/m^2],
    k2 - strength of sextupole lens in [1/m^3],
    tilt - tilt of lens in [rad].
    """

    def __init__(self, l=0., x_offs=0., y_offs=0., k1=0., tilt=0., eid=None):
        super().__init__(eid=eid)
        self.l = l
        self.k1 = k1
        self.x_offs = x_offs
        self.y_offs = y_offs
        self.tilt = tilt

    def __str__(self):
        s = 'XYQuadrupole('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'k1=%8.6e, ' % self.k1 if np.abs(self.k1) > 1e-15 else ""
        s += 'x_offs=%8.6e, ' % self.x_offs if np.abs(self.x_offs) > 1e-15 else ""
        s += 'y_offs=%8.6e, ' % self.y_offs if np.abs(self.y_offs) > 1e-15 else ""
        s += 'tilt=%8.6e, ' % self.tilt if np.abs(self.tilt) > 1e-15 else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s

    def create_r_matrix(self):
        k1 = self.k1

        if self.l == 0:
            hx = 0.
            hy = 0.
        else:
            hx = k1 * self.x_offs
            hy = -k1 * self.y_offs

        def r_mtx(z, k1, hx, hy, sum_tilts=0., energy=0.):
            # r = self.l/self.angle
            #  +K - focusing lens , -K - defoc
            gamma = energy / m_e_GeV

            kx2 = (k1 + hx * hx)
            ky2 = hy * hy - k1
            kx = np.sqrt(kx2 + 0.j)
            ky = np.sqrt(ky2 + 0.j)
            cx = np.cos(z * kx).real
            cy = np.cos(z * ky).real
            sy = (np.sin(ky * z) / ky).real if ky != 0 else z

            igamma2 = 0.

            if gamma != 0:
                igamma2 = 1. / (gamma * gamma)

            beta = np.sqrt(1. - igamma2)

            if kx != 0:
                sx = (np.sin(kx * z) / kx).real
                dx = hx / kx2 * (1. - cx)
                dy = hy / ky2 * (1. - cy)
                r56 = hx * hx * (z - sx) / kx2 / beta ** 2 + hy * hy * (z - sy) / ky2 / beta ** 2
            else:
                sx = z
                dx = z * z * hx / 2.
                dy = z * z * hy / 2.
                r56 = hx * hx * z ** 3 / 6. / beta ** 2 + hy * hy * z ** 3 / 6. / beta ** 2

            r56 -= z / (beta * beta) * igamma2

            u_matrix = np.array([[cx, sx, 0., 0., 0., dx / beta],
                                 [-kx2 * sx, cx, 0., 0., 0., sx * hx / beta],
                                 [0., 0., cy, sy, 0., dy / beta],
                                 [0., 0., -ky2 * sy, cy, 0., sy * hy / beta],
                                 [hx * sx / beta, dx / beta, hy * sy / beta, dy / beta, 1., r56],
                                 [0., 0., 0., 0., 0., 1.]])
            if sum_tilts != 0:
                u_matrix = np.dot(np.dot(rot_mtx(-sum_tilts), u_matrix), rot_mtx(sum_tilts))
            return u_matrix

        r_z_e = lambda z, energy: r_mtx(z, k1, hx=hx, hy=hy, sum_tilts=0, energy=energy)
        return r_z_e

    def get_T_z_e_func(self):
        return lambda z, energy: np.zeros((6, 6, 6))

    def create_first_order_main_params(self, energy: float, delta_length: float = None) -> FirstOrderParams:
        r_z_e = self.create_r_matrix()
        R = r_z_e(delta_length if delta_length is not None else self.l, energy)
        B = self._default_B(R)
        return FirstOrderParams(R, B, self.tilt)

    def create_second_order_main_params(self, energy: float, delta_length: float = 0.0) -> SecondOrderParams:
        T = t_nnn(delta_length if delta_length is not None else self.l, 0., 0., 0., energy)
        first_order_params = self.create_first_order_main_params(energy, delta_length)
        return SecondOrderParams(first_order_params.R, first_order_params.B, T, self.tilt, self.dx, self.dy)