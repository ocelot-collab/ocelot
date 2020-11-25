import numpy as np

from ocelot.cpbd.high_order import m_e_GeV
from ocelot.cpbd.r_matrix import rot_mtx
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.elements.sbend import SBend


class XYQuadrupole(SBend):
    """
    Quadrupole with offsets (linear element). The element is to test a transport feature and it is not tested.

    l - length of magnet in [m],
    k1 - strength of quadrupole lens in [1/m^2],
    x_offs - offset in horizontal direction in [m]
    y_offs - offset in vertical direction in [m]
    tilt - tilt of lens in [rad],
    """

    def __init__(self, l=0., x_offs=0.0, y_offs=0.0, k1=0.0, tilt=0.0, eid=None):
        Element.__init__(self, eid)
        self.l = l
        self.k1 = k1
        self.x_offs = x_offs
        self.y_offs = y_offs
        self.tilt = tilt

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