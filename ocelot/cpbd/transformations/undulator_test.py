from copy import copy

import numpy as np

from ocelot.cpbd.transformations.optics import m_e_GeV
from ocelot.cpbd.transformations.first_order import TransferMap


class UndulatorTestTM(TransferMap):
    def __init__(self, lperiod, Kx, ax=0, ndiv=10):
        TransferMap.__init__(self)
        self.lperiod = lperiod
        self.Kx = Kx
        self.ax = ax
        self.ndiv = ndiv
        self.map = lambda X, energy: self.map4undulator(X, self.length, self.lperiod, self.Kx, self.ax, energy,
                                                        self.ndiv)

    def map4undulator(self, u, z, lperiod, Kx, ax, energy, ndiv):
        kz = 2. * np.pi / lperiod
        if ax == 0:
            kx = 0
        else:
            kx = 2. * np.pi / ax
        zi = np.linspace(0., z, num=ndiv)
        h = zi[1] - zi[0]
        kx2 = kx * kx
        kz2 = kz * kz
        ky2 = kz * kz + kx * kx
        ky = np.sqrt(ky2)
        gamma = energy / m_e_GeV
        h0 = 0.
        if gamma != 0:
            h0 = 1. / (gamma / Kx / kz)
        h02 = h0 * h0
        h = h / (1. + u[5])
        x = u[0]
        y = u[2]
        for z in range(len(zi) - 1):
            chx = np.cosh(kx * x)
            chy = np.cosh(ky * y)
            shx = np.sinh(kx * x)
            shy = np.sinh(ky * y)
            u[1] -= h / 2. * chx * shx * (kx * ky2 * chy * chy + kx2 * kx * shy * shy) / (ky2 * kz2) * h02
            u[3] -= h / 2. * chy * shy * (ky2 * chx * chx + kx2 * shx * shx) / (ky * kz2) * h02
            u[4] -= h / 2. / (1. + u[5]) * ((u[1] * u[1] + u[3] * u[3]) + chx * chx * chy * chy / (
                    2. * kz2) * h02 + shx * shx * shy * shy * kx2 / (2. * ky2 * kz2) * h02)
            u[0] = x + h * u[1]
            u[2] = y + h * u[3]
        return u

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.delta_e = m.delta_e_z(s)
        m.map = lambda X, energy: m.map4undulator(X, m.length, m.lperiod, m.Kx, m.ax, energy, m.ndiv)
        return m

    @classmethod
    def create_from_element(cls, element, params=None):
        return cls(lperiod=element.lperiod, Kx=element.Kx, ax=element.ax,
                   ndiv=element.ndiv if hasattr(element, 'ndiv') else 5)
