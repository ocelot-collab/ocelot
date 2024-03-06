import numpy as np

from ocelot.cpbd.transformations.transformation import TMTypes
from ocelot.cpbd.high_order import m_e_GeV
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.transformations.transfer_map import TransferMap, TMTypes


class UndulatorTestTM(TransferMap):
    """[summary]
    Implementation of a Undulator Test Transformation.
    The concrete element atom have to implement: 
    create_undulator_test_tm_main_params(self) -> UndulatorTestParams
    """

    def __init__(self, create_tm_param_func, delta_e_func, tm_type: TMTypes, length: float, delta_length: float, **params) -> None:
        ndiv = params.get('ndiv')
        self.ndiv = ndiv if ndiv else 5
        super().__init__(create_tm_param_func, delta_e_func, tm_type, length, delta_length=delta_length)

    @classmethod
    def from_element(cls, element: Element, tm_type: TMTypes = TMTypes.MAIN, delta_l=None, **params):
        return cls.create(entrance_tm_params_func=None,
                          delta_e_func=element.create_delta_e,
                          main_tm_params_func=element.create_undulator_test_tm_main_params,
                          exit_tm_params_func=None,
                          has_params=False,
                          tm_type=tm_type, length=element.l, delta_length=delta_l, **params)

    def get_params(self):
        return self.create_tm_param_func()

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

    def map_function(self, X, energy: float):
        params = self.get_params()
        return self.map4undulator(X, self.delta_length if self.delta_length != None else self.length, params.lperiod, params.Kx, params.ax, energy, self.ndiv)
