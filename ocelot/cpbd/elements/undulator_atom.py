import logging

import numpy as np

from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.field_map import FieldMap
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.tm_params.runge_kutta_params import RungeKuttaParams
from ocelot.cpbd.tm_params.undulator_test_params import UndulatorTestParams
from ocelot.cpbd.high_order import m_e_GeV, m_e_eV
from ocelot.common.globals import speed_of_light, pi

_logger = logging.getLogger(__name__)

try:
    import numba as nb

    nb_flag = True
except:
    _logger.info("radiation_py.py: module NUMBA is not installed. Install it to speed up calculation")
    nb_flag = False


def und_field_py(x, y, z, lperiod, Kx, nperiods=None):
    kx = 0.
    kz = 2 * pi / lperiod
    ky = np.sqrt(kz * kz + kx * kx)
    c = speed_of_light
    m0 = m_e_eV
    B0 = Kx * m0 * kz / c
    k1 = -B0 * kx / ky
    k2 = -B0 * kz / ky

    kx_x = kx * x
    ky_y = ky * y
    kz_z = kz * z

    cosz = np.cos(kz_z)

    if nperiods is not None:
        ph_shift = np.pi / 2.
        def heaviside(x): return 0.5 * (np.sign(x) + 1)
        z_coef = (0.25 * heaviside(z) + 0.5 * heaviside(z - lperiod / 2.) + 0.25 * heaviside(z - lperiod)
                  - 0.25 * heaviside(z - (nperiods - 1) * lperiod) - 0.5 * heaviside(
            z - (nperiods - 0.5) * lperiod)
            - 0.25 * heaviside(z - nperiods * lperiod))
        cosz = np.cos(kz_z + ph_shift) * z_coef

    cosx = np.cos(kx_x)
    sinhy = np.sinh(ky_y)
    # cosz = np.cos(kz_z + ph_shift)*z_coef
    Bx = k1 * np.sin(kx_x) * sinhy * cosz  # // here kx is only real
    By = B0 * cosx * np.cosh(ky_y) * cosz
    Bz = k2 * cosx * sinhy * np.sin(kz_z)
    return (Bx, By, Bz)


und_field = und_field_py if not nb_flag else nb.jit(forceobj=False)(und_field_py)


class UndulatorAtom(Element):
    """
    Undulator
    lperiod - undulator period in [m];\n
    nperiod - number of periods;\n
    Kx - undulator paramenter for vertical field; \n
    Ky - undulator parameter for horizantal field;\n
    field_file - absolute path to magnetic field data;\n
    mag_field - None by default, the magnetic field map function - (Bx, By, Bz) = f(x, y, z)
    eid - id of undulator.
    """

    def __init__(self, lperiod=0., nperiods=0, Kx=0., Ky=0., field_file=None, eid=None):
        Element.__init__(self, eid)
        self.lperiod = lperiod
        self.nperiods = nperiods
        self.l = lperiod * nperiods
        self.Kx = Kx
        self.Ky = Ky
        self.solver = "linear"  # can be "lin" is linear matrix,  "sym" - symplectic method and "rk" is Runge-Kutta
        self.phase = 0.  # phase between Bx and By + pi/4 (spiral undulator)

        self.ax = -1  # width of undulator, when ax is negative undulator width is infinite
        # I need this for analytic description of undulator

        self.field_file = field_file
        self.field_map = FieldMap(self.field_file)
        self.mag_field = None  # the magnetic field map function - (Bx, By, Bz) = f(x, y, z)
        self.v_angle = 0.
        self.h_angle = 0.

    def __str__(self):
        s = 'Undulator('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'nperiods=%7.2f, ' % self.nperiods if np.abs(self.nperiods) > 1e-15 else ""
        s += 'lperiod=%7.4f, ' % self.lperiod if np.abs(self.lperiod) > 1e-15 else ""
        s += 'Kx=%7.3f, ' % self.Kx if np.abs(self.Kx) > 1e-15 else ""
        s += 'Ky=%7.3f, ' % self.Ky if np.abs(self.Ky) > 1e-15 else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s

    def R_main_matrix(self, energy, length):
        """
        in OCELOT coordinates:
        R56 = - Lu/(gamma**2 * beta**2) * (1 + 0.5 * K**2 * beta**2)
        S.Tomin, Varenna, 2017.
        """

        def undulator_r_z(z, lperiod, Kx, Ky, energy):
            gamma = energy / m_e_GeV
            r = np.eye(6)
            r[0, 1] = z
            if gamma != 0 and lperiod != 0 and Kx != 0:
                beta = 1 / np.sqrt(1.0 - 1.0 / (gamma * gamma))

                omega_x = np.sqrt(2.0) * np.pi * Kx / (lperiod * gamma * beta)
                omega_y = np.sqrt(2.0) * np.pi * Ky / (lperiod * gamma * beta)
                r[2, 2] = np.cos(omega_x * z)
                r[2, 3] = np.sin(omega_x * z) / omega_x
                r[3, 2] = -np.sin(omega_x * z) * omega_x
                r[3, 3] = np.cos(omega_x * z)

                r[4, 5] = - z / (gamma * beta) ** 2 * (1 + 0.5 * (Kx * beta) ** 2)

            else:
                r[2, 3] = z
            return r

        R = undulator_r_z(length, lperiod=self.lperiod, Kx=self.Kx, Ky=self.Ky, energy=energy)
        return R

    def create_first_order_main_params(self, energy: float, delta_length: float = None) -> FirstOrderParams:
        R = self.R_main_matrix(energy=energy, length=delta_length if delta_length != None else self.l)
        B = self._default_B(R)
        return FirstOrderParams(R, B, self.tilt)

    def create_runge_kutta_main_params(self):
        return RungeKuttaParams(mag_field=lambda x, y, z: und_field(x, y, z, self.lperiod, self.Kx))

    def und_field(self):
        return lambda x, y, z: und_field(x, y, z, self.lperiod, self.Kx)

    def create_undulator_test_tm_main_params(self) -> UndulatorTestParams:
        return UndulatorTestParams(self.lperiod, self.Kx, self.ax)
