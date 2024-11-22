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


def und_field_py(x, y, z, lperiod, Kx, nperiods=None, phase=0, end_poles='1'):
    """
    Parameters
    ----------
    x : 
        particle coordinate in x.
    y : 
        particle coordinate in y.
    z : 
        particle coordinate in z.
    lperiod : 
        undulator period length.
    Kx : 
        undulator K parameter in x direction.
    nperiods : 
        Number of undulator periods. The default is None.
    phase : optional
        Phase from which the magnetic field is stated. The default is 0, which is cos().
    end_poles : optional
        Correction poles to close magnetic field integrals. 
        Might be: 
            '0' - magnetic field starts from 0 value or sin()-like
            '1' - magnetic field starts from maximum value or cos()-like
            '3/4' - magnetic field starts with 0, 1/4, -3/4, +1, -1  (or 0, -1/4, +3/4, -1, +1) poles sequence and finishes with it (fraction is from maximum value of the field)
            '1/2' - magnetic field starts with 1/2 and finishes with -1/2 poles
        The default is '1'.

    Raises
    ------
    ValueError
        'end_poles' must be either '1', '0', '3/4' or '1/2'. 
        if '0' the field starts from 0, 
        if '1' the field starts from its maximum value,
        if 3/4' the following end poles sequence is added 0, 1/4, -3/4, +1, -1,
        if 1/2 the following end poles sequence is added 0, 1/2, -1, +1".

    Returns
    -------
    Bx : 
        Magnetic field in x direction.
    By : 
        Magnetic field in y direction.
    Bz : 
        Magnetic field in z direction.

    """
    kx = 0.
    kz = 2 * pi / lperiod
    ky = np.sqrt(kz * kz + kx * kx)
    c = speed_of_light
    m0 = m_e_eV
    B0 = Kx * m0 * kz / c
    k1 = B0 * kx / ky
    k2 = B0 * kz / ky

    kx_x = kx * x
    ky_y = ky * y
    kz_z = kz * z + phase

    if nperiods is not None:

        ph_shift = np.pi / 2.

        def heaviside(x):
            return 0.5 * (np.sign(x) + 1)

        if end_poles == '0':
            z_coef = 1
            cosz = np.cos(kz_z + ph_shift) * z_coef
            sinz = np.sin(kz_z + ph_shift) * z_coef

        elif end_poles == '1':
            z_coef = 1
            cosz = np.cos(kz_z) * z_coef
            sinz = np.sin(kz_z) * z_coef

        elif end_poles == '3/4':
            z_coef = (0.25 * heaviside(z) + 0.5 * heaviside(z - lperiod / 2.) + 0.25 * heaviside(z - lperiod)
                      - 0.25 * heaviside(z - (nperiods - 1) * lperiod) - 0.5 * heaviside(z - (nperiods - 0.5) * lperiod)
                      - 0.25 * heaviside(z - nperiods * lperiod))

            cosz = np.cos(kz_z + ph_shift) * z_coef
            sinz = np.sin(kz_z + ph_shift) * z_coef

        elif end_poles == '1/2':
            z_coef = (-0.5 * heaviside(z - lperiod * 0.5) - 0.5 * heaviside(z + lperiod)
                      + 0.5 * heaviside(z - (nperiods - 0.5) * lperiod)
                      + 0.5 * heaviside(z - nperiods * lperiod))

            cosz = np.cos(kz_z + ph_shift) * z_coef
            sinz = np.sin(kz_z + ph_shift) * z_coef

        else:
            raise ValueError("'end_poles' must be either '1', '0', '3/4' or '1/2';" +
                             "if '1' the field starts from max value;" +
                             "if '0' the field starts from 0;" +
                             "if 3/4' the following end poles sequence is added 1/4, -3/4, +1, -1 instead of existing +1, -1, +1, -1" +
                             "if 1/2 the following end poles sequence is added 1/2, -1, +1 instead of existing +1, -1, +1.")
    else:
        cosz = np.cos(kz_z)
        sinz = np.sin(kz_z)

    coshx = np.cosh(kx_x)
    sinhx = np.sinh(kx_x)
    coshy = np.cosh(ky_y)
    sinhy = np.sinh(ky_y)
    # cosz = np.cos(kz_z + ph_shift)*z_coef
    Bx =  k1 * sinhx * sinhy * cosz  # // here kx is only real
    By =  B0 * coshx * coshy * cosz
    Bz = -k2 * coshx * sinhy * sinz
    return (Bx, By, Bz)


und_field = und_field_py if not nb_flag else nb.jit(forceobj=False, nopython=True)(und_field_py)


class UndulatorAtom(Element):
    """
    Undulator
    lperiod - undulator period in [m];\n
    nperiod - number of periods;\n
    Kx - undulator paramenter for vertical field; \n
    Ky - undulator parameter for horizontal field;\n
    field_file - absolute path to magnetic field data;\n
    mag_field - None by default, the magnetic field map function - (Bx, By, Bz) = f(x, y, z)
    end_poles = "1"
                '0' - magnetic field starts from 0 value or sin()-like
                '1' - magnetic field starts from maximum value or cos()-like
                '3/4' - magnetic field starts with 0, 1/4, -3/4, +1, -1  (or 0, -1/4, +3/4, -1, +1) poles sequence
                        and finishes with it (fraction is from maximum value of the field)
                '1/2' - magnetic field starts with 1/2 and finishes with -1/2 poles
                The default is '1'
    phase : optional
            Phase from which the magnetic field is stated. The default is 0, which is cos().
    eid - id of undulator.
    """

    def __init__(self, lperiod=0., nperiods=0, Kx=0., Ky=0., phase=0, end_poles='1', field_file=None, eid=None):
        Element.__init__(self, eid)
        self.lperiod = lperiod
        self.nperiods = nperiods
        self.l = lperiod * nperiods
        self.Kx = Kx
        self.Ky = Ky
        self.solver = "linear"  # can be "lin" is linear matrix,  "sym" - symplectic method and "rk" is Runge-Kutta
        self.phase = phase  # phase between Bx and By + pi/4 (spiral undulator)
        self.end_poles = end_poles

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
            r[2, 3] = z

            if gamma != 0 and lperiod != 0:
                beta = np.sqrt(1.0 - 1.0 / (gamma * gamma))
                r[4, 5] = - z / (gamma * beta) ** 2 * (1 + 0.5 * (np.sqrt(Kx * Kx + Ky * Ky).real * beta) ** 2)

                if Kx.real > 0 and Kx.imag == 0:
                    omega_x = np.sqrt(2.0) * np.pi * Kx / (lperiod * gamma * beta)
                    r[2, 2] = np.cos(omega_x * z)
                    r[2, 3] = np.sin(omega_x * z) / omega_x
                    r[3, 2] = -np.sin(omega_x * z) * omega_x
                    r[3, 3] = np.cos(omega_x * z)

                elif Kx.real == 0 and Kx.imag != 0:
                    omega_x = np.sqrt(2.0) * np.pi * Kx / (lperiod * gamma * beta)
                    r[2, 2] = np.real(np.cosh(omega_x * z))
                    r[2, 3] = np.real(np.sinh(omega_x * z) / omega_x)
                    r[3, 2] = np.real(-np.sinh(omega_x * z) * omega_x)
                    r[3, 3] = np.real(np.cosh(omega_x * z))
                else:
                    r[2, 3] = z

                if Ky.real > 0 and Ky.imag == 0:
                    omega_y = np.sqrt(2.0) * np.pi * Ky / (lperiod * gamma * beta)
                    r[0, 0] = np.cos(omega_y * z)
                    r[0, 1] = np.sin(omega_y * z) / omega_y
                    r[1, 0] = -np.sin(omega_y * z) * omega_y
                    r[1, 1] = np.cos(omega_y * z)

                elif Ky.real == 0 and Ky.imag != 0:
                    omega_y = np.sqrt(2.0) * np.pi * Ky / (lperiod * gamma * beta)
                    r[0, 0] = np.real(np.cosh(omega_y * z))
                    r[0, 1] = np.real(np.sinh(omega_y * z) / omega_y)
                    r[1, 0] = np.real(-np.sinh(omega_y * z) * omega_y)
                    r[1, 1] = np.real(np.cosh(omega_y * z))
                else:
                    r[0, 1] = z

            return r

        R = undulator_r_z(length, lperiod=self.lperiod, Kx=self.Kx, Ky=self.Ky, energy=energy)
        return R

    def create_first_order_main_params(self, energy: float, delta_length: float = None) -> FirstOrderParams:
        R = self.R_main_matrix(energy=energy, length=delta_length if delta_length != None else self.l)
        B = self._default_B(R)
        return FirstOrderParams(R, B, self.tilt)

    def create_runge_kutta_main_params(self, energy):
        return RungeKuttaParams(mag_field=lambda x, y, z: und_field(x, y, z, self.lperiod, self.Kx))

    def und_field(self):
        return lambda x, y, z: und_field(x, y, z, self.lperiod, self.Kx)

    def create_undulator_test_tm_main_params(self) -> UndulatorTestParams:
        return UndulatorTestParams(self.lperiod, self.Kx, self.ax)
