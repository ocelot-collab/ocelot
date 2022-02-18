import logging

import numpy as np

from ocelot.cpbd.transformations.cavity import CavityTM
from ocelot.cpbd.high_order import m_e_GeV
from ocelot.common.globals import speed_of_light
from ocelot.cpbd.r_matrix import uni_matrix
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.tm_params.cavity_params import CavityParams

logger = logging.getLogger(__name__)


class CavityAtom(Element):
    """
    Standing wave RF cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    vx_{up/down}, vy_{up/down} - zero order kick of a {up/down}stream coupler
    vxx_{up/down}, vxy_{up/down} - first order kick  a {up/down}stream coupler
    """

    def __init__(self, l=0., v=0., phi=0., freq=0., vx_up=0, vy_up=0, vxx_up=0, vxy_up=0,
                 vx_down=0, vy_down=0, vxx_down=0, vxy_down=0, eid=None):
        super().__init__(eid=eid, has_edge=True)
        self.l = l
        self.v = v  # in GV
        self.freq = freq  # Hz
        self.phi = phi  # in grad
        self.E = 0
        self.vx_up = vx_up
        self.vy_up = vy_up
        self.vxx_up = vxx_up
        self.vxy_up = vxy_up
        self.vx_down = vx_down
        self.vy_down = vy_down
        self.vxx_down = vxx_down
        self.vxy_down = vxy_down

    def __str__(self):
        s = 'Cavity('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'v=%8.6e, ' % self.v if self.v != 0. else ""
        s += 'freq=%8.6e, ' % self.freq if np.abs(self.freq) > 1e-15 else ""
        s += 'phi=%8.6e, ' % self.phi if np.abs(self.phi) > 1e-15 else ""
        s += "vx_up=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vx_up) if self.vx_up != 0. else ""
        s += "vy_up=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vy_up) if self.vy_up != 0. else ""
        s += "vxx_up=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vxx_up) if self.vxx_up != 0. else ""
        s += "vxy_up=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vxy_up) if self.vxy_up != 0. else ""
        s += "vx_down=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vx_down) if self.vx_down != 0. else ""
        s += "vy_down=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vy_down) if self.vy_down != 0. else ""
        s += "vxx_down=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vxx_down) if self.vxx_down != 0. else ""
        s += "vxy_down=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vxy_down) if self.vxy_down != 0. else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s

    def _R_edge_matrix(self, energy: float, vxx: float, vxy: float):

        def ck_matrix(v, phi, vxx, vxy, energy):
            """
            matrix for coupler kick

            :param v: voltage of the cavity in GV
            :param phi: phase [deg] of the cavity
            :param vxx: first order coefficients of the coupler kicks
            :param vxy: first order coefficients of the coupler kicks
            :param energy: beam energy in GeV
            :return:
            """
            phi = phi * np.pi / 180.
            m21 = (vxx * v * np.exp(1j * phi)).real / energy
            m43 = - m21
            m23 = (vxy * v * np.exp(1j * phi)).real / energy

            coupl_kick = np.array([[1, 0., 0., 0., 0., 0.],
                                   [m21, 1, m23, 0., 0., 0.],
                                   [0., 0., 1, 0., 0., 0.],
                                   [m23, 0., m43, 1, 0., 0.],
                                   [0., 0., 0., 0., 1., 0.],
                                   [0., 0., 0., 0., 0., 1]])
            return coupl_kick

        r_z_e = ck_matrix(v=self.v, phi=self.phi,
                          vxx=vxx, vxy=vxy, energy=energy)
        return r_z_e

    def _R_main_matrix(self, energy: float, length: float):

        def cavity_R_z(z, V, E, freq, phi=0.):
            """
            :param z: length
            :param de: delta E
            :param freq: frequency
            :param E: initial energy
            :return: matrix
            """

            phi = phi * np.pi / 180.
            de = V * np.cos(phi)
            # pure pi-standing-wave case
            eta = 1
            # gamma = (E + 0.5 * de) / m_e_GeV
            Ei = E / m_e_GeV
            Ef = (E + de) / m_e_GeV
            Ep = (Ef - Ei) / z  # energy derivative
            if Ei == 0:
                logger.error("CAVITY: Initial energy is 0, check ParticleArray.E or Twiss.E OR cavity.v must be 0")

            cos_phi = np.cos(phi)
            alpha = np.sqrt(eta / 8.) / cos_phi * np.log(Ef / Ei)
            sin_alpha = np.sin(alpha)

            cos_alpha = np.cos(alpha)
            r11 = (cos_alpha - np.sqrt(2. / eta) * cos_phi * sin_alpha)

            if abs(Ep) > 1e-10:
                r12 = np.sqrt(8. / eta) * Ei / Ep * cos_phi * sin_alpha
            else:
                r12 = z
            r21 = -Ep / Ef * (cos_phi / np.sqrt(2. * eta) + np.sqrt(eta / 8.) / cos_phi) * sin_alpha

            r22 = Ei / Ef * (cos_alpha + np.sqrt(2. / eta) * cos_phi * sin_alpha)
            # print(f"z = {z}, V = {V}, E = {E}, phi = {phi}, alpha = {alpha}, r11 = {r11}")

            r56 = 0.
            beta0 = 1
            beta1 = 1

            k = 2. * np.pi * freq / speed_of_light
            r55_cor = 0.
            if V != 0 and E != 0:
                gamma2 = Ei * Ei
                beta0 = np.sqrt(1. - 1 / gamma2)
                gamma2 = Ef * Ef
                beta1 = np.sqrt(1. - 1 / gamma2)

                # r56 = (beta0 / beta1 - 1) * Ei / (Ef - Ei) * z
                r56 = - z / (Ef * Ef * Ei * beta1) * (Ef + Ei) / (beta1 + beta0)
                g0 = Ei
                g1 = Ef
                r55_cor = k * z * beta0 * V / m_e_GeV * np.sin(phi) * (g0 * g1 * (beta0 * beta1 - 1) + 1) / (
                    beta1 * g1 * (g0 - g1) ** 2)

            r66 = Ei / Ef * beta0 / beta1
            r65 = k * np.sin(phi) * V / (Ef * beta1 * m_e_GeV)
            cav_matrix = np.array([[r11, r12, 0., 0., 0., 0.],
                                   [r21, r22, 0., 0., 0., 0.],
                                   [0., 0., r11, r12, 0., 0.],
                                   [0., 0., r21, r22, 0., 0.],
                                   [0., 0., 0., 0., 1. + r55_cor, r56],
                                   [0., 0., 0., 0., r65, r66]]).real

            return cav_matrix

        if self.v == 0.:
            R = uni_matrix(length, 0., hx=0., sum_tilts=self.tilt, energy=energy)
        else:
            R = cavity_R_z(length, V=self.v * length / self.l, E=energy, freq=self.freq,
                           phi=self.phi)
        return R

    def kick_b(self, v, vx, vy, phi, energy):
        phi = phi * np.pi / 180.
        dxp = (vx * v * np.exp(1j * phi)).real / energy
        dyp = (vy * v * np.exp(1j * phi)).real / energy
        b = np.array([[0.], [dxp], [0.], [dyp], [0.], [0.]])
        return b

    def create_first_order_main_params(self, energy: float, delta_length: float) -> FirstOrderParams:
        R = self._R_main_matrix(energy=energy, length=delta_length if delta_length is not None else self.l)
        B = self._default_B(R)
        return FirstOrderParams(R, B, self.tilt)

    def create_first_order_entrance_params(self, energy: float, delta_length: float) -> FirstOrderParams:
        R = self._R_edge_matrix(energy=energy, vxx=self.vxx_up, vxy=self.vxy_up)
        B = self.kick_b(self.v, self.vx_up, self.vy_up, self.phi, energy)
        return FirstOrderParams(R, B, self.tilt)

    def create_first_order_exit_params(self, energy: float, delta_length: float) -> FirstOrderParams:
        R = self._R_edge_matrix(energy=energy, vxx=self.vxx_down, vxy=self.vxy_down)
        B = self.kick_b(self.v, self.vx_down, self.vy_down, self.phi, energy)
        return FirstOrderParams(R, B, self.tilt)

    def create_cavity_tm_main_params(self, energy: float, delta_length: float) -> CavityParams:
        fo_params = self.create_first_order_main_params(energy, delta_length)
        return CavityParams(R=fo_params.R, B=fo_params.B, tilt=self.tilt, v=self.v, freq=self.freq, phi=self.phi)

    def create_cavity_tm_entrance_params(self, energy: float, delta_length: float) -> CavityParams:
        fo_params = self.create_first_order_entrance_params(energy, delta_length)
        return CavityParams(R=fo_params.R, B=fo_params.B, tilt=self.tilt, v=self.v, freq=self.freq, phi=self.phi)

    def create_cavity_tm_exit_params(self, energy: float, delta_length: float) -> CavityParams:
        fo_params = self.create_first_order_exit_params(energy, delta_length)
        return CavityParams(R=fo_params.R, B=fo_params.B, tilt=self.tilt, v=self.v, freq=self.freq, phi=self.phi)

    def create_delta_e(self, total_length, delta_length=None) -> float:
        if delta_length is not None:
            phase_term = np.cos(self.phi * np.pi / 180.)
            return phase_term * self.v * delta_length / total_length if total_length != 0 else phase_term * self.v
        else:
            return self.v * np.cos(self.phi * np.pi / 180.)
