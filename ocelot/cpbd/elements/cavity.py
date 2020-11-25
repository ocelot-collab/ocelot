import logging

import numpy as np

from ocelot.cpbd.transformations.cavity import CavityTM
from ocelot.cpbd.high_order import m_e_GeV
from ocelot.common.globals import speed_of_light
from ocelot.cpbd.r_matrix import uni_matrix
from ocelot.cpbd.elements.element import Element

logger = logging.getLogger(__name__)


class Cavity(Element):
    """
    Standing wave RF cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    vx_{up/down}, vy_{up/down} - zero order kick of a {up/down}stream coupler
    vxx_{up/down}, vxy_{up/down} - first order kick  a {up/down}stream coupler
    """

    default_tm = CavityTM
    additional_tms = []

    def __init__(self, l=0., v=0., phi=0., freq=0., vx_up=0, vy_up=0, vxx_up=0, vxy_up=0,
                 vx_down=0, vy_down=0, vxx_down=0, vxy_down=0, eid=None):
        Element.__init__(self, eid)
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
        s = 'Cavity : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l    =%8.4f m\n' % self.l
        s += 'v    =%8.5f GV\n' % self.v
        s += 'freq =%8.1e Hz\n' % self.freq
        s += 'phi  =%8.2f deg\n' % self.phi
        s += "\nCoupler kick: \n"
        s += "vx_up    = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vx_up)
        s += "vy_up    = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vy_up)
        s += "vxx_up   = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vxx_up)
        s += "vxy_up   = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vxy_up)
        s += "vx_down  = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vx_down)
        s += "vy_down  = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vy_down)
        s += "vxx_down = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vxx_down)
        s += "vxy_down = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vxy_down)
        return s

    def create_r_matrix(self):

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
            r_z_e = lambda z, energy: uni_matrix(z, 0., hx=0., sum_tilts=self.tilt, energy=energy)
        else:
            r_z_e = lambda z, energy: cavity_R_z(z, V=self.v * z / self.l, E=energy, freq=self.freq,
                                                 phi=self.phi)
        return r_z_e
