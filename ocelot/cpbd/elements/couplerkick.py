import numpy as np

from ocelot.cpbd.transformations.coupler_kick import CouplerKickTM
from ocelot.cpbd.elements.element import Element


class CouplerKick(Element):
    """
    Coupler Kick element for Cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    vx, vy - zero order kick of a stream coupler
    vxx, vxy - first order kick  a stream coupler
    """

    default_tm = CouplerKickTM
    additional_tms = []

    def __init__(self, v=0., phi=0., freq=0., vx=0., vy=0., vxx=0., vxy=0., eid=None):
        Element.__init__(self, eid)
        self.l = 0.
        self.v = v  # in GV
        self.freq = freq  # Hz
        self.phi = phi  # in grad
        self.vx = vx
        self.vy = vy
        self.vxx = vxx
        self.vxy = vxy

    def __str__(self):
        s = 'CouplerKick : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'v    =%8.5f GV\n' % self.v
        s += 'freq =%8.1e Hz\n' % self.freq
        s += 'phi  =%8.2f deg\n' % self.phi
        s += "vx   = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vx)
        s += "vy   = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vy)
        s += "vxx  = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vxx)
        s += "vxy  = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vxy)
        return s

    def create_r_matrix(self):

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

        r_z_e = lambda z, energy: ck_matrix(v=self.v, phi=self.phi,
                                            vxx=self.vxx, vxy=self.vxy, energy=energy)
        return r_z_e