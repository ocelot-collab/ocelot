import numpy as np

from ocelot.cpbd.transformations.tw_cavity import TWCavityTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.common.globals import speed_of_light
from ocelot.cpbd.r_matrix import uni_matrix
from ocelot.cpbd.elements.element import Element


class TWCavity(Element):
    """
    Traveling wave cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    """

    default_tm = TWCavityTM
    additional_tms = [SecondTM]

    def __init__(self, l=0., v=0., phi=0., freq=0., eid=None):
        Element.__init__(self, eid)
        self.l = l
        self.v = v  # in GV
        self.freq = freq  # Hz
        self.phi = phi  # in grad
        self.E = 0

    def create_r_matrix(self):

        def tw_cavity_R_z(z, V, E, freq, phi=0.):
            """
            :param z: length
            :param de: delta E
            :param f: frequency
            :param E: initial energy
            :return: matrix
            """
            phi = phi * np.pi / 180.
            de = V * np.cos(phi)
            r12 = z * E / de * np.log(1. + de / E) if de != 0 else z
            r22 = E / (E + de)
            r65 = V * np.sin(phi) / (E + de) * (2 * np.pi / (speed_of_light / freq)) if freq != 0 else 0
            r66 = r22
            cav_matrix = np.array([[1, r12, 0., 0., 0., 0.],
                                   [0, r22, 0., 0., 0., 0.],
                                   [0., 0., 1, r12, 0., 0.],
                                   [0., 0., 0, r22, 0., 0.],
                                   [0., 0., 0., 0., 1., 0],
                                   [0., 0., 0., 0., r65, r66]]).real
            return cav_matrix

        def f_entrance(z, V, E, phi=0.):
            phi = phi * np.pi / 180.
            de = V * np.cos(phi)
            r = np.eye(6)
            r[1, 0] = -de / z / 2. / E
            r[3, 2] = r[1, 0]
            return r

        def f_exit(z, V, E, phi=0.):
            phi = phi * np.pi / 180.
            de = V * np.cos(phi)
            r = np.eye(6)
            r[1, 0] = +de / z / 2. / (E + de)
            r[3, 2] = r[1, 0]
            return r

        def cav(z, V, E, freq, phi):
            R_z = np.dot(tw_cavity_R_z(z, V, E, freq, phi), f_entrance(z, V, E, phi))
            R = np.dot(f_exit(z, V, E, phi), R_z)
            return R

        if self.v == 0.:
            r_z_e = lambda z, energy: uni_matrix(z, 0., hx=0., sum_tilts=self.tilt, energy=energy)
        else:
            r_z_e = lambda z, energy: cav(z, V=self.v * z / self.l, E=energy, freq=self.freq,
                                          phi=self.phi)
        return r_z_e