import numpy as np

from ocelot.cpbd.high_order import m_e_GeV
from ocelot.cpbd.elements.element import Element


class Solenoid(Element):
    """
    Solenoid
    l - length in m,
    k - strength B0/(2B*rho)
    """

    def __init__(self, l=0., k=0., eid=None):
        Element.__init__(self, eid)
        self.k = k  # B0/(2B*rho)
        self.l = l

    def __str__(self):
        s = 'Cavity : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l =%8.4f m\n' % self.l
        s += 'k =%8.3f 1/m\n' % self.k
        return s

    def create_r_matrix(self):

        def sol(l, k, energy):
            """
            K.Brown, A.Chao.
            :param l: effective length of solenoid
            :param k: B0/(2*Brho), B0 is field inside the solenoid, Brho is momentum of central trajectory
            :return: matrix
            """
            gamma = energy / m_e_GeV
            c = np.cos(l * k)
            s = np.sin(l * k)
            if k == 0:
                s_k = l
            else:
                s_k = s / k
            r56 = 0.
            if gamma != 0:
                gamma2 = gamma * gamma
                beta = np.sqrt(1. - 1. / gamma2)
                r56 -= l / (beta * beta * gamma2)
            sol_matrix = np.array([[c * c, c * s_k, s * c, s * s_k, 0., 0.],
                                   [-k * s * c, c * c, -k * s * s, s * c, 0., 0.],
                                   [-s * c, -s * s_k, c * c, c * s_k, 0., 0.],
                                   [k * s * s, -s * c, -k * s * c, c * c, 0., 0.],
                                   [0., 0., 0., 0., 1., r56],
                                   [0., 0., 0., 0., 0., 1.]]).real
            return sol_matrix

        r_z_e = lambda z, energy: sol(z, k=self.k, energy=energy)
        return r_z_e