import numpy as np

from ocelot.cpbd.transformations.multipole import MultipoleTM
from ocelot.cpbd.elements.element import Element


class Multipole(Element):
    """
    kn - list of strengths
    """

    default_tm = MultipoleTM
    additional_tms = []

    def __init__(self, kn=0., eid=None):
        Element.__init__(self, eid)
        kn = np.array([kn]).flatten()
        if len(kn) < 2:
            self.kn = np.append(kn, [0.])
        else:
            self.kn = kn
        self.n = len(self.kn)
        self.l = 0.

    def __str__(self):
        s = 'Multipole : '
        s += 'id = ' + str(self.id) + '\n'
        for i, k in enumerate(self.kn):
            s += 'k%i =%8.4f m\n' % (i, k)
        return s

    def create_r_matrix(self):
        r = np.eye(6)
        r[1, 0] = -self.kn[1]
        r[3, 2] = self.kn[1]
        r[1, 5] = self.kn[0]
        r_z_e = lambda z, energy: r
        return r_z_e