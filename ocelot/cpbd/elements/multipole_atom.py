import numpy as np

from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.tm_params.multipole_params import MultipoleParams


class MultipoleAtom(Element):
    """
    kn - list of strengths
    """

    dtilt = 0.  # Legacy pickles can carry the roll error on the wrapper.

    def __init__(self, kn=0., eid=None, tilt=0., dtilt=0., **kwargs):
        super().__init__(eid, **kwargs)
        self.tilt = tilt
        self.dtilt = dtilt
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

    def create_multipole_tm_main_params(self):
        return MultipoleParams(self.kn, self.dx, self.dy, self.tilt + self.dtilt)

    def create_first_order_main_params(self, energy: float, delta_length: float = None) -> FirstOrderParams:
        R = np.eye(6)
        R[1, 0] = -self.kn[1]
        R[3, 2] = self.kn[1]
        R[1, 5] = self.kn[0]
        return self._aligned_map_params(R, self.tilt + self.dtilt)
