import numpy as np

from ocelot.cpbd.elements.magnet import Magnet


class QuadrupoleAtom(Magnet):
    """
    quadrupole
    l - length of lens in [m],
    k1 - strength of quadrupole lens in [1/m^2],
    k2 - strength of sextupole lens in [1/m^3],
    tilt - tilt of lens in [rad].
    """

    def __init__(self, l=0., k1=0, k2=0., tilt=0., eid=None):
        # Element.__init__(self, eid)
        super().__init__(eid=eid)
        self.l = l
        self.k1 = k1
        self.k2 = k2
        self.tilt = tilt

    def __str__(self):
        s = 'Quadrupole('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'k1=%8.6e, ' % self.k1 if np.abs(self.k1) > 1e-15 else ""
        s += 'k2=%8.6e, ' % self.k2 if np.abs(self.k2) > 1e-15 else ""
        s += 'tilt=%8.6e, ' % self.tilt if np.abs(self.tilt) > 1e-15 else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s