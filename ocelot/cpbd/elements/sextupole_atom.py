import numpy as np

from ocelot.cpbd.elements.magnet import Magnet


class SextupoleAtom(Magnet):
    """
    sextupole
    l - length of lens in [m],
    k2 - strength of sextupole lens in [1/m^3].
    """

    def __init__(self, l=0., k2=0., tilt=0., eid=None):
        super().__init__(eid)
        self.l = l
        self.k2 = k2
        self.tilt = tilt

    def __str__(self):
        s = 'Sextupole : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l    =%8.4f m\n' % self.l
        s += 'k2   =%8.3f 1/m^3\n' % self.k2
        s += 'tilt =%8.2f deg\n' % (self.tilt * 180.0 / np.pi)
        return s
