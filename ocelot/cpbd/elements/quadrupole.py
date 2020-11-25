import numpy as np

from ocelot.cpbd.elements.element import Element


class Quadrupole(Element):
    """
    quadrupole
    l - length of lens in [m],
    k1 - strength of quadrupole lens in [1/m^2],
    k2 - strength of sextupole lens in [1/m^3],
    tilt - tilt of lens in [rad].
    """

    def __init__(self, l=0., k1=0, k2=0., tilt=0., eid=None):
        # Element.__init__(self, eid)
        super(Quadrupole, self).__init__(eid=eid)
        self.l = l
        self.k1 = k1
        self.k2 = k2
        self.tilt = tilt

    def __str__(self):
        s = 'Quadrupole : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l    =%8.4f m\n' % self.l
        s += 'k1   =%8.3f 1/m^2\n' % self.k1
        s += 'k2   =%8.3f 1/m^3\n' % self.k2
        s += 'tilt =%8.2f deg\n' % (self.tilt * 180.0 / np.pi)
        return s