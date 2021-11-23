import numpy as np

from ocelot.cpbd.elements.cor_atom import CorAtom


class HcorAtom(CorAtom):
    """
    horizontal corrector,
    l - length of magnet in [m],
    angle - angle of bend in [rad],
    """

    def __init__(self, l=0., angle=0., eid=None):
        super().__init__(eid=eid)
        self.has_edge = False
        self.l = l
        self.angle = angle
        self.tilt = 0.

    def __str__(self):
        s = 'Hcor : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l     =%8.4f m\n' % self.l
        s += 'angle =%8.3f deg\n' % (self.angle * 180.0 / np.pi)
        s += 'tilt  =%8.2f deg\n' % (self.tilt * 180.0 / np.pi)
        return s
