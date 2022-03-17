import numpy as np

from ocelot.cpbd.elements.cor_atom import CorAtom


class VcorAtom(CorAtom):
    """
    horizontal corrector,
    l - length of magnet in [m],
    angle - angle of bend in [rad],
    """

    def __init__(self, l=0., angle=0., eid=None):
        super().__init__(l=l, angle=angle, eid=eid)
        self.has_edge = False
        self.l = l
        self.angle = angle
        self.tilt = np.pi / 2.

    def __str__(self):
        s = 'Vcor('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'angle=%8.6e, ' % self.angle if np.abs(self.angle) > 1e-15 else ""
        s += 'tilt=%8.6e, ' % self.tilt if np.abs(self.tilt) > 1e-15 else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s