import numpy as np

from ocelot.cpbd.elements.element import Element


class Bend(Element):
    """
    bending magnet
    l - length of magnet in [m],
    angle - angle of bend in [rad],
    k1 - strength of quadrupole lens in [1/m^2],
    k2 - strength of sextupole lens in [1/m^3],
    tilt - tilt of lens in [rad],
    e1 - the angle of inclination of the entrance face [rad],
    e2 - the angle of inclination of the exit face [rad].
    fint - fringe field integral
    fintx - allows (fintx > 0) to set fint at the element exit different from its entry value.
    gap - the magnet gap [m], NOTE in MAD and ELEGANT: HGAP = gap/2
    h_pole1 - the curvature (1/r) of the entrance face
    h_pole1 - the curvature (1/r) of the exit face
    """

    def __init__(self, l=0., angle=0., k1=0., k2=0., e1=0., e2=0., tilt=0.0,
                 gap=0., h_pole1=0., h_pole2=0., fint=0., fintx=None, eid=None):
        Element.__init__(self, eid)
        self.l = l
        self.angle = angle
        self.k1 = k1
        self.k2 = k2
        self.e1 = e1
        self.e2 = e2
        self.gap = gap
        self.h_pole1 = h_pole1
        self.h_pole2 = h_pole2
        self.fint = fint
        self.fintx = fint
        if fintx is not None:
            self.fintx = fintx
        self.tilt = tilt

    def __str__(self):
        s = 'Bend : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l       =%8.4f m\n' % self.l
        s += 'angle   =%8.3f deg\n' % (self.angle * 180.0 / np.pi)
        s += 'e1      =%8.3f deg\n' % (self.e1 * 180.0 / np.pi)
        s += 'e2      =%8.3f deg\n' % (self.e2 * 180.0 / np.pi)
        s += 'tilt    =%8.3f deg\n' % (self.tilt * 180.0 / np.pi)
        s += 'fint    =%8.3f\n' % self.fint
        s += 'fintx   =%8.3f\n' % self.fintx
        s += 'gap     =%8.4f m\n' % self.gap
        s += 'h_pole1 =%8.4f 1/m\n' % self.h_pole1
        s += 'h_pole2 =%8.4f 1/m\n' % self.h_pole2
        return s