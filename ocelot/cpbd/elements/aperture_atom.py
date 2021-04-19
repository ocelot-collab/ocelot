import numpy as np

from ocelot.cpbd.elements.element import Element


class ApertureAtom(Element):
    """
    Aperture
    xmax - half size in horizontal plane in [m],
    ymax - half size in vertical plane in [m],
    type - "rect" or "elliptical".
    """

    def __init__(self, xmax=np.inf, ymax=np.inf, dx=0, dy=0, type="rect", eid=None):
        Element.__init__(self, eid)
        self.l = 0.
        self.xmax = xmax
        self.ymax = ymax
        self.dx = dx
        self.dy = dy
        self.type = type

    def __str__(self):
        s = 'Aperture : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l    =%8.4f m\n' % self.l
        s += 'xmax =%8.5f m\n' % self.xmax
        s += 'ymax =%8.5f m\n' % self.ymax
        s += 'dx   =%8.5f m\n' % self.dx
        s += 'dy   =%8.5f m\n' % self.dy
        s += 'type = %s \n' % self.type
        return s