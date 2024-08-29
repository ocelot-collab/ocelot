import numpy as np

from ocelot.cpbd.elements.element import Element


class ApertureAtom(Element):
    """
    Aperture
    xmax - half size in horizontal plane in [m],
    ymax - half size in vertical plane in [m],
    type - "rect" or "ellipt".
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
        s = 'Aperture('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'xmax=%8.6e, ' % self.xmax if self.xmax != np.inf else ""
        s += 'ymax=%8.6e, ' % self.ymax if self.ymax != np.inf else ""
        s += 'dx=%8.6e, ' % self.dx if np.abs(self.dx) > 1e-15 else ""
        s += 'dy=%8.6e, ' % self.dy if np.abs(self.dy) > 1e-15 else ""
        s += 'type=%s, ' % self.type if self.type != "rect" else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s