import numpy as np

from ocelot.cpbd.high_order import fringe_ent, fringe_ext
from ocelot.cpbd.elements.bend import Bend
from ocelot.cpbd.elements.element import Element


class Edge(Bend):
    def __init__(self, l=0., angle=0.0, k1=0., edge=0.,
                 tilt=0.0, dx=0.0, dy=0.0,
                 h_pole=0., gap=0., fint=0., pos=1, eid=None):
        Element.__init__(self, eid)
        if l != 0.:
            self.h = angle / l
        else:
            self.h = 0
        self.l = 0.
        self.k1 = k1
        self.h_pole = h_pole
        self.gap = gap
        self.fint = fint
        self.edge = edge
        self.dx = dx
        self.dy = dy
        self.tilt = tilt
        self.pos = pos

    def __str__(self):
        s = 'Edge : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'h      =%8.4f 1/m\n' % self.h
        s += 'fint   =%8.3f\n' % self.fint
        s += 'gap    =%8.4f m\n' % self.gap
        s += 'h_pole =%8.4f 1/m\n' % self.h_pole
        s += 'tilt   =%8.3f deg\n' % (self.tilt * 180.0 / np.pi)
        return s

    def create_r_matrix(self):
        sec_e = 1. / np.cos(self.edge)
        phi = self.fint * self.h * self.gap * sec_e * (1. + np.sin(self.edge) ** 2)
        # phi = self.fint * self.h * self.gap * sec_e * (1. + np.sin(2*self.edge) )
        r = np.eye(6)
        r[1, 0] = self.h * np.tan(self.edge)
        r[3, 2] = -self.h * np.tan(self.edge - phi)
        r_z_e = lambda z, energy: r
        return r_z_e

    def get_T_z_e_func(self):
        if self.pos == 1:
            _, T = fringe_ent(h=self.h, k1=self.k1, e=self.edge, h_pole=self.h_pole,
                              gap=self.gap, fint=self.fint)
        else:
            _, T = fringe_ext(h=self.h, k1=self.k1, e=self.edge, h_pole=self.h_pole,
                              gap=self.gap, fint=self.fint)
        return lambda z, energy: T