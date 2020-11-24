from copy import copy

from ocelot.cpbd.high_order import rk_field
from ocelot.cpbd.transformations.first_order import TransferMap


class RungeKuttaTM(TransferMap):
    def __init__(self, s_start=0, npoints=200):
        TransferMap.__init__(self)
        self.s_start = s_start
        self.npoints = npoints
        self.long_dynamics = True
        self.mag_field = lambda x, y, z: (0, 0, 0)
        self.map = lambda X, energy: rk_field(X, self.s_start, self.length, self.npoints, energy, self.mag_field,
                                              self.long_dynamics)

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.delta_e = m.delta_e_z(s)
        m.map = lambda X, energy: rk_field(X, m.s_start, s, m.npoints, energy, m.mag_field, m.long_dynamics)
        return m

    @classmethod
    def create_from_element(cls, element, params=None):
        tm = cls(s_start=element.s_start if hasattr(element, 's_start') else 0.,
                 npoints=element.npoints if hasattr(element, 'npoints') else 200)
        tm.mag_field = element.mag_field
        return tm
