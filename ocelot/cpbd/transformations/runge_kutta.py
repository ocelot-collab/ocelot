from copy import copy

from ocelot.cpbd.high_order import rk_field
from ocelot.cpbd.transformations.transfer_map import TransferMap


class RungeKuttaTM(TransferMap):
    def __init__(self, s_start=0, npoints=200):
        TransferMap.__init__(self)
        self.s_start = s_start
        self.npoints = npoints
        self.long_dynamics = True
        self.mag_field = lambda x, y, z: (0, 0, 0)

    def map_function(self, delta_length=None, length=None):
        return lambda X, energy: rk_field(X, self.s_start, delta_length if delta_length else self.length, self.npoints, energy, self.mag_field, self.long_dynamics)

    @classmethod
    def create_from_element(cls, element, params=None):
        tm = cls(s_start=element.s_start if hasattr(element, 's_start') else 0.,
                 npoints=element.npoints if hasattr(element, 'npoints') else 200)
        tm.mag_field = element.mag_field
        return tm
