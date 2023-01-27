import numpy as np

from ocelot.cpbd.elements.magnet import Magnet
from ocelot.cpbd.tm_params.kick_params import KickParams


class OctupoleAtom(Magnet):
    """
    octupole
    k3 - strength of octupole lens in [1/m^4],
    l - length of lens in [m].
    """

    def __init__(self, l=0., k3=0., tilt=0., eid=None):
        super().__init__(eid)
        self.l = l
        self.k3 = k3
        self.tilt = tilt

    def __str__(self):
        s = 'Octupole('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'k3=%8.6e, ' % self.k3 if np.abs(self.k3) > 1e-15 else ""
        s += 'tilt=%8.6e, ' % self.tilt if np.abs(self.tilt) > 1e-15 else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s

    def create_kick_main_params(self) -> KickParams:
        return KickParams(dx=self.dx, dy=self.dy, angle=self.angle, tilt=self.tilt, k1=self.k1, k2=self.k2, k3=self.k3)