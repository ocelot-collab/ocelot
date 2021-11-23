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
        s = 'Octupole : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l    =%8.4f m\n' % self.l
        s += 'k3   =%8.3f 1/m^4\n' % self.k3
        s += 'tilt =%8.2f deg\n' % (self.tilt * 180.0 / np.pi)
        return s

    def create_kick_main_params(self) -> KickParams:
        return KickParams(dx=self.dx, dy=self.dy, angle=self.angle, tilt=self.tilt, k1=self.k1, k2=self.k2, k3=self.k3)