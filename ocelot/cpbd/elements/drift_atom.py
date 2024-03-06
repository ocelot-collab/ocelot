from ocelot.cpbd.elements.magnet import Magnet
from ocelot.cpbd.tm_params.runge_kutta_params import RungeKuttaParams


class DriftAtom(Magnet):
    """
    drift - free space
    l - length of drift in [m]
    """

    def __init__(self, l=0., eid=None):
        super().__init__(eid)
        self.l = l

    def __str__(self):
        s = 'Drift('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s

    def create_runge_kutta_main_params(self, energy):
        return RungeKuttaParams(mag_field=lambda x, y, z: (0., 0., 0,))
