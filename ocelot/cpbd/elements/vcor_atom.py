import numpy as np

from ocelot.cpbd.r_matrix import uni_matrix
from ocelot.cpbd.elements.cor_atom import CorAtom
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams


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
        s = 'Vcor : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l     =%8.4f m\n' % self.l
        s += 'angle =%8.4f deg\n' % (self.angle * 180.0 / np.pi)
        s += 'tilt  =%8.2f deg\n' % (self.tilt * 180.0 / np.pi)
        return s

    def create_first_order_main_params(self, energy: float, delta_length: float = None) -> FirstOrderParams:
        z = delta_length if delta_length is not None else self.l
        R = uni_matrix(z, 0, hx=0, sum_tilts=0, energy=energy)
        B = self.kick_b(z=z, l=self.l, angle_x=0., angle_y=self.angle)
        return FirstOrderParams(R, B, self.tilt)
