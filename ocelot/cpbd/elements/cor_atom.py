import numpy as np
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.tm_params.second_order_params import SecondOrderParams
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.high_order import t_nnn
from ocelot.cpbd.r_matrix import uni_matrix


class CorAtom(Element):
    def __init__(self, l=0., angle=0., eid=None):
        super().__init__(eid=eid)
        self.l = l
        self.angle = angle

    def kick_b(self, z, l, angle, tilt):
        angle_x = angle * np.cos(tilt)
        angle_y = angle * np.sin(tilt)
        if l == 0:
            hx = 0.
            hy = 0.
        else:
            hx = angle_x / l
            hy = angle_y / l

        dx = hx * z * z / 2.
        dy = hy * z * z / 2.
        dx1 = hx * z if l != 0 else angle_x
        dy1 = hy * z if l != 0 else angle_y
        b = np.array([[dx], [dx1], [dy], [dy1], [0.], [0.]])
        return b

    def create_first_order_main_params(self, energy: float, delta_length: float = None) -> FirstOrderParams:
        z = delta_length if delta_length is not None else self.l
        R = uni_matrix(z, 0, hx=0, sum_tilts=0, energy=energy)
        B = self.kick_b(z=z, l=self.l, angle=self.angle, tilt=self.tilt)
        return FirstOrderParams(R, B, self.tilt)

    def create_second_order_main_params(self, energy: float, delta_length: float = 0) -> SecondOrderParams:
        params = self.create_first_order_main_params(energy, delta_length)
        T = t_nnn(delta_length if delta_length is not None else self.l, 0, 0, 0, energy)
        return SecondOrderParams(params.R, params.B, T, self.tilt, self.dx, self.dy)