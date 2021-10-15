import numpy as np
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.tm_params.second_order_params import SecondOrderParams
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.high_order import t_nnn


class CorAtom(Element):
    def __init__(self, l=0., angle=0., eid=None):
        super().__init__(eid=eid)
        self.l = l
        self.angle = angle

    def kick_b(self, z, l, angle_x, angle_y):
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
        raise NotImplementedError()

    def create_second_order_main_params(self, energy: float, delta_length: float) -> SecondOrderParams:
        params = self.create_first_order_main_params(energy, delta_length)
        T = t_nnn(delta_length if delta_length is not None else self.l, 0, 0, 0, energy)
        # TODO: Question: Why we have to set tilt = 0 ? Vcor is set to pi/2 in __init__.
        return SecondOrderParams(params.R, params.B, T, 0., self.dx, self.dy)
