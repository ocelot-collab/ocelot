import numpy as np
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.tm_params.kick_params import KickParams
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.tm_params.second_order_params import SecondOrderParams
from ocelot.cpbd.high_order import t_nnn
from ocelot.cpbd.r_matrix import uni_matrix
from ocelot.cpbd.tm_utils import map_transform_with_offsets


class Magnet(Element):
    def __init__(self, eid=None, has_edge=False):
        super().__init__(eid=eid, has_edge=has_edge)
        self.angle = 0.  # Magnets Drift, Bend, Correctors (just angle)
        self.k1 = 0.  # Magnets quadropole
        self.k2 = 0.  # Magnets Sextupole

    def create_first_order_main_params(self, energy: float, delta_length: float = None) -> FirstOrderParams:
        k1 = self.k1
        if self.l == 0:
            hx = 0.
        else:
            hx = self.angle / self.l
        if delta_length is not None:
            R = uni_matrix(delta_length, k1, hx=hx, sum_tilts=0, energy=energy)
        else:
            R = uni_matrix(self.l, k1, hx=hx, sum_tilts=0, energy=energy)
        B = self._default_B(R)
        return FirstOrderParams(R, B, self.tilt)

    def create_second_order_main_params(self, energy: float, delta_length: float = 0.0) -> SecondOrderParams:
        T = t_nnn(delta_length if delta_length is not None else self.l, 0. if self.l == 0 else self.angle / self.l,
                  self.k1, self.k2, energy)
        first_order_params = self.create_first_order_main_params(energy, delta_length)
        m_off = np.array([[-self.dx], [0], [-self.dy], [0], [0], [0]])
        # alternative way
        # B = first_order_params.B + np.dot(np.dot(m_off.T, T), m_off)[0]
        B, R = map_transform_with_offsets(m_off, first_order_params.R, T)
        return SecondOrderParams(R, B, T, self.tilt, self.dx, self.dy)

    def create_kick_entrance_params(self) -> KickParams:
        return KickParams(dx=self.dx, dy=self.dy, angle=self.angle, tilt=self.tilt, k1=self.k1, k2=self.k2)

    def create_kick_main_params(self) -> KickParams:
        return KickParams(dx=self.dx, dy=self.dy, angle=self.angle, tilt=self.tilt, k1=self.k1, k2=self.k2)

    def create_kick_exit_params(self) -> KickParams:
        return KickParams(dx=self.dx, dy=self.dy, angle=self.angle, tilt=self.tilt, k1=self.k1, k2=self.k2)
