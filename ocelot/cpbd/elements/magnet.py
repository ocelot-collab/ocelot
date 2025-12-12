import numpy as np
from ocelot.common.math_op import get_tilt_matrix
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.tm_params.kick_params import KickParams
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.tm_params.second_order_params import SecondOrderParams
from ocelot.cpbd.high_order import t_nnn
from ocelot.cpbd.r_matrix import uni_matrix
from ocelot.cpbd.tm_utils import map_transform_with_offsets


class Magnet(Element):
    def __init__(self, eid=None, has_edge=False, **kwargs):
        super().__init__(eid=eid, has_edge=has_edge, **kwargs)
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

    def get_transfer_geometry(self):
        # 1. Fallback for straight elements
        if self.angle == 0:
            return super().get_transfer_geometry()

        # 2. Curved Geometry Logic
        rho = 0.0
        if self.l != 0:
            rho = self.l / self.angle

        # --- A. End Point (Full Angle) ---
        ca = np.cos(self.angle)
        sa = np.sin(self.angle)

        R_end = np.array([
            rho * (ca - 1),
            0,
            rho * sa
        ])

        S_end = np.array([
            [ca, 0, -sa],
            [0, 1, 0],
            [sa, 0, ca]
        ])

        # --- B. Midpoint (Half Angle) - CRITICAL SECTION ---
        # Make sure you are NOT doing: R_mid = np.array([0, 0, self.l / 2.0])

        half_angle = self.angle / 2.0
        ca_mid = np.cos(half_angle)
        sa_mid = np.sin(half_angle)

        R_mid = np.array([
            rho * (ca_mid - 1),  # <--- This gives the X offset (Sagitta)
            0,
            rho * sa_mid
        ])

        # Rotation at midpoint
        S_mid = np.array([
            [ca_mid, 0, -sa_mid],
            [0, 1, 0],
            [sa_mid, 0, ca_mid]
        ])

        # --- C. Apply Tilt ---
        if self.tilt != 0:
            from ocelot.common.math_op import get_tilt_matrix
            T = get_tilt_matrix(self.tilt)

            R_end = T @ R_end
            R_mid = T @ R_mid

            S_end = T @ S_end @ T.T
            S_mid = T @ S_mid @ T.T

        return R_end, S_end, R_mid, S_mid