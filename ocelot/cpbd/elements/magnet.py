import numpy as np
from ocelot.common.math_op import get_tilt_matrix
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.tm_params.kick_params import KickParams
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.tm_params.second_order_params import SecondOrderParams
from ocelot.cpbd.high_order import t_nnn
from ocelot.cpbd.r_matrix import linear_magnet_matrix
from ocelot.cpbd.tm_utils import map_transform_with_offsets


class Magnet(Element):
    def __init__(self, eid=None, has_edge=False, **kwargs):
        super().__init__(eid=eid, has_edge=has_edge, **kwargs)
        self.angle = 0.  # Magnets Drift, Bend, Correctors (just angle)
        self.k1 = 0.  # Magnets quadropole
        self.k2 = 0.  # Magnets Sextupole

    def _curvature(self) -> float:
        return 0.0 if self.l == 0 else self.angle / self.l

    def linear_r_main(self, energy: float, delta_length: float = None, xp=np):
        length = self.l if delta_length is None else delta_length
        return linear_magnet_matrix(length, self.k1, hx=self._curvature(), sum_tilts=0.0, energy=energy, xp=xp)

    def create_first_order_main_params(self, energy: float, delta_length: float = None) -> FirstOrderParams:
        R = self.linear_r_main(energy, delta_length, xp=np)
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
        """
        Geometry for magnets.

        - If angle == 0: falls back to straight-element geometry (Element.get_transfer_geometry()).
        - If angle != 0: uses a sector-like arc in the local x-s plane (y is vertical).
          Then applies TILT as a rotation around local s, which rotates the bending plane:
            * tilt = 0      -> horizontal bend (x-s plane)
            * tilt = pi/2   -> vertical bend (y-s plane)

        Returns:
            (R_end, S_end, R_mid, S_mid) in local coordinates.

        Important:
            - For dipoles, TILT *does* affect the reference trajectory (bending plane rotation),
              hence both R_* and S_* are rotated by the roll matrix.
            - For quadrupoles/sextupoles (angle==0), TILT should NOT change trajectory;
              it is handled in the map/field and/or drawing.
        """
        # 1) Straight / non-bending magnet
        if getattr(self, "angle", 0.0) == 0.0:
            return super().get_transfer_geometry()

        L = float(getattr(self, "l", 0.0))
        ang = float(getattr(self, "angle", 0.0))

        # Guard against division by zero; a "bend" with L==0 is degenerate
        if abs(ang) < 1e-15 or abs(L) < 1e-15:
            return super().get_transfer_geometry()

        # Curvature radius (sector-like)
        rho = L / ang

        # --- A) End point (full angle) in the *un-tilted* x-s plane ---
        ca = np.cos(ang)
        sa = np.sin(ang)

        R_end = np.array([rho * (ca - 1.0), 0.0, rho * sa], dtype=float)

        # Rotation of the local frame at the end of a bend (about local y)
        S_end = np.array([
            [ca, 0.0, -sa],
            [0.0, 1.0, 0.0],
            [sa, 0.0, ca]
        ], dtype=float)

        # --- B) Midpoint (half angle) ---
        half = 0.5 * ang
        ca2 = np.cos(half)
        sa2 = np.sin(half)

        R_mid = np.array([rho * (ca2 - 1.0), 0.0, rho * sa2], dtype=float)

        S_mid = np.array([
            [ca2, 0.0, -sa2],
            [0.0, 1.0, 0.0],
            [sa2, 0.0, ca2]
        ], dtype=float)

        # --- C) Apply TILT (roll about local s-axis) to rotate bending plane ---
        tilt = float(getattr(self, "tilt", 0.0))
        if tilt != 0.0:
            from ocelot.common.math_op import get_tilt_matrix
            T = get_tilt_matrix(tilt)  # roll about s

            # Rotate the arc displacement into the tilted plane
            R_end = T @ R_end
            R_mid = T @ R_mid

            # Similarity transform rotates the bend rotation into the tilted plane
            S_end = T @ S_end @ T.T
            S_mid = T @ S_mid @ T.T

        return R_end, S_end, R_mid, S_mid
