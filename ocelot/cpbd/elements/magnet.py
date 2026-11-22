from __future__ import annotations

import numpy as np
from ocelot.common.globals import speed_of_light, m_e_GeV
from ocelot.common.geometry import get_tilt_matrix
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.r_matrix import uni_matrix


class Magnet(Element):
    """
    Atom base class for magnetic elements (quadrupoles, bends, solenoids, etc.).

    Physics Parameters
    ===================
    - angle: Bending angle or deflection [rad] (zero for quadrupoles, non-zero for bends)
    - k1: Quadrupole strength [1/m²] (focusing/defocusing)
    - k2: Sextupole strength [1/m³] (chromatic, nonlinear effects)
    - dtilt: Mechanical roll error about the displaced entrance axis [rad]

    Available Hook Methods
    =======================
    Inherited from Element with magnetic-specific implementations:

    - create_first_order_main_params(energy) → FirstOrderParams
    - create_second_order_main_params(energy) → SecondOrderParams
    - create_kick_entrance/main/exit_params() → KickParams

    Subclasses like BendAtom add edge hooks:
    - create_first_order_entrance/exit_params(energy)
    - create_second_order_entrance/exit_params(energy)

    See Also
    --------
    https://ocelot-collab.github.io/docs/docu/elements/magnet/
    Element : Base atom class
    """
    dtilt = 0.  # Default for standalone magnet atoms pickled before this field existed.

    def __init__(self, eid=None, has_edge=False, **kwargs):
        self.dtilt = kwargs.pop("dtilt", 0.)
        super().__init__(eid=eid, has_edge=has_edge, **kwargs)
        self.angle = 0.  # Magnets Drift, Bend, Correctors (just angle)
        self.k1 = 0.  # Magnets quadropole
        self.k2 = 0.  # Magnets Sextupole
        self.k3 = 0.  # Magnets Octupole

    @property
    def _map_tilt(self):
        # A curved reference also needs the mechanical change of normal planes.
        return self.tilt if self.angle else self.tilt + self.dtilt

    @staticmethod
    def _brho(energy: float) -> float:
        gamma = energy / m_e_GeV
        beta = np.sqrt(1. - 1. / (gamma * gamma)) if gamma != 0 else 1.
        return beta * energy * 1e9 / speed_of_light

    def _hard_edge_field_components(self, x, y, energy: float):
        brho = self._brho(energy)
        h = self.angle / self.l if self.l != 0. else 0.
        k1 = getattr(self, "k1", 0.)
        k2 = getattr(self, "k2", 0.)
        k3 = getattr(self, "k3", 0.)

        x2 = x * x
        y2 = y * y
        by = brho * (
            h
            + k1 * x
            + 0.5 * k2 * (x2 - y2)
            + k3 * (x * x2 - 3. * x * y2) / 6.
        )
        bx = brho * (
            k1 * y
            + k2 * x * y
            + k3 * (3. * x2 * y - y * y2) / 6.
        )
        return bx, by

    def _design_mag_field(self, energy: float):
        """
        Hard-edge magnetic field from the element strengths.

        The field is longitudinally uniform and follows the normal multipole
        convention used by the kick maps.
        """
        dx = self.dx
        dy = self.dy
        tilt = self.tilt
        cos_t = np.cos(tilt)
        sin_t = np.sin(tilt)

        def field(x, y, z):
            x_local = cos_t * (x - dx) + sin_t * (y - dy)
            y_local = -sin_t * (x - dx) + cos_t * (y - dy)
            bx_local, by_local = self._hard_edge_field_components(x_local, y_local, energy)

            bx = cos_t * bx_local - sin_t * by_local
            by = sin_t * bx_local + cos_t * by_local
            return bx, by, 0.

        return field

    def _roll_field(self, field):
        if not self.dtilt:
            return field
        dx, dy = self.dx, self.dy
        c, s = np.cos(self.dtilt), np.sin(self.dtilt)

        def rolled_field(x, y, z):
            bx, by, bz = field(dx + c*(x-dx) + s*(y-dy), dy - s*(x-dx) + c*(y-dy), z)
            return c*bx - s*by, s*bx + c*by, bz

        return rolled_field

    def default_mag_field(self, energy: float):
        return self._roll_field(self._design_mag_field(energy))

    def get_mag_field(self, energy: float):
        return self.default_mag_field(energy) if self.mag_field is None else self._roll_field(self.mag_field)

    def create_runge_kutta_main_params(self, energy):
        from ocelot.cpbd.tm_params.runge_kutta_params import RungeKuttaParams
        if not self.dtilt:
            return RungeKuttaParams(self.get_mag_field(energy))
        design_field = self._design_mag_field(energy) if self.mag_field is None else self.mag_field
        return RungeKuttaParams(self._roll_field(design_field), design_field=design_field,
                                alignment=(self.dx, self.dy, self.dtilt))

    def _body_matrix(self, length, energy):
        hx = self.angle / self.l if self.l else 0.
        return uni_matrix(length, self.k1, hx=hx, energy=energy)

    def create_first_order_main_params(self, energy, delta_length=None, *, start_length=0.):
        length = self.l if delta_length is None else delta_length
        params = self._aligned_map_params(self._body_matrix(length, energy), self._map_tilt)
        return self._body_roll(params, start_length, start_length + length, energy)

    def create_second_order_main_params(self, energy, delta_length=0.0, *, start_length=0.):
        from ocelot.cpbd.high_order import t_nnn
        length = self.l if delta_length is None else delta_length
        T = t_nnn(length, self.angle / self.l if self.l else 0., self.k1, self.k2, energy)
        params = self._aligned_map_params(self._body_matrix(length, energy), self._map_tilt, T)
        return self._body_roll(params, start_length, start_length + length, energy)

    def _body_roll(self, params, start, stop, energy):
        if self.angle and self.dtilt:
            from ocelot.cpbd.magnet_alignment import roll_map
            return roll_map(params, self._rolled_face(start), self._rolled_face(stop), energy)
        return params

    def _rolled_face(self, s):
        from ocelot.cpbd.magnet_alignment import rolled_face
        position, frame = self._design_face(s)
        return rolled_face(position, frame, self.dx, self.dy, self.dtilt)

    def create_kick_main_params(self) -> KickParams:
        from ocelot.cpbd.tm_params.kick_params import KickParams
        return KickParams(dx=self.dx, dy=self.dy, angle=self.angle, tilt=self._map_tilt,
                          k1=self.k1, k2=self.k2, k3=self.k3,
                          rolled_face=self._rolled_face if self.angle and self.dtilt else None)

    def create_kick_entrance_params(self) -> KickParams:
        return self.create_kick_main_params()

    def create_kick_exit_params(self) -> KickParams:
        return self.create_kick_main_params()

    def _design_face(self, s):
        """Survey origin and axes at arc length s, in the nominal entrance frame."""
        if self.angle == 0. or self.l == 0.:
            return super()._design_face(s)
        a = self.angle * (s / self.l)
        ca, sa = np.cos(a), np.sin(a)
        # Arc displacement, using sinc to retain precision at small angles.
        position = s * np.array([-np.sin(a / 2) * np.sinc(a / (2 * np.pi)), 0., np.sinc(a / np.pi)])
        axes = np.array([[ca, 0., -sa], [0., 1., 0.], [sa, 0., ca]])
        roll = get_tilt_matrix(self.tilt)
        return roll @ position, roll @ axes @ roll.T

    def _rk_nominal_face(self, s):
        # Even if the exit points forward again, a sector that has already
        # turned through 90 degrees cannot be integrated using Cartesian Z.
        if self.l and abs(self.angle * (s / self.l)) >= np.pi / 2:
            raise ValueError("Nominal RK reference must move forward in Cartesian Z")
        return self._design_face(s)

    def get_transfer_geometry(self):
        """Design survey at exit and midpoint, independent of installation errors."""
        return (*self._design_face(self.l), *self._design_face(self.l / 2.))
