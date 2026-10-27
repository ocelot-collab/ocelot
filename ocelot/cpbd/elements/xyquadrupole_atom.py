import numpy as np

from ocelot.common.globals import m_e_GeV
from ocelot.cpbd.elements.magnet import Magnet
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams


def _quadrupole_functions(k, length):
    """Cosine- and sine-like solutions of u'' + k*u = 0, including k=0."""
    phase = np.sqrt(abs(k)) * length
    if k > 0:
        return np.cos(phase), length * np.sinc(phase / np.pi)
    if k < 0:
        return np.cosh(phase), np.sinh(phase) / np.sqrt(-k)
    return 1.0, length


def _arrival_integral(k, length, c, s):
    """Integral of z*C(z)*S(z), evaluated regularly near k=0."""
    v = k * length**2
    if abs(v) < 0.1:
        term = 1.0 / 3.0
        result = term
        for n in range(1, 9):
            term *= -4.0 * v / (2 * n * (2 * n + 3))
            result += term
        return length**3 * result
    return length * s**2 / 2.0 - (length - c * s) / (4.0 * k)


def _s_minus_lc(k, length, c, s):
    """S(L)-L*C(L), avoiding cancellation for small k*L**2."""
    v = k * length**2
    if abs(v) < 0.1:
        term = result = v / 3.0
        for n in range(2, 10):
            term *= -v / (2 * (n - 1) * (2 * n + 1))
            result += term
        return length * result
    return s - length * c


def _angle_arrival_integral(k, length, c, s):
    """Integral of S(z)*(S(z)-z*C(z))/k for the angle-squared timing term.

    Its regular limit is L**5/15 - 4*k*L**7/315 + ... . The closed
    expression subtracts terms of order k to obtain a term of order k**2.
    """
    v = k * length**2
    if abs(v) < 0.1:
        term = result = 1.0 / 15.0
        for n in range(1, 9):
            term *= -4 * v * (n + 1) / (n * (2 * n + 4) * (2 * n + 5))
            result += term
        return length**5 * result
    return (3 * (length - c * s) - 2 * k * length * s**2) / (4 * k**2)


class XYQuadrupoleAtom(Magnet):
    """Weak-curvature, first-order map about an offset and tilted reference.

    See ``XYQuadrupole`` for the public coordinate and validity contract.
    """

    def __init__(self, l=0., x_offs=0., y_offs=0., k1=0., tilt=0., eid=None,
                 x_angle=0., y_angle=0., **kwargs):
        super().__init__(eid=eid, **kwargs)
        self.l = l
        self.k1 = k1
        self.x_offs = x_offs
        self.y_offs = y_offs
        self.x_angle = x_angle
        self.y_angle = y_angle
        self._validate_angles()
        self.tilt = tilt

    def _validate_angles(self):
        norm2 = self.x_angle**2 + self.y_angle**2
        if not np.isfinite(norm2) or norm2 >= 1:
            raise ValueError("XYQuadrupole requires finite x_angle, y_angle with x_angle**2 + y_angle**2 < 1")

    def _entrance_frame(self):
        """Minimal rotation from the magnet axis to the entrance unit tangent."""
        self._validate_angles()
        tx, ty = self.x_angle, self.y_angle
        tz = np.sqrt(1 - tx**2 - ty**2)
        # Rodrigues rotation for E_Z -> (tx, ty, tz), with no additional roll.
        return np.array([[1 - tx**2 / (1 + tz), -tx * ty / (1 + tz), tx],
                         [-tx * ty / (1 + tz), 1 - ty**2 / (1 + tz), ty],
                         [-tx, -ty, tz]])

    def __str__(self):
        s = 'XYQuadrupole('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'k1=%8.6e, ' % self.k1 if np.abs(self.k1) > 1e-15 else ""
        s += 'x_offs=%8.6e, ' % self.x_offs if np.abs(self.x_offs) > 1e-15 else ""
        s += 'y_offs=%8.6e, ' % self.y_offs if np.abs(self.y_offs) > 1e-15 else ""
        s += 'x_angle=%8.6e, ' % self.x_angle if np.abs(self.x_angle) > 1e-15 else ""
        s += 'y_angle=%8.6e, ' % self.y_angle if np.abs(self.y_angle) > 1e-15 else ""
        s += 'tilt=%8.6e, ' % self.tilt if np.abs(self.tilt) > 1e-15 else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s

    def _linear_map(self, length, energy):
        """Entrance-to-length map, revised v4 Eqs. (37)-(42), in tau coordinates."""
        self._validate_angles()
        # As for other Ocelot magnetic maps, energy=0 denotes beta=1.
        igamma2 = (m_e_GeV / energy)**2 if energy != 0 else 0.0
        beta2 = 1.0 - igamma2
        beta = np.sqrt(beta2)
        k = self.k1
        cx, sx = _quadrupole_functions(k, length)
        cy, sy = _quadrupole_functions(-k, length)
        r = np.eye(6)
        geometric_timing = 0.0
        # Each plane obeys u'' + K*u = 0, with h=K*(offset*C+angle*S).
        for i, K, offset, angle, c, s in (
                (0, k, self.x_offs, self.x_angle, cx, sx),
                (2, -k, self.y_offs, self.y_angle, cy, sy)):
            difference = _s_minus_lc(K, length, c, s)
            dispersion = (K * offset * length * s + angle * difference) / (2 * beta)
            dispersion_prime = K * (offset * (s + length * c) + angle * length * s) / (2 * beta)
            r[i:i+2, i:i+2] = [[c, s], [-K * s, c]]
            r[i:i+2, 5] = [dispersion, dispersion_prime]
            # Symplectic timing row, Eq. (39), with Ocelot's tau=-sigma.
            r[4, i:i+2] = [K * s * dispersion + c * dispersion_prime,
                           s * dispersion_prime - c * dispersion]
            # Appendix A, Eq. (A18): offset^2, offset*angle and angle^2 terms.
            geometric_timing += 0.5 * (
                K**2 * offset**2 * _arrival_integral(K, length, c, s)
                + offset * angle * K * s * difference
                + K**2 * angle**2 * _angle_arrival_integral(K, length, c, s))
        r[4, 5] = (geometric_timing - length * igamma2) / beta2
        return r

    def create_r_matrix(self):
        return self._linear_map

    def get_T_z_e_func(self):
        return lambda z, energy: np.zeros((6, 6, 6))

    def create_first_order_main_params(self, energy: float, delta_length=None, *, start_length=0.0):
        length = self.l if delta_length is None else delta_length
        if length == 0:
            r = np.eye(6)
        elif start_length == 0:
            r = self._linear_map(length, energy)
        else:
            entrance = self._linear_map(start_length, energy)
            exit_map = self._linear_map(start_length + length, energy)
            r = np.linalg.solve(entrance.T, exit_map.T).T
        params = FirstOrderParams(r, np.zeros((6, 1)), self.tilt)
        # Alignment translations are in beam coordinates, after rolling R.
        params.B = self._default_B(params.get_rotated_R())
        return params

    def create_second_order_main_params(self, energy: float, delta_length=None):
        raise NotImplementedError(
            "XYQuadrupole implements only the analytical first-order map; "
            "a second-order map requires the cubic Hamiltonian."
        )

    def _reference_geometry(self, s_values):
        """Sample displacement and frame relative to the entrance reference.

        Appendix B reports d=W_in.T*(r-r_in), Q=W_in.T*W. Integrating
        these directly starts at d=0, Q=I: the constant entrance rotation
        cancels. The entrance angles still determine the varying curvatures.
        """
        from scipy.integrate import solve_ivp
        from ocelot.common.geometry import get_tilt_matrix

        self._validate_angles()
        if self.l == 0:
            return np.zeros((len(s_values), 3)), np.repeat(np.eye(3)[None], len(s_values), axis=0)

        def rhs(s, state):
            frame = state[3:].reshape(3, 3)
            cx, sx = _quadrupole_functions(self.k1, s)
            cy, sy = _quadrupole_functions(-self.k1, s)
            hx = self.k1 * (self.x_offs * cx + self.x_angle * sx)
            hy = -self.k1 * (self.y_offs * cy + self.y_angle * sy)
            omega = np.array([[0.0, 0.0, -hx], [0.0, 0.0, -hy], [hx, hy, 0.0]])
            return np.r_[frame[:, 2], (frame @ omega).ravel()]

        solution = solve_ivp(rhs, (0.0, self.l), np.r_[np.zeros(3), np.eye(3).ravel()],
                             t_eval=s_values, method="DOP853", rtol=1e-11, atol=1e-13)
        if not solution.success:
            raise RuntimeError(f"XYQuadrupole reference survey failed: {solution.message}")
        roll = get_tilt_matrix(self.tilt)
        positions = solution.y[:3].T @ roll.T
        frames = roll @ solution.y[3:].T.reshape(-1, 3, 3) @ roll.T
        return positions, frames

    def default_mag_field(self, energy):
        """Ideal quadrupole field in the reference entrance frame, for CSR RK.

        Transform positions into magnet coordinates and rotate B back. The
        incoming reference has offsets (x_offs,y_offs) and the minimal-rotation
        entrance frame of Appendix B. A tilted axis gives a nonzero Bz in
        reference entrance coordinates, even though B_Z=0 in the magnet.
        This field does not change the wrapper's analytical-only TM contract.
        """
        from ocelot.common.geometry import get_tilt_matrix

        to_magnet = self._entrance_frame() @ get_tilt_matrix(self.tilt).T
        brho_k = self._brho(energy) * self.k1

        def field(x, y, z):
            x_local = (to_magnet[0, 0] * (x - self.dx) + to_magnet[0, 1] * (y - self.dy)
                       + to_magnet[0, 2] * z + self.x_offs)
            y_local = (to_magnet[1, 0] * (x - self.dx) + to_magnet[1, 1] * (y - self.dy)
                       + to_magnet[1, 2] * z + self.y_offs)
            bx, by = brho_k * y_local, brho_k * x_local
            return tuple(to_magnet[0, i] * bx + to_magnet[1, i] * by for i in range(3))

        return field

    def get_transfer_geometry(self):
        """Survey the same varying-curvature Bishop reference used by the map."""
        self._validate_angles()
        if self.l == 0 or self.k1 == 0 or (self.x_offs == self.y_offs == self.x_angle == self.y_angle == 0):
            return super().get_transfer_geometry()
        positions, frames = self._reference_geometry([self.l / 2, self.l])
        return positions[1], frames[1], positions[0], frames[0]
