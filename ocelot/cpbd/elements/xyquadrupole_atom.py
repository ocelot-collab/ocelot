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
    """Integral of z*C(z)*S(z), with the regular series from Appendix C."""
    v = k * length**2
    if abs(v) < 0.1:
        term = 1.0 / 3.0
        result = term
        for n in range(1, 9):
            term *= -4.0 * v / (2 * n * (2 * n + 3))
            result += term
        return length**3 * result
    return length * s**2 / 2.0 - (length - c * s) / (4.0 * k)


class XYQuadrupoleAtom(Magnet):
    """Weak-curvature, first-order map about a parallel-entry reference.

    See ``XYQuadrupole`` for the public coordinate and validity contract.
    The reference curvatures are k1*x_offs*Cx(s) and -k1*y_offs*Cy(s).
    """

    def __init__(self, l=0., x_offs=0., y_offs=0., k1=0., tilt=0., eid=None, **kwargs):
        super().__init__(eid=eid, **kwargs)
        self.l = l
        self.k1 = k1
        self.x_offs = x_offs
        self.y_offs = y_offs
        self.tilt = tilt

    def __str__(self):
        s = 'XYQuadrupole('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'k1=%8.6e, ' % self.k1 if np.abs(self.k1) > 1e-15 else ""
        s += 'x_offs=%8.6e, ' % self.x_offs if np.abs(self.x_offs) > 1e-15 else ""
        s += 'y_offs=%8.6e, ' % self.y_offs if np.abs(self.y_offs) > 1e-15 else ""
        s += 'tilt=%8.6e, ' % self.tilt if np.abs(self.tilt) > 1e-15 else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s

    def _linear_map(self, length, energy):
        """Entrance-to-length map, Eqs. (37), (41), (43)-(44), in tau coordinates."""
        # As for other Ocelot magnetic maps, energy=0 denotes beta=1.
        igamma2 = (m_e_GeV / energy)**2 if energy != 0 else 0.0
        beta2 = 1.0 - igamma2
        beta = np.sqrt(beta2)
        k = self.k1
        cx, sx = _quadrupole_functions(k, length)
        cy, sy = _quadrupole_functions(-k, length)
        hx = k * self.x_offs
        hy = -k * self.y_offs

        r = np.eye(6)
        r[:2, :2] = [[cx, sx], [-k * sx, cx]]
        r[2:4, 2:4] = [[cy, sy], [k * sy, cy]]
        r[:4, 5] = [hx * length * sx, hx * (sx + length * cx),
                     hy * length * sy, hy * (sy + length * cy)]
        r[:4, 5] /= 2.0 * beta
        # tau=-sigma reverses the paper's timing row (R55 remains one).
        r[4, :4] = [hx * (length + cx * sx), hx * sx**2,
                     hy * (length + cy * sy), hy * sy**2]
        r[4, :4] /= 2.0 * beta
        ix = _arrival_integral(k, length, cx, sx)
        iy = _arrival_integral(-k, length, cy, sy)
        r[4, 5] = (0.5 * (hx**2 * ix + hy**2 * iy) - length * igamma2) / beta2
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

    def get_transfer_geometry(self):
        """Survey the model's curved reference and Bishop frame (Appendix D).

        Displacements start at the reference entrance, not the magnet axis.
        Integrate the varying rotations in order: their generators generally
        do not commute when both entrance offsets are nonzero.
        """
        if self.l == 0 or self.k1 == 0 or (self.x_offs == 0 and self.y_offs == 0):
            return super().get_transfer_geometry()

        from scipy.integrate import solve_ivp
        from ocelot.common.geometry import get_tilt_matrix

        def rhs(s, state):
            frame = state[3:].reshape(3, 3)
            cx, _ = _quadrupole_functions(self.k1, s)
            cy, _ = _quadrupole_functions(-self.k1, s)
            hx, hy = self.k1 * self.x_offs * cx, -self.k1 * self.y_offs * cy
            omega = np.array([[0.0, 0.0, -hx], [0.0, 0.0, -hy], [hx, hy, 0.0]])
            return np.r_[frame[:, 2], (frame @ omega).ravel()]

        solution = solve_ivp(rhs, (0.0, self.l), np.r_[np.zeros(3), np.eye(3).ravel()],
                             t_eval=[self.l / 2, self.l], method="DOP853", rtol=1e-11, atol=1e-13)
        if not solution.success:
            raise RuntimeError(f"XYQuadrupole reference survey failed: {solution.message}")
        roll = get_tilt_matrix(self.tilt)
        mid, end = solution.y.T
        return (roll @ end[:3], roll @ end[3:].reshape(3, 3) @ roll.T,
                roll @ mid[:3], roll @ mid[3:].reshape(3, 3) @ roll.T)
