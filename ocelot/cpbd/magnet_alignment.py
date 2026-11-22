"""Rigid mechanical roll between the normal planes of two reference frames."""

import numpy as np

from ocelot.common.geometry import get_tilt_matrix
from ocelot.common.globals import m_e_GeV


class FaceTransform:
    """Source origin and axes expressed in destination coordinates."""

    def __init__(self, shift, axes):
        self.shift = shift
        self.axes = axes

    def inverse(self):
        return FaceTransform(-self.axes.T @ self.shift, self.axes.T)

    def _crossing(self, X, energy):
        # Straight flight from the source face to z=0 of the destination.
        inv_beta = 1. / np.sqrt(1. - (m_e_GeV / energy)**2) if energy else 1.
        px, py, delta = X[1], X[3], X[5]
        ps2 = 1. + delta * (2. * inv_beta + delta) - px**2 - py**2
        if np.any(ps2 <= 0):
            raise ValueError("Magnet alignment requires forward particle momentum")
        ps = np.sqrt(ps2)
        p = self.axes @ np.array([px, py, ps])
        if np.any(p[2] <= 0):
            raise ValueError("Particle must cross the rolled magnet face forwards")
        position = self.shift[:, None] + self.axes[:, :2] @ X[[0, 2]]
        flight = -position[2] / p[2]
        relative_energy = inv_beta + delta
        out = X.copy()
        out[[0, 2]] = position[:2] + p[:2] * flight
        out[[1, 3]] = p[:2]
        out[4] += relative_energy * flight
        return out, ps, p, relative_energy, flight

    def apply(self, X, energy):
        X[:] = self._crossing(X, energy)[0]

    def taylor(self, orbit, energy, second_order):
        """Value, Jacobian and quadratic coefficients about one incoming orbit."""
        value, ps, p, e, flight = self._crossing(orbit, energy)
        ps, e, flight = ps[0], e[0], flight[0]
        p = p[:, 0]
        # Differentiate ps^2 = 1 + 2 delta/beta + delta^2 - px^2 - py^2.
        dps = np.array([0., -orbit[1, 0], 0., -orbit[3, 0], 0., e]) / ps
        dp = np.outer(self.axes[:, 2], dps)
        dp[:, [1, 3]] += self.axes[:, :2]
        dr = np.zeros((3, 6))
        dr[:, [0, 2]] = self.axes[:, :2]
        df = -(dr[2] + flight * dp[2]) / p[2]
        # Position and arrival time both advance by flight times (px, py, E/p0).
        velocity = np.array([p[0], p[1], e])
        dv = np.vstack((dp[:2], [0., 0., 0., 0., 0., 1.]))
        R = np.eye(6)
        R[[0, 2]] = dr[:2]
        R[[0, 2, 4]] += flight * dv + np.outer(velocity, df)
        R[[1, 3]] = dp[:2]
        if not second_order:
            return value, R, None

        ps2_term = (np.diag([0., -1., 0., -1., 0., 1.]) - np.outer(dps, dps)) / (2. * ps)
        f2 = -(flight * self.axes[2, 2] * ps2_term + np.outer(df, dp[2])) / p[2]
        T = np.zeros((6, 6, 6))
        T[[1, 3]] = self.axes[:2, 2, None, None] * ps2_term
        T[[0, 2, 4]] = velocity[:, None, None] * f2 + dv[:, :, None] * df
        T[[0, 2]] += flight * T[[1, 3]]
        return value, R, (T + T.swapaxes(1, 2)) / 2.


def rolled_face(position, frame, dx, dy, dtilt):
    """Physical face in nominal face coordinates after a roll at the entrance."""
    roll = get_tilt_matrix(dtilt)
    # Rotate the arm from the entrance pivot, then express both its change
    # and the new axes in the local survey frame at this arc length.
    arm = position - np.array([dx, dy, 0.])
    return FaceTransform(frame.T @ ((roll - np.eye(3)) @ arm), frame.T @ roll @ frame)


def roll_map(params, entrance, exit_face, energy):
    """Express a native magnet map and its two boundary crossings in survey axes."""
    from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
    from ocelot.cpbd.tm_params.second_order_params import SecondOrderParams
    from ocelot.cpbd.tm_utils import transfer_maps_mult

    second = isinstance(params, SecondOrderParams)
    zero = np.zeros((6, 1))
    b, r, t = entrance.inverse().taylor(zero, energy, second)
    body_r = params.get_rotated_R()
    if not second:
        out, jacobian, _ = exit_face.taylor(params.B + body_r @ b, energy, False)
        return FirstOrderParams(jacobian @ body_r @ r, out, 0.)

    b, r, t = transfer_maps_mult(b, r, t, params.B, body_r, params.get_rotated_T())
    # Expand the exit transformation at the orbit delivered by the body.
    out, jacobian, quadratic = exit_face.taylor(b, energy, True)
    b, r, t = transfer_maps_mult(zero, r, t, out, jacobian, quadratic)
    return SecondOrderParams(r, b, t, 0., params.dx, params.dy)
