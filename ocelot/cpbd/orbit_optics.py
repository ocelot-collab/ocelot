"""Twiss parameters and dispersion about the orbit a particle actually follows.

Ordinary Twiss propagation uses the first-order design maps. A displaced
quadrupole can steer the beam into a downstream sextupole, where the field
gradient depends on the beam's horizontal and vertical position. To include
this focusing, track a reference particle through the active maps and
linearize each map where the particle enters it. The selected tracking
method must include the field term responsible for the effect.

Use ``twiss(lattice, tws0, reference_orbit=particle)`` for a given launch.
Here ``tws0`` gives the entrance beam ellipse and dispersion. The particle
gives the coordinates of the trajectory about which the optics are calculated.
In a fixed-momentum ring, ``periodic_twiss(lattice, reference_orbit="closed")``
finds the closed orbit and matches the transverse optics about it.
"""

from copy import deepcopy
from dataclasses import dataclass
from typing import Iterable, Literal, Self

import numpy as np

from ocelot.cpbd.beam import Particle, Twiss
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.transformations.linearization import map_derivatives
from ocelot.cpbd.transformations.transformation import Transformation


@dataclass
class _OrbitMap:
    """Reference particle and linear transport from the entrance to its current face.

    ``orbit`` holds (x, px, y, py, tau, delta) at ``s`` [m] and reference
    energy ``energy`` [GeV]. ``matrix`` maps a small entrance perturbation
    to its six-dimensional separation from this particle at the current face.
    Both ``trace_orbit_optics`` and the one-turn search in ``closed_orbit_optics``
    create this internal state at the lattice entrance with an identity matrix
    and advance it through the tracking maps.
    """

    orbit: np.ndarray
    matrix: np.ndarray
    energy: float
    s: float

    @classmethod
    def entrance(cls, tws0: Twiss) -> Self:
        return cls(np.array([tws0.x, tws0.xp, tws0.y, tws0.yp, tws0.tau, tws0.p], dtype=float),
                   np.eye(6), tws0.E, tws0.s)

    def advance(self, tm: Transformation) -> "_OrbitMap":
        jacobian, _ = map_derivatives(tm, self.energy, self.orbit)
        orbit = tm.map_function(self.orbit[:, None].copy(), self.energy)[:, 0]
        length = tm.length if tm.delta_length is None else tm.delta_length
        return _OrbitMap(orbit, jacobian @ self.matrix, self.energy + tm.get_delta_e(), self.s + length)


def twiss_at_launch(tws: Twiss | None, particle: Particle) -> Twiss:
    """Copy the entrance beam optics and set the reference orbit from the particle.

    The Twiss parameters describe the beam around this trajectory. The particle
    supplies its six launch coordinates and can also supply the reference energy.
    """
    if not isinstance(particle, Particle):
        raise TypeError("For twiss(), pass reference_orbit=Particle() for a zero launch "
                        "or supply a Particle with entrance coordinates")
    tws0 = deepcopy(tws) if tws is not None else Twiss()
    if particle.E and tws0.E and not np.isclose(particle.E, tws0.E, rtol=1e-12, atol=0.):
        raise ValueError("Particle and Twiss specify different reference energies")
    tws0.E = particle.E or tws0.E
    tws0.x, tws0.xp, tws0.y, tws0.yp = particle.x, particle.px, particle.y, particle.py
    tws0.tau, tws0.p = particle.tau, particle.p
    return tws0


def _matched_modes(matrix: np.ndarray) -> np.ndarray:
    """Find stable modes, with Re(v)ᵀ J Im(v) = 1 in (x, px, y, py)."""
    from ocelot.cpbd.optics import UnstableLatticeError

    J2 = np.array([[0., 1.], [-1., 0.]])
    J = np.kron(np.eye(2), J2)
    # An off-axis truncated SecondTM can have a small area defect (the
    # modulus tolerance below is 1e-7). Unit eigenvalue moduli alone do
    # not make the eigenvectors a canonical change of coordinates.
    if not np.allclose(matrix.T @ J @ matrix, J, rtol=0., atol=2e-7):
        raise ValueError("Closed-orbit normal modes require a symplectic transverse one-turn map")

    resolution = 64 * np.finfo(float).eps * max(1., np.linalg.norm(matrix))
    if np.max(abs(matrix[:2, 2:])) <= resolution and np.max(abs(matrix[2:, :2])) <= resolution:
        # An exactly decoupled resonant ring has two known physical planes.
        # Roundoff from a roll by pi/2 must not select arbitrary mixtures.
        matrix = matrix.copy()
        matrix[:2, 2:] = matrix[2:, :2] = 0.
    values, vectors = np.linalg.eig(matrix)
    upper = values.imag > resolution
    half_trace_x = np.trace(matrix[:2, :2]) / 2
    half_trace_y = np.trace(matrix[2:, 2:]) / 2
    if np.count_nonzero(upper) != 2 or not np.allclose(abs(values), 1., atol=1e-7, rtol=0.):
        raise UnstableLatticeError(half_trace_x, half_trace_y, detail=f"transverse one-turn eigenvalues {values}")

    modes = vectors[:, upper]
    separation = abs(values[upper][0] - values[upper][1])
    if separation <= np.sqrt(resolution):
        # A defective repeated eigenvalue can split by O(sqrt(eps)) in eig().
        # At its midpoint a Jordan block has only one independent eigenvector.
        midpoint = np.mean(values[upper])
        singular = np.linalg.svd(matrix - midpoint * np.eye(4), compute_uv=False)
        if singular[-1] <= resolution and singular[-2] > np.sqrt(resolution):
            raise UnstableLatticeError(half_trace_x, half_trace_y, detail="defective transverse one-turn eigenspace")
    # At coincident tunes eig() may return two nonorthogonal vectors in
    # the same eigenspace. Choose its horizontal- and vertical-like axes
    # by extremizing the horizontal fraction of the symplectic action.
    if separation <= resolution:
        from scipy.linalg import eigh

        action = modes.conj().T @ J @ modes / (2j)
        horizontal = modes[:2].conj().T @ J2 @ modes[:2] / (2j)
        signs = np.linalg.eigvalsh(action)
        if signs[0] * signs[1] > 0:
            sign = np.sign(signs[0])
            _, axes = eigh(sign * horizontal, sign * action)
            modes = modes @ axes
        elif matrix[:2, 2:].any() or matrix[2:, :2].any():
            raise ValueError("Coincident coupled tunes of opposite action do not define unique mode emittances")

    # Conjugate a negative-action eigenvector to retain the tune above pi.
    horizontal = np.imag(modes[0].conj() * modes[1])
    action = horizontal + np.imag(modes[2].conj() * modes[3])
    if np.any(abs(action) <= np.finfo(float).eps):
        raise UnstableLatticeError(half_trace_x, half_trace_y, detail="vanishing transverse mode action")
    modes[:, action < 0.] = modes[:, action < 0.].conj()
    modes = modes / np.sqrt(abs(action))
    fraction = horizontal / action
    if abs(fraction[0] - fraction[1]) <= 1e-12:
        # Equal plane participation leaves no horizontal-like mode. Put
        # the smaller fractional tune first to make unequal weights reproducible.
        phases = (np.sign(action) * np.angle(values[upper])) % (2 * np.pi)
        order = np.argsort(phases)
    else:
        order = np.argsort(fraction)[::-1]
    modes = modes[:, order]
    # Roundoff mixes nearly resonant eigenvectors. Remove the component in
    # the first mode's canonical plane before assigning unequal emittances.
    q, p = modes[:, 0].real, modes[:, 0].imag
    v = modes[:, 1]
    v = v + q * (p @ J @ v) - p * (q @ J @ v)
    norm = v.real @ J @ v.imag
    if norm <= np.finfo(float).eps:
        raise ValueError("Transverse normal modes are numerically indistinguishable")
    modes[:, 1] = v / np.sqrt(norm)
    return modes


def closed_orbit_optics(
    lattice: MagneticLattice, tws: Twiss | None, reference: Particle | Literal["closed"]
) -> tuple[Twiss, np.ndarray]:
    """Find the entrance orbit, matched optics and two transverse normal modes.

    The particle's four transverse coordinates give the initial guess. The
    solver varies them at fixed delta until they repeat after one turn.
    Arrival time may slip in a magnetic ring.
    The one-turn Jacobian supplies the matched dispersion and mode vectors.
    """
    from scipy.optimize import root

    if isinstance(reference, str) and reference == "closed":
        reference = Particle()
    elif not isinstance(reference, Particle):
        raise TypeError('For periodic_twiss(), pass reference_orbit="closed" or supply a Particle as the initial guess')
    tws0 = twiss_at_launch(tws, reference)
    orbit_guess = np.array([reference.x, reference.px, reference.y, reference.py], dtype=float)

    def one_turn(x: np.ndarray, tau: float = tws0.tau) -> _OrbitMap:
        state = _OrbitMap(np.r_[x, tau, tws0.p], np.eye(6), tws0.E, tws0.s)
        for element in lattice.sequence:
            for tm in element.tms:
                state = state.advance(tm)
        return state

    def residual(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        turn = one_turn(x)
        return turn.orbit[:4] - x, turn.matrix[:4, :4] - np.eye(4)

    solution = root(residual, orbit_guess, jac=True, tol=1e-11)
    turn = one_turn(solution.x)
    if not np.allclose(turn.orbit[:4], solution.x, atol=1e-11, rtol=0.):
        raise ValueError(f"Transverse closed orbit did not converge: {solution.message}")
    matrix = turn.matrix
    if (not np.isclose(turn.energy, tws0.E, atol=1e-13, rtol=1e-12)
            or not np.isclose(turn.orbit[5], tws0.p, atol=1e-12, rtol=0.)
            or not np.allclose(matrix[5], np.eye(6)[5], atol=1e-9, rtol=0.)
            or not np.allclose(matrix[:4, 4], 0., atol=1e-9, rtol=0.)):
        raise ValueError("Closed-orbit Twiss requires a fixed-momentum ring with time-independent transverse motion")

    tau_slip = turn.orbit[4] - tws0.tau
    if abs(tau_slip) > 1e-12:
        # A magnetic ring may slip in arrival time. A time-dependent kick at
        # a stationary phase can nevertheless vanish on the first turn and
        # change the closed orbit or its focusing at the next arrival.
        following = one_turn(solution.x, turn.orbit[4])
        expected = turn.orbit.copy()
        expected[4] += tau_slip
        if (not np.allclose(following.orbit, expected, atol=1e-11, rtol=0.)
                or not np.allclose(following.matrix, matrix, atol=1e-9, rtol=1e-9)):
            raise ValueError("Closed orbit or its linear optics changes with the arrival time on the next turn")

    modes = _matched_modes(matrix[:4, :4])
    # Map two unit-beta ellipses with the specified mode emittances into beam coordinates.
    basis = np.eye(6)
    basis[:4, :4] = np.column_stack((modes[:, 0].real, modes[:, 0].imag,
                                      modes[:, 1].real, modes[:, 1].imag))
    unit = Twiss(E=tws0.E, beta_x=1., beta_y=1., emit_x=tws0.emit_x, emit_y=tws0.emit_y)
    unit._pending_emit_xn = tws0._pending_emit_xn
    unit._pending_emit_yn = tws0._pending_emit_yn
    matched = Twiss.track(basis, unit)
    matched.Dx, matched.Dxp, matched.Dy, matched.Dyp = np.linalg.solve(
        np.eye(4) - matrix[:4, :4], matrix[:4, 5])
    matched.x, matched.xp, matched.y, matched.yp = solution.x
    matched.tau, matched.p, matched.s = tws0.tau, tws0.p, tws0.s
    matched.mux, matched.muy = tws0.mux, tws0.muy
    return matched, modes


def _projected_mode_advance(old: np.ndarray, new: np.ndarray, plane: int) -> float:
    """Phase change of a transported eigenvector viewed in x (0) or y (1)."""
    row = 2 * plane
    if abs(new[row]) <= 64 * np.finfo(float).eps * np.linalg.norm(new):
        return np.nan  # Phase is undefined where the position projection vanishes.
    ratio = new[row] / old[row]
    advance = np.angle(ratio)
    other = slice(2, 4) if plane == 0 else slice(0, 2)
    in_one_plane = not (new[other].any() or old[other].any())
    # A decoupled plane or complete mode rotation has a forward eigenphase
    # in [0, 2 pi). An arbitrary coupled projection can rotate backwards.
    # Follow its nearest branch between observations.
    if advance < -1e-13 and (in_one_plane or np.allclose(new, ratio * old, rtol=1e-12, atol=1e-14)):
        advance += 2 * np.pi
    return advance


def trace_orbit_optics(
    lattice: MagneticLattice,
    tws0: Twiss,
    n_points: int | None,
    targets: Iterable[tuple[OpticElement, int]],
    modes: np.ndarray | None = None,
) -> list[Twiss]:
    """Track the reference particle and report linear optics along its path.

    Each sample combines the reference orbit with the beam ellipse propagated
    by the cumulative Jacobian. Dispersion gives the changes in transverse
    position and momentum per unit change in the local relative energy coordinate.
    """
    dispersion = np.array([tws0.Dx, tws0.Dxp, tws0.Dy, tws0.Dyp, 0., 1.])
    last_phase = np.array([tws0.mux, tws0.muy], dtype=float)
    last_modes = modes.copy() if modes is not None else None

    def sample(state: _OrbitMap, eid: str) -> Twiss:
        # The cumulative Jacobian retains longitudinal feedback from intervening maps.
        out = Twiss.track(state.matrix, tws0)
        out.x, out.xp, out.y, out.yp, out.tau, out.p = state.orbit
        out.s, out.E, out.id = state.s, state.energy, eid
        derivative = state.matrix @ dispersion
        # Express dispersion per unit of local relative momentum deviation.
        out.Dx, out.Dxp, out.Dy, out.Dyp = (derivative[:4] / derivative[5]
                                             if derivative[5] != 0. else np.full(4, np.nan))
        if modes is None:
            advance = (np.array([out.mux, out.muy]) - last_phase) % (2 * np.pi)
            advance[np.isclose(advance, 2 * np.pi, atol=1e-13, rtol=0.)] = 0.
            last_phase[:] += advance
            out.mux, out.muy = last_phase
        else:
            transported = state.matrix[:4, :4] @ modes
            phase = np.full(2, np.nan)
            for i in range(2):
                advance = _projected_mode_advance(last_modes[:, i], transported[:, i], i)
                if np.isfinite(advance):
                    last_phase[i] += advance
                    phase[i] = last_phase[i]
                    last_modes[:, i] = transported[:, i]
            out.mux, out.muy = phase
        return out

    grid = None if n_points is None else tws0.s + np.linspace(0., lattice.totalLen, n_points)
    samples = [tws0] if grid is None or len(grid) else []
    point = 1
    previous = tws0
    state = _OrbitMap.entrance(tws0)
    exits = []
    for elem in lattice.sequence:
        entrance = state
        if grid is not None and elem.l > 0 and point < len(grid) and abs(grid[point] - entrance.s) <= 1e-13:
            samples.append(previous)
            point += 1
        for tm in elem.tms:
            # Sample a body in beam order, after its entrance face and before
            # its exit. Each partial map starts at the original element
            # entrance so the plotted grid cannot change the tracked orbit.
            if grid is not None and tm.length:
                while point < len(grid) and grid[point] < entrance.s + elem.l - 1e-13:
                    length = grid[point] - entrance.s
                    if length > 1e-13:
                        partial = entrance
                        for section in elem.get_section_tms(length):
                            partial = partial.advance(section)
                        previous = sample(partial, elem.id)
                    samples.append(previous)
                    point += 1
            state = state.advance(tm)
            previous = sample(state, elem.id)
            if grid is None:
                samples.append(previous)
        exits.append(previous)
    if grid is not None:
        samples.extend([previous] * max(0, len(grid) - point))
    for elem, index in targets:
        elem.tws = exits[index]
    return samples
