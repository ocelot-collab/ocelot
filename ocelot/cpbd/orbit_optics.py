"""Linear optics along a tracked reference trajectory."""

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
    orbit: np.ndarray
    matrix: np.ndarray
    energy: float
    s: float

    @classmethod
    def entrance(cls, seed: Twiss) -> Self:
        return cls(np.array([seed.x, seed.xp, seed.y, seed.yp, seed.tau, seed.p], dtype=float),
                   np.eye(6), seed.E, seed.s)

    def advance(self, tm: Transformation) -> "_OrbitMap":
        jacobian, _ = map_derivatives(tm, self.energy, self.orbit)
        orbit = tm.map_function(self.orbit[:, None].copy(), self.energy)[:, 0]
        length = tm.length if tm.delta_length is None else tm.delta_length
        return _OrbitMap(orbit, jacobian @ self.matrix, self.energy + tm.get_delta_e(), self.s + length)


def orbit_seed(tws: Twiss | None, particle: Particle) -> Twiss:
    """Copy the optics and set the six entrance coordinates in the tracking frame."""
    if not isinstance(particle, Particle):
        raise TypeError("For twiss(), pass reference_orbit=Particle() for a zero launch "
                        "or supply a Particle with entrance coordinates")
    seed = deepcopy(tws) if tws is not None else Twiss()
    if particle.E and seed.E and not np.isclose(particle.E, seed.E, rtol=1e-12, atol=0.):
        raise ValueError("Particle and Twiss specify different reference energies")
    seed.E = particle.E or seed.E
    seed.x, seed.xp, seed.y, seed.yp = particle.x, particle.px, particle.y, particle.py
    seed.tau, seed.p = particle.tau, particle.p
    return seed


def _matched_modes(matrix: np.ndarray) -> np.ndarray:
    """Return the two positive-action transverse modes with unit symplectic norm."""
    from ocelot.cpbd.optics import UnstableLatticeError

    values, vectors = np.linalg.eig(matrix)
    sx = np.imag(vectors[0].conj() * vectors[1])
    sy = np.imag(vectors[2].conj() * vectors[3])
    action = sx + sy
    keep = action > 1e-14
    if np.count_nonzero(keep) != 2 or not np.allclose(abs(values), 1., atol=1e-7, rtol=0.):
        raise UnstableLatticeError(np.trace(matrix[:2, :2]) / 2, np.trace(matrix[2:, 2:]) / 2,
                                   detail=f"transverse one-turn eigenvalues {values}")
    # Label the mode with the larger fractional horizontal action as mode x.
    order = np.argsort((sx[keep] / action[keep]))[::-1]
    return (vectors[:, keep] / np.sqrt(action[keep]))[:, order]


def closed_orbit_optics(
    lattice: MagneticLattice, tws: Twiss | None, reference: Particle | Literal["closed"]
) -> tuple[Twiss, np.ndarray]:
    """Solve the fixed-momentum transverse fixed point of the active one-turn map."""
    from scipy.optimize import root

    if isinstance(reference, str) and reference == "closed":
        reference = Particle()
    elif not isinstance(reference, Particle):
        raise TypeError('For periodic_twiss(), pass reference_orbit="closed" or supply a Particle as the initial guess')
    seed = orbit_seed(tws, reference)
    entrance = _OrbitMap.entrance(seed)

    def one_turn(x: np.ndarray) -> _OrbitMap:
        state = _OrbitMap(np.r_[x, entrance.orbit[4:]], np.eye(6), seed.E, seed.s)
        for element in lattice.sequence:
            for tm in element.tms:
                state = state.advance(tm)
        return state

    def residual(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        turn = one_turn(x)
        return turn.orbit[:4] - x, turn.matrix[:4, :4] - np.eye(4)

    solution = root(residual, entrance.orbit[:4], jac=True, tol=1e-11)
    turn = one_turn(solution.x)
    if not np.allclose(turn.orbit[:4], solution.x, atol=1e-11, rtol=0.):
        raise ValueError(f"Transverse closed orbit did not converge: {solution.message}")
    matrix = turn.matrix
    if (not np.isclose(turn.energy, seed.E, atol=1e-13, rtol=1e-12)
            or not np.isclose(turn.orbit[5], seed.p, atol=1e-12, rtol=0.)
            or not np.allclose(matrix[5], np.eye(6)[5], atol=1e-9, rtol=0.)
            or not np.allclose(matrix[:4, 4], 0., atol=1e-9, rtol=0.)):
        raise ValueError("Closed-orbit Twiss requires a fixed-momentum ring with time-independent transverse motion")

    modes = _matched_modes(matrix[:4, :4])
    # Map two unit-beta ellipses with the specified mode emittances into beam coordinates.
    basis = np.eye(6)
    basis[:4, :4] = np.column_stack((modes[:, 0].real, modes[:, 0].imag,
                                      modes[:, 1].real, modes[:, 1].imag))
    unit = Twiss(E=seed.E, beta_x=1., beta_y=1., emit_x=seed.emit_x, emit_y=seed.emit_y)
    unit._pending_emit_xn = seed._pending_emit_xn
    unit._pending_emit_yn = seed._pending_emit_yn
    matched = Twiss.track(basis, unit)
    matched.Dx, matched.Dxp, matched.Dy, matched.Dyp = np.linalg.solve(
        np.eye(4) - matrix[:4, :4], matrix[:4, 5])
    matched.x, matched.xp, matched.y, matched.yp = solution.x
    matched.tau, matched.p, matched.s = seed.tau, seed.p, seed.s
    matched.mux, matched.muy = seed.mux, seed.muy
    return matched, modes


def trace_orbit_optics(
    lattice: MagneticLattice,
    seed: Twiss,
    n_points: int | None,
    targets: Iterable[tuple[OpticElement, int]],
    modes: np.ndarray | None = None,
) -> list[Twiss]:
    """Transport the orbit and the full entrance-to-observation-point Jacobian."""
    dispersion = np.array([seed.Dx, seed.Dxp, seed.Dy, seed.Dyp, 0., 1.])
    initial_phase = np.array([seed.mux, seed.muy])

    def sample(state: _OrbitMap, previous: Twiss, eid: str) -> Twiss:
        # The cumulative Jacobian retains longitudinal feedback from intervening maps.
        out = Twiss.track(state.matrix, seed)
        out.x, out.xp, out.y, out.yp, out.tau, out.p = state.orbit
        out.s, out.E, out.id = state.s, state.energy, eid
        derivative = state.matrix @ dispersion
        # Express dispersion per unit of local relative momentum deviation.
        out.Dx, out.Dxp, out.Dy, out.Dyp = (derivative[:4] / derivative[5]
                                             if derivative[5] != 0. else np.full(4, np.nan))
        phase = np.array([out.mux, out.muy])
        if modes is not None:
            projected = (state.matrix[:4, :4] @ modes)[[0, 2], [0, 1]]
            phase = initial_phase + np.angle(projected / modes[[0, 2], [0, 1]])
        # Twiss.track returns a phase modulo one turn for a single cumulative matrix.
        last = np.array([previous.mux, previous.muy])
        advance = (phase - last) % (2 * np.pi)
        advance[np.isclose(advance, 2 * np.pi, atol=1e-13, rtol=0.)] = 0.
        out.mux, out.muy = last + advance
        return out

    grid = None if n_points is None else seed.s + np.linspace(0., lattice.totalLen, n_points)
    samples = [seed] if grid is None or len(grid) else []
    point = 1
    previous = seed
    state = _OrbitMap.entrance(seed)
    exits = []
    for elem in lattice.sequence:
        # Each interior sample starts from the same entrance. Slicing a nonlinear
        # map successively would make the exit trajectory depend on the plot grid.
        if grid is not None:
            end = state.s + elem.l
            while point < len(grid) and grid[point] < end - 1e-13:
                length = grid[point] - state.s
                partial = state
                if length > 1e-13:
                    for tm in elem.get_section_tms(length):
                        partial = partial.advance(tm)
                samples.append(sample(partial, previous, elem.id))
                point += 1
        for tm in elem.tms:
            state = state.advance(tm)
            previous = sample(state, previous, elem.id)
            if grid is None:
                samples.append(previous)
        exits.append(previous)
    if grid is not None:
        samples.extend([previous] * max(0, len(grid) - point))
    for elem, index in targets:
        elem.tws = exits[index]
    return samples
