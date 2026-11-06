"""Boundary, residual, and compatibility contracts for matcher phase targets."""

import numpy as np
import pytest

from ocelot.cpbd.beam import Twiss
from ocelot.cpbd.elements import Cavity, Drift, Marker, Quadrupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.matcher import FloquetPhaseTarget, MatchProblem, PhaseAdvanceTarget
from ocelot.cpbd.optics import twiss


def _seed():
    tw = Twiss()
    tw.beta_x, tw.beta_y = 3.0, 5.0
    tw.alpha_x, tw.alpha_y = 0.2, -0.3
    tw.mux, tw.muy = 0.25, 0.4
    tw.E = 0.1
    return tw


@pytest.mark.parametrize("plane", ["x", "y"])
@pytest.mark.parametrize("internal", [False, True])
@pytest.mark.parametrize("same_element", [False, True])
def test_phase_target_includes_finite_boundaries(plane, internal, same_element):
    start = Quadrupole(l=0.7, k1=0.8, eid="Q")
    end = start if same_element else Drift(l=0.4, eid="D")
    prefix = [Cavity(l=0.3, v=0.02, phi=0.0, freq=1.3e9)] if internal else []
    section = [start] if same_element else [start, end]
    lat = MagneticLattice(prefix + section)
    seed = _seed()
    entrance = twiss(MagneticLattice(prefix), seed)[-1] if prefix else seed
    _, matrix, _ = lat.transfer_maps(energy=entrance.E, start=start, stop=end)
    i = 0 if plane == "x" else 2
    # Independent Courant-Snyder phase from the inclusive section matrix.
    expected = np.arctan2(
        matrix[i, i + 1],
        getattr(entrance, f"beta_{plane}") * matrix[i, i]
        - getattr(entrance, f"alpha_{plane}") * matrix[i, i + 1],
    ) % (2.0 * np.pi)
    problem = MatchProblem(lat, seed)
    problem.target_phase_advance(start, end, plane, expected, tol=1e-12)
    merit, reports, _, _ = problem.evaluate()
    assert merit == 0.0
    assert reports[0].met
    assert reports[0].details["actual"] == pytest.approx(expected)
    assert expected > 0.0


def test_phase_target_preserves_positional_options_and_normalizes_plane():
    q = Quadrupole(l=1.0, k1=1.0)
    problem = MatchProblem(MagneticLattice([q]), _seed())
    target = problem.target_phase_advance(
        q, q, " X ", 0.5, "phase", 3.0, 0.1, False, "tag", periodicity=np.pi
    )
    assert isinstance(target, PhaseAdvanceTarget)
    assert (target.name, target.weight, target.tol, target.active, target.tag) == (
        "phase", 3.0, 0.1, False, "tag"
    )
    assert target.plane == "x"
    assert target.periodicity == np.pi


@pytest.mark.parametrize("kwargs, message", [
    ({"periodicity": 0.0}, "finite and positive"),
    ({"periodicity": -1.0}, "finite and positive"),
    ({"periodicity": np.nan}, "finite and positive"),
    ({"periodicity": np.inf}, "finite and positive"),
    ({"periodicity": np.complex128(1 + 1j)}, "finite and positive real"),
    ({"value": np.nan}, "finite"),
    ({"value": np.inf}, "finite"),
    ({"value": 1j}, "target_floquet_phase"),
    ({"plane": "z"}, "Unsupported plane"),
])
def test_phase_target_rejects_invalid_inputs(kwargs, message):
    q = Quadrupole(l=1.0, k1=1.0)
    problem = MatchProblem(MagneticLattice([q]), _seed())
    args = dict(start=q, end=q, plane="x", value=1.0)
    args.update(kwargs)
    with pytest.raises(ValueError, match=message):
        problem.target_phase_advance(**args)


@pytest.mark.parametrize("periodicity", [np.pi, 2.0 * np.pi])
@pytest.mark.parametrize(
    "error, expected", [(0.05, 0.0), (0.100001, 1e-6), (-0.100001, -1e-6)]
)
def test_phase_target_wraps_before_applying_tolerance(periodicity, error, expected):
    q = Quadrupole(l=1.0, k1=1.0)
    problem = MatchProblem(MagneticLattice([q]), _seed())
    _, _, _, state = problem.evaluate()
    actual = state.twiss_end.mux - state.twiss_start.mux
    target = problem.target_phase_advance(
        q, q, "x", actual + 3 * periodicity - error, tol=0.1, periodicity=periodicity
    )
    assert target.residuals(state)[0] == pytest.approx(expected, abs=1e-14)
    report = target.report(state)
    assert report.tol == 0.1
    assert report.met == (expected == 0.0)
    assert report.details["actual_mod_periodicity"] == pytest.approx(
        (actual + periodicity / 2) % periodicity - periodicity / 2
    )
    assert report.details["diff_mod_periodicity"] == pytest.approx(error)


def test_phase_target_validates_occurrences_and_order():
    shared = Drift(l=0.4, eid="D")
    start, end = Marker(eid="S"), Marker(eid="E")
    lat = MagneticLattice([start, shared, shared, end])
    problem = MatchProblem(lat, _seed())
    _, _, _, state = problem.evaluate()
    target = PhaseAdvanceTarget(start, end, "x", 0.0)
    assert np.isfinite(target.residuals(state)).all()
    for a, b, message in [
        (shared, end, "indices"),
        (start, shared, "indices"),
        (Drift(l=0.4, eid="D"), end, "not present"),
        (end, start, "before start"),
    ]:
        with pytest.raises(ValueError, match=message):
            PhaseAdvanceTarget(a, b, "x", 0.0).residuals(state)


def test_phase_target_default_boundaries_and_twiss_delta_convention():
    q = Quadrupole(l=1.0, k1=1.0)
    problem = MatchProblem(MagneticLattice([q]), _seed())
    _, _, _, state = problem.evaluate()
    expected = state.twiss_end.mux - state.twiss_start.mux
    problem.target_phase_advance(None, None, "x", expected, tol=1e-12)
    # A Twiss difference remains exit-to-exit, unlike the inclusive target.
    problem.target_twiss_delta(q, q, "mux", 0.0, wrap_phase=True)
    merit, reports, _, _ = problem.evaluate()
    assert merit == 0.0
    assert all(report.met for report in reports)
    assert expected > 0.0


def test_floquet_class_remains_public_and_includes_single_element():
    q = Quadrupole(l=1.0, k1=1.0)
    problem = MatchProblem(MagneticLattice([q]), _seed())
    target = problem.target_floquet_phase(q, q, "x", 1.0, tol=1e-12)
    assert isinstance(target, FloquetPhaseTarget)
    problem.add_target(FloquetPhaseTarget(q, q, "y", 1j, tol=1e-12))
    merit, reports, _, _ = problem.evaluate()
    assert merit == 0.0
    assert all(report.met for report in reports)
