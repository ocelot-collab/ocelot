import importlib.util
import json
import os
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Dict, List

import numpy as np
import pytest

from ocelot import MagneticLattice, twiss
from ocelot.cpbd.beam import Twiss
from ocelot.cpbd.elements import Bend, Cavity, Drift, Marker, Quadrupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice as SimpleMagneticLattice
from ocelot.cpbd.matcher import MatchProblem, Objective


BENCH_BASELINE_FILE = Path(__file__).with_name("matcher_benchmark_baseline.json")

LENSES_VAR4B = [
    0.214716163625,
    -0.142348719548,
    0.24641054782,
    -0.258945660977,
    0.0,
    -0.0,
    0.140650315717,
    -0.28678140370010374,
    0.163714891881,
    -0.272437886893,
]


@dataclass
class RunResult:
    duration_s: float
    success: bool
    merit: float
    end_beta_x: float
    end_beta_y: float
    end_Dx: float
    max_beta_x: float
    max_beta_y: float
    r23: float
    variables: Dict[str, float]


def _import_module_from_path(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot create module spec for {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def t3_module():
    t3_file = Path(__file__).with_name("t3.py")
    if not t3_file.exists():
        pytest.skip(f"Local lattice file not found: {t3_file}")
    return _import_module_from_path("matcher_test_t3", t3_file)


def _twiss_seed():
    tw0 = Twiss()
    tw0.beta_x = 10.0
    tw0.beta_y = 10.0
    tw0.alpha_x = 0.0
    tw0.alpha_y = 0.0
    tw0.E = 1.0
    return tw0


def _simple_drift_lattice(length=1.0):
    start = Marker(eid="START")
    dvar = Drift(l=length, eid="DVAR")
    end = Marker(eid="END")
    lat = SimpleMagneticLattice((start, dvar, end))
    return lat, start, dvar, end


def _simple_cavity_lattice(v=0.02, phi=0.0):
    start = Marker(eid="START")
    cav = Cavity(l=0.5, v=v, phi=phi, freq=1.3e9, eid="CAV")
    end = Marker(eid="END")
    lat = SimpleMagneticLattice((start, cav, end))
    return lat, start, cav, end


def test_twiss_target_on_element():
    lat, _start, dvar, end = _simple_drift_lattice(0.7)
    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(dvar, "l", limits=(0.0, None), name="L_DVAR")
    problem.target_twiss(end, "s", value=2.2, weight=1.0e6, tol=0.0)

    result = problem.solve(solver="ls_trf", max_iter=80, tol=1.0e-10)
    assert result.success
    assert np.isclose(dvar.l, 2.2, atol=1.0e-4)


def test_vary_drift_length_with_finite_limits():
    lat, _start, dvar, end = _simple_drift_lattice(3.0)
    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(dvar, quantity="l", limits=(2.0, 5.0), name="D1_l")
    problem.target_twiss(end, "s", value=4.2, weight=1.0e6, tol=0.0)

    result = problem.solve(solver="ls_trf", max_iter=120, tol=1.0e-12)
    assert result.success
    assert np.isclose(dvar.l, 4.2, atol=1.0e-4)
    assert 2.0 <= result.variables["D1_l"] <= 5.0


def test_vary_bend_angle_with_limits():
    start = Marker(eid="START")
    bvar = Bend(l=1.0, angle=0.20, eid="BVAR")
    end = Marker(eid="END")
    lat = SimpleMagneticLattice((start, bvar, end))

    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(bvar, "angle", limits=(0.15, 0.30), name="BVAR_angle")
    problem.target_twiss(end, "Dx", value=0.12, weight=1.0e6, tol=0.0)

    result = problem.solve(solver="ls_trf", max_iter=120, tol=1.0e-10)
    assert result.success
    assert np.isclose(bvar.angle, 0.2411665855, atol=1.0e-4)
    assert 0.15 <= result.variables["BVAR_angle"] <= 0.30


def test_bounds_with_unsupported_solver_raise_error():
    lat, _start, dvar, end = _simple_drift_lattice(1.0)
    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(dvar, "l", limits=(0.0, 5.0))
    problem.target_twiss(end, "s", value=2.0, weight=1.0e6)

    with pytest.raises(ValueError, match="does not support bounds"):
        problem.solve(solver="bfgs", max_iter=40, tol=1.0e-10)


def test_initial_guess_outside_bounds_reports_variable_name():
    lat, _start, dvar, end = _simple_drift_lattice(1.0)
    dvar.l = -0.5

    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(dvar, "l", limits=(0.0, 5.0), name="DVAR_l")
    problem.target_twiss(end, "s", value=2.0, weight=1.0e6)

    with pytest.raises(ValueError, match="DVAR_l"):
        problem.solve(solver="ls_trf", max_iter=40, tol=1.0e-10)


def test_vary_cavity_voltage_to_match_end_energy():
    lat, _start, cav, end = _simple_cavity_lattice(v=0.02, phi=0.0)
    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(cav, "v", limits=(0.0, 0.05), name="CAV_v")
    problem.target_twiss(end, "E", value=1.03, weight=1.0e6, tol=0.0)

    result = problem.solve(solver="ls_trf", max_iter=120, tol=1.0e-12)
    assert result.success
    assert np.isclose(cav.v, 0.03, atol=1.0e-6)
    assert 0.0 <= result.variables["CAV_v"] <= 0.05


def test_vary_cavity_phase_to_match_end_energy():
    lat, _start, cav, end = _simple_cavity_lattice(v=0.03, phi=20.0)
    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(cav, "phi", limits=(0.0, 60.0), name="CAV_phi")
    problem.target_twiss(end, "E", value=1.0 + 0.03 * np.cos(np.deg2rad(30.0)), weight=1.0e6, tol=0.0)

    result = problem.solve(solver="ls_trf", max_iter=120, tol=1.0e-12)
    assert result.success
    assert np.isclose(cav.phi, 30.0, atol=1.0e-4)
    assert 0.0 <= result.variables["CAV_phi"] <= 60.0


def test_twiss_delta_target():
    lat, start, dvar, end = _simple_drift_lattice(0.9)
    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(dvar, "l", limits=(0.0, None))
    problem.target_twiss_delta(start, end, "s", value=1.8, weight=1.0e6)

    result = problem.solve(solver="ls_trf", max_iter=80, tol=1.0e-10)
    assert result.success
    assert np.isclose(dvar.l, 1.8, atol=1.0e-4)


def test_twiss_delta_wrap_phase_equivalence():
    lat, start, _dvar, end = _simple_drift_lattice(1.2)

    # Build one reference state to get the actual phase delta.
    probe = MatchProblem(lat, _twiss_seed())
    _m, _t, _o, state = probe.evaluate()
    dmu_actual = state.twiss_at(end).mux - state.twiss_at(start).mux

    # A phase-equivalent target shifted by -2*pi.
    target_equiv = dmu_actual - 2.0 * np.pi

    problem_wrapped = MatchProblem(lat, _twiss_seed())
    problem_wrapped.target_twiss_delta(
        start,
        end,
        "mux",
        value=target_equiv,
        relation="==",
        wrap_phase=True,
        weight=1.0,
    )
    _m_wrap, reports_wrap, _o_wrap, _s_wrap = problem_wrapped.evaluate()

    problem_raw = MatchProblem(lat, _twiss_seed())
    problem_raw.target_twiss_delta(
        start,
        end,
        "mux",
        value=target_equiv,
        relation="==",
        wrap_phase=False,
        weight=1.0,
    )
    _m_raw, reports_raw, _o_raw, _s_raw = problem_raw.evaluate()

    assert np.isclose(reports_wrap[0].residual_norm, 0.0, atol=1.0e-12)
    assert reports_raw[0].residual_norm > np.pi


def test_periodic_twiss_target_between_elements():
    lat, start, _dvar, end = _simple_drift_lattice(1.0)
    tw0 = _twiss_seed()

    problem = MatchProblem(lat, tw0)
    problem.vary_twiss("alpha_x", limits=(-1.0, 1.0), name="alpha_x0")
    problem.target_periodic_twiss("beta_x", start, end, weight=1.0e6, tol=1.0e-10)

    result = problem.solve(solver="ls_trf", max_iter=120, tol=1.0e-12)
    assert result.success

    _merit, reports, _objectives, state = problem.evaluate()
    assert reports[0].details["type"] == "twiss_periodic"
    assert np.isclose(
        state.twiss_at(end).beta_x,
        state.twiss_at(start).beta_x,
        atol=1.0e-8,
    )


def test_target_periodic_twiss_defaults_to_lattice_boundaries():
    lat, _start, _dvar, _end = _simple_drift_lattice(1.0)
    tw0 = _twiss_seed()

    problem = MatchProblem(lat, tw0)
    problem.vary_twiss("alpha_y", limits=(-1.0, 1.0), name="alpha_y0")
    problem.target_periodic_twiss("beta_y", weight=1.0e6, tol=1.0e-10)

    result = problem.solve(solver="ls_trf", max_iter=120, tol=1.0e-12)
    assert result.success

    _merit, reports, _objectives, state = problem.evaluate()
    assert reports[0].details["start"] == "twiss_start"
    assert reports[0].details["end"] == "twiss_end"
    assert np.isclose(state.twiss_end.beta_y, state.twiss_start.beta_y, atol=1.0e-8)


def test_partial_periodic_twiss_can_mix_with_strict_targets():
    start = Marker(eid="START")
    d0 = Drift(l=0.5, eid="D0")
    q1 = Quadrupole(l=0.2, k1=0.5, eid="Q1")
    d1 = Drift(l=0.5, eid="D1")
    q2 = Quadrupole(l=0.2, k1=-0.5, eid="Q2")
    d2 = Drift(l=0.5, eid="D2")
    end = Marker(eid="END")
    lat = SimpleMagneticLattice((start, d0, q1, d1, q2, d2, end))

    tw0 = Twiss()
    tw0.beta_x = 9.0
    tw0.beta_y = 10.0
    tw0.E = 1.0

    problem = MatchProblem(lat, tw0, periodic=False)
    problem.vary_element(q1, quantity="k1", limits=(-5.0, 5.0), name="Q1.k1")
    problem.vary_element(q2, quantity="k1", limits=(-5.0, 5.0), name="Q2.k1")
    problem.vary_twiss("alpha_x", limits=(-5.0, 5.0), name="twiss0.alpha_x")
    problem.vary_twiss("alpha_y", limits=(-5.0, 5.0), name="twiss0.alpha_y")

    problem.target_twiss(end, "alpha_x", 0.0, weight=1.0e6, tol=1.0e-10)
    problem.target_twiss(end, "beta_y", 9.0, weight=1.0e6, tol=1.0e-10)
    problem.target_periodic_twiss("beta_x", start, end, weight=1.0e6, tol=1.0e-10)
    problem.target_periodic_twiss("alpha_y", start, end, weight=1.0e6, tol=1.0e-10)

    result = problem.solve(solver="ls_trf", max_iter=300, tol=1.0e-12)
    assert result.success
    assert result.merit < 1.0e-12
    assert not np.isclose(problem.twiss0.alpha_x, tw0.alpha_x)
    assert not np.isclose(problem.twiss0.alpha_y, tw0.alpha_y)

    _merit, reports, _objectives, state = problem.evaluate()
    assert all(report.met for report in reports)
    assert np.isclose(state.twiss_at(end).alpha_x, 0.0, atol=1.0e-8)
    assert np.isclose(state.twiss_at(end).beta_y, 9.0, atol=1.0e-8)
    assert np.isclose(state.twiss_at(end).beta_x, state.twiss_at(start).beta_x, atol=1.0e-8)
    assert np.isclose(state.twiss_at(end).alpha_y, state.twiss_at(start).alpha_y, atol=1.0e-8)


def test_rmatrix_entry_target():
    lat, start, dvar, end = _simple_drift_lattice(1.1)
    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(dvar, "l", limits=(0.0, None))
    problem.target_rmatrix(start, end, i=0, j=1, value=3.5, weight=1.0e6)

    result = problem.solve(solver="ls_trf", max_iter=80, tol=1.0e-10)
    assert result.success
    _b, r_seg, _t = lat.transfer_maps(energy=1.0, start=start, stop=end)
    assert np.isclose(r_seg[0, 1], 3.5, atol=1.0e-4)


def test_rmatrix_entry_target_between_internal_elements():
    start = Marker(eid="START")
    d0 = Drift(l=0.4, eid="D0")
    m1 = Marker(eid="M1")
    dvar = Drift(l=1.1, eid="DVAR")
    m2 = Marker(eid="M2")
    d2 = Drift(l=0.7, eid="D2")
    end = Marker(eid="END")
    lat = SimpleMagneticLattice((start, d0, m1, dvar, m2, d2, end))

    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(dvar, "l", limits=(0.0, 5.0), name="DVAR_l")
    problem.target_rmatrix(m1, m2, i=0, j=1, value=2.7, weight=1.0e6)

    result = problem.solve(solver="ls_trf", max_iter=100, tol=1.0e-10)
    assert result.success
    assert np.isclose(dvar.l, 2.7, atol=1.0e-4)

    _b, r_seg, _t = lat.transfer_maps(energy=1.0, start=m1, stop=m2)
    assert np.isclose(r_seg[0, 1], 2.7, atol=1.0e-4)


def test_rmatrix_block_target():
    lat, start, dvar, end = _simple_drift_lattice(0.8)
    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(dvar, "l", limits=(0.0, None))
    problem.target_rmatrix_block(
        start,
        end,
        target_matrix=np.array([[1.0, 2.4], [0.0, 1.0]]),
        rows=[0, 1],
        cols=[0, 1],
        weight=1.0e6,
    )

    result = problem.solve(solver="ls_trf", max_iter=100, tol=1.0e-10)
    assert result.success
    assert np.isclose(dvar.l, 2.4, atol=1.0e-4)


def test_power_supply_linked_variable():
    start = Marker(eid="S")
    q1 = Quadrupole(l=0.2, k1=0.7, eid="Q1")
    q2 = Quadrupole(l=0.2, k1=-0.7, eid="Q2")
    end = Marker(eid="E")
    lat = SimpleMagneticLattice((start, q1, q2, end))

    problem = MatchProblem(lat, _twiss_seed())
    ps = problem.vary_linked_elements([q1, q2], scales=[1.0, -1.0], quantity="k1", name="PS_Q12")

    ps.set(2.0)
    assert np.isclose(q1.k1, 2.0, atol=1.0e-12)
    assert np.isclose(q2.k1, -2.0, atol=1.0e-12)

    ps.set(-1.3)
    assert np.isclose(q1.k1, -1.3, atol=1.0e-12)
    assert np.isclose(q2.k1, 1.3, atol=1.0e-12)


def test_report_generation():
    lat, _start, dvar, end = _simple_drift_lattice(1.0)
    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(dvar, "l", limits=(0.0, None))
    problem.target_twiss(end, "s", value=1.6, weight=1.0e6)

    result = problem.solve(solver="ls_trf", max_iter=80, tol=1.0e-10)
    assert result.success
    assert len(result.target_reports) == 1
    assert "actual" in result.target_reports[0].details


def test_generic_function_objective_minimization():
    lat, _start, dvar, _end = _simple_drift_lattice(1.3)
    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(dvar, "l", limits=(0.0, 5.0))
    problem.minimize_function(lambda state: state.twiss_end.s, name="min_total_s", weight=1.0e4)

    result = problem.solve(solver="ls_trf", max_iter=80, tol=1.0e-10)
    assert result.success
    assert dvar.l < 1.0e-4
    assert len(result.objective_reports) == 1
    assert result.objective_reports[0].name == "min_total_s"


def test_i5_objective_is_regular_objective_term():
    lat, _start, dvar, end = _simple_drift_lattice(1.0)
    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(dvar, "l", limits=(0.0, 5.0))
    problem.target_twiss(end, "s", value=1.2, weight=1.0e6)
    problem.minimize_i5_integral(weight=1.0)

    merit, _targets, objectives, _state = problem.evaluate()
    assert np.isfinite(merit)
    assert len(objectives) == 1
    assert objectives[0].details["type"] == "i5"


def test_custom_objective_class_can_be_added():
    lat, _start, dvar, _end = _simple_drift_lattice(0.5)

    class EndSObjective(Objective):
        def __init__(self, target_s: float, **kwargs):
            super().__init__(**kwargs)
            self.target_s = float(target_s)

        def residuals(self, state):
            return np.array([state.twiss_end.s - self.target_s], dtype=float)

    problem = MatchProblem(lat, _twiss_seed())
    problem.vary_element(dvar, "l", limits=(0.0, None), name="L_DVAR")
    problem.add_objective(EndSObjective(target_s=2.0, name="end_s_obj", weight=1.0e6))

    result = problem.solve(solver="ls_trf", max_iter=120, tol=1.0e-10)
    assert result.success
    assert np.isclose(dvar.l, 2.0, atol=1.0e-4)
    assert any(rep.name == "end_s_obj" for rep in result.objective_reports)


def test_match_end_twiss_by_varying_initial_twiss():
    lat, _start, _dvar, end = _simple_drift_lattice(1.0)
    tw0 = _twiss_seed()
    tw0.beta_x = 10.0
    tw0.alpha_x = 0.0

    problem = MatchProblem(lat, tw0, periodic=False)
    problem.vary_twiss(quantity="beta_x", limits=(1.0, 50.0), name="bx0")
    problem.target_twiss(end, "alpha_x", value=-0.5, weight=1.0e6)

    result = problem.solve(solver="ls_trf", max_iter=120, tol=1.0e-12)
    assert result.success
    assert np.isclose(problem.twiss0.beta_x, 2.0, atol=1.0e-4)
    assert np.isclose(result.variables["bx0"], 2.0, atol=1.0e-4)


def _build_lat(t3):
    return MagneticLattice(t3.cell, stop=t3.otrc_2560_t3)


def _variant_4b_lenses(t3) -> List:
    return [
        t3.qe_2480_t3,
        t3.qe_2495_t3,
        t3.qe_2506_t3,
        t3.qe_2518_t3,
        t3.qh_2529_t3,
        t3.qh_2532_t3,
        t3.qh_2544_t3,
        t3.qm_2549_t3,
        t3.qm_2554_t3,
        t3.qm_2559_t3,
    ]


def _set_variant_4b(t3):
    for lens, val in zip(_variant_4b_lenses(t3), LENSES_VAR4B):
        lens.k1 = val


def _extract_state(lat, t3):
    tws = twiss(lat, t3.twiss0, attach2elem=True)
    end = tws[-1]
    max_bx = max(t.beta_x for t in tws)
    max_by = max(t.beta_y for t in tws)

    seq = lat.get_sequence_part(t3.vcchirp_2490_t3, t3.otrc_2560_t3)
    seg_lat = MagneticLattice(seq)
    _, rmat, _ = seg_lat.transfer_maps(energy=14.0)
    r23 = float(rmat[2, 3])
    return end, max_bx, max_by, r23


def _run_case1(t3) -> RunResult:
    lat = _build_lat(t3)
    _set_variant_4b(t3)

    problem = MatchProblem(lat, t3.twiss0, periodic=False)
    problem.vary_element(t3.qm_2554_t3, "k1", limits=(-1, 1))
    problem.target_twiss(t3.otrc_2560_t3, "Dx", 0.0, weight=1e6)
    problem.target_global("beta_x", 300.0, relation="<=", weight=1e6, tol=0.0, name="beta_x_max")
    problem.target_global("beta_y", 200.0, relation="<=", weight=1e6, name="beta_y_max")

    t0 = perf_counter()
    result = problem.solve(solver="ls_trf", max_iter=3000, tol=1.0e-9)
    dt = perf_counter() - t0

    end, max_bx, max_by, r23 = _extract_state(lat, t3)
    return RunResult(
        duration_s=dt,
        success=result.success,
        merit=result.merit,
        end_beta_x=float(end.beta_x),
        end_beta_y=float(end.beta_y),
        end_Dx=float(end.Dx),
        max_beta_x=float(max_bx),
        max_beta_y=float(max_by),
        r23=r23,
        variables=result.variables,
    )


def _run_case2(t3) -> RunResult:
    lat = _build_lat(t3)
    _set_variant_4b(t3)

    problem = MatchProblem(lat, t3.twiss0, periodic=False)
    problem.vary_element(t3.qe_2506_t3, "k1", limits=(-1, 1))
    problem.vary_element(t3.qe_2518_t3, "k1", limits=(-1, 1))
    problem.vary_element(t3.qh_2529_t3, "k1", limits=(-1, 1))
    problem.vary_element(t3.qh_2532_t3, "k1", limits=(-1, 1))
    problem.vary_element(t3.qh_2544_t3, "k1", limits=(-1, 1))

    problem.target_twiss(t3.otrc_2560_t3, "beta_x", 6.7, weight=1e6)
    problem.target_twiss(t3.otrc_2560_t3, "beta_y", 50.0, relation="<=", weight=1e6)
    problem.target_rmatrix(
        t3.vcchirp_2490_t3,
        t3.otrc_2560_t3,
        i=2,
        j=3,
        value=-36.0,
        relation="==",
        weight=1e6,
        name="R23_target",
    )
    problem.target_global("beta_x", 150.0, relation="<=", weight=1e6, tol=0.0, name="beta_x_max")
    problem.target_global("beta_y", 150.0, relation="<=", weight=1e6, name="beta_y_max")

    t0 = perf_counter()
    result = problem.solve(solver="ls_trf", max_iter=3000, tol=1.0e-9)
    dt = perf_counter() - t0

    end, max_bx, max_by, r23 = _extract_state(lat, t3)
    return RunResult(
        duration_s=dt,
        success=result.success,
        merit=result.merit,
        end_beta_x=float(end.beta_x),
        end_beta_y=float(end.beta_y),
        end_Dx=float(end.Dx),
        max_beta_x=float(max_bx),
        max_beta_y=float(max_by),
        r23=r23,
        variables=result.variables,
    )


def test_match_optics_case1_regression(t3_module):
    out = _run_case1(t3_module)
    assert out.success
    assert abs(out.end_Dx) < 1.0e-6
    assert out.max_beta_x <= 300.0 + 1.0e-9
    assert out.max_beta_y <= 200.0 + 1.0e-9
    assert -1.0 <= out.variables["QM.2554.T3.k1"] <= 1.0


def test_match_optics_case2_regression(t3_module):
    out = _run_case2(t3_module)
    assert out.success
    assert np.isclose(out.end_beta_x, 6.7, atol=1.0e-3)
    assert out.end_beta_y <= 50.0 + 1.0e-6
    assert np.isclose(out.r23, -36.0, atol=1.0e-3)
    assert out.max_beta_x <= 150.0 + 1.0e-9
    assert out.max_beta_y <= 150.0 + 1.0e-9


def _bench_env_flag(primary: str, legacy: str) -> bool:
    return os.environ.get(primary, os.environ.get(legacy, "0")) == "1"


def test_match_optics_benchmark(t3_module):
    case1 = _run_case1(t3_module)
    case2 = _run_case2(t3_module)

    bench = {
        "case1_sec": case1.duration_s,
        "case2_sec": case2.duration_s,
        "total_sec": case1.duration_s + case2.duration_s,
    }
    print(
        f"\nbenchmark case1={bench['case1_sec']:.6f}s "
        f"case2={bench['case2_sec']:.6f}s total={bench['total_sec']:.6f}s"
    )

    if _bench_env_flag("MATCHER_BENCH_UPDATE", "NEW_MATCH_BENCH_UPDATE"):
        BENCH_BASELINE_FILE.write_text(json.dumps(bench, indent=2))

    enforce = _bench_env_flag("MATCHER_BENCH_ENFORCE", "NEW_MATCH_BENCH_ENFORCE")
    if enforce and BENCH_BASELINE_FILE.exists():
        baseline = json.loads(BENCH_BASELINE_FILE.read_text())
        max_slowdown = float(os.environ.get("MATCHER_BENCH_MAX_SLOWDOWN", "1.25"))
        assert bench["case1_sec"] <= baseline["case1_sec"] * max_slowdown
        assert bench["case2_sec"] <= baseline["case2_sec"] * max_slowdown
