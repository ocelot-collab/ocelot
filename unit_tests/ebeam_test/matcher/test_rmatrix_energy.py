"""Segment matrices must use the reference energy at the segment entrance."""

import numpy as np
import pytest

from ocelot.common.globals import m_e_GeV
from ocelot.cpbd.beam import Twiss
from ocelot.cpbd.elements import Cavity, Drift, Marker
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.matcher import MatchProblem


@pytest.mark.parametrize("cavity_occurrences", [1, 2])
def test_rmatrix_after_acceleration_uses_local_energy(cavity_occurrences):
    cavity = Cavity(l=0.5, v=0.1, phi=0.0, freq=1.3e9)
    start, end = Marker(), Marker()
    drift = Drift(l=1.0)
    # Repeated upstream cavities are valid; only the segment boundaries
    # need to identify unique lattice occurrences.
    lattice = MagneticLattice([cavity] * cavity_occurrences + [start, drift, end])
    seed = Twiss(beta_x=10.0, beta_y=10.0, E=0.1)
    segment_energy = seed.E + cavity_occurrences * cavity.v
    expected_r56 = -drift.l / ((segment_energy / m_e_GeV) ** 2 - 1.0)

    problem = MatchProblem(lattice, seed)
    problem.target_rmatrix(start, end, i=4, j=5, value=expected_r56, tol=1e-14)
    merit, reports, _objectives, state = problem.evaluate()

    assert not state.failed, state.failure_reason
    assert state.twiss_at(start).E == pytest.approx(segment_energy)
    np.testing.assert_allclose(state.r_matrix(start, end)[4, 5], expected_r56, rtol=1e-12)
    assert reports[0].met
    assert merit == 0.0


@pytest.mark.parametrize("upstream_cavity", [False, True])
def test_rmatrix_starting_with_cavity_uses_its_entrance_energy(upstream_cavity):
    prefix = [Cavity(l=0.5, v=0.1, phi=0.0, freq=1.3e9)] if upstream_cavity else []
    start = Cavity(l=0.5, v=0.05, phi=20.0, freq=1.3e9)
    end = Marker()
    lattice = MagneticLattice(prefix + [start, Drift(l=1.0), end])
    seed = Twiss(beta_x=10.0, beta_y=10.0, E=0.1)
    segment_energy = 0.2 if upstream_cavity else 0.1
    expected = lattice.transfer_maps(energy=segment_energy, start=start, stop=end)[1]

    state = MatchProblem(lattice, seed).evaluate()[3]

    assert not state.failed, state.failure_reason
    assert state.twiss_at(start).E > segment_energy
    np.testing.assert_allclose(state.r_matrix(start, end), expected, rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(
        state.r_matrix(None, None), lattice.transfer_maps(energy=seed.E)[1],
        rtol=1e-12, atol=1e-15,
    )


def test_rmatrix_matching_responds_to_upstream_cavity_voltage():
    cavity = Cavity(l=0.5, v=0.03, phi=0.0, freq=1.3e9)
    start, end = Marker(), Marker()
    drift = Drift(l=1.0)
    lattice = MagneticLattice([cavity, start, drift, end])
    seed = Twiss(beta_x=10.0, beta_y=10.0, E=0.1)
    target_energy = 0.2
    target_r56 = -drift.l / ((target_energy / m_e_GeV) ** 2 - 1.0)
    problem = MatchProblem(lattice, seed)
    problem.vary_element(cavity, "v", limits=(0.0, 0.2))
    problem.target_rmatrix(start, end, i=4, j=5, value=target_r56, weight=1e6, tol=1e-14)

    result = problem.solve(max_iter=80, tol=1e-10)

    assert result.success, result.message
    assert cavity.v == pytest.approx(target_energy - seed.E, abs=1e-8)
    assert result.target_reports[0].met
    actual = lattice.transfer_maps(energy=seed.E + cavity.v, start=start, stop=end)[1]
    np.testing.assert_allclose(actual[4, 5], target_r56, rtol=1e-7)
