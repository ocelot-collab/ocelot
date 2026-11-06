"""CSR Cartesian geometry, local frames and physical arc-length mapping."""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from ocelot.common.globals import m_e_GeV, speed_of_light
from ocelot.cpbd.beam import ParticleArray
from ocelot.cpbd.csr import CSR, CSRConfigurationError, SaferCSR
from ocelot.cpbd.elements import Bend, Drift, Marker, Undulator
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.navi import Navigator


def prepared(elements, rk=False, step=.001, cls=CSR, **kwargs):
    lattice = MagneticLattice([*elements, Marker(eid='csr_end')])
    process = cls(rk_traj=rk, energy=.13, traj_step=step, **kwargs)
    Navigator(lattice).add_physics_proc(process, lattice.sequence[0], lattice.sequence[-1])
    return lattice, process


@pytest.mark.parametrize('rk', [False, True])
@pytest.mark.parametrize('angles,tilts', [((.1,), (0.,)), ((.1, -.2), (0., 0.)),
                                        ((.1, .1, -.05), (0., np.pi/2, -.4)),
                                        ((1.8, .2), (.3, 1.2))])
def test_bends_follow_survey_in_transported_frames(rk, angles, tilts):
    elements = [Bend(l=.5, angle=a, tilt=t, eid=f'b{i}') for i, (a, t) in enumerate(zip(angles, tilts))]
    lattice, process = prepared(elements, rk)
    _, survey = lattice.survey()
    traj = process.csr_traj
    np.testing.assert_allclose(traj[1:4, -1], survey[-1]['r_end'], atol=2e-12, rtol=0)
    np.testing.assert_allclose(traj[4:, -1], survey[-1]['W'][:, 2], atol=2e-12, rtol=0)
    np.testing.assert_allclose(np.linalg.norm(traj[4:], axis=0), 1., atol=2e-14, rtol=0)
    assert np.all(np.diff(traj[0]) > 0)
    assert np.all(np.linalg.norm(np.diff(traj[1:4], axis=1), axis=0) > 0)


def test_rk_opposite_bends_have_finite_exit_kernel_and_kick():
    _, process = prepared([Bend(l=.5, angle=.1, eid='b1'), Bend(l=.5, angle=-.2, eid='b2')], True)
    p = ParticleArray(1000)
    p.E = .13
    p.q_array[:] = 1e-12
    p.rparticles[4] = np.random.default_rng(1).normal(0, 1e-4, p.n)
    process.z0 = 1.
    process.apply(p, .01)
    assert np.all(np.isfinite(p.p()))
    assert np.linalg.norm(p.p()) > 0


def test_mixed_plane_rk_and_analytic_kernels_agree():
    elements = [Bend(l=.3, angle=.12, eid='bx'), Bend(l=.4, angle=-.08, tilt=1.2, eid='by')]
    _, exact = prepared(elements)
    _, rk = prepared(elements, True)
    for index in (350, 500, 700):
        expected = exact.CSR_K1(index, exact.csr_traj, [100, 2e-5], gamma=.13/m_e_GeV)
        actual = rk.CSR_K1(index, rk.csr_traj, [100, 2e-5], gamma=.13/m_e_GeV)
        # The finite-energy formula subtracts nearly equal terms; assess the
        # error relative to the kernel norm, including its zero crossings.
        assert np.linalg.norm(actual - expected) / np.linalg.norm(expected) < 1e-8


def test_stop_before_undulator_does_not_require_rk_energy():
    bend = Bend(l=.5, angle=.1, eid='b')
    und = Undulator(lperiod=.04, nperiods=2, Kx=3, eid='u')
    lattice = MagneticLattice([bend, und])
    process = CSR()
    Navigator(lattice).add_physics_proc(process, bend, und)
    assert not process.rk_traj
    assert process.csr_traj[0, -1] == pytest.approx(bend.l)


@pytest.mark.parametrize('energy', [None, 0., np.nan, m_e_GeV])
def test_rk_requires_valid_reference_energy(energy):
    lattice = MagneticLattice([Bend(l=.1, angle=.1, eid='b'), Marker(eid='end')])
    with pytest.raises(CSRConfigurationError, match='energy'):
        Navigator(lattice).add_physics_proc(CSR(rk_traj=True, energy=energy), *lattice.sequence)


@pytest.fixture(scope='module')
def strong_undulator():
    return prepared([Undulator(lperiod=.4, nperiods=9, Kx=44.81, eid='wiggler'),
                     Drift(l=1., eid='out')], cls=SaferCSR, step=.0005)[1]


def test_wiggler_keeps_physical_arc_length_and_maps_navigator_exit(strong_undulator):
    from scipy.special import ellipk

    process = strong_undulator
    arc_exit = process._lattice_to_trajectory_s(3.6)
    index = np.argmin(abs(process.csr_traj[0] - arc_exit))
    assert arc_exit > 3.62
    # For B_y=-B0*cos(kZ), n_x=K/(gamma*beta)*sin(kZ). Integrating
    # ds/dZ=1/sqrt(1-n_x**2) over whole periods gives an elliptic integral.
    amplitude = 44.81 / np.sqrt((.13/m_e_GeV)**2 - 1.)
    expected_arc = 3.6 * 2/np.pi * ellipk(amplitude**2)
    assert arc_exit == pytest.approx(expected_arc, abs=2e-11)
    assert process.csr_traj[3, index] == pytest.approx(3.6, abs=1e-12)
    assert process._lattice_to_trajectory_s(4.6) == pytest.approx(arc_exit + 1.)
    np.testing.assert_allclose(np.linalg.norm(process.csr_traj[4:], axis=0), 1., atol=2e-14, rtol=0)


def test_wiggler_kick_selects_correct_position_and_integrates_physical_length(strong_undulator, monkeypatch):
    from ocelot.cpbd import csr as csr_module

    process = strong_undulator
    indices = []
    def kernel(index, traj, ndw, gamma):
        indices.append(index)
        return np.ones(ndw[0] + 1)
    monkeypatch.setattr(process, 'CSR_K1', kernel)
    monkeypatch.setattr(process.bin_smoth, 'Q2EQUI', lambda *args: (-.01, .01, 100, np.ones(100)))
    monkeypatch.setattr(csr_module, 'csr_convolution', lambda a, b: np.ones(len(a) + len(b) - 1))
    p = ParticleArray(20)
    p.E = .13
    p.q_array[:] = 1e-12
    p.tau()[:] = np.linspace(-1e-4, 1e-4, p.n)
    process.z0 = 3.6
    process.apply(p, .1)
    assert process.csr_traj[3, indices[-1]] == pytest.approx(3.6, abs=1e-12)
    path_length = process._lattice_to_trajectory_s(3.6) - process._lattice_to_trajectory_s(3.5)
    expected = path_length / .0002 * 1e-9 / np.sqrt(p.E**2 - m_e_GeV**2)
    np.testing.assert_allclose(p.p(), expected, rtol=2e-13)


def test_safer_csr_uses_same_wiggler_position(strong_undulator, monkeypatch):
    process = strong_undulator
    radii = []
    monkeypatch.setattr(process, '_derbenev_criterion', lambda x, z, r: radii.append(r) or 0.)
    process.z0 = 3.61  # Already in the drift, despite extra wiggler arc length.
    process.query_derbenev_criterion(1e-4, 1e-4, 1e-4)
    assert radii == []


@pytest.mark.parametrize('end_poles', [False, True])
def test_undulator_phase_and_vertical_strength_affect_orbit(end_poles):
    _, horizontal = prepared([Undulator(lperiod=.04, nperiods=3, Kx=3, eid='h')], end_poles=end_poles)
    _, shifted = prepared([Undulator(lperiod=.04, nperiods=3, Kx=3, phase=.4, eid='p')], end_poles=end_poles)
    _, vertical = prepared([Undulator(lperiod=.04, nperiods=3, Ky=3, eid='v')], end_poles=end_poles)
    # Compare at a fixed physical Z rather than assuming equal arc-grid sizes.
    def position_at(process, row):
        return np.interp(.03, process.csr_traj[3], process.csr_traj[row])
    assert abs(position_at(shifted, 1) - position_at(horizontal, 1)) > 1e-6
    assert np.max(abs(vertical.csr_traj[2])) > 1e-6
    np.testing.assert_allclose(position_at(vertical, 2), -position_at(horizontal, 1), atol=1e-12)


def test_custom_3d_wiggler_matches_independent_lorentz_solution_after_bend():
    period, length = .08, .24
    def field(x, y, z):
        return .12*np.sin(2*np.pi*z/period), .2*np.cos(2*np.pi*z/period), .03
    bend = Bend(l=.2, angle=.25, tilt=.7, eid='entrance')
    und = Undulator(lperiod=period, nperiods=3, mag_field=field, eid='custom')
    _, process = prepared([bend, und], True, step=.0005)
    beta = np.sqrt(1 - (m_e_GeV/.13)**2)
    scale = speed_of_light / (.13e9 * beta)
    def rhs(arc, state):
        return np.r_[state[3:], scale*np.cross(state[3:], field(*state[:3]))]
    def exit_plane(arc, state):
        return state[2] - length
    exit_plane.terminal = True
    solution = solve_ivp(rhs, (0., length*2), [0., 0., 0., 0., 0., 1.],
                         events=exit_plane, dense_output=True, method='DOP853',
                         rtol=2e-12, atol=2e-14, max_step=.001)
    assert solution.success and len(solution.t_events[0]) == 1
    _, survey = MagneticLattice([bend]).survey()
    frame, origin = survey[-1]['W'], survey[-1]['r_end']
    segment = process.csr_traj[:, process.csr_traj[0] >= bend.l]
    expected = solution.sol(segment[0] - bend.l)
    np.testing.assert_allclose(segment[1:4], origin[:, None] + frame @ expected[:3], atol=2e-10, rtol=0)
    np.testing.assert_allclose(segment[4:], frame @ expected[3:], atol=5e-9, rtol=0)
    assert segment[0, -1] - bend.l == pytest.approx(solution.t_events[0][0], abs=2e-11)
    reference_traj = process.csr_traj.copy()
    reference_traj[1:4, -segment.shape[1]:] = origin[:, None] + frame @ expected[:3]
    reference_traj[4:, -segment.shape[1]:] = frame @ expected[3:]
    index = process.csr_traj.shape[1] - 1
    reference_kernel = process.CSR_K1(index, reference_traj, [100, 2e-5], gamma=.13/m_e_GeV)
    kernel = process.CSR_K1(index, process.csr_traj, [100, 2e-5], gamma=.13/m_e_GeV)
    assert np.linalg.norm(kernel - reference_kernel) / np.linalg.norm(reference_kernel) < 1e-5


def test_undulator_map_uses_all_components_and_actual_length(tmp_path):
    path = tmp_path / 'field.txt'
    # A line map with both transverse fields, starting at a nonzero map Z.
    np.savetxt(path, np.column_stack((np.linspace(10., 110., 11), np.full(11, .03), np.full(11, .1))))
    lattice, mapped = prepared([Undulator(field_file=str(path), eid='map')])
    _, explicit = prepared([Undulator(lperiod=.1, nperiods=1,
                                     mag_field=lambda x, y, z: (.03, .1, 0.), eid='func')])
    assert lattice.totalLen == pytest.approx(.1)
    assert mapped.csr_traj[3, -1] == pytest.approx(.1, abs=1e-12)
    np.testing.assert_allclose(mapped.csr_traj, explicit.csr_traj, atol=2e-13, rtol=0)
