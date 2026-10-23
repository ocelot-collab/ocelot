"""Physics checks independent of the nominal element transfer matrices."""
import copy

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from ocelot.common.globals import m_e_GeV, speed_of_light
from ocelot.cpbd.beam import ParticleArray
from ocelot.cpbd.elements import Drift, Quadrupole
from ocelot.cpbd.elements.undulator_atom import und_field
from ocelot.cpbd.high_order import rk_field, rk_field_to_ocelot, rk_track_in_field
from ocelot.cpbd.transformations.exact_drift import ExactDriftTM
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaGlobalTM, RungeKuttaOcelotTM


def zero_field(x, y, z):
    return 0. * x, 0. * x, 0. * x


@pytest.mark.parametrize('energy', [0.001, 1.])
@pytest.mark.parametrize('tracker', [rk_field, rk_field_to_ocelot])
def test_rk_matches_exact_drift_with_energy_angles_and_arrival_offsets(energy, tracker):
    particles = np.array([[0.001, 0.02, -0.003, -0.01, 0.004, 0.02],
                          [0., 0., 0., 0., -0.003, -0.01]]).T
    expected = particles.copy()
    ExactDriftTM._drift(expected, energy, 0.7)
    tracker(particles, 0.25, 0.95, 41, energy, zero_field)
    np.testing.assert_allclose(particles, expected, atol=2e-14, rtol=0)


@pytest.mark.parametrize('beta0', [0.8, 0.999])
@pytest.mark.parametrize('delta', [0., 0.04])
@pytest.mark.parametrize('sign', [-1., 1.])
def test_raw_rk_circle_uses_particle_rigidity_and_cartesian_z(beta0, delta, sign):
    energy = m_e_GeV / np.sqrt(1 - beta0**2)
    h = sign * 0.3
    brho = energy * 1e9 * beta0 / speed_of_light
    initial = np.zeros((6, 1))
    initial[5] = delta
    trajectory, path = rk_track_in_field(initial, 0.7, 81, energy,
                                          lambda x, y, z: (0.*x, brho*h + 0.*x, 0.*x),
                                          s_start=0.2, return_path=True)
    z = trajectory[4::9, 0] - 0.2
    curvature = h / np.sqrt(1 + 2*delta/beta0 + delta**2)
    cosine = np.sqrt(1 - (curvature*z)**2)
    expected_x = -curvature*z*z / (1 + cosine)
    np.testing.assert_allclose(trajectory[0::9, 0], expected_x, atol=3e-12, rtol=0)
    np.testing.assert_allclose(trajectory[1::9, 0], -curvature*z/cosine, atol=3e-14, rtol=0)
    np.testing.assert_allclose(path[:, 0], np.arcsin(curvature*z)/curvature, atol=3e-12, rtol=0)


@pytest.mark.parametrize('npoints', [81, 321])
def test_offset_field_full_map_matches_slices_in_all_six_coordinates(npoints):
    quad = Quadrupole(l=1., k1=1., dx=-0.02, dy=-0.03,
                      tilt=0.2, tm=RungeKuttaOcelotTM, npoints=npoints)
    initial = ParticleArray(n=3)
    initial.E = 0.001
    initial.rparticles[:] = np.array([[0., 0., 0., 0., 0., 0.],
                                     [0.001, 0.002, -0.001, -0.003, 0.004, 0.02],
                                     [-0.002, 0.001, 0.002, 0.001, -0.003, -0.01]]).T
    full, sliced = copy.deepcopy(initial), copy.deepcopy(initial)
    for tm in quad.tms:
        tm.apply(full)
    boundaries = [0., 0.13, 0.51, 0.76, 1.]
    for start, stop in zip(boundaries[:-1], boundaries[1:]):
        for tm in quad.get_section_tms(stop-start, start_l=start):
            tm.apply(sliced)
    np.testing.assert_allclose(full.rparticles, sliced.rparticles, atol=2e-10, rtol=0)
    np.testing.assert_allclose(full.rparticles[:, 0], 0, atol=2e-10)
    # No cached reference depends on which occurrence was tracked before.
    repeated = copy.deepcopy(initial)
    for tm in quad.tms:
        tm.apply(repeated)
    np.testing.assert_array_equal(full.rparticles, repeated.rparticles)


def test_curved_plane_conversion_matches_independent_lorentz_solver():
    """Integrate in particle arc length, with a plane event (not production Z RK4)."""
    energy = 0.001
    beta0 = np.sqrt(1 - (m_e_GeV/energy)**2)
    brho = energy*1e9*beta0/speed_of_light
    length = 0.8

    def field(x, y, z):
        # Offset quadrupole plus a longitudinal field exercises frame transport.
        return brho*(y + .03), brho*(x + .02), brho*(.06 + .02*z)

    def reference_rhs(s, w):
        r, tangent, ex = w.reshape(3, 3)
        dt = np.cross(tangent, np.asarray(field(*r))/brho)
        return np.r_[tangent, dt, -tangent*np.dot(ex, dt)]

    ref = solve_ivp(reference_rhs, [0, length], [0,0,0, 0,0,1, 1,0,0],
                    method='DOP853', rtol=2e-13, atol=2e-15).y[:, -1]
    origin, tangent, ex = ref.reshape(3, 3)
    frame = np.column_stack((ex, np.cross(tangent, ex), tangent))
    q = np.array([.003, .01, -.002, -.015, .007, .04])
    magnitude = np.sqrt(1 + 2*q[5]/beta0 + q[5]**2)
    beta = magnitude/(1/beta0 + q[5])
    momentum = np.array([q[1], q[3], np.sqrt(magnitude**2-q[1]**2-q[3]**2)])

    def rhs(path, state):
        r, direction = state[:3], state[3:6]
        return np.r_[direction, np.cross(direction, np.asarray(field(*r))/(brho*magnitude)), 1/beta]

    def exit_plane(path, state):
        return np.dot(state[:3]-origin, tangent)

    exit_plane.terminal = True
    solution = solve_ivp(rhs, [0, 2*length], np.r_[[q[0],q[2],0], momentum/magnitude, 0.],
                         events=exit_plane, method='DOP853', rtol=2e-13, atol=2e-15)
    assert len(solution.t_events[0]) == 1
    final = solution.y[:, -1]
    position = frame.T @ (final[:3]-origin)
    momentum = frame.T @ (final[3:6]*magnitude)
    expected = [position[0], momentum[0], position[1], momentum[1],
                q[4]+final[6]-length/beta0, q[5]]
    tracked = q[:, None].copy()
    rk_field_to_ocelot(tracked, 0, length, 321, energy, field)
    np.testing.assert_allclose(tracked[:, 0], expected, atol=2e-12, rtol=0)
    sliced = q[:, None].copy()
    for start, stop in zip([0., .23, .51], [.23, .51, length]):
        rk_field_to_ocelot(sliced, start, stop, 321, energy, field, reference_length=length)
    np.testing.assert_allclose(sliced[:, 0], expected, atol=2e-12, rtol=0)


@pytest.mark.parametrize('npoints', [1, 2.5, True])
def test_invalid_grid_is_rejected(npoints):
    with pytest.raises(ValueError, match='npoints'):
        rk_track_in_field(np.zeros((6, 1)), 1, npoints, 1, zero_field)


def test_nonforward_particle_is_rejected():
    particle = np.zeros((6, 1))
    particle[1] = 1.1
    with pytest.raises(ValueError, match='forward'):
        rk_field(particle, 0, 1, 20, 1, zero_field)


@pytest.mark.parametrize('long_dynamics', [False, True])
@pytest.mark.parametrize('start', [0., .35])
def test_global_endpoint_matches_stored_trajectory(start, long_dynamics):
    energy = .001
    beta0 = np.sqrt(1 - (m_e_GeV/energy)**2)
    brho = energy*1e9*beta0/speed_of_light
    particles = np.array([[.001, .02, -.003, -.01, .004, .02],
                          [-.002, -.01, .001, .02, -.003, -.01]]).T
    expected = particles.copy()
    magnitude2 = 1 + 2*particles[5]/beta0 + particles[5]**2
    beta = np.sqrt(magnitude2)/(1/beta0 + particles[5])
    slopes = particles.copy()
    pz = np.sqrt(magnitude2 - particles[1]**2 - particles[3]**2)
    slopes[[1, 3]] /= pz

    def field(x, y, z):
        return brho*(y + .03), brho*(x + .02*np.cos(3*z)), brho*(.06 + 0*x)

    trajectory, paths = rk_track_in_field(slopes, start+.7, 81, energy, field,
                                          s_start=start, return_path=True)
    exit_slopes = trajectory[-9:-3]
    exit_pz = np.sqrt(magnitude2/(1 + exit_slopes[1]**2 + exit_slopes[3]**2))
    expected[:4] = exit_slopes[:4]
    expected[[1, 3]] *= exit_pz
    if long_dynamics:
        expected[4] += paths[-1]/beta - .7/beta0
    result = rk_field(particles, start, start+.7, 81, energy, field, long_dynamics)
    assert result is particles
    np.testing.assert_allclose(particles, expected, atol=2e-14, rtol=0)


@pytest.mark.parametrize('npoints', [1, 2.5, True])
def test_global_invalid_grid_is_rejected(npoints):
    with pytest.raises(ValueError, match='npoints'):
        rk_field(np.zeros((6, 1)), 0., 1., npoints, 1., zero_field)


def test_global_reversed_interval_is_rejected():
    with pytest.raises(ValueError, match='s_stop'):
        rk_field(np.zeros((6, 1)), .5, .2, 20, 1., zero_field)


@pytest.mark.parametrize('count', [0, 2])
def test_global_zero_length_preserves_particles(count):
    particles = np.full((6, count), .001)
    expected = particles.copy()
    rk_field(particles, .3, .3, 20, .001, zero_field)
    np.testing.assert_array_equal(particles, expected)


def test_global_memory_does_not_grow_with_step_count():
    """Storing every sample would add over 5 MB for this small bunch."""
    import tracemalloc

    particles = np.zeros((6, 128))
    peaks = []
    # Warm up before measuring allocations from the integration itself.
    rk_field(particles, 0., 1., 2, 1., zero_field)
    if tracemalloc.is_tracing():
        pytest.skip('requires its own allocation measurement')
    for npoints in [33, 513]:
        tracemalloc.start()
        try:
            rk_field(particles, 0., 1., npoints, 1., zero_field)
            peaks.append(tracemalloc.get_traced_memory()[1])
        finally:
            tracemalloc.stop()
    assert peaks[1] < 2*peaks[0] + 64_000


def test_sinusoidal_field_has_fourth_order_grid_convergence():
    energy, delta, curvature, wave_number, length = .001, .03, .2, 30., 1.1
    beta0 = np.sqrt(1-(m_e_GeV/energy)**2)
    brho = energy*1e9*beta0/speed_of_light
    momentum2 = 1+2*delta/beta0+delta**2
    initial = np.zeros((6,1))
    initial[5] = delta
    field = lambda x,y,z: (0*x, 0*x+brho*curvature*np.cos(wave_number*z), 0*x)
    errors = []
    for count in [101, 201]:
        trajectory = rk_track_in_field(initial, length, count, energy, field)
        z = trajectory[4::9, 0]
        px = -curvature/wave_number*np.sin(wave_number*z)
        slope = px/np.sqrt(momentum2-px**2)
        errors.append(np.max(abs(trajectory[1::9, 0]-slope)))
    assert errors[1] < errors[0]/12
    assert errors[1] < 3e-9


@pytest.mark.parametrize('delta', [-.01, 0., .01])
def test_global_undulator_arrival_time_matches_exact_planar_orbit(delta):
    """Global tau includes the extra flight time along the undulator orbit."""
    energy, period, periods, K = .6, .4, 9, 44.81
    length = period * periods
    beta0 = np.sqrt(1 - (m_e_GeV / energy)**2)
    magnitude = np.sqrt(1 + 2*delta/beta0 + delta**2)
    beta = magnitude / (1/beta0 + delta)
    amplitude = K * m_e_GeV / (energy * beta0)

    # In the symmetry plane PX=-amplitude*sin(2*pi*Z/period).
    # Integrate its exact path over one period, independent of the RK solver.
    path = periods * quad(
        lambda z: magnitude / np.sqrt(
            magnitude**2 - (amplitude*np.sin(2*np.pi*z/period))**2),
        0., period, epsabs=1e-14, epsrel=1e-13)[0]
    particles = np.zeros((6, 1))
    particles[4, 0] = .003
    particles[5, 0] = delta
    field = lambda x, y, z: und_field(x, y, z, period, K)
    rk_field(particles, 0., length, 6000, energy, field)
    expected = [0., 0., 0., 0., .003 + path/beta - length/beta0, delta]
    np.testing.assert_allclose(particles[:, 0], expected, atol=2e-12, rtol=0)


def test_navigator_slices_with_repeated_element_occurrences():
    from ocelot.cpbd.magnetic_lattice import MagneticLattice
    from ocelot.cpbd.navi import Navigator
    from ocelot.cpbd.track import tracking_step

    quad = Quadrupole(l=.4, k1=.8, dx=-.02, dy=-.03,
                      tm=RungeKuttaOcelotTM, npoints=201)
    lattice = MagneticLattice([quad, Drift(l=.1), quad], method={Quadrupole: RungeKuttaOcelotTM})
    full = ParticleArray(n=1)
    full.E = .001
    full.rparticles[:,0] = [.001,.002,-.001,-.003,.004,.02]
    sliced = copy.deepcopy(full)
    for element in lattice.sequence:
        for tm in element.tms:
            tm.apply(full)
    navigator = Navigator(lattice)
    for ds in np.diff(np.linspace(0,lattice.totalLen,32)):
        tracking_step(lattice,sliced,dz=ds,navi=navigator)
    np.testing.assert_allclose(sliced.rparticles, full.rparticles, atol=2e-12, rtol=0)


@pytest.mark.parametrize('method', [RungeKuttaGlobalTM, RungeKuttaOcelotTM])
def test_empty_beam_can_advance_after_particle_losses(method):
    particles = ParticleArray(n=0)
    particles.E = 1.
    drift = Drift(l=.5, tm=method)
    drift.tms[0].apply(particles)
    assert particles.rparticles.shape == (6, 0)
    assert particles.s == pytest.approx(.5)
