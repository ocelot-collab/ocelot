import importlib

import numpy as np
import pytest

from ocelot.cpbd.beam import ParticleArray
from ocelot.cpbd.elements import Drift, Undulator
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaTM
from ocelot.rad.radiation_py import BeamTraject, integ_beta2, x2xgaus
from ocelot.rad.screen import Screen
from ocelot.rad.spline_py import cspline_coef


radiation_module = importlib.import_module("ocelot.rad.radiation_py")


class CustomUndulator(Undulator):
    pass


def _integ_beta2_scalar_reference(x, y):
    """Original scalar accumulation used as a numerical reference."""
    a, b, c, d, _ = cspline_coef(x, y)
    integral = 0.0
    values = [0.0]

    for i, step in enumerate(np.diff(x)):
        integral += step * (
            d[i] * d[i]
            + step
            * (
                c[i] * d[i]
                + step
                * (
                    (c[i] * c[i] + 2.0 * b[i] * d[i]) / 3.0
                    + step
                    * (
                        (b[i] * c[i] + a[i] * d[i]) / 2.0
                        + step
                        * (
                            (b[i] * b[i] + 2.0 * a[i] * c[i]) / 5.0
                            + step
                            * (
                                a[i] * b[i] / 3.0
                                + a[i] * a[i] * step / 7.0
                            )
                        )
                    )
                )
            )
        )
        values.append(integral)

    return np.asarray(values)


def test_x2xgaus_matches_three_point_gauss_nodes():
    coordinates = np.array([0.0, 2.0, 4.0, 6.0])
    sqrt_3_over_5 = np.sqrt(3.0 / 5.0)
    local_nodes = np.array(
        [
            1.0 - sqrt_3_over_5,
            1.0,
            1.0 + sqrt_3_over_5,
        ]
    )
    expected = np.concatenate(
        (
            coordinates[:1],
            (coordinates[:-1, np.newaxis] + local_nodes).reshape(-1),
            coordinates[-1:],
        )
    )

    np.testing.assert_array_equal(x2xgaus(coordinates), expected)


def test_integ_beta2_matches_scalar_reference():
    coordinates = np.linspace(0.0, 2.0, 17)
    beta = 0.03 * np.sin(1.7 * coordinates) + 0.01 * coordinates

    expected = _integ_beta2_scalar_reference(coordinates, beta)

    np.testing.assert_allclose(
        integ_beta2(coordinates, beta),
        expected,
        rtol=5.0e-16,
        atol=0.0,
    )


def test_integ_beta2_numba_and_numpy_fallback_are_equivalent(monkeypatch):
    coordinates = np.linspace(0.0, 3.6, 257)
    beta = 0.02 * np.sin(2.3 * coordinates) - 0.005 * coordinates
    expected = radiation_module._integ_beta2_numpy(coordinates, beta)

    if radiation_module._integ_beta2_compiled is not None:
        np.testing.assert_array_equal(
            radiation_module._integ_beta2_compiled(coordinates, beta),
            expected,
        )

    monkeypatch.setattr(radiation_module, "_integ_beta2_compiled", None)
    np.testing.assert_array_equal(
        radiation_module.integ_beta2(coordinates, beta),
        expected,
    )


@pytest.mark.skipif(
    radiation_module.gintegrator_over_spectrum is None,
    reason="parallel spectrum kernel requires Numba",
)
def test_parallel_spectrum_kernel_matches_general_kernel():
    n_motion = 4
    n_gauss_points = 3 * (n_motion - 1) + 2
    z = np.linspace(0.0, 3000.0, n_gauss_points)
    x = 0.02 * np.sin(z / 500.0)
    y = 0.01 * np.cos(z / 700.0)
    bx = 2.0e-5 * np.cos(z / 500.0)
    by = -1.0e-5 * np.sin(z / 700.0)
    field_x = 0.05 * np.sin(z / 300.0)
    field_y = 0.7 * np.cos(z / 400.0)
    beta_integral_x = np.linspace(0.0, 2.0e-6, n_gauss_points)
    beta_integral_y = np.linspace(0.0, 1.0e-6, n_gauss_points)

    energies = np.linspace(0.001, 0.003, 65)
    initial_phase = np.linspace(-0.2, 0.3, len(energies))
    general_outputs = [
        np.zeros((len(energies), 1)),
        np.zeros((len(energies), 1)),
        np.zeros((len(energies), 1)),
        np.zeros((len(energies), 1)),
        initial_phase.reshape(-1, 1).copy(),
    ]
    spectrum_outputs = [
        np.zeros(len(energies)),
        np.zeros(len(energies)),
        np.zeros(len(energies)),
        np.zeros(len(energies)),
        initial_phase.copy(),
    ]

    half_step = (z[-1] - z[0]) / 2.0 / (n_motion - 1)
    n_end = n_gauss_points - 2
    gamma = 1200.0
    distance = 1.0e6
    x_screen = 0.1
    y_screen = -0.2

    radiation_module.gintegrator_over_traj(
        n_motion,
        np.array([x_screen]),
        np.array([[y_screen]]),
        energies.reshape(-1, 1),
        n_end,
        gamma,
        half_step,
        distance,
        x,
        y,
        z,
        bx,
        by,
        beta_integral_x,
        beta_integral_y,
        field_x,
        field_y,
        *general_outputs,
    )
    radiation_module.gintegrator_over_spectrum(
        n_motion,
        x_screen,
        y_screen,
        energies,
        n_end,
        gamma,
        half_step,
        distance,
        x,
        y,
        z,
        bx,
        by,
        beta_integral_x,
        beta_integral_y,
        field_x,
        field_y,
        *spectrum_outputs,
    )

    for general, spectrum in zip(general_outputs, spectrum_outputs):
        np.testing.assert_allclose(
            spectrum,
            general.reshape(-1),
            rtol=2.0e-14,
            atol=2.0e-15,
        )


def test_beam_trajectory_accessors_concatenate_segments_in_order():
    first = np.zeros((18, 2))
    second = np.zeros((27, 2))

    for segment_index, segment in enumerate((first, second), start=1):
        steps = segment.shape[0] // 9
        for coordinate_index in range(6):
            values = (
                1000.0 * segment_index
                + 100.0 * coordinate_index
                + np.arange(steps)[:, np.newaxis] * 10.0
                + np.arange(2)
            )
            segment[coordinate_index::9, :] = values

    trajectory = BeamTraject([first, second])

    for accessor, coordinate_index in (
        (trajectory.x, 0),
        (trajectory.xp, 1),
        (trajectory.y, 2),
        (trajectory.yp, 3),
        (trajectory.z, 4),
        (trajectory.p, 5),
    ):
        expected = np.concatenate(
            (
                first[coordinate_index::9, 1],
                second[coordinate_index::9, 1],
            )
        )
        np.testing.assert_array_equal(accessor(1), expected)


def test_screen_rebuild_preserves_memory_screen_views():
    screen = Screen()
    screen.nx = 2
    screen.ny = 2
    screen.num_energy = 3
    screen.update()
    screen.arReEx[:] = np.linspace(1.0, 2.0, len(screen.arReEx))
    screen.arImEx[:] = np.linspace(0.1, 0.2, len(screen.arImEx))
    screen.arReEy[:] = np.linspace(2.0, 3.0, len(screen.arReEy))
    screen.arImEy[:] = np.linspace(0.2, 0.3, len(screen.arImEy))

    screen.rebuild_efields()

    for field in (
        screen.arReEx,
        screen.arImEx,
        screen.arReEy,
        screen.arImEy,
        screen.arPhase,
    ):
        assert np.shares_memory(field, screen.memory_screen)


def test_track4rad_beam_accepts_undulator_subclass_and_honors_npoints(
    monkeypatch,
):
    npoints = 37
    undulator = CustomUndulator(
        lperiod=0.02,
        nperiods=5,
        Kx=1.0,
        npoints=npoints,
    )
    lattice = MagneticLattice((undulator,))
    particles = ParticleArray(n=1)
    particles.E = 1.0
    calls = []

    def fake_rk_track_in_field(y0, length, point_count, energy, mag_field, *, return_path=False):
        calls.append(point_count)
        trajectory = np.zeros((point_count * 9, y0.shape[1]))
        trajectory[4::9, :] = np.linspace(0.0, length, point_count)[:, None]
        paths = np.linspace(0., length, point_count)[:, None] * np.sqrt(1 + y0[1]**2 + y0[3]**2)
        return (trajectory, paths) if return_path else trajectory

    monkeypatch.setattr(
        radiation_module,
        "rk_track_in_field",
        fake_rk_track_in_field,
    )

    trajectories, energies = radiation_module.track4rad_beam(
        particles,
        lattice,
        accuracy=5,
    )

    assert calls == [npoints]
    assert trajectories[0].shape == (npoints * 9, 1)
    assert energies == [particles.E]


def test_radiation_npoints_prefers_active_runge_kutta_override():
    undulator = Undulator(
        lperiod=0.02,
        nperiods=5,
        Kx=1.0,
        npoints=37,
    )
    undulator.set_tm(RungeKuttaTM, npoints=53)

    assert radiation_module._undulator_trajectory_points(
        undulator,
        accuracy=5,
    ) == 53


def test_track4rad_beam_includes_trailing_non_undulator_section(monkeypatch):
    undulator = Undulator(
        lperiod=0.02,
        nperiods=5,
        Kx=1.0,
        npoints=4,
    )
    drift = Drift(l=0.05)
    lattice = MagneticLattice((undulator, drift))
    particles = ParticleArray(n=1)
    particles.E = 1.0
    particles.px()[0] = 1.0e-3

    def fake_rk_track_in_field(y0, length, point_count, energy, mag_field, *, return_path=False):
        trajectory = np.zeros((point_count * 9, y0.shape[1]))
        trajectory[0::9, :] = y0[0]
        trajectory[1::9, :] = y0[1]
        trajectory[2::9, :] = y0[2]
        trajectory[3::9, :] = y0[3]
        trajectory[4::9, :] = np.linspace(0.0, length, point_count)[:, None]
        trajectory[5::9, :] = y0[5]
        paths = np.linspace(0., length, point_count)[:, None] * np.sqrt(1 + y0[1]**2 + y0[3]**2)
        return (trajectory, paths) if return_path else trajectory

    monkeypatch.setattr(
        radiation_module,
        "rk_track_in_field",
        fake_rk_track_in_field,
    )

    trajectories, energies = radiation_module.track4rad_beam(
        particles,
        lattice,
    )

    assert len(trajectories) == 2
    assert energies == [particles.E, particles.E]
    assert trajectories[-1][4::9, 0][-1] == pytest.approx(lattice.totalLen)
    assert particles.x()[0] == pytest.approx(
        particles.px()[0] * drift.l,
        rel=1.0e-6,
    )


@pytest.mark.parametrize(
    "npoints,exception",
    [
        (3, ValueError),
        (12.5, TypeError),
        (True, TypeError),
    ],
)
def test_track4rad_beam_rejects_invalid_undulator_npoints(npoints, exception):
    undulator = Undulator(
        lperiod=0.02,
        nperiods=5,
        Kx=1.0,
        npoints=npoints,
    )
    lattice = MagneticLattice((undulator,))
    particles = ParticleArray(n=1)
    particles.E = 1.0

    with pytest.raises(exception, match="npoints"):
        radiation_module.track4rad_beam(particles, lattice)


def test_track4rad_beam_rejects_lattice_without_trackable_elements():
    particles = ParticleArray(n=1)
    particles.E = 1.0
    lattice = MagneticLattice((Drift(l=0.0),))

    with pytest.raises(ValueError, match="no trackable nonzero-length"):
        radiation_module.track4rad_beam(particles, lattice)


@pytest.mark.parametrize('energy', [0.001, 1.])
def test_radiation_endpoint_matches_fixed_frame_rk_with_arrival_offsets(energy):
    import copy
    from ocelot.cpbd.high_order import rk_field

    field = lambda x, y, z: (0.*x, 0.0001 + 0.*x, 0.*x)
    undulator = Undulator(lperiod=.02, nperiods=5, Kx=0., npoints=81, mag_field=field)
    particles = ParticleArray(n=2)
    particles.E = energy
    particles.rparticles[:] = np.array([[.001, .01, -.002, -.003, .004, .02],
                                       [0., 0., 0., 0., -.002, -.01]]).T
    expected = particles.rparticles.copy()
    rk_field(expected, 0, undulator.l, 81, energy, field)
    trajectories, _ = radiation_module.track4rad_beam(particles, MagneticLattice([undulator]))
    np.testing.assert_allclose(particles.rparticles, expected, atol=2e-14, rtol=0)
    assert particles.s == pytest.approx(undulator.l)
    # Restoration must reproduce the saved SR endpoint exactly. Its slope
    # conversions can differ by roundoff from direct momentum tracking above.
    saved_endpoint = particles.rparticles.copy()
    saved = BeamTraject(trajectories, end_particle_array=particles)
    restored = copy.deepcopy(particles)
    restored.rparticles[:] = 0
    saved.p_array_end(restored)
    np.testing.assert_array_equal(restored.rparticles, saved_endpoint)
    # Repeated restoration does not accumulate path length again.
    saved.p_array_end(restored)
    np.testing.assert_array_equal(restored.rparticles, saved_endpoint)


def test_standalone_trajectory_conversion_preserves_timing_and_momenta():
    from ocelot.cpbd.high_order import rk_track_in_field
    from ocelot.cpbd.transformations.exact_drift import ExactDriftTM

    particles = ParticleArray(n=2)
    particles.E = .001
    particles.rparticles[:] = np.array([[0., .01, 0., -.02, .003, .02],
                                       [0., 0., 0., 0., -.001, -.01]]).T
    expected = particles.rparticles.copy()
    ExactDriftTM._drift(expected, particles.E, .3)
    beta0 = np.sqrt(1-(radiation_module.m_e_GeV/particles.E)**2)
    slopes = particles.rparticles.copy()
    pz = np.sqrt(1+2*slopes[5]/beta0+slopes[5]**2-slopes[1]**2-slopes[3]**2)
    slopes[1] /= pz
    slopes[3] /= pz
    trajectory = rk_track_in_field(slopes, .3, 31, particles.E, lambda x,y,z: (0*x,0*x,0*x))
    BeamTraject([trajectory]).p_array_end(particles)
    np.testing.assert_allclose(particles.rparticles, expected, atol=2e-15, rtol=0)


def test_non_undulator_samples_include_true_entrance_and_do_not_use_tau_as_z():
    particles = ParticleArray(n=1)
    particles.E = 1.
    particles.x()[:] = .001
    particles.px()[:] = .002
    particles.tau()[:] = .03
    trajectory, _ = radiation_module._track_non_undulator_segment(particles, [Drift(l=.2)], .4, 1., 1.)
    assert trajectory[0, 0] == pytest.approx(.001)
    assert trajectory[4, 0] == pytest.approx(.4)
    assert trajectory[-5, 0] == pytest.approx(.6)
    assert particles.x()[0] == pytest.approx(.0014)


def test_radiation_loss_changes_energy_without_changing_direction(monkeypatch):
    from ocelot.common.globals import m_e_GeV

    monkeypatch.setattr(radiation_module, 'energy_loss_und', lambda *args: .01)
    undulator = Undulator(lperiod=.02, nperiods=5, Kx=0., npoints=11)
    particles = ParticleArray(n=1)
    particles.E = 1.
    particles.px()[:] = .02
    particles.p()[:] = .04
    p0_before = np.sqrt(particles.E**2-m_e_GeV**2)
    energy_before = particles.E + particles.p()[0]*p0_before
    magnitude_before = np.sqrt(energy_before**2-m_e_GeV**2)/p0_before
    direction_x = particles.px()[0]/magnitude_before
    radiation_module.track4rad_beam(particles, MagneticLattice([undulator]), energy_loss=True)
    p0_after = np.sqrt(particles.E**2-m_e_GeV**2)
    energy_after = particles.E + particles.p()[0]*p0_after
    magnitude_after = np.sqrt(energy_after**2-m_e_GeV**2)/p0_after
    assert particles.E == pytest.approx(.99)
    assert energy_after == pytest.approx(energy_before-.01)
    assert particles.px()[0]/magnitude_after == pytest.approx(direction_x, abs=1e-15)


def test_radiation_segments_use_a_common_field_normalization(monkeypatch):
    from ocelot.cpbd.beam import Beam

    beam = Beam()
    beam.E, beam.I = 1., .1
    screen = Screen()
    screen.z = 10
    screen.nx = screen.ny = screen.num_energy = 1
    screen.start_energy = screen.end_energy = 100

    def trajectories(particles, lattice, **kwargs):
        parts = [np.zeros((36, 1)), np.zeros((36, 1))]
        parts[0][4::9, 0] = np.linspace(0, .1, 4)
        parts[1][4::9, 0] = np.linspace(.1, .2, 4)
        particles.E = .8  # after the second aggregate energy loss
        return parts, [1., .9]

    def equal_normalized_fields(gamma, trajectory, output):
        output.arReEx[:] += 1.

    monkeypatch.setattr(radiation_module, 'track4rad_beam', trajectories)
    monkeypatch.setattr(radiation_module, 'radiation_py', equal_normalized_fields)
    radiation_module.calculate_radiation(MagneticLattice([Drift(l=.2)]), screen, beam, energy_loss=True)
    # Physical field amplitudes scale with gamma. Relative to the entrance
    # normalization their coherent sum is 1 + 0.9, irrespective of final energy.
    np.testing.assert_allclose(np.hypot(screen.arReEx, screen.arImEx), 1.9, atol=1e-14, rtol=0)


def test_coherent_radiation_keeps_initial_arrival_offsets(monkeypatch):
    from ocelot.cpbd.high_order import rk_field

    field = lambda x,y,z: (0*x, 0*x, 0*x)
    undulator = Undulator(lperiod=.02, nperiods=5, Kx=0, npoints=21, mag_field=field)
    particles = ParticleArray(n=2)
    particles.E = .001
    particles.rparticles[4] = [.002, -.003]
    particles.rparticles[5] = [.02, -.01]
    expected = particles.rparticles.copy()
    rk_field(expected, 0, .1, 21, particles.E, field)
    screen = Screen()
    screen.z = 10.
    screen.nx = screen.ny = screen.num_energy = 1
    screen.start_energy = screen.end_energy = 100.
    monkeypatch.setattr(radiation_module, 'radiation_py', lambda *args: None)
    radiation_module.coherent_radiation(MagneticLattice([undulator]), screen, particles, verbose=False)
    np.testing.assert_allclose(particles.rparticles, expected, atol=2e-14, rtol=0)
