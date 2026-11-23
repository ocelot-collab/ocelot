import numpy as np
import pytest

from ocelot.cpbd.beam import Particle, Twiss
from ocelot.cpbd.elements import Cavity, Drift, Matrix, Multipole, Quadrupole, SBend, Sextupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.optics import UnstableLatticeError, periodic_twiss, twiss
from ocelot.cpbd.transformations import KickTM, SecondTM


def coordinates(tws):
    return np.array([tws.x, tws.xp, tws.y, tws.yp, tws.tau, tws.p])


def track_probes(lattice, entrance, energy):
    # Finite difference the entire line independently of its local Jacobians.
    step = 1e-7
    points = entrance[:, None] + np.column_stack((np.zeros(6), step * np.eye(6), -step * np.eye(6)))
    for elem in lattice.sequence:
        for tm in elem.tms:
            points = tm.map_function(points, energy)
            energy += tm.get_delta_e()
    return points[:, 0], (points[:, 1:7] - points[:, 7:]) / (2 * step)


@pytest.mark.parametrize("method", [SecondTM, KickTM])
def test_sextupole_feeddown_matches_tracked_beam(method):
    sext = Sextupole(l=.3, k2=95., dx=.0002, tilt=.14)
    line = MagneticLattice([Quadrupole(l=.2, k1=2., dx=.001), Drift(l=.8), sext,
                            Drift(l=.4), sext], method={Sextupole: method})
    entrance = Particle(x=.001, y=-.0007, p=.002, E=.3)
    seed = Twiss(beta_x=3., beta_y=2., emit_x=1e-6, emit_y=2e-6, Dx=.02)
    sample = twiss(line, seed, reference_orbit=entrance, attach2elem=[line.sequence[0]])
    centre, matrix = track_probes(line, coordinates(sample[0]), entrance.E)
    np.testing.assert_allclose(coordinates(sample[-1]), centre, atol=1e-14)
    cov = matrix[:4, :4] @ seed._transverse_covariance() @ matrix[:4, :4].T
    np.testing.assert_allclose(sample[-1]._transverse_covariance(), cov, rtol=1e-8, atol=1e-17)
    disp = matrix @ np.array([.02, 0., 0., 0., 0., 1.])
    out = sample[-1]
    np.testing.assert_allclose([out.Dx, out.Dxp, out.Dy, out.Dyp], disp[:4] / disp[5], atol=2e-11)
    assert line.sequence[0].tws is sample[1]
    assert abs(out.beta_x - twiss(line, seed)[-1].beta_x) > .01
    assert seed.E == seed.x == 0. and entrance.x == .001
    grid = twiss(line, seed, nPoints=29, reference_orbit=entrance, return_df=True)
    assert len(grid) == 29
    np.testing.assert_allclose(grid.s, np.linspace(0., line.totalLen, 29))
    np.testing.assert_allclose(grid.beta_x.iloc[-1], out.beta_x)


def test_affine_design_optics_and_orbit_sampling_at_bend_and_thin_kick():
    bend = SBend(l=.6, angle=.08, e1=.02, e2=.03, dx=.001)
    line = MagneticLattice([bend, Multipole(kn=[0., .7])])
    seed = Twiss(E=1., beta_x=3., beta_y=2., mux=.5, s=2.)
    particle = Particle(x=.002)
    samples = twiss(line, seed, reference_orbit=particle)
    grid = twiss(line, seed, reference_orbit=particle, nPoints=3)
    for sample, expected in zip(samples[:-1], twiss(line, seed)[:-1]):
        assert sample.beta_x == pytest.approx(expected.beta_x)
        assert sample.mux == pytest.approx(expected.mux)
    half = np.array([[.002, 0., 0., 0., 0., 0.]]).T
    for tm in bend.get_section_tms(.3):
        half = tm.map_function(half, seed.E)
    np.testing.assert_allclose(coordinates(grid[1]), half[:, 0], atol=1e-14)
    np.testing.assert_allclose(coordinates(grid[-1]), coordinates(samples[-1]), atol=1e-14)
    assert [sample.s for sample in grid] == pytest.approx([2., 2.3, 2.6])


@pytest.mark.parametrize("offset", [0., .001])
def test_thin_sextupole_measures_offset_from_magnetic_axis(offset):
    sext = Multipole(kn=[0., 0., 100.], dx=.002)
    line = MagneticLattice([sext])
    seed = Twiss(beta_x=4., beta_y=2.)
    out = twiss(line, seed, reference_orbit=Particle(x=sext.dx + offset))[-1]
    assert out.xp == pytest.approx(-50. * offset**2, abs=1e-15)
    assert out.alpha_x == pytest.approx(seed.beta_x * 100. * offset, abs=1e-13)
    assert out.alpha_y == pytest.approx(-seed.beta_y * 100. * offset, abs=1e-13)


def test_accelerating_line_retains_longitudinal_feedback():
    line = MagneticLattice([SBend(l=.6, angle=.12), Cavity(l=.2, v=.015, freq=1.3e9, phi=30.),
                            SBend(l=.6, angle=-.12)], method={SBend: SecondTM})
    seed = Twiss(E=.08, beta_x=2., beta_y=3., emit_x=1e-6, emit_y=2e-6, Dx=.2)
    result = twiss(line, seed, reference_orbit=Particle(x=.001, tau=.0004, p=.002))
    orbit, matrix = track_probes(line, coordinates(result[0]), seed.E)
    cov = matrix[:4, :4] @ seed._transverse_covariance() @ matrix[:4, :4].T
    np.testing.assert_allclose(result[-1]._transverse_covariance(), cov, atol=3e-16)
    np.testing.assert_allclose(coordinates(result[-1]), orbit, atol=1e-14)
    disp = matrix @ [seed.Dx, 0., 0., 0., 0., 1.]
    assert result[-1].Dx == pytest.approx(disp[0] / disp[5], abs=2e-10)
    assert result[-1].E == pytest.approx(.08 + .015 * np.cos(np.pi / 6))


@pytest.mark.parametrize("method,dy", [(SecondTM, 0.), (KickTM, .0002)])
def test_ring_closes_on_actual_orbit_and_tunes_match_tracking(method, dy):
    qf = Quadrupole(l=.25, k1=1.8, dx=.0003, dy=dy)
    drift = Drift(l=.7)
    ring = MagneticLattice([qf, drift, Sextupole(l=.15, k2=40.), Quadrupole(l=.25, k1=-1.7), drift],
                           method={Sextupole: method})
    seed = Twiss(E=.3, emit_x=1e-6, emit_y=2e-6)
    samples = periodic_twiss(ring, seed, reference_orbit="closed")
    start, end = samples[0], samples[-1]
    orbit, one_turn = track_probes(ring, coordinates(start), seed.E)
    assert abs(start.x) > 1e-3
    np.testing.assert_allclose(coordinates(start)[:4], orbit[:4], atol=2e-12)
    # A truncated thick SecondTM can have a small area defect off axis.
    # Its matched ellipse shape is periodic and its transported area follows tracking.
    np.testing.assert_allclose([end.beta_x, end.alpha_x, end.beta_y, end.alpha_y],
                               [start.beta_x, start.alpha_x, start.beta_y, start.alpha_y], atol=1e-12)
    transverse = one_turn[:4, :4]
    np.testing.assert_allclose(end._transverse_covariance(),
                               transverse @ start._transverse_covariance() @ transverse.T, atol=2e-16)
    if method is KickTM:
        np.testing.assert_allclose(start._transverse_covariance(), end._transverse_covariance(), atol=2e-16)
    eigenvalues = np.linalg.eigvals(one_turn[:4, :4])
    phases = np.sort(np.angle(eigenvalues[np.imag(eigenvalues) > 0]))
    np.testing.assert_allclose(np.sort([end.mux, end.muy]), phases, atol=2e-10)
    guess = Particle(x=-.0002, E=seed.E)
    np.testing.assert_allclose(coordinates(ring.periodic_twiss(seed, reference_orbit=guess)),
                               coordinates(start), atol=1e-12)
    plus, minus = [ring.periodic_twiss(seed, reference_orbit=Particle(p=dp)) for dp in (1e-6, -1e-6)]
    np.testing.assert_allclose([start.Dx, start.Dxp, start.Dy, start.Dyp],
                               (coordinates(plus) - coordinates(minus))[:4] / 2e-6, atol=2e-10)


def test_nonclosing_or_unstable_ring_does_not_return_optics():
    constant_kick = Matrix(b1=.001)
    constant_kick.r = np.eye(6)
    with pytest.raises(ValueError, match="closed orbit did not converge"):
        periodic_twiss(MagneticLattice([constant_kick]), reference_orbit="closed")
    with pytest.raises(UnstableLatticeError):
        periodic_twiss(MagneticLattice([Quadrupole(l=.3, k1=1.)]), reference_orbit="closed")
    with pytest.raises(ValueError, match="fixed-momentum"):
        periodic_twiss(MagneticLattice([Cavity(l=.2, v=.01)]), Twiss(E=.1), reference_orbit="closed")
