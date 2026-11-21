import numpy as np
import pytest

from ocelot.cpbd.beam import ParticleArray, Twiss, generate_parray, optics_from_moments
from ocelot.cpbd.elements import Cavity, Drift, Quadrupole, Solenoid
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.matcher import MatchProblem
from ocelot.cpbd.optics import twiss


def beam_covariance():
    # Unequal projected emittances and a correlated entrance beam.
    factor = 1e-3 * np.array([[2.0, 0.0, 0.0, 0.0],
                               [0.1, 0.5, 0.0, 0.0],
                               [0.6, 0.2, 1.0, 0.0],
                               [-0.2, 0.1, 0.4, 0.3]])
    return factor @ factor.T


def seed_from_covariance(cov):
    ex, ey = (np.sqrt(np.linalg.det(block)) for block in (cov[:2, :2], cov[2:, 2:]))
    return Twiss(E=0.1, emit_x=ex, emit_y=ey,
                 beta_x=cov[0, 0] / ex, alpha_x=-cov[0, 1] / ex,
                 beta_y=cov[2, 2] / ey, alpha_y=-cov[2, 3] / ey,
                 xy=cov[0, 2], xpy=cov[0, 3], ypx=cov[1, 2], pxpy=cov[1, 3])


def assert_projected(optics, cov):
    for plane, block in (("x", cov[:2, :2]), ("y", cov[2:, 2:])):
        emit = np.sqrt(np.linalg.det(block))
        assert getattr(optics, f"emit_{plane}") == pytest.approx(emit, rel=1e-11, abs=0.)
        assert getattr(optics, f"beta_{plane}") == pytest.approx(block[0, 0] / emit)
        assert getattr(optics, f"alpha_{plane}") == pytest.approx(-block[0, 1] / emit)
        assert getattr(optics, f"sigma_{plane}") == pytest.approx(np.sqrt(block[0, 0]))
    np.testing.assert_allclose([[optics.xy, optics.xpy], [optics.ypx, optics.pxpy]],
                               cov[:2, 2:], atol=1e-18)


@pytest.mark.parametrize("correlated", [False, True])
def test_open_twiss_matches_particle_covariance_at_every_exit(correlated):
    cov = beam_covariance()
    if not correlated:
        cov[:2, 2:] = cov[2:, :2] = 0.
    seed = seed_from_covariance(cov)
    lat = MagneticLattice([Drift(l=0.15), Quadrupole(l=0.3, k1=2., tilt=np.pi / 4),
                           Drift(l=1.), Cavity(l=0.3, v=0.025, freq=1.3e9),
                           Solenoid(l=0.2, k=1.1)])
    # Eight symmetric particles have exactly the requested second moments.
    particles = ParticleArray(8)
    particles.E = seed.E
    particles.rparticles[:4] = 2 * np.linalg.cholesky(cov) @ np.c_[np.eye(4), -np.eye(4)]
    samples = iter(twiss(lat, seed)[1:])
    for element in lat.sequence:
        for tm in element.first_order_tms:
            tm.apply(particles)
            optics = next(samples)
            expected = np.cov(particles.rparticles[:4], bias=True)
            assert_projected(optics, expected)
            np.testing.assert_allclose([optics.xx, optics.xpx, optics.pxpx,
                                        optics.yy, optics.ypy, optics.pypy],
                                       expected[[0, 0, 1, 2, 2, 3], [0, 1, 1, 2, 3, 3]], atol=1e-18)

    expected = np.cov(particles.rparticles[:4], bias=True)
    assert_projected(twiss(lat, seed, nPoints=7)[-1], expected)
    problem = MatchProblem(lat, seed)
    problem.target_twiss(lat.sequence[-1], "sigma_y", np.sqrt(expected[2, 2]), tol=1e-14)
    assert problem.evaluate()[1][0].met
    assert_projected(seed, cov)


def test_dispersion_transfers_between_planes():
    lat = MagneticLattice([Quadrupole(l=0.3, k1=2., tilt=np.pi / 4), Drift(l=1.)])
    seed = Twiss(beta_x=4., beta_y=1., Dx=0.1, Dxp=-0.01)
    _, matrix, _ = lat.transfer_maps(energy=0.)
    expected = (matrix @ [seed.Dx, seed.Dxp, 0., 0., 0., 1.])[:4]
    assert abs(expected[2]) > 0.01
    out = twiss(lat, seed)[-1]
    np.testing.assert_allclose([out.Dx, out.Dxp, out.Dy, out.Dyp], expected)


def test_emittance_free_twiss_has_equal_emittance_optics_without_a_beam_size():
    lat = MagneticLattice([Quadrupole(l=0.3, k1=2., tilt=np.pi / 4), Drift(l=1.),
                           Solenoid(l=0.2, k=1.1)])
    seed = Twiss(beta_x=4., beta_y=1.)
    physical = Twiss(seed)
    physical.emit_x = physical.emit_y = 1e-6
    expected = twiss(lat, physical)[-1]
    after_quad = twiss(MagneticLattice(lat.sequence[:1]), seed)[-1]
    out = twiss(MagneticLattice(lat.sequence[1:]), Twiss(after_quad))[-1]
    np.testing.assert_allclose([out.beta_x, out.alpha_x, out.beta_y, out.alpha_y],
                               [expected.beta_x, expected.alpha_x, expected.beta_y, expected.alpha_y])
    assert out.emit_x == out.emit_y == out.sigma_x == out.sigma_y == 0.
    restored = Twiss.from_series(after_quad.to_series())
    assert twiss(MagneticLattice(lat.sequence[1:]), restored)[-1].beta_y == pytest.approx(out.beta_y)


def test_optics_from_moments_preserves_betatron_cross_correlations():
    betatron = beam_covariance()
    dispersion = np.array([0.2, -0.06, 0.11, 0.03])
    cov = np.zeros((6, 6))
    cov[5, 5] = 1e-4
    cov[:4, 5] = cov[5, :4] = dispersion * cov[5, 5]
    cov[:4, :4] = betatron + np.outer(dispersion, dispersion) * cov[5, 5]
    seed = optics_from_moments(np.zeros(6), cov, energy=0.1)
    assert_projected(seed, betatron)
    matrix = MagneticLattice([Solenoid(l=0.3, k=1.5)]).transfer_maps(seed.E)[1]
    assert_projected(Twiss.track(matrix, seed), matrix[:4, :4] @ betatron @ matrix[:4, :4].T)


def test_normalized_emittances_can_be_resolved_after_coupled_transport():
    line = MagneticLattice([Quadrupole(l=0.3, k1=2., tilt=np.pi / 4), Drift(l=1.)])
    with pytest.warns(UserWarning, match="Normalized emittance"):
        pending = Twiss(beta_x=4., beta_y=1., emit_xn=1e-6, emit_yn=3e-6)
    physical = Twiss(pending)
    physical.E = 1.
    expected = twiss(line, physical)[-1]
    result = twiss(line, pending)[-1]
    result.E = physical.E
    np.testing.assert_allclose([result.emit_x, result.emit_y, result.beta_x, result.beta_y,
                                result.xy, result.xpy, result.ypx, result.pxpy],
                               [expected.emit_x, expected.emit_y, expected.beta_x, expected.beta_y,
                                expected.xy, expected.xpy, expected.ypx, expected.pxpy], atol=1e-20)


def test_generate_bunch_from_coupled_twiss(monkeypatch):
    cov = beam_covariance()

    def symmetric_particles(mean, covariance, count):
        assert count == 8
        return (2 * np.linalg.cholesky(covariance) @ np.c_[np.eye(4), -np.eye(4)]).T + mean

    monkeypatch.setattr(np.random, "multivariate_normal", symmetric_particles)
    particles = generate_parray(tws=seed_from_covariance(cov), nparticles=8, sigma_p=0., chirp=0.)
    np.testing.assert_allclose(np.cov(particles.rparticles[:4], bias=True), cov, atol=1e-18)


def test_match_skew_quadrupole_to_tracked_beam_size():
    quad = Quadrupole(l=0.3, k1=2., tilt=np.pi / 4)
    end = Drift(l=1.)
    lat = MagneticLattice([quad, end])
    cov = np.diag([4., 0.25, 1., 1.]) * 1e-6
    particles = ParticleArray(8)
    particles.E = 0.1
    particles.rparticles[:4] = 2 * np.linalg.cholesky(cov) @ np.c_[np.eye(4), -np.eye(4)]
    for element in lat.sequence:
        for tm in element.tms:
            tm.apply(particles)
    target = np.std(particles.y())
    assert target == pytest.approx(2.15420249e-3)
    quad.k1 = 1.5
    problem = MatchProblem(lat, seed_from_covariance(cov))
    problem.vary_element(quad, "k1", limits=(1., 3.))
    problem.target_twiss(end, "sigma_y", target, tol=1e-12, weight=1e12)
    assert problem.solve(tol=1e-12).success
    assert quad.k1 == pytest.approx(2., abs=1e-8)


def test_coupled_matrix_target_accepts_energy_only_seed():
    quad = Quadrupole(l=0.3, k1=2., tilt=np.pi / 4)
    lat = MagneticLattice([quad])
    expected = lat.transfer_maps(energy=0.1)[1][1, 2]
    problem = MatchProblem(lat, Twiss(E=0.1))
    problem.target_rmatrix(quad, quad, i=1, j=2, value=expected, tol=1e-14)
    assert problem.evaluate()[1][0].met
