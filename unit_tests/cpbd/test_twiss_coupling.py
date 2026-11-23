import numpy as np
import pytest

from ocelot.cpbd.beam import ParticleArray, Twiss, generate_parray, optics_from_moments
from ocelot.cpbd.elements import Cavity, Drift, Matrix, Quadrupole, Solenoid
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.matcher import MatchProblem
from ocelot.cpbd.optics import UnstableLatticeError, periodic_twiss, twiss


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


def test_reference_beam_keeps_relative_plane_areas_before_coupling():
    # Unequal transverse damping (e.g. an effective Matrix) changes the
    # relative emittances before a subsequent rotation mixes the planes.
    damping = np.diag([1., .4, 1., .8, 1., 1.])
    rotation, _ = rotated_oscillators(.6, [.2, 1.5])
    seed = Twiss(beta_x=4., beta_y=1.)
    whole = Twiss.track(rotation.r @ damping @ damping, seed)
    after_damping = Twiss.track(damping, Twiss.track(damping, seed))
    split = Twiss.track(rotation.r, after_damping)
    np.testing.assert_allclose([split.beta_x, split.alpha_x, split.beta_y, split.alpha_y],
                               [whole.beta_x, whole.alpha_x, whole.beta_y, whole.alpha_y])


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


@pytest.mark.parametrize("assign", ["geometric", "energy_first", "energy_last"])
def test_assign_beam_sizes_to_coupled_reference_ellipse(assign):
    seed = Twiss(beta_x=4., beta_y=1.)
    skew = np.eye(6)
    skew[1, 2] = skew[3, 0] = -.8
    reference = Twiss.track(skew, seed)
    # The skew kick gives px -> px - 0.8 y, py -> py - 0.8 x.
    unit = np.array([[4., 0., 0., -3.2], [0., .89, -.8, 0.],
                     [0., -.8, 1., 0.], [-3.2, 0., 0., 3.56]])
    beam = Twiss(E=1., emit_x=1e-6, emit_y=3e-6)
    scale_x = np.sqrt(beam.emit_x / np.sqrt(np.linalg.det(unit[:2, :2])))
    scale_y = np.sqrt(beam.emit_y / np.sqrt(np.linalg.det(unit[2:, 2:])))
    scale = np.diag([scale_x, scale_x, scale_y, scale_y])
    if assign == "geometric":
        reference.emit_x, reference.emit_y = beam.emit_x, beam.emit_y
    elif assign == "energy_first":
        reference.E = beam.E
        reference.emit_xn, reference.emit_yn = beam.emit_xn, beam.emit_yn
    else:
        with pytest.warns(UserWarning, match="Normalized emittance"):
            reference.emit_xn, reference.emit_yn = beam.emit_xn, beam.emit_yn
        reference = Twiss.track(np.eye(6), reference)
        reference.E = beam.E
    assert_projected(reference, scale @ unit @ scale)


def test_generate_bunch_from_coupled_twiss(monkeypatch):
    cov = beam_covariance()

    def symmetric_particles(mean, covariance, count):
        assert count == 8
        return (2 * np.linalg.cholesky(covariance) @ np.c_[np.eye(4), -np.eye(4)]).T + mean

    monkeypatch.setattr(np.random, "multivariate_normal", symmetric_particles)
    seed = seed_from_covariance(cov)
    seed.mux, seed.muy = np.pi / 2, .17
    particles = generate_parray(tws=seed, nparticles=8, sigma_p=0., chirp=0.)
    np.testing.assert_allclose(np.cov(particles.rparticles[:4], bias=True), cov, atol=1e-18)
    assert (seed.mux, seed.muy) == (np.pi / 2, .17)


def test_geometric_emittances_survive_coupled_dataframe_restart_without_energy():
    cov = beam_covariance()
    seed = seed_from_covariance(cov)
    seed.E = 0.
    line = MagneticLattice([Drift(l=.4)])
    row = twiss(line, seed, return_df=True).iloc[-1]
    restarted = Twiss.from_series(row)
    assert_projected(restarted, twiss(line, seed)[-1]._transverse_covariance())
    expected = twiss(MagneticLattice([Drift(l=.8)]), seed)[-1]._transverse_covariance()
    assert_projected(twiss(line, restarted)[-1], expected)


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


def rotated_oscillators(roll, phases):
    """Two independent oscillators viewed in laboratory axes rotated by roll."""
    c, s = np.cos(roll), np.sin(roll)
    lab = np.kron([[c, -s], [s, c]], np.eye(2))
    # (beta, alpha) = (4, 0.6) and (0.25, -0.3) in the oscillator frame.
    axes = np.array([[2., 0., 0., 0.], [-.3, .5, 0., 0.],
                     [0., 0., .5, 0.], [0., 0., .6, 2.]])
    axes = lab @ axes
    rotation = np.zeros((4, 4))
    for i, mu in zip((0, 2), phases):
        rotation[i:i+2, i:i+2] = [[np.cos(mu), np.sin(mu)], [-np.sin(mu), np.cos(mu)]]
    turn = Matrix()
    turn.r = np.eye(6)
    turn.r[:4, :4] = axes @ rotation @ np.linalg.inv(axes)
    return turn, axes


@pytest.mark.parametrize("roll", [0., 1e-8, .68])
def test_matched_rotated_oscillators(roll):
    phases = [.73, 4.2]
    turn, axes = rotated_oscillators(roll, phases)
    line = MagneticLattice([turn])
    seed = Twiss(emit_x=1e-6, emit_y=4e-6)
    # Uniform populations at the four cardinal points of each oscillator
    # have the known diagonal second moments in their own frame.
    covariance = (axes * [seed.emit_x, seed.emit_x, seed.emit_y, seed.emit_y]) @ axes.T
    start, end = periodic_twiss(line, seed, reference_orbit="closed")
    assert_projected(start, covariance)
    assert_projected(end, covariance)
    np.testing.assert_allclose([end.mux, end.muy], phases, atol=1e-13)
    if roll == 0.:
        design_start, design_end = periodic_twiss(line, seed)
        assert_projected(design_start, covariance)
        np.testing.assert_allclose([design_end.mux, design_end.muy], phases, atol=1e-13)
    else:
        with pytest.raises(ValueError, match="reference_orbit='closed'"):
            periodic_twiss(line, seed)


def test_equal_plane_participation_labels_lower_fractional_tune_first():
    turn, axes = rotated_oscillators(np.pi / 4, [4.2, .73])
    start, end = periodic_twiss(MagneticLattice([turn]), Twiss(emit_x=1e-6, emit_y=4e-6),
                                reference_orbit="closed")
    assert_projected(start, (axes * [4e-6, 4e-6, 1e-6, 1e-6]) @ axes.T)
    np.testing.assert_allclose([end.mux, end.muy], [.73, 4.2], atol=1e-13)


@pytest.mark.parametrize("split", [1e-7, 1e-11, 0.])
def test_nearly_coincident_modes_preserve_eigenemittances(split):
    turn, axes = rotated_oscillators(.5, [.8, .8 + split])
    seed = Twiss(emit_x=1e-6, emit_y=4e-6)
    start, end = periodic_twiss(MagneticLattice([turn]), seed, reference_orbit="closed")
    sigma = start._transverse_covariance()
    J = np.kron(np.eye(2), [[0, 1], [-1, 0]])
    emittances = np.sort(np.linalg.eigvals(J @ sigma).imag)[2:]
    np.testing.assert_allclose(emittances, [seed.emit_x, seed.emit_y], rtol=2e-12, atol=0.)
    np.testing.assert_allclose(end._transverse_covariance(), sigma, atol=2e-17)
    np.testing.assert_allclose([end.mux, end.muy], [.8, .8 + split], atol=2e-12)
    if split == 0.:
        # For equal weights the covariance is independent of the choice of
        # axes within a degenerate positive-action eigenspace.
        equal = Twiss(emit_x=1e-6, emit_y=1e-6)
        matched = MagneticLattice([turn]).periodic_twiss(equal, reference_orbit="closed")
        assert_projected(matched, 1e-6 * axes @ axes.T)


def test_unit_eigenvalue_moduli_insufficient_for_canonical_modes():
    turn, _ = rotated_oscillators(0., [.5, 1.2])
    # A noncanonical change of variables keeps both eigenvalue pairs exactly
    # on the unit circle, but destroys their symplectic orthogonality.
    shear = np.eye(4)
    shear[0, 2] = .3
    turn.r[:4, :4] = shear @ turn.r[:4, :4] @ np.linalg.inv(shear)
    with pytest.raises(ValueError, match="symplectic"):
        periodic_twiss(MagneticLattice([turn]), reference_orbit="closed")


@pytest.mark.parametrize("phase", [np.pi, 2*np.pi])
def test_integer_and_half_integer_tunes_have_no_unique_matched_ellipse(phase):
    turn, _ = rotated_oscillators(0., [phase, .6])
    with pytest.raises(UnstableLatticeError):
        periodic_twiss(MagneticLattice([turn]), reference_orbit="closed")


def test_opposite_action_resonance():
    phases = [.8, 2*np.pi - .8]
    turn, _ = rotated_oscillators(.5, phases)
    with pytest.raises(ValueError, match="opposite action"):
        periodic_twiss(MagneticLattice([turn]), reference_orbit="closed")
    # Roundoff coupling in an uncoupled ring still admits the physical axes.
    turn, _ = rotated_oscillators(1e-17, phases)
    end = periodic_twiss(MagneticLattice([turn]), reference_orbit="closed")[-1]
    np.testing.assert_allclose([end.mux, end.muy], phases, atol=1e-13)


def test_symplectic_jordan_block_has_no_periodic_covariance():
    # A transverse roll commutes with a drift of equal length in x and y.
    # After n turns the drift length is n, despite unit-modulus eigenvalues.
    c, s = np.cos(.6), np.sin(.6)
    turn = Matrix()
    turn.r = np.eye(6)
    turn.r[:4, :4] = np.kron([[c, s], [-s, c]], [[1., 1.], [0., 1.]])
    with pytest.raises(UnstableLatticeError, match="defective"):
        periodic_twiss(MagneticLattice([turn]), reference_orbit="closed")


@pytest.mark.parametrize("roll, tunes, projected_betas", [
    (.025, [.0726104573653, .0450359581345], [5.66070449324, 5.64640212588]),
    (.080, [.0850539268260, .00574121670147], [19.9273465560, 22.7826309803]),
])
def test_rolled_fodo_matches_madx_tunes_and_projected_beam(roll, tunes, projected_betas):
    # MAD-X 5.09.03 TWISS, RIPKEN for a 1 GeV positron, EX=1e-6, EY=3e-6.
    # Derive the projected betas from MAD-X SIGij at the entrance. Its
    # coupled BETX and BETY columns use a different parameterization.
    drift = Drift(l=.7)
    ring = MagneticLattice([Quadrupole(l=.25, k1=1.8, dtilt=roll), drift,
                            Quadrupole(l=.25, k1=-1.7), drift])
    seed = Twiss(E=1., emit_x=1e-6, emit_y=3e-6)
    start, *_, end = periodic_twiss(ring, seed, reference_orbit="closed")
    np.testing.assert_allclose(np.array([end.mux, end.muy]) / (2 * np.pi), tunes, atol=5e-13, rtol=0.)
    np.testing.assert_allclose([start.beta_x, start.beta_y], projected_betas, rtol=2e-11)
