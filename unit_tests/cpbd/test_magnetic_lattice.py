import pytest
import numpy as np

from ocelot.cpbd.magnetic_lattice import MagneticLattice, insert_markers_by_predicate, merger
from ocelot.cpbd.elements import Drift, SBend, Marker, Matrix, Quadrupole, Cavity, Sextupole
from ocelot.cpbd.beam import ParticleArray, Twiss
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


def test_transfer_maps_tracks_offset_orbit_through_sextupole():
    lattice = MagneticLattice(
        [Quadrupole(l=0.4, k1=2.0, dx=1e-3), Drift(l=1.0), Sextupole(l=0.2, k2=20.0)],
        method={"global": SecondTM},
    )
    B, R, T = lattice.transfer_maps(energy=1.0)

    # Sample the orbit and its first two derivatives along an entrance perturbation.
    direction = np.array([1., 0., 0., 0., 0., 1.])
    step = 1e-5
    particles = ParticleArray(n=3)
    particles.E = 1.0
    particles.rparticles[:] = direction[:, None] * np.array([-step, 0., step])
    for elem in lattice.sequence:
        for tm in elem.tms:
            tm.apply(particles)
    minus, orbit, plus = particles.rparticles.T

    np.testing.assert_allclose(B[:, 0], orbit, rtol=0, atol=1e-14)
    np.testing.assert_allclose(R @ direction, (plus - minus) / (2 * step), rtol=0, atol=1e-8)
    np.testing.assert_allclose(np.einsum("ijk,j,k->i", T, direction, direction),
                               (plus - 2 * orbit + minus) / (2 * step**2), rtol=0, atol=1e-7)


def test_magnetic_lattice_find_indices():
    d1 = Drift(l=0.1)
    sb = SBend(l=0.2, angle=3)
    d2 = Drift(l=0.2)
    sequence = [d1, sb, d2, sb]
    mlat = MagneticLattice(sequence)
    indices = mlat.find_indices(SBend)
    assert indices == [1, 3]


def test_magnetic_lattice_find_indices_by_predicate():
    d1 = Drift(l=0.1)
    sb = SBend(l=0.2, angle=3)
    d2 = Drift(l=0.2)

    sequence = [d1, sb, d2, sb]
    mlat = MagneticLattice(sequence)
    def pred(element):
        return element.l == 0.2 and isinstance(element, SBend)

    indices = mlat.find_indices_by_predicate(pred)
    assert indices == [1, 3]


def test_magnetic_lattice_resolves_element_occurrences_by_identity():
    shared = Drift(l=0.1, eid="D")
    same_id_but_distinct = Drift(l=0.2, eid="D")
    lattice = MagneticLattice((shared, same_id_but_distinct, shared))

    assert lattice.find_element_indices(shared) == [0, 2]
    assert lattice.find_element_indices(same_id_but_distinct) == [1]
    assert lattice.resolve_element_index(shared, occurrence=0) == 0
    assert lattice.resolve_element_index(shared, occurrence=1) == 2
    assert lattice.resolve_element_index(same_id_but_distinct) == 1

    with pytest.raises(ValueError, match=r"indices \[0, 2\]"):
        lattice.resolve_element_index(shared)
    with pytest.raises(ValueError, match="occurrence 2 is out of range"):
        lattice.resolve_element_index(shared, occurrence=2)


def test_get_sequence_part_rejects_an_ambiguous_element_boundary():
    shared = Drift(l=0.1, eid="D")
    marker = Marker(eid="END")
    lattice = MagneticLattice((shared, marker, shared))

    with pytest.raises(ValueError, match=r"indices \[0, 2\]"):
        lattice.get_sequence_part(start=shared, stop=marker)

    with pytest.raises(ValueError, match=r"indices \[0, 2\]"):
        lattice.transfer_maps(energy=1.0, start=shared, stop=marker)


def test_magnetic_lattice_insert_markers_by_predicate():
    d1 = Drift(l=0.1)
    sb = SBend(l=0.2, angle=3)
    d2 = Drift(l=0.2)

    sequence = [d1, sb, d2, sb]
    mlat = MagneticLattice(sequence)
    def pred(element):
        return element.l == 0.2 and isinstance(element, SBend)

    before_suffix = "before-suffix"
    after_suffix = "after-suffix"    
    insert_markers_by_predicate(sequence, pred,
                                before_suffix=before_suffix,
                                after_suffix=after_suffix)
    isb = sequence.index(sb)
    before = sequence[isb - 1]
    after = sequence[isb + 1]    



    assert isinstance(before, Marker)
    assert isinstance(after, Marker)    
    assert before.id == f"{sb.id}{before_suffix}"
    assert after.id == f"{sb.id}{after_suffix}"    

    
def test_magnetic_lattice_periodic_twiss_with_cav():
    d2 = Drift(l=0.5)
    qf = Quadrupole(l=0.2, k1=0.3)
    qdh = Quadrupole(l=0.2 / 2, k1=-0.3)
    c = Cavity(l=0.5, v=0.01, freq=1.3e9, phi=0)
    cell = (qdh, d2, c, d2, qf, d2, c, d2, qdh)
    lat = MagneticLattice(cell)

    tws0 = Twiss()
    tws0.E = 0.005
    tws0.beta_x = 10
    tws0.beta_y = 10
    tw_p = lat.periodic_twiss(tws=Twiss(E=0.005))
    beta_x = 1.60031832321874
    beta_y = 1.8790818475215212
    alpha_x = -0.8882452561429366
    alpha_y = -0.9590585725276709
    gamma_x = 1.1178898654751601
    gamma_y = 1.0216656331766445
    E = 0.005
    s = 0.0
    ref = np.array([beta_x, beta_y, alpha_x, alpha_y, gamma_x, gamma_y, E, s])
    res = np.array([tw_p.beta_x, tw_p.beta_y, tw_p.alpha_x, tw_p.alpha_y, tw_p.gamma_x, tw_p.gamma_y, tw_p.E, tw_p.s])
    np.testing.assert_allclose(res, ref, rtol=1.e-12, atol=1.e-14)

def test_magnetic_lattice_periodic_twiss():
    d2 = Drift(l=0.5)
    qf = Quadrupole(l=0.2, k1=0.3)
    qdh = Quadrupole(l=0.2 / 2, k1=-0.3)
    c = Drift(l=0.5)
    cell = (qdh, d2, c, d2, qf, d2, c, d2, qdh)
    lat = MagneticLattice(cell)

    tws0 = Twiss()
    tws0.E = 0.005
    tws0.beta_x = 10
    tws0.beta_y = 10
    tw_p = lat.periodic_twiss(tws=Twiss(E=0.005))
    beta_x = 33.09417032395688
    beta_y = 36.43155958891432
    alpha_x = 0.0
    alpha_y = 5.675938272265563e-16
    gamma_x = 0.030216802240728775
    gamma_y = 0.027448728829722893
    E = 0.005
    s = 0.0
    ref = np.array([beta_x, beta_y, alpha_x, alpha_y, gamma_x, gamma_y, E, s])
    res = np.array([tw_p.beta_x, tw_p.beta_y, tw_p.alpha_x, tw_p.alpha_y, tw_p.gamma_x, tw_p.gamma_y, tw_p.E, tw_p.s])
    np.testing.assert_allclose(res, ref, rtol=1.e-12, atol=1.e-14)


def test_merger_preserves_local_sextupole_feeddown_and_retained_linear_matrix():
    quad = Quadrupole(l=.2, k1=2., dx=.003)
    linear = Matrix(t211=900.)
    linear.r = np.eye(6)
    source = MagneticLattice([quad, Drift(l=.5), Sextupole(l=.3, k2=100.), Drift(l=.1), linear],
                             method={"global": TransferMap, Sextupole: SecondTM, Matrix: SecondTM})
    linear.set_tm(TransferMap)  # Stored T is deliberately inactive in this element.
    probes = np.array([[0., .001, -.002], [0., .0002, -.0001], [0., -.001, .002],
                       [0., -.0001, .0003], [0., 0., 0.], [0., -.003, .001]])
    expected = probes.copy()
    for elem in source.sequence:
        elem.apply(expected, 1.)

    merged = merger(source, remaining_elems=[quad, linear], init_energy=1.)
    merged.update_transfer_maps()
    assert merged.sequence[0] is quad and merged.sequence[-1] is linear
    assert type(linear.tms[0]) is TransferMap
    actual = probes.copy()
    for elem in merged.sequence:
        elem.apply(actual, 1.)
    # With just one nonlinear element, direct particle tracking is exactly quadratic.
    np.testing.assert_allclose(actual, expected, atol=2e-16, rtol=0.)


@pytest.mark.parametrize("second_method", [TransferMap, SecondTM])
def test_merger_combines_drift_sections_with_different_orders(second_method):
    source = MagneticLattice([Drift(l=.2, tm=SecondTM), Marker(), Drift(l=.5, tm=second_method)], method={})
    merged = merger(source, init_energy=1.)
    merged.update_transfer_maps()
    expected = np.array([[0., .003, 0., -.002, 0., .007]]).T
    actual = expected.copy()
    for elem in source.sequence:
        elem.apply(expected, 1.)
    assert len(merged.sequence) == 1
    assert isinstance(merged.sequence[0], Drift if second_method is SecondTM else Matrix)
    merged.sequence[0].apply(actual, 1.)
    np.testing.assert_allclose(actual, expected, atol=2e-16, rtol=0.)


@pytest.mark.parametrize("kick_terms", [{}, {"t255": 10., "t615": -20.}, {"t215": 10., "t611": -5.}],
                         ids=["static", "time-dependent-orbit", "time-dependent-focusing"])
def test_closed_orbit_with_arrival_time_slip(kick_terms):
    # Canonical thin kicks generated by -10*x*tau**2 or -5*x**2*tau.
    # Both have the identity linearization at x=tau=0. At the next arrival
    # the former moves the orbit and the latter changes its linear focusing.
    kick = Matrix(**kick_terms)
    kick.r = np.eye(6)
    ring_map = Matrix(b5=.02)
    ring_map.r = np.eye(6)
    for plane, phase in ((slice(0, 2), .6), (slice(2, 4), .9)):
        ring_map.r[plane, plane] = [[np.cos(phase), np.sin(phase)], [-np.sin(phase), np.cos(phase)]]
    lattice = MagneticLattice([kick, ring_map], method={Matrix: SecondTM})
    if kick_terms:
        with pytest.raises(ValueError, match="next turn"):
            lattice.periodic_twiss(reference_orbit="closed")
    else:
        matched = lattice.periodic_twiss(reference_orbit="closed")
        np.testing.assert_allclose([matched.beta_x, matched.beta_y], [1., 1.], atol=1e-14)
