import copy

import numpy as np
import pytest

from ocelot.cpbd import errors, track
from ocelot.cpbd.elements import Drift, Marker, Quadrupole, SBend, Sextupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.transformations import SecondTM


def test_zero_rms_preserves_other_families_and_random_state():
    quad = Quadrupole(l=.4, k1=1., dx=.01, dy=.02, dtilt=.03, tilt=.4)
    bend = SBend(l=.2, angle=.02, dx=.04, dy=.05, dtilt=.06)
    marker = Marker()
    lattice = MagneticLattice([quad, bend, marker])
    bend.set_tm(SecondTM)
    bend_map = bend.tms[1]
    before = np.random.get_state()
    seeded, records = errors.errors_seed(lattice, {Quadrupole: {"offset": 0., "dtilt": 0.}})
    after = np.random.get_state()
    assert seeded is lattice
    assert records == [dict(dx=0., dy=0., dtilt=0.), None, None]
    assert (quad.dx, quad.dy, quad.dtilt, quad.tilt) == (0., 0., 0., .4)
    assert (bend.dx, bend.dy, bend.dtilt) == (.04, .05, .06)
    assert bend.tms[1] is bend_map
    assert not hasattr(marker, "dtilt")
    np.testing.assert_array_equal(before[1], after[1])
    assert before[2:] == after[2:]


@pytest.mark.parametrize("family,strength", [(Quadrupole, {"k1": .5}), (SBend, {"angle": .1})])
def test_independent_occurrences_shared_slices_and_map_replay(monkeypatch, family, strength):
    magnet = family(l=.4, eid="M", **strength)
    same_id = family(l=.4, eid="M", **strength)
    lattice = MagneticLattice([magnet, Marker(), magnet, same_id, magnet], method={"global": SecondTM})
    replay = copy.deepcopy(lattice)
    lattice.transfer_maps(1.)
    draws = iter(np.arange(1, 10) * 1e-5)
    monkeypatch.setattr(errors, "tgauss", lambda sigma: next(draws))
    _, records = errors.errors_seed(lattice, {family: {"offset": .001, "dtilt": .002}},
                                    shared_groups=[(4, 3)])
    assert records[1] is None
    assert records[3] == records[4]
    assert len({records[i]["dtilt"] for i in (0, 2, 3)}) == 3
    assert lattice.find_element_indices(magnet) == [0]
    for i in (0, 2, 3, 4):
        elem = lattice.sequence[i]
        assert lattice[elem] is elem
        assert records[i] == dict(dx=elem.dx, dy=elem.dy, dtilt=elem.dtilt)
    errors.apply_alignment_errors(replay, records)
    for actual, expected in zip(replay.transfer_maps(1.), lattice.transfer_maps(1.)):
        np.testing.assert_array_equal(actual, expected)


def test_record_equal_to_existing_alignment_preserves_shared_strength():
    quad = Quadrupole(l=.4, k1=2., dx=.001)
    line = MagneticLattice([quad, quad, quad])
    errors.apply_alignment_errors(line, [None, dict(dx=.001, dy=0., dtilt=0.), dict(dx=.002, dy=0., dtilt=0.)])
    assert line.find_element_indices(quad) == [0, 1]
    quad.k1 = 3.
    assert [element.k1 for element in line.sequence] == [3., 3., 2.]


def test_replay_preserves_untouched_occurrence_of_a_reused_element():
    shared = Sextupole(l=.2, dx=.001, dy=.002, dtilt=.003)
    lattice = MagneticLattice([shared, shared, shared])
    record = dict(dx=.1, dy=.2, dtilt=.3)
    errors.apply_alignment_errors(lattice, [record, None, record])
    first, middle, last = lattice.sequence
    assert middle is shared
    assert first is last and first is not shared
    assert (shared.dx, shared.dy, shared.dtilt) == (.001, .002, .003)
    assert (first.dx, first.dy, first.dtilt) == (.1, .2, .3)


@pytest.mark.parametrize("groups", [[(0, 1)], [(0, 0)], [(0, 2)], [(-1, 0)]])
def test_invalid_shared_groups_leave_lattice_intact(groups):
    quad = Quadrupole(dx=.1)
    lattice = MagneticLattice([quad, Marker()])
    with pytest.raises(ValueError, match="[Ss]lice index"):
        errors.errors_seed(lattice, {Quadrupole: {"offset": 0., "dtilt": 0.}}, shared_groups=groups)
    assert quad.dx == .1


def test_mpi_workers_replay_root_draws_on_expanded_lattice(monkeypatch):
    # Run the two rank paths in sequence, copying the broadcast as MPI does.
    class Comm:
        rank = 0
        payload = None

        def Get_size(self):
            return 2

        def Get_rank(self):
            return self.rank

        def bcast(self, data, root):
            if self.rank == root:
                self.payload = copy.deepcopy(data)
            return copy.deepcopy(self.payload)

        def scatter(self, chunks, root):
            return [self.rank]

        def gather(self, particles, root):
            return [particles, []] if self.rank == root else None

    tracked = []

    def capture(lattice, nturns, particles, nsuperperiods, save_track):
        assert nsuperperiods == 1
        tracked.append([(e.dx, e.dy, e.dtilt) for e in lattice.sequence])
        return particles

    comm = Comm()
    monkeypatch.setattr(track, "track_nturns", capture)
    lattice = MagneticLattice([Quadrupole(), SBend(dtilt=.03)])
    config = {Quadrupole: {"offset": .001, "dtilt": .002}}
    track.track_nturns_mpi(comm, lattice, 1, [0, 1], errors=config, nsuperperiods=2,
                           shared_groups=[(0, 2)])
    records = comm.payload
    assert records[0] == records[2] and records[1] is records[3] is None
    comm.rank = 1
    track.track_nturns_mpi(comm, lattice, 1, None, nsuperperiods=2)
    comm.rank = 0
    track.track_nturns_mpi(comm, lattice, 1, [0, 1], errors=records, nsuperperiods=2)
    expected = (records[0]["dx"], records[0]["dy"], records[0]["dtilt"])
    assert tracked == [[expected, (0., 0., .03)] * 2] * 3
    assert lattice.sequence[0].dtilt == 0.
    with pytest.raises(ValueError, match="per lattice occurrence"):
        track.track_nturns_mpi(comm, lattice, 1, [0, 1], errors=records)


def test_mpi_output_grids_follow_launch_order(monkeypatch):
    import sys
    from types import SimpleNamespace

    comm = SimpleNamespace(Get_rank=lambda: 0)
    monkeypatch.setitem(sys.modules, "mpi4py", SimpleNamespace(MPI=SimpleNamespace(COMM_WORLD=comm)))

    def tracked(comm, lattice, turns, particles, **kwargs):
        for i, particle in enumerate(particles):
            particle.turn = i
            particle.mux, particle.muy = particle.x, particle.y
        return particles

    def tunes(particles, lattice, turns, harm, nsuperperiods):
        assert nsuperperiods == 3
        return particles

    monkeypatch.setattr(track, "track_nturns_mpi", tracked)
    monkeypatch.setattr(track, "freq_analysis", tunes)
    xs, ys = [-1., 1., 2.], [-2., 2.]
    np.testing.assert_array_equal(track.da_mpi(None, 2, xs, ys), [[0, 1, 2], [3, 4, 5]])
    _, qx, qy = track.fma(None, 2, xs, ys, nsuperperiods=3)
    np.testing.assert_array_equal(qx, [xs, xs])
    np.testing.assert_array_equal(qy, [[-2.] * 3, [2.] * 3])


def test_empty_tracking_chunk_returns_no_particles():
    assert track.track_nturns(MagneticLattice([Drift(l=1.)]), 3, []).size == 0


def test_turn_tracking_with_rolled_quadrupole():
    drift = Drift(l=.5)
    line = MagneticLattice([Quadrupole(l=.2, k1=1.3, dtilt=.005), drift,
                            Quadrupole(l=.2, k1=-1.3), drift])
    particles = track.create_track_list([1e-4], [2e-4], [0.])
    tracked = track.track_nturns(line, 3, particles, print_progress=False)
    assert tracked[0].turn == 2
    assert np.all(np.isfinite(tracked[0].p_list))
