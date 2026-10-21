import numpy as np
import pytest

from ocelot.cpbd.beam import Twiss
from ocelot.cpbd.elements import Quadrupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.optics import twiss


def _normalization_matrix(beta, alpha):
    """Map normalized phase-space coordinates to physical coordinates."""
    return np.array([[np.sqrt(beta), 0.0], [-alpha / np.sqrt(beta), 1 / np.sqrt(beta)]])


@pytest.mark.parametrize("plane, offset", [("x", 0), ("y", 2)])
@pytest.mark.parametrize("phase", [
    0.0, 1e-9, np.pi / 2, np.pi - 1e-9, np.pi, np.pi + 1e-9,
    4.0, 3 * np.pi / 2, 2 * np.pi - 1e-9,
])
def test_twiss_track_recovers_phase_in_all_quadrants(plane, offset, phase):
    seed = Twiss(beta_x=2.0, alpha_x=0.3, beta_y=3.0, alpha_y=-0.4,
                 mux=7.0, muy=8.0, E=1.0)
    entrance = _normalization_matrix(getattr(seed, f"beta_{plane}"),
                                     getattr(seed, f"alpha_{plane}"))
    exit_matrix = _normalization_matrix(beta=5.0, alpha=-0.7)
    rotation = np.array([[np.cos(phase), np.sin(phase)],
                         [-np.sin(phase), np.cos(phase)]])
    matrix = np.eye(6)
    matrix[offset:offset + 2, offset:offset + 2] = exit_matrix @ rotation @ np.linalg.inv(entrance)

    result = Twiss.track(matrix, seed)

    # The known rotation must survive the change of entrance/exit Twiss basis.
    assert getattr(result, f"mu{plane}") == pytest.approx(
        getattr(seed, f"mu{plane}") + phase, rel=0, abs=1e-12,
    )
    assert getattr(result, f"beta_{plane}") == pytest.approx(5.0)
    assert getattr(result, f"alpha_{plane}") == pytest.approx(-0.7)
    other_plane = "y" if plane == "x" else "x"
    assert getattr(result, f"mu{other_plane}") == getattr(seed, f"mu{other_plane}")
    assert seed.mux == 7.0
    assert seed.muy == 8.0


def test_twiss_track_identity_adds_no_phase():
    seed = Twiss(beta_x=2.0, alpha_x=0.3, beta_y=3.0, alpha_y=-0.4,
                 mux=7.0, muy=8.0, E=1.0)

    result = Twiss.track(np.eye(6), seed)

    assert result.mux == seed.mux
    assert result.muy == seed.muy


@pytest.mark.parametrize("plane, k1", [("x", 1.0), ("y", -1.0)])
@pytest.mark.parametrize("n_points", [None, 17])
def test_twiss_long_quadrupole_phase_matches_split_elements(plane, k1, n_points):
    # In the focusing plane beta=1/sqrt(abs(k1))=1 m is matched, so mu(s)=s.
    seed = Twiss(beta_x=1.0, beta_y=1.0, mux=0.3, muy=0.7, E=1.0)
    whole = MagneticLattice([Quadrupole(l=4.0, k1=k1)])
    split = MagneticLattice([Quadrupole(l=2.0, k1=k1), Quadrupole(l=2.0, k1=k1)])
    initial_phase = getattr(seed, f"mu{plane}")

    for lattice in (whole, split):
        samples = twiss(lattice, seed, nPoints=n_points)

        np.testing.assert_allclose(
            [getattr(sample, f"mu{plane}") for sample in samples],
            [initial_phase + sample.s for sample in samples], rtol=0, atol=1e-12,
        )
        assert getattr(samples[-1], f"beta_{plane}") == pytest.approx(1.0)
        assert getattr(samples[-1], f"alpha_{plane}") == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("plane, k1", [("x", 1.0), ("y", -1.0)])
def test_twiss_accumulates_phase_beyond_one_turn(plane, k1):
    seed = Twiss(beta_x=1.0, beta_y=1.0, E=1.0)
    lattice = MagneticLattice([Quadrupole(l=4.0, k1=k1), Quadrupole(l=4.0, k1=k1)])

    samples = twiss(lattice, seed)

    # Each map advances by less than 2*pi; their sum must not be wrapped.
    assert getattr(samples[-1], f"mu{plane}") == pytest.approx(8.0, rel=0, abs=1e-12)
