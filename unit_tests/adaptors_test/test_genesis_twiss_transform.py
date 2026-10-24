"""Real beam rematching must remain usable by CPBD Twiss propagation."""
import numpy as np
import pytest

from ocelot.adaptors.genesis import find_transform, transform_beam_twiss
from ocelot.cpbd.beam import Twiss, generate_beam
from ocelot.cpbd.elements import Drift
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.optics import twiss


pytestmark = pytest.mark.filterwarnings('error::PendingDeprecationWarning')


@pytest.mark.parametrize('complex_eigenpairs', [False, True])
def test_rematched_beam_has_real_twiss_and_can_propagate(monkeypatch, complex_eigenpairs):
    if complex_eigenpairs:
        eig = np.linalg.eig
        def complex_eig(matrix):
            values, vectors = eig(matrix)
            return values.astype(complex), vectors.astype(complex)
        monkeypatch.setattr(np.linalg, 'eig', complex_eig)
    beam = generate_beam(E=14., I=7000., l_beam=5e-6, emit_n=.5e-6,
                         beta=20., l_window=6e-6, shape='flattop')
    target = Twiss(beam.pk())
    target.beta_x, target.alpha_x = 12., .7
    target.beta_y, target.alpha_y = 28., -.4
    transform_beam_twiss(beam, target)
    for name in ('beta_x', 'alpha_x', 'beta_y', 'alpha_y'):
        assert np.isrealobj(getattr(beam, name))
        np.testing.assert_allclose(getattr(beam, name), getattr(target, name), atol=1e-12)
    final = twiss(MagneticLattice([Drift(l=.3)]), Twiss(beam.pk()))[-1]
    assert np.isfinite(final.mux) and np.isfinite(final.muy)


def test_transform_preserves_twiss_matrix_and_area():
    g1 = np.array([[20., .4], [.4, (1+.4**2)/20.]])
    g2 = np.array([[12., -.7], [-.7, (1+.7**2)/12.]])
    inverse, forward = find_transform(g1, g2)
    assert type(inverse) is np.ndarray
    assert type(forward) is np.ndarray
    assert np.isrealobj(inverse) and np.isrealobj(forward)
    np.testing.assert_allclose(inverse.T @ g1 @ inverse, g2, atol=1e-12)
    assert np.linalg.det(inverse) == pytest.approx(1.)
