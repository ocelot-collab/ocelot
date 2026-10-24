"""The iris propagator must preserve complex lineouts on current NumPy."""

import numpy as np
import pytest

from ocelot.optics import wave


@pytest.mark.parametrize("account_first_loss", [False, True])
def test_iris_complex_lineouts(monkeypatch, account_first_loss):
    monkeypatch.setattr(wave, "nthread", 1)
    field = wave.RadiationField((1, 9, 11))
    field.dx = field.dy = 0.002
    field.dz = field.xlamds = 1e-4
    field.domain_z = 'f'
    field.fld[:] = 1 + 2j
    original = field.fld.copy()

    result, x, y, remaining, losses = wave.dfl_prop_iris(
        field, N=2, a=0.05, b=0.3, n_iter_per_iris=2,
        acount_first_cell_loss=account_first_loss,
    )

    assert x.shape == (11, 4)
    assert y.shape == (9, 4)
    assert x.dtype == y.dtype == np.dtype(np.complex128)
    np.testing.assert_array_equal(x[:, 0], original[0, 4, :])
    np.testing.assert_array_equal(y[:, 0], original[0, :, 5])
    np.testing.assert_array_equal(field.fld, original)
    for array in (result.fld, x, y, remaining, losses):
        assert np.isfinite(array).all()
