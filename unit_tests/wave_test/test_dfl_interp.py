import numpy as np
import pytest

from ocelot.optics.wave import RadiationField, dfl_interp


@pytest.mark.parametrize("method, degree", [("linear", 1), ("cubic", 3), ("quintic", 5)])
@pytest.mark.parametrize("new_length", [(4.0, 2.0), (6.0, 2.0), (4.0, 3.0), (6.0, 3.0)])
def test_dfl_interp_complex_polynomial_and_zero_fill(method, degree, new_length, monkeypatch):
    import scipy.interpolate

    # Also exercise the migration on SciPy versions that still expose interp2d.
    monkeypatch.delattr(scipy.interpolate, "interp2d", raising=False)

    def field(x, y):
        return 2 + x + 0.2 * y ** degree + 1j * (1 - 0.3 * y + 0.05 * x ** degree)

    dfl = RadiationField((2, 7, 9))
    dfl.dx, dfl.dy, dfl.dz = 4.0 / 9, 2.0 / 7, 1e-6
    dfl.xlamds = 1e-9
    x, y = np.meshgrid(np.linspace(-2, 2, 9), np.linspace(-1, 1, 7))
    dfl.fld[0] = field(x, y)
    dfl.fld[1] = (2 - 1j) * dfl.fld[0]
    original = dfl.fld.copy()

    result = dfl_interp(dfl, newN=(13, 11), newL=new_length, method=method)

    xnew, ynew = np.meshgrid(
        np.linspace(-new_length[0] / 2, new_length[0] / 2, 13),
        np.linspace(-new_length[1] / 2, new_length[1] / 2, 11),
    )
    outside = (abs(xnew) > 2) | (abs(ynew) > 1)
    expected = field(xnew, ynew)
    expected[outside] = 0
    expected /= np.linalg.norm(expected)
    for index, phase in enumerate([1, (2 - 1j) / abs(2 - 1j)]):
        np.testing.assert_allclose(
            result.fld[index] / np.linalg.norm(result.fld[index]), phase * expected,
            rtol=1e-12, atol=1e-14,
        )
    np.testing.assert_array_equal(result.fld[:, outside], 0)
    # Preserve the legacy normalization, which excludes the final grid row/column.
    np.testing.assert_allclose(
        np.sum(abs(result.fld) ** 2, axis=(1, 2)),
        np.sum(abs(original[:, :-1, :-1]) ** 2, axis=(1, 2)),
        rtol=1e-12,
    )
    assert result.fld.shape == (2, 11, 13)
    assert result.dx == pytest.approx(new_length[0] / 13)
    assert result.dy == pytest.approx(new_length[1] / 11)
    assert result.dz == dfl.dz
    assert result.xlamds == dfl.xlamds
    np.testing.assert_array_equal(dfl.fld, original)


def test_dfl_interp_cropped_grid_uses_existing_power_normalization():
    dfl = RadiationField((1, 7, 9))
    dfl.dx, dfl.dy, dfl.dz = 4.0 / 9, 2.0 / 7, 1e-6
    dfl.fld[:] = 1 + 2j

    result = dfl_interp(dfl, newN=(5, 3), newL=(2.0, 1.0))

    # The old crop selects four x samples and two y samples for normalization.
    np.testing.assert_allclose(result.fld, (1 + 2j) * np.sqrt(8 / 15))


def test_dfl_interp_zero_slice_stays_zero():
    dfl = RadiationField((2, 7, 9))
    dfl.dx = dfl.dy = dfl.dz = 1e-6
    dfl.fld[0] = 1 + 1j

    result = dfl_interp(dfl, newN=(13, 11))

    np.testing.assert_array_equal(result.fld[1], 0)
    assert np.isfinite(result.fld).all()


def test_dfl_interp_rejects_unsupported_method():
    dfl = RadiationField((1, 7, 9))
    dfl.dx = dfl.dy = dfl.dz = 1e-6
    dfl.fld[:] = 1

    with pytest.raises(ValueError, match="linear.*cubic.*quintic"):
        dfl_interp(dfl, newN=(13, 11), method="nearest")
