import numpy as np
import pytest

from ocelot.cpbd.elements.bend import Bend
from ocelot.cpbd.elements.drift import Drift
from ocelot.cpbd.elements.quadrupole import Quadrupole
from ocelot.cpbd.elements.sextupole import Sextupole
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


@pytest.mark.parametrize(
    "element",
    [
        Drift(l=0.2, dx=1.e-3, dy=-2.e-3),
        Bend(l=0.3, angle=0.01, dx=1.e-3, dy=-2.e-3),
        Quadrupole(l=0.4, k1=2.0, dx=1.e-3, dy=-2.e-3),
    ],
)
def test_element_constructor_accepts_offsets(element):
    assert element.dx == pytest.approx(1.e-3)
    assert element.dy == pytest.approx(-2.e-3)


def test_element_constructor_accepts_mag_field():
    def field(x, y, z):
        return x, y, z

    drift = Drift(l=0.2, mag_field=field)

    assert drift.mag_field is field


@pytest.mark.parametrize("method", [TransferMap, SecondTM])
def test_rolled_quadrupole_tracks_its_displaced_axis(method):
    quad = Quadrupole(l=0.4, k1=2.0, tilt=np.pi / 4, dx=1e-3, tm=method)
    axis = np.array([[quad.dx], [0.], [quad.dy], [0.], [0.], [0.]])
    tracked = quad.tms[0].map_function(axis.copy(), energy=1.0)
    np.testing.assert_allclose(tracked, axis, rtol=0, atol=1e-16)


@pytest.mark.parametrize("family,strengths", [(Quadrupole, {"k1": 2., "k2": 12.}),
                                              (Sextupole, {"k2": 20.})])
@pytest.mark.parametrize("length", [0.4, 0.17])
def test_offset_translates_rolled_second_order_map(family, strengths, length):
    energy = 1.
    displaced = family(l=0.4, tilt=-0.37, dx=1e-3, dy=-2e-3, tm=SecondTM, **strengths)
    centered = family(l=length, tilt=displaced.tilt, tm=SecondTM, **strengths)
    r, t = centered.R(energy)[0], centered.T(energy)[0]
    d = np.array([displaced.dx, 0., displaced.dy, 0., 0., 0.])

    # Expand d + F(X - d) using the centered map in beam coordinates.
    expected_b = d - r @ d + np.einsum("ijk,j,k->i", t, d, d)
    expected_r = r - (t + t.swapaxes(1, 2)) @ d
    params = displaced.get_section_tms(length)[0].get_params(energy)
    np.testing.assert_allclose(params.B[:, 0], expected_b, rtol=0, atol=1e-16)
    np.testing.assert_allclose(params.get_rotated_R(), expected_r, rtol=0, atol=1e-15)
    np.testing.assert_allclose(params.get_rotated_T(), t, rtol=0, atol=1e-15)
