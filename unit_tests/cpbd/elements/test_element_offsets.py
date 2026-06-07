import pytest

from ocelot.cpbd.elements.bend import Bend
from ocelot.cpbd.elements.drift import Drift
from ocelot.cpbd.elements.quadrupole import Quadrupole


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
