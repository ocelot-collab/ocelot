import pytest
import numpy as np

from ocelot.cpbd.elements import Bend, Drift, Octupole, Quadrupole, Sextupole, Solenoid, Undulator
from ocelot.cpbd.beam import Particle
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaGlobalTM, RungeKuttaOcelotTM


@pytest.mark.parametrize(
    "element,family",
    [
        (Drift(l=0.2, npoints=37), Drift),
        (Quadrupole(l=0.2, k1=1.0, npoints=37), Quadrupole),
        (Sextupole(l=0.2, k2=2.0, npoints=37), Sextupole),
        (Octupole(l=0.2, k3=3.0, npoints=37), Octupole),
        (Solenoid(l=0.2, k=0.1, npoints=37), Solenoid),
        (Undulator(lperiod=0.02, nperiods=5, Kx=1.0, npoints=37), Undulator),
    ],
)
def test_runge_kutta_npoints_constructor_param_used_by_lattice_method(element, family):
    MagneticLattice((element,), method={family: RungeKuttaOcelotTM})

    assert all(isinstance(tm, RungeKuttaOcelotTM) for tm in element.tms)
    assert all(tm.npoints == 37 for tm in element.tms)


def test_runge_kutta_npoints_constructor_param_used_for_edge_aware_bend():
    bend = Bend(l=0.3, angle=0.01, npoints=41)

    MagneticLattice((bend,), method={Bend: RungeKuttaGlobalTM})

    assert all(isinstance(tm, RungeKuttaGlobalTM) for tm in bend.tms)
    assert len(bend.tms) == 3
    assert all(tm.npoints == 41 for tm in bend.tms)


def test_runge_kutta_npoints_constructor_param_can_be_overridden_by_set_tm():
    quad = Quadrupole(l=0.2, k1=1.0, npoints=37)

    quad.set_tm(RungeKuttaOcelotTM, npoints=53)

    assert quad.tms[0].npoints == 53


def test_runge_kutta_section_maps_preserve_npoints_and_start_position():
    bend = Bend(l=1.0, angle=0.01, npoints=73)

    bend.set_tm(RungeKuttaOcelotTM)
    section = bend.get_section_tms(start_l=0.25, delta_l=0.2)

    assert len(section) == 1
    assert isinstance(section[0], RungeKuttaOcelotTM)
    assert section[0].npoints == 73
    assert section[0].s_start == pytest.approx(0.25)


def test_runge_kutta_section_map_tracks_absolute_field_interval():
    z_seen = []

    def field(x, y, z):
        z_seen.extend(np.ravel(np.asarray(z, dtype=float)).tolist())
        return 0.0 * x, 0.0 * x, 0.0 * x

    bend = Bend(l=1.0, angle=0.0, npoints=20, mag_field=field)
    bend.set_tm(RungeKuttaGlobalTM)

    section = bend.get_section_tms(start_l=0.25, delta_l=0.2)
    section[0].apply(Particle(E=1.0))

    assert min(z_seen) >= 0.25 - 1.0e-12
    assert max(z_seen) <= 0.45 + 1.0e-12
