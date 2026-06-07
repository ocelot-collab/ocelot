__author__ = 'Sergey Tomin'
import sys
sys.path.append("../..")

import numpy as np

from ocelot import *
from ocelot.cpbd.beam import generate_parray
from ocelot.cpbd.elements.undulator_atom import und_field
from ocelot.gui.accelerator import *
import copy


d1 = Drift(l=0.1, eid="D1")
d2 = Drift(l=0.4, eid="D2")

# A Drift can be used as a generic hard-edge field region when a user-defined
# magnetic field callback is attached.
field_region = Drift(l=0.2, eid="FIELD_REGION")
field_region.mag_field = lambda x, y, z: (
    0.0 * x,
    0.03 * np.cos(2.0 * np.pi * z / field_region.l),
    0.0 * x,
)

# Bends and multipoles have hard-edge default fields for Runge-Kutta tracking.
bend = Bend(l=0.3, angle=0.01, eid="B_RK")
quad = Quadrupole(l=0.2, k1=1.2, eid="Q_RK")
sext = Sextupole(l=0.15, k2=15.0, eid="SX_RK")
octu = Octupole(l=0.15, k3=150.0, eid="OC_RK")

# Undulators can still use either the analytic default field or a user callback.
und = Undulator(lperiod=0.4, nperiods=9, Kx=44.81, npoints=2000, eid="U_RK")
und.mag_field = lambda x, y, z: und_field(x, y, z, und.lperiod, und.Kx)


cell = (d1, field_region, bend, quad, sext, octu, und, d2)
method = {
    "global": SecondTM,
    Drift: RungeKuttaTM,
    Bend: RungeKuttaTM,
    Quadrupole: RungeKuttaTM,
    Sextupole: RungeKuttaTM,
    Octupole: RungeKuttaTM,
    Undulator: RungeKuttaTM,
}

lat = MagneticLattice(cell, method=method)

np.random.seed(33)

p_array = generate_parray(
    sigma_x=0,
    sigma_px=0,
    sigma_y=None,
    sigma_py=None,
    sigma_tau=100e-6 / 2.36,
    sigma_p=0,
    chirp=0.01 / 2.36,
    charge=0.5e-9,
    nparticles=5000,
    energy=0.6,
)

p_array_1 = copy.deepcopy(p_array)

navi = Navigator(lat)
tws_track, p_array_1 = track(lat, p_array_1, navi, calc_tws=False)

plt.plot(p_array.tau() * 1000, p_array.p(), ".", label="initial distrib.")
plt.plot(p_array_1.tau() * 1000, p_array_1.p(), ".", label="after RK elements")
plt.legend()
plt.grid(True)
plt.show()
