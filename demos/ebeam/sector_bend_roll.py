"""Compare a design roll and an installation roll of a sector dipole."""

import matplotlib.pyplot as plt
import numpy as np

import ocelot as ocl
from ocelot.cpbd.optics import trace_z


roll = 3e-3  # rad
angle = 0.15  # rad of horizontal design deflection
length = 0.6  # m
drift = ocl.Drift(l=1.4)
s = np.linspace(0., length + drift.l, 101)
fig, ax = plt.subplots(layout="constrained")

for tilt, dtilt, label in [(roll, 0., "Design roll"), (0., roll, "Installation roll")]:
    bend = ocl.SBend(l=length, angle=angle, tilt=tilt, dtilt=dtilt)
    line = ocl.MagneticLattice([bend, drift])
    orbit = trace_z(line, ocl.Particle(E=1.), s)
    # Each sample is measured from the design orbit of its own lattice.
    ax.plot(s, [1e3 * p.y for p in orbit], label=label)

# Rotating the design bend also rotates its surveyed orbit. An installation
# roll keeps the survey fixed and gives py/p0 = -sin(angle)*sin(roll) at exit.
ax.axvspan(0., length, color="0.9", zorder=0, label="Sector dipole")
ax.set(xlabel="Distance along the design orbit [m]", ylabel="Vertical orbit offset [mm]",
       title="A 3 mrad roll of a sector dipole")
ax.legend()
plt.show()
