"""Compare transverse RF focusing at zero net acceleration."""

import matplotlib.pyplot as plt
import numpy as np

import ocelot as ocl
from ocelot.cpbd.optics import trace_z


# Both structures have zero reference energy gain at 90 degrees. Alternating
# fields in the standing-wave structure still focus an off-axis particle.
length = 0.8  # m
fig, ax = plt.subplots(layout="constrained")
for cavity_type, label in [(ocl.Cavity, "Standing wave"), (ocl.TWCavity, "Travelling wave")]:
    cavity = cavity_type(l=length, v=0.030, phi=90., freq=1.3e9)
    line = ocl.MagneticLattice([cavity, ocl.Drift(l=3.)])
    particle = ocl.Particle(E=0.050, x=1e-3)  # 50 MeV, parallel entrance
    s = np.linspace(0., line.totalLen, 77)
    trajectory = trace_z(line, particle, s)
    ax.plot(s, [1e3 * p.x for p in trajectory], label=label)

ax.axvspan(0., length, color="0.9", zorder=0, label="RF structure")
ax.set(xlabel="Distance along the beamline [m]", ylabel="Horizontal position [mm]",
       title="RF focusing at the zero crossing")
ax.legend()
plt.show()
