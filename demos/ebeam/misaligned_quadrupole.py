"""Track an on-axis particle through a transversely displaced quadrupole."""

import matplotlib.pyplot as plt
import numpy as np

import ocelot as ocl


upstream = ocl.Drift(l=0.3)
downstream = ocl.Drift(l=1.5)
fig, (ax_x, ax_y) = plt.subplots(2, 1, sharex=True, layout="constrained")

for dx, dy, label in [(0., 0., "Aligned"), (0.5e-3, 0.3e-3, "Displaced")]:
    quad = ocl.Quadrupole(l=0.4, k1=2., dx=dx, dy=dy)  # offsets in metres
    line = ocl.MagneticLattice([upstream, quad, downstream])
    particle = ocl.Particle(E=1.)  # 1 GeV, x = px = y = py = 0
    navigator = ocl.Navigator(line)

    s = np.linspace(0., line.totalLen, 111)
    x, y = [particle.x], [particle.y]
    for ds in np.diff(s):
        ocl.tracking_step(line, [particle], ds, navigator)
        x.append(particle.x)
        y.append(particle.y)

    ax_x.plot(s, 1e3 * np.array(x), label=label)
    ax_y.plot(s, 1e3 * np.array(y))

# With positive k1 the particle bends towards the shifted axis in x and away in y.
for ax in (ax_x, ax_y):
    ax.axvspan(upstream.l, upstream.l + quad.l, color="0.85", label="Quadrupole", zorder=0)
ax_x.set(ylabel="x [mm]", title="Steering by a displaced quadrupole")
ax_y.set(xlabel="s [m]", ylabel="y [mm]")
ax_x.legend()
plt.show()
