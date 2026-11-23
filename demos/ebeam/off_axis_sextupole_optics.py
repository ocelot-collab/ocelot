"""A shifted quadrupole steers the orbit into a sextupole, changing its focusing."""

import matplotlib.pyplot as plt

import ocelot as ocl


quad = ocl.Quadrupole(l=0.4, k1=2., dx=1.5e-3)
upstream = ocl.Drift(l=1.)
downstream = ocl.Drift(l=1.2)
tws0 = ocl.Twiss(E=1., beta_x=4., beta_y=3.)  # Twiss parameters at the beamline entrance
fig, (ax_orbit, ax_beta) = plt.subplots(2, 1, sharex=True, layout="constrained")

for k2l, label in [(0., "Sextupole off"), (150., "Sextupole on")]:
    sext = ocl.Multipole(kn=[0., 0., k2l])  # thin sextupole, integrated K2 in m^-2
    line = ocl.MagneticLattice([quad, upstream, sext, downstream])
    # The design map uses the sextupole gradient at its centre, where it is zero.
    design = ocl.twiss(line, tws0, nPoints=131, return_df=True)
    # Track the on-axis entrance particle and propagate deviations around it.
    tracked = ocl.twiss(line, tws0, reference_orbit=ocl.Particle(), nPoints=131, return_df=True)
    ax_orbit.plot(tracked.s, 1e3 * tracked.x, label=label)
    ax_beta.plot(tracked.s, tracked.beta_x - design.beta_x, label=label)

# The quadrupole offset steers both orbits equally before the sextupole.
# There the local horizontal kick gradient is -K2L * x, so beta changes downstream.
for ax in (ax_orbit, ax_beta):
    ax.axvline(quad.l + upstream.l, color="0.4", linestyle=":", label="Thin sextupole")
ax_orbit.set(ylabel="x [mm]", title="Sextupole feed-down on a steered orbit")
ax_beta.set(xlabel="s [m]", ylabel=r"$\beta_x$ (tracked $-$ design) [m]")
ax_orbit.legend()
plt.show()
