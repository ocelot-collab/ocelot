"""A rolled focusing quadrupole couples the two transverse planes of a FODO ring."""

import matplotlib.pyplot as plt
import numpy as np

import ocelot as ocl


drift = ocl.Drift(l=0.7)
qd = ocl.Quadrupole(l=0.25, k1=-1.7)
# emit_x and emit_y give the horizontal-like and vertical-like mode emittances.
# The one-turn map determines the matched optics.
ring_beam = ocl.Twiss(E=1., emit_x=1e-6, emit_y=3e-6)
fig, (ax_beta, ax_corr) = plt.subplots(2, 1, sharex=True, layout="constrained")

for roll, label, style in [(0., "Aligned", "--"), (0.025, "25 mrad roll", "-")]:
    qf = ocl.Quadrupole(l=0.25, k1=1.8, dtilt=roll)
    ring = ocl.MagneticLattice([qf, drift, qd, drift])
    optics = ocl.periodic_twiss(ring, ring_beam, reference_orbit="closed", nPoints=77, return_df=True)

    # Divide the one-turn phase advances by 2π to obtain the two mode tunes.
    q1 = (optics.mux.iloc[-1] - optics.mux.iloc[0]) / (2 * np.pi)
    q2 = (optics.muy.iloc[-1] - optics.muy.iloc[0]) / (2 * np.pi)
    print(f"{label}: Q1 = {q1:.6f}, Q2 = {q2:.6f}")

    # Each projected beta is the position variance divided by that plane's emittance.
    ax_beta.plot(optics.s, optics.beta_x, style, color="C0", label=r"$\beta_x$" if roll else None)
    ax_beta.plot(optics.s, optics.beta_y, style, color="C1", label=r"$\beta_y$" if roll else None)
    correlation = optics.xy / np.sqrt(optics.xx * optics.yy)
    ax_corr.plot(optics.s, correlation, style, color="C2", label=label)

ax_beta.set(ylabel=r"Projected $\beta$ [m]", title="One period of a coupled FODO ring")
ax_corr.set(xlabel="s [m]", ylabel=r"$\langle xy\rangle/(\sigma_x\sigma_y)$")
ax_beta.legend(loc="center", ncol=2)
ax_corr.legend()
plt.show()
