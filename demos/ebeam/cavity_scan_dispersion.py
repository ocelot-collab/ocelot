"""Measure dispersion by changing the RF voltage with the dipole field fixed."""

import matplotlib.pyplot as plt
import numpy as np

import ocelot as ocl
from ocelot.cpbd.response_matrix import LinacCavityScanRM


energy = 0.15  # entrance total energy [GeV]
rf = ocl.Cavity(l=0.5, v=0.05, phi=0.)
bend = ocl.SBend(l=0.4, angle=0.1)
bpms = [ocl.Monitor() for _ in range(6)]
cell = [rf, bpms[0], bend, bpms[1]]
for bpm in bpms[2:]:
    cell.extend([ocl.Drift(l=0.5), bpm])
line = ocl.MagneticLattice(cell, method={"global": ocl.SecondTM})

# A 100 kV on-crest increment raises the energy by 0.05% of its nominal exit value.
step = 0.0001  # GV
exit_energy = energy + rf.v
scan = LinacCavityScanRM(line, [], [], bpms, cavity_weights={rf: 1.},
                         voltage_step=step, normalization_energy=exit_energy)
dx, _ = scan.read_virtual_dispersion(energy)
s = rf.l + np.r_[0., bend.l + np.arange(5) * 0.5]
# At high energy D at a screen after a sector is rho*(1-cos(theta)) + L*sin(theta).
theta = bend.angle
sector_dx = bend.l / theta * (1 - np.cos(theta)) + (s[1:] - s[1]) * np.sin(theta)

fig, ax = plt.subplots(layout="constrained")
ax.plot(s[1:], 1e3 * sector_dx, label="First-order sector prediction")
ax.plot(s, 1e3 * dx, "o", label="RF scan at BPMs")
ax.axvspan(s[0], s[1], color="0.9", label="Sector bend", zorder=-1)
ax.set(xlabel="Distance from entrance [m]", ylabel="Horizontal dispersion [mm]",
       title="An RF voltage change moves the orbit after the bend")
ax.legend()
plt.show()
