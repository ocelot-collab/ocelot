"""Interference from two transversely separated electrons in one undulator."""

import matplotlib.pyplot as plt

import ocelot as ocl
from ocelot.common.globals import h_eV_s, q_e
from ocelot.rad import Screen, coherent_radiation


energy = 0.25  # GeV
undulator = ocl.Undulator(lperiod=0.03, nperiods=12, Kx=0.6, npoints=301)
line = ocl.MagneticLattice([undulator])
gamma = energy / ocl.m_e_GeV
wavelength = undulator.lperiod * (1 + undulator.Kx**2 / 2) / (2 * gamma**2)

screen = Screen()
screen.z, screen.size_x, screen.nx = 1000., 0.4, 161  # m, half-width [m], pixels
screen.ny = screen.num_energy = 1
screen.start_energy = screen.end_energy = h_eV_s * ocl.speed_of_light / wavelength

fig, ax = plt.subplots(layout="constrained")
for spacing, label in [(0., "Coincident"), (0.4e-3, "Separated by 0.4 mm")]:
    pair = ocl.ParticleArray(n=2)
    pair.E = energy
    pair.q_array[:] = q_e
    pair.x()[:] = [-spacing / 2, spacing / 2]  # simultaneous, parallel electrons
    coherent_radiation(line, screen, pair, verbose=False)
    # Screen.Xph is in mm. Dividing by the distance in m gives mrad.
    ax.plot(screen.Xph / screen.z, screen.Total, label=label)

# The path difference in the far field is spacing * theta. The first
# destructive interference appears at theta = wavelength / (2 * spacing).
ax.set(xlabel="Horizontal observation angle [mrad]",
       ylabel="Photons per bunch / mm² / 0.1% BW",
       title="Two-electron interference at the undulator fundamental")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncols=2)
plt.show()
