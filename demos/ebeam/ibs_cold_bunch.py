"""How much energy spread does a cold, round bunch acquire in a short drift?"""

from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np

import ocelot as ocl
from ocelot.common.globals import m_e_eV, q_e, ro_e
from ocelot.cpbd.physics_proc import IBS

np.random.seed(17)
sigma_x = 100e-6  # round beam at a waist [m]
sigma_tau = 100e-6  # rms bunch duration times c [m]
emit_n = 1e-6  # normalized rms emittance [m rad]
energy = 0.25  # total energy [GeV]
beta_gamma = np.sqrt((energy / ocl.m_e_GeV)**2 - 1)
beam = ocl.generate_parray(sigma_x=sigma_x, sigma_px=emit_n / (beta_gamma * sigma_x),
                           sigma_tau=sigma_tau, sigma_p=0., chirp=0.,
                           energy=energy, nparticles=30000)

# The drift is short compared with the waist beta function (about 5 m).
drift = ocl.Drift(l=0.1)
end = ocl.Marker()
line = ocl.MagneticLattice([drift, end])
charges = np.linspace(0., 250., 9)  # pC
spreads = []
for charge in charges:
    bunch = deepcopy(beam)
    bunch.q_array[:] = charge * 1e-12 / bunch.n
    nav = ocl.Navigator(line, unit_step=drift.l)
    ibs = IBS(method="Nagaitsev", update_Clog=False, Clog=8., bounds=None)
    nav.add_physics_proc(ibs, drift, end)
    ocl.track(line, bunch, nav, calc_tws=False, print_progress=False)
    spreads.append(bunch.p0c * 1e9 * np.std(bunch.p()))

# Cold-beam variance grows in proportion to the number of electrons.
charge_curve = np.linspace(0., charges[-1], 200)
variance_gamma = 8. * ro_e**2 * (charge_curve * 1e-12 / q_e) * drift.l / (
    4 * sigma_x * emit_n * sigma_tau)
fig, ax = plt.subplots(layout="constrained")
ax.plot(charge_curve, m_e_eV * np.sqrt(variance_gamma), label="Cold-beam prediction")
ax.plot(charges, spreads, "o", label="Tracked bunch")
ax.set(xlabel="Bunch charge [pC]", ylabel="Final rms energy spread [eV]",
       title="Longitudinal IBS heating in a 10 cm drift", xlim=(0, 260), ylim=(0, None))
ax.legend()
plt.show()
