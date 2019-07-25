__author__ = 'Sergey Tomin'

import matplotlib
from ocelot import *
from ocelot.rad import *
from ocelot.gui import *
import numpy as np


chirp_coeff = 0.01/2.36
sigma_tau = 100e-6/2.36

tau = np.array([-0.7, -0.3, 0, 0.3, 0.7])*100e-6/2.36*0

p_array = ParticleArray(n=5)
p_array.E = 0.6
p_array.rparticles[4, :] = tau[:]
p_array.rparticles[5, :] = -chirp_coeff*tau/sigma_tau
p_array.q_array[:] = 1e-10

und = Undulator(lperiod=0.4, nperiods=9, Kx=44.81)
lat = MagneticLattice( ( und, ))


screen = Screen()
screen.z = 1000.0
screen.size_x = 15
screen.size_y = 15
screen.nx = 1001
screen.ny = 1
screen.start_energy = 0.00850446  # eV
screen.end_energy = 15e-3  # eV
screen.num_energy = 1


screen_i = coherent_radiation(lat, screen, p_array, accuracy=1)

show_flux(screen_i, unit="mm", title="R56 != 0")
plt.show()
