__author__ = 'Sergey Tomin'
import sys
import matplotlib
from ocelot.rad import *

from ocelot.gui import *
from ocelot import *

from ocelot.rad.radiation_py import *
from ocelot.rad.undulator_params import *
import copy


sigma_tau = 100e-6/2.36
tau_p_cor = 0.013/2.36
tau = np.array([-1, 0, 1])*sigma_tau
phi = tau/1.45859E-04*360



font = {'size'   : 14}
matplotlib.rc('font', **font)

p_array_init = ParticleArray(n=3)
p_array_init.tau()[:] = np.array([-1, 0, 1])*sigma_tau
p_array_init.p()[:] = tau_p_cor*tau/sigma_tau
p_array_init.E = 0.6
p_array_init.q_array[:] = 1e-10


p_array = copy.deepcopy(p_array_init)

screen = Screen()
screen.z = 1000.0
screen.size_x = 15
screen.size_y = 15
screen.nx = 2000
screen.ny = 1
screen.start_energy = 0.00850
screen.end_energy = 15e-3
screen.num_energy = 1
screen.update()

und = Undulator(lperiod=0.4, nperiods=9, Kx=44.821)
lat = MagneticLattice((und,))

screen = coherent_radiation(lat, screen, p_array, accuracy=2, end_poles=False)

show_flux(screen)
plt.show()


