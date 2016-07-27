__author__ = 'Sergey Tomin'

from ocelot.rad import *
from ocelot import *
from ocelot.gui import *

font = {'size'   : 14}
matplotlib.rc('font', **font)

beam = Beam()
beam.E = 2.5

beam.I = 0.1

beam.beta_x = 12.84
beam.beta_y = 6.11
beam.Dx = 0.526

und = Undulator(Kx = 0.43, nperiods = 280, lperiod=0.007, eid= "und")

lat = MagneticLattice((und))

screen = Screen()
screen.z = 100.0
screen.size_x = 0.0
screen.size_y = 0.0
screen.nx = 1
screen.ny = 1


screen.start_energy = 7300 #eV
screen.end_energy = 7900 #eV
screen.num_energy = 1000

#traj, screen = generaSR.calculateSR_py(lat, beam, screen, runParameters = None)
screen = calculate_radiation(lat, screen, beam )
show_flux(screen, unit="mrad")