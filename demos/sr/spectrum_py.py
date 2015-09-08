__author__ = 'Sergey Tomin'


from ocelot.cpbd.elements import *
from ocelot.lib.genera.src.python.radiation.radiation_py import *
from ocelot.lib.genera.src.python.radiation.em_screen import *
from ocelot.common.screen import *
from ocelot.cpbd.optics import *
import matplotlib
from ocelot.cpbd.beam import *

font = {'size'   : 14}
matplotlib.rc('font', **font)

beam = Beam()
beam.E = 2.5

beam.I = 0.1

beam.beta_x = 12.84
beam.beta_y = 6.11
beam.Dx = 0.526

und = Undulator(Kx = 0.43, nperiods = 280, lperiod=0.007, id = "und")

lat = MagneticLattice((und))

screen = Screen()
screen.z = 100.0
screen.size_x = 0.0
screen.size_y = 0.0
screen.nx = 1
screen.ny = 1


screen.start_energy = 7000 #eV
screen.end_energy = 7900 #eV
screen.num_energy = 1000
em_screen = EMScreen(screen)
#traj, screen = generaSR.calculateSR_py(lat, beam, screen, runParameters = None)
em_screen = calculate_radiation(lat, em_screen, beam )
show_flux(em_screen, unit="mrad")