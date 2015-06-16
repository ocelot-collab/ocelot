__author__ = 'Sergey Tomin'

from ocelot.lib.genera.src.python.radiation.radiation_py import *
from ocelot.lib.genera.src.python.radiation.em_screen import *
import matplotlib
font = {'size'   : 14}
matplotlib.rc('font', **font)

beam = Beam()
beam.E = 2.5

beam.I = 0.1

beam.beta_x = 12.84
beam.beta_y = 6.11
beam.Dx = 0.526

und = Undulator (Kx = 0.43, nperiods = 500, lperiod=0.007, id = "und")

lat = MagneticLattice((und), energy = beam.E)

screen = Screen()
screen.z = 100.0
screen.size_x = 0.005
screen.size_y = 0.0
screen.nx = 100
screen.ny = 1


screen.start_energy = 7761.2 #eV
screen.end_energy = 7900 #eV
screen.num_energy = 1
em_screen = EMScreen(screen)

em_screen = calculate_radiation(lat, em_screen, beam )
show_flux(em_screen, unit="mrad")



screen = Screen()
screen.z = 100.0
screen.size_x = 0.002 # m
screen.size_y = 0.002 # m
screen.nx = 1
screen.ny = 1


screen.start_energy = 7700 #eV
screen.end_energy = 7800 #eV
screen.num_energy = 100
em_screen = EMScreen(screen)

em_screen = calculate_radiation(lat, em_screen, beam )
show_flux(em_screen, unit="mrad")