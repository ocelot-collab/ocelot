__author__ = 'Sergey Tomin'


from ocelot.cpbd.elements import *
from ocelot.lib.genera.src.python.radiation import generaSR
from ocelot.common.screen import *
from ocelot.cpbd.optics import *
import matplotlib
from  ocelot.lib.genera.src.python.radiation.em_screen import show_flux
from ocelot.cpbd.beam import *
font = {'size'   : 14}
matplotlib.rc('font', **font)

beam = Beam()
beam.E = 2.5

beam.I = 0.1

beam.beta_x = 12.84
beam.beta_y = 6.11
beam.Dx = 0.526

und = Undulator (Kx = 0.49, nperiods = 280, lperiod=0.007, id="und")

lat = MagneticLattice((und))

screen = Screen()
screen.z = 100.0
screen.size_x = 0.0
screen.size_y = 0.0
screen.nx = 1
screen.ny = 1


screen.start_energy = 7400 #eV
screen.end_energy = 7900 #eV
screen.num_energy = 1000

traj, screen = generaSR.calculateSR_py(lat, beam, screen, runParameters = None)

show_flux(screen, unit="mrad")