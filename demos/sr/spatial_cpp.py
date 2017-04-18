__author__ = 'Sergey Tomin'


from ocelot.lib.genera.src.python.radiation import generaSR
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
#beam.emit_x = 90e-9
#beam.emit_y = 1e-9

und = Undulator (Kx = 0.43, nperiods = 500, lperiod=0.007, eid= "und")

lat = MagneticLattice((und))

screen = Screen()
screen.z = 100.0
screen.size_x = 0.002
screen.size_y = 0.0
screen.nx = 100
screen.ny = 1


screen.start_energy = 7761.2 #eV
screen.end_energy = 7900 #eV
screen.num_energy = 1

traj, screen = generaSR.calculateSR_py(lat, beam, screen, runParameters = None)

show_flux(screen, unit="mrad")



screen = Screen()
screen.z = 100.0
screen.size_x = 0.002 # m
screen.size_y = 0.002 # m
screen.nx = 1
screen.ny = 1


screen.start_energy = 7700 #eV
screen.end_energy = 7800 #eV
screen.num_energy = 100

traj, screen = generaSR.calculateSR_py(lat, beam, screen, runParameters = None)

show_flux(screen, unit="mrad")


# 2D
screen = Screen()
screen.z = 100.0
screen.size_x = 0.002
screen.size_y = 0.002
screen.nx = 51
screen.ny = 51


screen.start_energy = 7761.2 #eV
screen.end_energy = 7900 #eV
screen.num_energy = 1

traj, screen = generaSR.calculateSR_py(lat, beam, screen, runParameters = None)
show_flux(screen, unit="mrad")
