__author__ = 'Sergey Tomin'
"""
can read different types of files. By default, mag_file is in [mm] .
for python version, only vertical component of magnetic field (By) is taken into account.
In order to overcome this limitation, someone have to change function radiation_py.field_map2field_func(z, By).
Sergey Tomin 04.11.2016.
"""
from ocelot.rad import *
from ocelot import *
from ocelot.gui import *

font = {'size'   : 14}
matplotlib.rc('font', **font)

beam = Beam()
beam.E = 1.25

beam.I = 0.1

beam.beta_x = 12.84
beam.beta_y = 6.11
beam.Dx = 0.526

und = Undulator(field_file="mag_file.txt", eid="und")

lat = MagneticLattice((und))

screen = Screen()
screen.z = 50.0
screen.size_x = 0.0
screen.size_y = 0.0
screen.nx = 1
screen.ny = 1


screen.start_energy = 0.0001 #eV
screen.end_energy = 0.1 #eV
screen.num_energy = 1000

screen = calculate_radiation(lat, screen, beam, accuracy=2)
show_flux(screen, unit="mrad")