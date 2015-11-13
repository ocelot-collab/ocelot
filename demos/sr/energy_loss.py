__author__ = 'Sergey Tomin'

from time import time
from pylab import *
from ocelot.rad.radiation_py import *
from ocelot.rad.screen import *
from ocelot.cpbd.elements import *
#from ocelot.cpbd.optics import *
#from ocelot.cpbd.e_beam_params import *
from ocelot.cpbd.beam import *
from copy import deepcopy

font = {'size'   : 20}
matplotlib.rc('font', **font)


beam = Beam()
beam.E = 17.5
beam.I = 0.1 #A
screen = Screen()
screen.z = 5000.0
screen.size_x = 0.00
screen.size_y = 0.00
screen.nx = 1
screen.ny = 1

screen.start_energy = 8030 #eV
screen.end_energy = 8090 #eV
screen.num_energy = 500

start = time()

Nunduls = 15
Nperiods = 125
nperiods = 1
U40_short = Undulator(nperiods = nperiods, lperiod=0.040, Kx = 4, id = "und")

seg = (U40_short,)*int(Nperiods*Nunduls/nperiods)

lat = MagneticLattice(seg)


screen_no = calculate_radiation(lat, deepcopy(screen), beam, energy_loss = False)
screen = calculate_radiation(lat, deepcopy(screen), beam, energy_loss = True)

E_no = screen_no.Eph
t_no = screen_no.Total
max_I = max(t_no)
print "time = ", time() - start


plot(E_no, t_no/max_I,"k",  screen.Eph, screen.Total/max_I, "r", lw = 2)
grid(True)
legend(["No e.loss", "e.loss", "no e.loss"], loc = 2)
xlabel(r"$E_{ph}$")
#ylabel("Flux, ph/s/mm^2/(0.1%BW)")
ylabel("Normalized intens")
show()