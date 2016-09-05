__author__ = 'Sergey Tomin'

from copy import deepcopy
from ocelot.gui import *
from ocelot import *
from ocelot.rad import *
from time import time

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
U40_short = Undulator(nperiods = nperiods, lperiod=0.040, Kx = 4, eid= "und")

seg = (U40_short,)*int(Nperiods*Nunduls/nperiods)

lat = MagneticLattice(seg)


screen_no = calculate_radiation(lat, deepcopy(screen), beam, energy_loss = False)
screen = calculate_radiation(lat, deepcopy(screen), beam, energy_loss = True)

E_no = screen_no.Eph
t_no = screen_no.Total
max_I = max(t_no)
print("time = ", time() - start)


plt.plot(E_no, t_no/max_I,"k",  screen.Eph, screen.Total/max_I, "r", lw = 2)
plt.grid(True)
plt.legend(["No e.loss", "e.loss", "no e.loss"], loc = 2)
plt.xlabel(r"$E_{ph}$")
#plt.ylabel("Flux, ph/s/mm^2/(0.1%BW)")
plt.ylabel("Normalized intens")
plt.show()