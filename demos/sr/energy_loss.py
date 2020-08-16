__author__ = 'Sergey Tomin'

from copy import deepcopy
from ocelot.gui import *
from ocelot import *
from ocelot.rad import *
from ocelot.rad.radiation_py import *
from time import time
from ocelot.rad.undulator_params import *
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
nperiods = 5

U40_short = Undulator(nperiods=nperiods, lperiod=0.040, Kx=4, eid="und")
seq = (U40_short,)*int(Nperiods*Nunduls/nperiods)

lat = MagneticLattice(seq)
E = energy_loss_und(energy=beam.E, Kx=4, lperiod=0.04, L=0.04*Nperiods, energy_loss=True)

print("energy loss per undulator = ", E*1e3, "MeV")

screen_no = calculate_radiation(lat, deepcopy(screen), beam, energy_loss=False)
screen_loss = calculate_radiation(lat, deepcopy(screen), beam, energy_loss=True)
E_no = screen_no.Eph
t_no = screen_no.Total
max_I = max(t_no)
print("time = ", time() - start)


plt.plot(E_no, t_no/max_I, lw=2, label="no e.loss")
plt.plot(screen_loss.Eph, screen_loss.Total/max_I, lw=2, label="e.loss")


# tapering
seq = []
energy = beam.E
K = 4
for i in range(Nunduls):
    U40 = Undulator(nperiods=125, lperiod=0.040, Kx=K, eid="und")
    seq.append(U40)
    E_loss = energy_loss_und(energy=energy, Kx=K, lperiod=U40.lperiod, L=U40.lperiod * Nperiods, energy_loss=True)
    K = np.sqrt(2 * (((energy - E_loss) / beam.E) ** 2 * (1 + 4 ** 2 / 2) - 1))
    energy -= E_loss



lat = MagneticLattice(seq)

screen_comp = calculate_radiation(lat, deepcopy(screen), beam, energy_loss=True)
plt.plot(screen_comp.Eph, screen_comp.Total/max_I, "--", lw=2, label="tapering")

plt.grid(True)
plt.legend()
plt.xlabel(r"$E_{ph}$")
#plt.ylabel("Flux, ph/s/mm^2/(0.1%BW)")
plt.ylabel("Normalized intens")
plt.show()