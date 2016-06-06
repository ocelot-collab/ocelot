__author__ = 'Sergey Tomin'

from ocelot.rad import *
from ocelot import *
from ocelot.gui import *
import numpy as np
import time
font = {'size'   : 14}
matplotlib.rc('font', **font)
#from scipy.optimize import curve_fit
from ocelot.demos.sr.k_analysis import *
#from ocelot.lib.genera.src.python.radiation import generaSR

font = {'size'   : 14}
matplotlib.rc('font', **font)

beam = Beam()
beam.E = 17.5
beam.I = 0.1


und = Undulator(Kx = 4., nperiods = 125, lperiod=0.04, eid= "und")
lat = MagneticLattice((und))

screen = Screen()
screen.z = 500.0
screen.size_x = 0.
screen.size_y = 0.
screen.nx = 1
screen.ny = 1


screen.start_energy = 7950 #eV
screen.end_energy = 8200 #eV
screen.num_energy = 1000

screen = calculate_radiation(lat, screen, beam)
show_flux(screen, unit="mrad")

# K-mono scan

beam_energy = 17.5  # GeV
b_energy_jit = 1e-4 # dE/E

screen = Screen()
screen.z = 500.0
screen.size_x = 0.01
screen.size_y = 0.01
screen.nx = 51
screen.ny = 51


ds = screen.size_x/screen.nx*screen.size_y/screen.ny

n_scan_points = 30
n_shots = 5
scan_Kmono_energy = np.linspace(start=8000, stop=8150, num=n_scan_points)


start = time.time()
flux = []
Etotal = []

for n, eph in enumerate(scan_Kmono_energy):
    print(n, "/", n_scan_points)
    for i in range(n_shots):
        beam.E = np.random.normal(beam_energy, beam_energy*b_energy_jit, 1)
        print("beam energy: ", beam.E)
        screen.start_energy = eph # 8078.2 - 50 + i*100/30. #eV
        screen.num_energy = 1
        screen = calculate_radiation(lat, screen, beam)
        flux.append(sum(screen.Total)*ds)
        Etotal.append(eph)

print("time cpp = ", start - time.time())


e_fin, polynom = data_analysis(Etotal, flux=flux, method="least")



print("Eph_fin = ", e_fin)
x = np.linspace(Etotal[0], Etotal[-1], num=100)
plt.plot(Etotal, flux, "r.", lw =2, label="exp data")
plt.plot(x, polynom(x), "b", label="fit func")
plt.plot(e_fin, polynom(e_fin), "go", lw = 3, label=r"$E_{ph}=$" + str(np.around(e_fin, decimals=2)))
plt.xlabel(r"$E_{ph}$, eV")
plt.grid(True)
plt.legend()
plt.show()


