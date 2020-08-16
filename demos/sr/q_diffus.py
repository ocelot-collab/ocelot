__author__ = 'Sergey Tomin'

from pylab import *
from ocelot.rad import *
from ocelot import *
from copy import deepcopy
from time import time
import scipy.stats as stats
font = {'size': 20}
matplotlib.rc('font', **font)


beam = Beam()
beam.E = 17.5
beam.I = 0.1    #A

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

nund = 1               # Number of undulators
nperiods = 125          # Period numbers per undulator
kp = 1                  # applying energy diffusion "kick" every N period
npls = 150              # particles number

U40_short = Undulator(nperiods=kp, lperiod=0.040, Kx=4, eid="und")

seg = (U40_short,)*int(nperiods*nund/kp)

lat = MagneticLattice(seg)


screen_no = calculate_radiation(lat, deepcopy(screen), beam, energy_loss=False, quantum_diff=False)
total = np.zeros(screen.num_energy)
Uq = []
for i in range(npls):
    print(i)
    screen = calculate_radiation(lat, deepcopy(screen), beam, energy_loss=False, quantum_diff=True)
    total += screen.Total/npls
    Uq.append(screen.Ef_electron)

sigma = sigma_gamma_quat(beam.E, U40_short.Kx, U40_short.lperiod, U40_short.lperiod*nund*nperiods)


num_bins = 50
# the histogram of the data
n, bins, patches = plt.hist((np.array(Uq) - beam.E)/beam.E, num_bins, normed=1, facecolor='green', alpha=0.5)

print("sigma = ", sigma)
print(sigma/sqrt(nund*nperiods))
y = stats.norm.pdf(bins, 0, sigma)
plt.plot(bins, y, 'r--', lw=2)
plt.xlabel(r'$\delta E/ E $')
#plt.ylabel('particle number')
#plt.title(r'Histogram of IQ: $\mu=100$, $\sigma=$')

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.show()



E_no = screen_no.Eph
t_no = screen_no.Total
max_I = max(t_no)
print("time = ", time() - start)


plt.plot(E_no, t_no/max_I, lw=2, label="no quantum fluct.")
plt.plot(screen.Eph, total/max_I, lw=2, label="quantum fluct.")
plt.grid(True)
legend()
xlabel(r"$E_{ph}$")
#ylabel("Flux, ph/s/mm^2/(0.1%BW)")
ylabel("Normalized intens")
show()