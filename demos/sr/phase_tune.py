__author__ = 'Sergey Tomin'

from ocelot.rad import *
from ocelot import *
from ocelot.gui import *

font = {'size'   : 14}
matplotlib.rc('font', **font)

beam = Beam()
beam.E = 17.5

beam.I = 0.1

beam.beta_x = 12.84
beam.beta_y = 6.11
beam.Dx = 0.526

und = Undulator(Kx = 4., nperiods=125, lperiod=0.04, eid= "und")
D = Drift(l=0.5, eid="D")
b1 = Hcor(l=0.1, angle = 5*-0.00001, eid="b1")
b2 = Hcor(l=0.2, angle = 5*0.00002, eid="b2")
b3 = Hcor(l=0.1, angle = 5*-0.00001, eid="b3")
phase_shift =  (b1, b2, b3)
cell = (und, D, phase_shift, D, und)
lat = MagneticLattice(cell)

screen = Screen()
screen.z = 100.0
screen.size_x = 0.0
screen.size_y = 0.0
screen.nx = 1
screen.ny = 1


screen.start_energy = 7900 #eV
screen.end_energy = 8200 #eV
screen.num_energy = 1000

print_rad_props(beam, K=und.Kx, lu=und.lperiod, L=und.l, distance=screen.z)


screen = calculate_radiation(lat, screen, beam)

# trajectory
for u in screen.motion:
    plt.plot(u[4::9], u[0::9], "r")
plt.show()

show_flux(screen, unit="mrad")

und = Undulator(Kx = 4., nperiods=125, lperiod=0.04, eid= "und")
D = Drift(l=0.5, eid="D")
b1 = Hcor(l=0.1, angle = 10*-0.00001, eid="b1")
b2 = Hcor(l=0.2, angle = 10*0.00002, eid="b2")
b3 = Hcor(l=0.1, angle = 10*-0.00001, eid="b3")
phase_shift =  (b1, b2, b3)
cell = (und, D, phase_shift, D, und)
lat = MagneticLattice(cell)

screen = calculate_radiation(lat, screen, beam)

# trajectory
for u in screen.motion:
    plt.plot(u[4::9], u[0::9], "r")
plt.show()

show_flux(screen, unit="mrad")