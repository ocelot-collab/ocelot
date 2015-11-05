__author__ = 'Sergey Tomin'


from ocelot.cpbd.elements import *
from ocelot.lib.genera.src.python.radiation.radiation_py import *
from ocelot.lib.genera.src.python.radiation.em_screen import *
from ocelot.common.screen import *
from ocelot.cpbd.optics import *
from ocelot.lib.genera.unit2unit import *
import matplotlib
from ocelot.cpbd.beam import *

font = {'size'   : 14}
matplotlib.rc('font', **font)

beam = Beam()
beam.E = 17.5

beam.I = 0.1

beam.beta_x = 12.84
beam.beta_y = 6.11
beam.Dx = 0.526

und = Undulator(Kx = 4., nperiods=125, lperiod=0.04, id = "und")
D = Drift(l=0.5, id="D")
b1 = Hcor(l=0.1, angle = 5*-0.00001, id="b1")
b2 = Hcor(l=0.2, angle = 5*0.00002, id="b2")
b3 = Hcor(l=0.1, angle = 5*-0.00001, id="b3")
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

em_screen = EMScreen(screen)
#traj, screen = generaSR.calculateSR_py(lat, beam, screen, runParameters = None)
em_screen = calculate_radiation(lat, em_screen, beam)

# trajectory
for u in em_screen.motion:
    plt.plot(u[4::9], u[0::9], "r")
plt.show()

show_flux(em_screen, unit="mrad")