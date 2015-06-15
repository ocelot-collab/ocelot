__author__ = 'Sergey Tomin'

from ocelot.cpbd.optics import *
from ocelot.cpbd.elements import *
from ocelot.gui.accelerator import *
import pylab as plt

und = Undulator (Kx = 0.49, nperiods=200, lperiod=0.007, id = "und")
D1 = Drift (l = 0.5, id = "D1")
Q1 = Quadrupole (l = 0.3, k1 = 3., id = "Q1")
Q2 = Quadrupole (l = 0.3, k1 = -3, id = "Q2")

line  = (D1, Q1, D1, und, D1, Q2, D1)

beam = Beam()
beam.E = 2.5
beam.I = 0.1

tw0 = Twiss(beam)


lat = MagneticLattice(line, energy = 2.5)

tws = twiss(lat, tw0, nPoints = 1000)

plot_opt_func(lat, tws)

plt.show()

p1 = Particle(x = 0.001, y = 0.001)

navi = Navigator()
dz = 0.01
P1 = []
for i in range(int(lat.totalLen/dz)):
    track(lat, [p1], dz = dz, navi = navi)
    P1.append(copy(p1))

s = [f.s for f in P1]
x = [f.x for f in P1]
y = [f.y for f in P1]


plt.plot(s, x, "r", label = "X")
plt.plot(s, y, "b", label = "Y")
#plt.plot(s, [f.px for f in P1], "b", label = "Y")
plt.title("Energy is not zero. Nonlinear transfer map")
plt.legend()
plt.xlabel("S, m")
plt.ylabel("X/Y, m")
plt.grid(True)
plt.show()


lat = MagneticLattice(line, energy = 0)


p1 = Particle(x = 0.001, y = 0.001)

navi = Navigator()
dz = 0.01
P1 = []
for i in range(int(lat.totalLen/dz)):
    track(lat, [p1], dz = dz, navi = navi)
    P1.append(copy(p1))

s = [f.s for f in P1]
x = [f.x for f in P1]
y = [f.y for f in P1]


plt.plot(s, x, "r", label = "X")
plt.plot(s, y, "b", label = "Y")
plt.title("Beam energy is zero")
plt.legend()
plt.xlabel("S, m")
plt.ylabel("X/Y, m")
plt.grid(True)
plt.show()




