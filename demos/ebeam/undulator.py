__author__ = 'Sergey Tomin'

from ocelot import *
from ocelot.gui import *
import copy

und = Undulator(Kx=2, nperiods=100, lperiod=0.01, eid="und")
#und.solver = "sym"
und.ax = 0.1
D1 = Drift(l=0.5, eid="D1")
Q1 = Quadrupole(l=0.3, k1=3., eid="Q1")
Q2 = Quadrupole(l=0.3, k1=-3, eid="Q2")

line = (D1, Q1, D1, und, D1, Q2, D1)

beam = Beam()
beam.E = 2.5
beam.I = 0.1

tw0 = Twiss(beam)

method = MethodTM()
method.params[Undulator] = UndulatorTestTM
lat = MagneticLattice(line, method=method)

tws = twiss(lat, tw0, nPoints=100)

plot_opt_func(lat, tws)

plt.show()

p1 = Particle(x=0.001, y=0.002)
p1.E = beam.E
navi = Navigator()
dz = 0.01
P1 = []
for i in range(int(lat.totalLen/dz)):
    tracking_step(lat, [p1], dz=dz, navi=navi)
    P1.append(copy.copy(p1))

s = [f.s for f in P1]
x = [f.x for f in P1]
y = [f.y for f in P1]

plt.figure(1)
plt.plot(s, x, "g", label="X sym")
plt.plot(s, y, "r", label="Y sym")
plt.grid(True)

#und.solver = "lin"
lat = MagneticLattice(line)

p1 = Particle(x=0.001, y=0.002)
p1.E = beam.E
navi = Navigator()
dz = 0.01
P1 = []
for i in range(int(lat.totalLen/dz)):
    tracking_step(lat, [p1], dz=dz, navi=navi)
    P1.append(copy.copy(p1))

s = [f.s for f in P1]
x = [f.x for f in P1]
y = [f.y for f in P1]

#plt.figure(2)
plt.plot(s, x, "y", label="X lin")
plt.plot(s, y, "b", label="Y lin")
# plt.plot(s, [f.px for f in P1], "b", label = "Y")
#plt.title("E != 0. Linear transfer map")
plt.legend()
plt.xlabel("S, m")
plt.ylabel("X/Y, m")
plt.grid(True)
plt.show()