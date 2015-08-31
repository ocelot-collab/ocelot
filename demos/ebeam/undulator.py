__author__ = 'Sergey Tomin'

from ocelot.cpbd.elements import *
from ocelot.gui.accelerator import *
from ocelot.cpbd.beam import *
from ocelot.cpbd.track import *

# TODO: check the result

und = Undulator (Kx=2, nperiods=200, lperiod=0.007, id="und")
und.solver = "sym"
D1 = Drift(l=0.5, id="D1")
Q1 = Quadrupole(l=0.3, k1=3., id="Q1")
Q2 = Quadrupole(l=0.3, k1=-3, id="Q2")

line = (D1, Q1, D1, und, D1, Q2, D1)

beam = Beam()
beam.E = 2.5
beam.I = 0.1

tw0 = Twiss(beam)


lat = MagneticLattice(line)

tws = twiss(lat, tw0, nPoints=100)

plot_opt_func(lat, tws)

plt.show()

p1 = Particle(x=0.001, y=0.002)
p1.E = beam.E
navi = Navigator()
dz = 0.01
P1 = []
for i in range(int(lat.totalLen/dz)):
    step(lat, [p1], dz=dz, navi=navi)
    P1.append(copy(p1))

s = [f.s for f in P1]
x = [f.x for f in P1]
y = [f.y for f in P1]


plt.plot(s, x, "r", label="X")
plt.plot(s, y, "b", label="Y")
# plt.plot(s, [f.px for f in P1], "b", label = "Y")
plt.title("E != 0. Symplectic transfer map")
plt.legend()
plt.xlabel("S, m")
plt.ylabel("X/Y, m")
plt.grid(True)
plt.show()

lat = MagneticLattice(line)
und.solver = "lin"
p1 = Particle(x=0.001, y=0.002)
p1.E = beam.E
navi = Navigator()
dz = 0.01
P1 = []
for i in range(int(lat.totalLen/dz)):
    step(lat, [p1], dz=dz, navi=navi)
    P1.append(copy(p1))

s = [f.s for f in P1]
x = [f.x for f in P1]
y = [f.y for f in P1]


plt.plot(s, x, "r", label="X")
plt.plot(s, y, "b", label="Y")
# plt.plot(s, [f.px for f in P1], "b", label = "Y")
plt.title("E != 0. Linear transfer map")
plt.legend()
plt.xlabel("S, m")
plt.ylabel("X/Y, m")
plt.grid(True)
plt.show()


lat = MagneticLattice(line)


p1 = Particle(x=0.001, y=0.002)
p1.E = 0.0
navi = Navigator()
dz = 0.01
P1 = []
for i in range(int(lat.totalLen/dz)):
    step(lat, [p1], dz=dz, navi=navi)
    P1.append(copy(p1))

s = [f.s for f in P1]
x = [f.x for f in P1]
y = [f.y for f in P1]


plt.plot(s, x, "r", label="X")
plt.plot(s, y, "b", label="Y")
plt.title("E = 0")
plt.legend()
plt.xlabel("S, m")
plt.ylabel("X/Y, m")
plt.grid(True)
plt.show()