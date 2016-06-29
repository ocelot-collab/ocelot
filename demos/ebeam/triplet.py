__author__ = 'Sergey Tomin'

from ocelot import *
from ocelot.gui import *
import copy
import numpy as np


Q1 = Quadrupole(l = 0.3, k1 = 5)
Q2 = Quadrupole(l = 0.3, k1 = -5)
D = Drift(l = 0.5)
De = Drift(l = 0.5)
lattice = [D, Q1, D, Q2, D, Q1, De]

method1 = MethodTM()
method1.global_method = TransferMap
lat1 = MagneticLattice(copy.deepcopy(lattice), method=method1)

method2 = MethodTM()
method2.global_method = SecondTM
lat2 = MagneticLattice(copy.deepcopy(lattice), method=method2)

tw0 = Twiss()
tw0.beta_x = 5.
tw0.alpha_x = -0.87
tw0.beta_y = 2.1
tw0.alpha_y = 0.96

tws = twiss(lat1, tw0, nPoints=None)

plot_opt_func(lat1, tws)
plt.show()

# track ellipse
t = np.linspace(0, 2*pi, num = 100)
x, px = 0.1*np.cos(t), 0.1*np.sin(t)
plist = []
for xi, pxi in zip(x, px):
    plist.append(Particle(x = xi, px= pxi))

plist_1 = copy.deepcopy(plist)

navi = Navigator()
dz = lat1.totalLen
tracking_step(lat2, plist, dz = dz, navi = navi)  # R + T
#for e in lat1.sequence:
#    print("check", e.transfer_map.__class__)
tracking_step(lat1, plist_1, dz = dz, navi = navi)  # R only

tau2 = [f.tau for f in plist]
x2 = [f.x for f in plist]

tau1 = [f.tau for f in plist_1]
x1 = [f.x for f in plist_1]

plt.suptitle("Tracking w/o 2d order matrices")
plt.subplot(121)
plt.title("$(x, x')$ plane. S = 0 m")
plt.plot(x, px, "r.-", label = "X")
plt.xlabel("X, m")
plt.ylabel("Xp, rad")
plt.grid(True)

plt.subplot(122)
plt.title("$(\Delta l, \Delta p/p)$ plane. End of lattice")
plt.plot(x2, tau2, "r.-", label="R+T")
plt.plot(x1, tau1, "b.-", label="R")
plt.legend()
plt.xlabel("X, m")
plt.ylabel("$\delta l$, m")
plt.grid(True)
plt.show()

# trajectory with energy offset

p1 = Particle(x = 0.01, p = 0.02)
p2 = Particle(x = 0.01, p = 0.02)
P1 = [copy.copy(p1)]
P2 = [copy.copy(p2)]
navi1 = Navigator()
navi2 = Navigator()
dz = 0.01

for i in range(int(lat1.totalLen/dz)):
    tracking_step(lat1, [p1], dz = dz, navi = navi1)  # R only
    tracking_step(lat2, [p2], dz = dz, navi = navi2)  # R + T
    P1.append(copy.copy(p1))
    P2.append(copy.copy(p2))

s = [p.s for p in P1]
x1 = [p.x for p in P1]
x2 = [p.x for p in P2]
plt.title("Trajectories with $\Delta p/p = 0.02$")
plt.xlabel("S, m")
plt.ylabel("X, m")
plt.plot(s, x1,'r', s, x2,'b')
plt.legend(["R", "R+T"])
plt.show()