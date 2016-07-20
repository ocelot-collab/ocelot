__author__ = 'Sergey Tomin'

from ocelot import *
from ocelot.gui import *
import copy as cp
import numpy as np
"""
ro*tan(phi)+
"""

Q1 = Quadrupole(l=0.2, k1=5)
Q2 = Quadrupole(l=0.2, k1=-5)
Db = Drift(l = 2)
Dc = Drift(l = 3/2.)

# to get non dispersive section
angle = 45.*pi/180.
phi = Q1.l*np.sqrt(Q1.k1)
Lc = 2.*Dc.l + Q2.l
ro = (1./np.sqrt(Q1.k1)*(Lc*np.sqrt(Q1.k1)*np.cos(phi) + 2.*np.sin(phi))/(Lc*np.sqrt(Q1.k1)*np.sin(phi) - 2.*np.cos(phi)) - Db.l)/np.tan(angle/2.)


B1 = SBend(l = ro*angle, angle=-angle)
B2 = SBend(l = ro*angle, angle=angle)
lattice = [B1, Db, Q1, Dc, Q2, Dc, Q1, Db, B2]
method1 = MethodTM()
#method1.global_method = TransferMap
method2 = MethodTM()
method2.global_method = SecondTM
lat1 = MagneticLattice(cp.deepcopy(lattice), method=method1)
lat2 = MagneticLattice(cp.deepcopy(lattice), method=method2)
tw0 = Twiss()
tw0.beta_x = 5.
tw0.alpha_x = 1.4
tw0.beta_y = 16
tw0.alpha_y = 5.1

tws = twiss(lat1, tw0, nPoints=500)

plot_opt_func(lat1,tws)
plt.show()

# trajectory with energy offset

p1 = Particle(x=0.00, p=0.02)
p2 = Particle(x=0.00, p=0.02)
P1 = [cp.copy(p1)]
P2 = [cp.copy(p2)]
navi1 = Navigator()
navi2 = Navigator()
dz = 0.01
for i in range(int(lat1.totalLen/dz)):
    tracking_step(lat1, [p1], dz = dz, navi = navi1)  # R only
    tracking_step(lat2, [p2], dz = dz, navi = navi2)  # R + T
    P1.append(cp.copy(p1))
    P2.append(cp.copy(p2))

s = [p.s for p in P1]
x1 = [p.x for p in P1]
x2 = [p.x for p in P2]
plt.title("Trajectories with $\Delta p/p = 0.02$")
plt.xlabel("S, m")
plt.ylabel("X, m")
plt.plot(s, x1,'r', s, x2,'b')
plt.legend(["R", "R+T"], loc=2)
plt.show()

# phase trajectory

t = np.linspace(0, 2*pi, num = 100)
x, xp = 0.1*np.cos(t), 0.1*np.sin(t)
plist = []
for xi, xpi in zip(x,xp):
    plist.append(Particle(x = xi, px= xpi))

plist_1 = cp.deepcopy(plist)
navi = Navigator()

tracking_step(lat1, plist_1, dz=lat1.totalLen, navi=cp.copy(navi))
tracking_step(lat2, plist, dz=lat2.totalLen, navi=cp.copy(navi))
x2 = [f.x for f in plist]
xp2 = [f.px for f in plist]

x1 = [f.x for f in plist_1]
xp1 = [f.px for f in plist_1]

plt.suptitle("Tracking")
plt.subplot(121)
plt.title("S = 0 m")
plt.plot(x, xp, "r.-", label = "X")
plt.xlabel("X, m")
plt.ylabel("Xp, rad")
plt.grid(True)

plt.subplot(122)
plt.title("end of section")
plt.plot(x2, xp2, "r.-",label="R + T")
plt.plot(x1, xp1, 'b.-', label="R")
plt.legend(loc=2)
plt.xlabel("X, m")
plt.ylabel("Xp, rad")
plt.grid(True)
plt.show()