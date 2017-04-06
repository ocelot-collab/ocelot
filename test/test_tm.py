__author__ = 'Sergey Tomin'

from ocelot import *
from copy import deepcopy
import numpy as np


d = Drift(l=0.3)
b1 = Bend(l=0.4, k1=1.3, angle=0.4)
b2 = RBend(l=0.4, k1=1.3, angle=0.4)
b3 = SBend(l=0.4, k1=1.3, angle=0.4)
s = Sextupole(l=0.2, k2=11)
oq = Octupole(l=0.1, k3=6)
c = Cavity(l=0.6, v=13e6, phi=12*pi/180., freq=1e6)
u = Undulator(lperiod=0.02, nperiods=50, Kx=1.2)
sol = Solenoid(l=1.2, k=2.3)
m = Multipole(kn=[0.1, 0.2, 0.3, 4])
hc = Hcor(l=0.2, angle=1e-3)
vc = Vcor(l=0.2, angle=1e-3)

cell = [d, b1, b2, b3, s, oq, c, u, sol, m, hc, vc]

method = MethodTM()
method.global_method = TransferMap
lat = MagneticLattice(cell, method=method)

p0 = Particle(x=0.0002, px=-0.0012, y=-0.0005, py=0.0004, tau=0.00003, p=0.001, E=1)

for elem in lat.sequence:
    p1 = deepcopy(p0)
    p2 = deepcopy(p0)
    x0 = np.array([p1.x, p1.px, p1.y, p1.py, p1.tau, p1.p])
    p1 = elem.transfer_map*p1
    elem.transfer_map.apply(p2)

    R = elem.transfer_map.R(p1.E)
    B = elem.transfer_map.B(p1.E)
    x = np.array([p1.x, p1.px, p1.y, p1.py, p1.tau, p1.p])
    a = np.add(np.transpose(np.dot(R, np.transpose(x0.reshape(1, 6)))), B).reshape(6)
    #print(a, x)
    print(elem.__class__.__name__, p1.__str__() == p2.__str__(), np.array_equal(x, a))