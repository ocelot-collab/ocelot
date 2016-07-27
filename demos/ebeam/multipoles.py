__author__ = 'Sergey Tomin'

from copy import copy
from ocelot import *
from ocelot.gui import *
import numpy as np

C = 1000.
Ncells = 16
Nbends = 32
D = Drift(l=C/Ncells/4, eid="D")
Qf = Multipole(kn=[0., 0.021/2.], eid="Qf")
Qd = Multipole(kn=[0., -0.02], eid="Qd")
B = Multipole(kn=2.*pi/Nbends)
Sf = Multipole(kn=(0., 0., 0.0), eid="Sf")
Sd = Multipole(kn=(0., 0., -0.0), eid="Sd")
F = Multipole(kn=[0., 0., 0., 0., 0.1])
cell = (Qf,Sf, D,F,B, D, Qd, Sd, D, B, D, Sf, Qf)

lat = MagneticLattice(Ncells*cell)
tws = twiss(lat)
plot_opt_func(lat, tws)
plt.show()

compensate_chromaticity(lat, ksi_x_comp=0, ksi_y_comp=0)
# Single particle tracking with energy shift
p1 = Particle(x = -0.0001)
p2 = Particle(x = 0.001)

navi = Navigator()
dz = 1.
P1 = []
P2 = []
for i in range(int(lat.totalLen/dz)):
    tracking_step(lat, [p1, p2], dz=dz, navi=navi)
    P1.append(copy(p1))
    P2.append(copy(p2))

s = [f.s for f in P1]#map(lambda f: f.s, P1)
x = [f.x for f in P1]#map(lambda f: f.x, P1)
y = [f.y for f in P1]#map(lambda f: f.y, P1)
tau = [f.tau for f in P1]# map(lambda f: f.tau, P1)
#s = map(lambda f: f.s, P1)
x2 = [f.x for f in P2]#map(lambda f: f.x, P2)
y2 = [f.y for f in P2]#map(lambda f: f.y, P2)
tau2 = [f.tau for f in P2]#map(lambda f: f.tau, P1)




plt.plot(s, x, "r", label = "x1")
plt.plot(s, x2, "b", label = "x2")
plt.xlabel("S, m")
plt.ylabel("x1/x2, m")
plt.legend()
plt.show()
nturns = 1000
x_array = np.linspace(0.4, 0.60, num=100)
track_list = create_track_list(x_array, [0.], [0.], energy=0.)
track_list = track_nturns(lat, nturns, track_list, save_track=True)
track_list = stable_particles(track_list, nturns)
print ("number stable particles = ", len(track_list))
for xyp in track_list:
    x = xyp.get_x()
    px = xyp.get_xp()
    plt.plot(x, px, ".")
plt.show()
