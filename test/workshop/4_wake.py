"""
This script was created by Sergey Tomin for Workshop: Designing future X-ray FELs. Source and license info is on GitHub.
August 2016.
"""
# this python library provides generic shallow (copy) and deep copy (deepcopy) operations
from copy import deepcopy

# import from Ocelot main modules and functions
from ocelot import *

# import from Ocelot graphical modules
from ocelot.gui.accelerator import *

from ocelot.adaptors.astra2ocelot import *

# LATTICE
D00m25 = Drift(l = 0.25)
D01m = Drift(l = 1)
D02m = Drift(l = 2)

# Create markers for defining places of the wakes applying
w1_start = Marker()
w1_stop = Marker()

w2_start = Marker()
w2_stop = Marker()

w3_start = Marker()
w3_stop = Marker()

w4_start = Marker()
w4_stop = Marker()

w5_start = Marker()
w5_stop = Marker()

w6_start = Marker()
w6_stop = Marker()
# quadrupoles
Q1 = Quadrupole(l = 0.5, k1 = 0.215)

# lattice
lattice = (D01m, w1_start, D02m, w1_stop, w2_start, D02m, w2_stop, w3_start, D02m, w3_stop, D00m25, Q1,
           D00m25, w4_start, D02m, w4_stop, w5_start, D02m, w5_stop, w6_start, D02m, w6_stop, D01m)

# creation MagneticLattice
method = MethodTM()
method.global_method = SecondTM
lat = MagneticLattice(lattice, method=method)

# Load beam file
p_array_init = astraBeam2particleArray(filename='beam_chirper.ast')

# Initialization of the wakes and the places of their applying

from ocelot.cpbd.wake3D import *

# load wake tables of corrugated structures
wk_vert = WakeTable('wake_vert_1m.txt')
wk_hor = WakeTable('wake_hor_1m.txt')

# creation of wake object with parameters
wake_v1 = Wake()
wake_v1.w_sampling = 500
wake_v1.wake_table = wk_vert
wake_v1.step = 1 # step in Navigator.unit_step, dz = Navigator.unit_step * wake.step [m]

wake_h1 = Wake()
wake_h1.w_sampling = 500
wake_h1.wake_table = wk_hor
wake_h1.step = 1

wake_v2 = deepcopy(wake_v1)

wake_h2 = deepcopy(wake_h1)

wake_v3 = deepcopy(wake_v1)

wake_h3 = deepcopy(wake_h1)

navi = Navigator(lat)

# Add the wakes in the lattice
# add physics proccesses
navi.add_physics_proc(wake_v1, w1_start, w1_stop)
navi.add_physics_proc(wake_h1, w2_start, w2_stop)
navi.add_physics_proc(wake_v2, w3_start, w3_stop)
navi.add_physics_proc(wake_h2, w4_start, w4_stop)
navi.add_physics_proc(wake_v3, w5_start, w5_stop)
navi.add_physics_proc(wake_h3, w6_start, w6_stop)

# definiing unit step in [m]
navi.unit_step = 0.2

# deep copy of the initial beam distribution
p_array = deepcopy(p_array_init)
tws_track, p_array = track(lat, p_array, navi)

# Longitudinal beam distribution
tau0 = p_array_init.tau()
p0 = p_array_init.p()

tau1 = p_array.tau()
p1 = p_array.p()

plt.figure(1)
plt.plot(-tau0*1000, p0, "r.", -tau1*1000, p1, "b.")
plt.legend(["before", "after"], loc=4)
plt.grid(True)
plt.xlabel(r"$\tau$, mm")
plt.ylabel(r"$\frac{\Delta E}{E}$")

# Beam distribution
tau = np.array([p.tau for p in p_array])
dp = np.array([p.p for p in p_array])
x = np.array([p.x for p in p_array])
y = np.array([p.y for p in p_array])

plt.figure(2)
ax1 = plt.subplot(311)
ax1.plot(-tau*1000, x*1000, 'r.')
plt.setp(ax1.get_xticklabels(), visible=False)
plt.ylabel("x, mm")
plt.grid(True)

ax2 = plt.subplot(312, sharex=ax1)
ax2.plot(-tau*1000, y*1000, 'r.')
plt.setp(ax2.get_xticklabels(), visible=False)
plt.ylabel("y, mm")
plt.grid(True)

ax3 = plt.subplot(313, sharex=ax1)
ax3.plot(-tau*1000, dp, 'r.')
plt.ylabel("dp/p")
plt.xlabel("s, mm")
plt.grid(True)
plt.show()