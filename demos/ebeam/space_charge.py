"""
This script was created by Sergey Tomin for Workshop: Designing future X-ray FELs. Source and license info is on GitHub.
August 2016.

Modified in 2017,
S.Tomin
"""
from time import time

# this python library provides generic shallow (copy) and deep copy (deepcopy) operations
from copy import deepcopy

# import from Ocelot main modules and functions
from ocelot import *

# import from Ocelot graphical modules
from ocelot.gui.accelerator import *

from ocelot.adaptors.astra2ocelot import *


# *********************************** LATTICE: START ********************************
D_14 = Drift(l = 0.2216 + 0.0996, eid = 'D_14')
D_15 = Drift(l = 0.3459, eid = 'D_15')
D_22 = Drift(l = 0.2043, eid = 'D_22')
D_23 = Drift(l = 0.085+0.4579 + 0.2211 +0.085 , eid = 'D_23')

phi1=18.7268
V1=18.50662e-3/np.cos(phi1*pi/180)

C_A1_1_1_I1 = Cavity(l = 1.0377, v = V1, freq = 1.3e9, phi = phi1, eid = 'C.A1.1.1.I1')
C_A1_1_2_I1 = Cavity(l = 1.0377, v = V1, freq = 1.3e9, phi = phi1, eid = 'C.A1.1.2.I1')
C_A1_1_3_I1 = Cavity(l = 1.0377, v = V1, freq = 1.3e9, phi = phi1, eid = 'C.A1.1.3.I1')
C_A1_1_4_I1 = Cavity(l = 1.0377, v = V1, freq = 1.3e9, phi = phi1, eid = 'C.A1.1.4.I1')
C_A1_1_5_I1 = Cavity(l = 1.0377, v = V1, freq = 1.3e9, phi = phi1, eid = 'C.A1.1.5.I1')
C_A1_1_6_I1 = Cavity(l = 1.0377, v = V1, freq = 1.3e9, phi = phi1, eid = 'C.A1.1.6.I1')
C_A1_1_7_I1 = Cavity(l = 1.0377, v = V1, freq = 1.3e9, phi = phi1, eid = 'C.A1.1.7.I1')
C_A1_1_8_I1 = Cavity(l = 1.0377, v = V1, freq = 1.3e9, phi = phi1, eid = 'C.A1.1.8.I1')

phi13=180
V13=-20.2E-3/8/np.cos(phi13*pi/180)
C3_AH1_1_1_I1 = Cavity(l = 0.346, v = V13, freq = 3.9e9, phi = phi13, eid = 'C3.AH1.1.1.I1')
C3_AH1_1_2_I1 = Cavity(l = 0.346, v = V13, freq = 3.9e9, phi = phi13, eid = 'C3.AH1.1.2.I1')
C3_AH1_1_3_I1 = Cavity(l = 0.346, v = V13, freq = 3.9e9, phi = phi13, eid = 'C3.AH1.1.3.I1')
C3_AH1_1_4_I1 = Cavity(l = 0.346, v = V13, freq = 3.9e9, phi = phi13, eid = 'C3.AH1.1.4.I1')
C3_AH1_1_5_I1 = Cavity(l = 0.346, v = V13, freq = 3.9e9, phi = phi13, eid = 'C3.AH1.1.5.I1')
C3_AH1_1_6_I1 = Cavity(l = 0.346, v = V13, freq = 3.9e9, phi = phi13, eid = 'C3.AH1.1.6.I1')
C3_AH1_1_7_I1 = Cavity(l = 0.346, v = V13, freq = 3.9e9, phi = phi13, eid = 'C3.AH1.1.7.I1')
C3_AH1_1_8_I1 = Cavity(l = 0.346, v = V13, freq = 3.9e9, phi = phi13, eid = 'C3.AH1.1.8.I1')

Q_37_I1 = Quadrupole(l = 0.3, k1 = -1.537886, tilt = 0.0, eid = 'Q.37.I1')
Q_38_I1 = Quadrupole(l = 0.3, k1 = 1.435078, tilt = 0.0, eid = 'Q.38.I1')

start_sim = Marker()

cell = (start_sim, D_14, C_A1_1_1_I1, D_15, C_A1_1_2_I1,
D_15, C_A1_1_3_I1, D_15, C_A1_1_4_I1, D_15, C_A1_1_5_I1, D_15, C_A1_1_6_I1,
D_15, C_A1_1_7_I1, D_15, C_A1_1_8_I1, D_22, Q_37_I1, D_23, Q_38_I1)

# *********************************** LATTICE: STOP ********************************



# ******************************* GENERATE BEAM DISTRIBUTION: START **************************

# load astra file
# >> p_array_init = astraBeam2particleArray(filename='beam_6MeV.ast')
# >> tws = get_envelope(p_array_init)
# calculation std from astra beam file
# >> print(np.std(p_array_init.x()))
# >> print(np.std(p_array_init.px()))
# >> print(np.std(p_array_init.y()))
# >> print(np.std(p_array_init.py()))

sigma_x = 0.000231507245956
sigma_px =0.000204206874319
sigma_y = 0.000231583942392
sigma_py =0.000204272734636

x = np.random.randn(200000)*sigma_x
px = np.random.randn(200000)*sigma_px
y = np.random.randn(200000)*sigma_y
py = np.random.randn(200000)*sigma_py

# To reproduce the simplest tau-dp/p dependency the covariance matrix was calculated from astra beam distribution. See below:
# >> indx = np.argsort(p_array_init.tau())
# >> t = np.atleast_2d(p_array_init.tau()[indx[50000:150000]]).T
# >> ps = np.atleast_2d(p_array_init.p()[indx[50000:150000]]).T
# >> c = np.cov(np.hstack((t, ps)).T)
# >> print("cov=", c)

cov_t_p =  [[  6.89508231e-07,  -2.98688604e-07],
 [ -2.98688604e-07,   1.87434257e-07]]

long_dist = np.random.multivariate_normal((0, 0), cov_t_p, 200000)
tau = long_dist[:, 0]
dp = long_dist[:, 1]

p_array_init = ParticleArray(n=200000)
p_array_init.E = 0.0065 # GeV
p_array_init.rparticles[0] = x
p_array_init.rparticles[1] = px
p_array_init.rparticles[2] = y
p_array_init.rparticles[3] = py
p_array_init.rparticles[4] = tau
p_array_init.rparticles[5] = dp

# creating charge array
Q = 5e-10
p_array_init.q_array = np.ones(200000)*Q/200000

# beam transformation
beta_x  = 1.59966676201
beta_y  = 1.60018325757
alpha_x = -0.995487979563
alpha_y = -0.996116091572
mux = 0
muy = 0
bt = BeamTransform([alpha_x, beta_x, mux], [alpha_y, beta_y, muy])
bt.apply(p_array_init,dz=0)

sI1, I1 = get_current(p_array_init, charge=p_array_init.q_array[0], num_bins=200)

#plt.figure(1)
#plt.title("energy distribution: start")
#plt.plot(-p_array_init.tau()*1000, p_array_init.p(), 'r.')
#plt.xlabel("s, mm")
#plt.ylabel("dp/p")
#plt.grid(True)

#plt.figure(2)
#plt.title("current: start")
#plt.plot(sI1*1000, I1)
#plt.xlabel("s, mm")
#plt.ylabel("I, A")
#plt.grid(True)

# ******************************* GENERATE BEAM DISTRIBUTION: END **************************


# initialization of tracking method
method = MethodTM()

# for second order tracking we have to choose SecondTM
method.global_method = SecondTM

# for first order tracking uncomment next line
# method.global_method = TransferMap

# we will start simulation from point 3.2 from the gun. For this purpose  marker was created (start_sim=Marker())
# and placed in 3.2 m after gun
# Q_38_I1 is quadrupole between RF cavities 1.3 GHz and 3.9 GHz
# C3_AH1_1_8_I1 is the last section of the 3.9 GHz cavity
lat = MagneticLattice(cell, method=method)

sc1 = SpaceCharge()
sc1.nmesh_xyz = [63, 63, 63]
#sc1.low_order_kick = False
sc1.step = 1

sc5 = SpaceCharge()
sc5.nmesh_xyz = [63, 63, 63]
sc5.step = 5
#sc5.low_order_kick = False


navi = Navigator(lat)

# add physics processes from the first element to the last of the lattice
navi.add_physics_proc(sc1, lat.sequence[0], C_A1_1_2_I1)
navi.add_physics_proc(sc5, C_A1_1_2_I1, lat.sequence[-1])

# definiing of unit step in [m]
navi.unit_step = 0.02

# deep copy of the initial beam distribution
p_array = deepcopy(p_array_init)
start = time()
tws_track, p_array = track(lat, p_array, navi)
print("time exec: ", time() - start, "sec")

# you can change top_plot argument, for example top_plot=["alpha_x", "alpha_y"]
plot_opt_func(lat, tws_track, top_plot=["E"], fig_name=0, legend=False)


sI1_e, I1_e = get_current(p_array, charge=p_array.q_array[0], num_bins=200)

plt.figure(1)
plt.title("energy distribution: Start-End")
plt.plot(-p_array_init.tau()*1000, p_array_init.p(), 'r.', label="start")
plt.plot(-p_array.tau()*1000, p_array.p(), 'b.', label="end")
plt.legend()
plt.xlabel("s, mm")
plt.ylabel("dp/p")
plt.grid(True)

plt.figure(2)
plt.title("current: Start-End")
plt.plot(sI1*1000, I1, "r", label="start")
plt.plot(sI1_e*1000, I1_e, "b", label="end")
plt.legend()
plt.xlabel("s, mm")
plt.ylabel("I, A")
plt.grid(True)

plt.show()
