"""
This script was created by Sergey Tomin for Workshop: Designing future X-ray FELs. Source and license info is on GitHub.
August 2016.
"""
from time import time

# this python library provides generic shallow (copy) and deep copy (deepcopy) operations
from copy import deepcopy

# import from Ocelot main modules and functions
from ocelot import *

# import from Ocelot graphical modules
from ocelot.gui.accelerator import *

from ocelot.adaptors.astra2ocelot import *

# import injector lattice
from ocelot.test.workshop.injector_lattice import *


phi1=18.7268
V1=18.50662e-3/np.cos(phi1*pi/180)

C_A1_1_1_I1.v = V1; C_A1_1_1_I1.phi = phi1
C_A1_1_2_I1.v = V1; C_A1_1_2_I1.phi = phi1
C_A1_1_3_I1.v = V1; C_A1_1_3_I1.phi = phi1
C_A1_1_4_I1.v = V1; C_A1_1_4_I1.phi = phi1
C_A1_1_5_I1.v = V1; C_A1_1_5_I1.phi = phi1
C_A1_1_6_I1.v = V1; C_A1_1_6_I1.phi = phi1
C_A1_1_7_I1.v = V1; C_A1_1_7_I1.phi = phi1
C_A1_1_8_I1.v = V1; C_A1_1_8_I1.phi = phi1

phi13=180
V13=-20.2E-3/8/np.cos(phi13*pi/180)

C3_AH1_1_1_I1.v=V13; C3_AH1_1_1_I1.phi=phi13
C3_AH1_1_2_I1.v=V13; C3_AH1_1_2_I1.phi=phi13
C3_AH1_1_3_I1.v=V13; C3_AH1_1_3_I1.phi=phi13
C3_AH1_1_4_I1.v=V13; C3_AH1_1_4_I1.phi=phi13
C3_AH1_1_5_I1.v=V13; C3_AH1_1_5_I1.phi=phi13
C3_AH1_1_6_I1.v=V13; C3_AH1_1_6_I1.phi=phi13
C3_AH1_1_7_I1.v=V13; C3_AH1_1_7_I1.phi=phi13
C3_AH1_1_8_I1.v=V13; C3_AH1_1_8_I1.phi=phi13


p_array_init = astraBeam2particleArray(filename='beam_6MeV.ast')

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
lat = MagneticLattice(cell, start=start_sim, stop=Q_38_I1, method=method)

sc1 = SpaceCharge()
sc1.nmesh_xyz = [63, 63, 63]
sc1.low_order_kick = False
sc1.step = 1

sc5 = SpaceCharge()
sc5.nmesh_xyz = [63, 63, 63]
sc5.step = 5
sc5.low_order_kick = False


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
plt.show()

#particleArray2astraBeam(p_array, "beam_sc_after_RF.ast")
# Comparison with ASTRA
# Beam tracking with ASTRA was performed by Igor Zagorodnov (DESY).

sa, bx_sc, by_sc, bx_wo_sc, by_wo_sc = np.loadtxt("astra_sim.txt", usecols=(0, 1, 2, 3, 4), unpack=True)

s = [tw.s for tw in tws_track]
bx = [tw.beta_x for tw in tws_track]
by = [tw.beta_y for tw in tws_track]
ax = plot_API(lat, legend=False)
ax.plot(s, bx, "r", label="Ocelot, bx")
ax.plot(sa-3.2, bx_sc, "b-",label="ASTRA, bx")
ax.plot(s, by, "r", label="Ocelot, by")
ax.plot(sa-3.2, by_sc, "b-",label="ASTRA, by")
ax.legend()
plt.show()