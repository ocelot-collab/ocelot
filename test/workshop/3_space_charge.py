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

# import injector lattice
from ocelot.test.workshop.injector_lattice import *

p_array_init = astraBeam2particleArray(filename='beam_distrib.ast')

# initialization of tracking method
method = MethodTM()

# for second order tracking we have to choose SecondTM
method.global_method = SecondTM

# for first order tracking uncomment next line
# method.global_method = TransferMap

# we will start simulation from point 3.2 from the gun. For this purpose  marker was created (start_sim=Marker())
# and placed in 3.2 m after gun
# C3_AH1_1_8_I1 is the last section of the 3.9 GHz cavity
#lat = MagneticLattice(cell, start=start_sim, stop=C3_AH1_1_8_I1, method=method)
lat = MagneticLattice(cell, start=start_sim, stop=D_38, method=method)

sc = SpaceCharge()
sc.nmesh_xyz = [31, 31, 31]
sc.step = 1

navi = Navigator(lat)

# add physics processes from the first element to the last of the lattice
navi.add_physics_proc(sc, lat.sequence[0], lat.sequence[-1])

# definiing unit step in [m]
navi.unit_step = 0.1

# deep copy of the initial beam distribution
p_array = deepcopy(p_array_init)
tws_track, p_array = track(lat, p_array, navi)


# you can change top_plot argument, for example top_plot=["alpha_x", "alpha_y"]
plot_opt_func(lat, tws_track, top_plot=["E"], fig_name=0, legend=False)
plt.show()