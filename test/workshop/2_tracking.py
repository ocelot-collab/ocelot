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

lat = MagneticLattice(cell)

# initialization of Twiss object
tws0 = Twiss()
# defining initial twiss parameters
tws0.beta_x = 29.171
tws0.beta_y = 29.171
tws0.alpha_x = 10.955
tws0.alpha_y = 10.955
# defining initial electron energy in GeV
tws0.E = 0.005

# calculate twiss functions with initial twiss parameters
tws = twiss(lat, tws0, nPoints=None)
tws1 = tws[-1]
print(tws[-1])
# ploting twiss paramentrs.
plot_opt_func(lat, tws, top_plot=["Dx"], fig_name="i1", legend=False)

# Loading of beam distribution
p_array_init = astraBeam2particleArray(filename='beam_130MeV.ast')


# initialization of tracking method
method = MethodTM()

# for second order tracking we have to choose SecondTM
method.global_method = SecondTM

# for first order tracking uncomment next line
# method.global_method = TransferMap

# we will start simulation from the first quadrupole (QI.46.I1) after RF section.
# you can change stop element (and the start element, as well)
# START_73_I1 - marker before Dog leg
# START_96_I1 - marker before Bunch Compresion

lat = MagneticLattice(cell, start=QI_46_I1, stop=None, method=method)


navi = Navigator(lat)
p_array = deepcopy(p_array_init)
tws_track, p_array = track(lat, p_array, navi)

# you can change top_plot argument, for example top_plot=["alpha_x", "alpha_y"]
plot_opt_func(lat, tws_track, top_plot=["E"], fig_name=0, legend=False)
plt.show()

# Current profile
bins_start, hist_start = get_current(p_array, charge=p_array.q_array[0], num_bins=200)

plt.figure(4)
plt.title("current: end")
plt.plot(bins_start*1000, hist_start)
plt.xlabel("s, mm")
plt.ylabel("I, A")
plt.grid(True)
plt.show()

# Beam distribution

tau = np.array([p.tau for p in p_array])
dp = np.array([p.p for p in p_array])
x = np.array([p.x for p in p_array])
y = np.array([p.y for p in p_array])

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