__author__ = 'Sergey Tomin'

"""
tracking of the 200k particles. 2nd order + Space Change effect + Wakes of the ACC1 cavity. Till the end of ACC39.
"""

from ocelot import *
from ocelot.gui.accelerator import *
from ocelot.adaptors.astra2ocelot import *
from ocelot.test.workshop.injector_sim import *

tws0 = Twiss()
tws0.beta_x = 29.171
tws0.beta_y = 29.171
tws0.alpha_x = 10.955
tws0.alpha_y = 10.955
tws0.E = 0.005


method = MethodTM()
method.global_method = SecondTM

lat = MagneticLattice(cell, start=start_sim, stop=C3_AH1_1_8_I1, method=method)

p_array = astraBeam2particleArray(filename='Exfel.0320.ast')

navi = Navigator(lat)
navi.unit_step = 0.1

sc = SpaceCharge()
sc.step = 1

wake1 = Wake()
wake1.wake_table = WakeTable('TESLA_MODULE_WAKE_TAYLOR.dat')
wake1.step = 12
wake1.factor = 1./11.0688


navi.add_physics_proc(sc, lat.sequence[0], lat.sequence[-1])
navi.add_physics_proc(wake1, C_A1_1_1_I1, D_9)
tws_track, p_array = track(lat, p_array, navi)

# plotting

plot_opt_func(lat, tws_track, top_plot=["E"], fig_name=0, legend=False)

bins_start, hist_start = get_current(p_array, charge=p_array.q_array[0], num_bins=200)

tau = np.array([p.tau for p in p_array])
dp = np.array([p.p for p in p_array])
x = np.array([p.x for p in p_array])
y = np.array([p.y for p in p_array])

plt.figure(1)
plt.plot(tau*1000, x*1000, 'r.')
plt.xlabel("s, mm")
plt.ylabel("x, mm")
plt.grid(True)

plt.figure(2)
plt.plot(tau*1000, y*1000, 'r.')
plt.xlabel("s, mm")
plt.ylabel("y, mm")
plt.grid(True)

plt.figure(3)
plt.plot(tau*1000, dp, 'r.')
plt.xlabel("s, mm")
plt.ylabel("dp/p")
plt.grid(True)

plt.figure(4)
plt.title("current: end")
plt.plot(bins_start*1000, hist_start)
plt.xlabel("s, mm")
plt.ylabel("I, A")
plt.grid(True)
plt.show()
