ocelot_dir = "/Users/tomins/Nextcloud/DESY/repository/ocelot"
import sys

sys.path.append(ocelot_dir)

from accelerator.s2e_sections.sections import *
from ocelot.utils.section_track import *
from ocelot.gui.accelerator import *
import time
from ocelot.common.globals import *
import injector_lattice as i1
from ocelot.cpbd.match import *


tws0 = Twiss()
tws0.E = 0.005
tws0.beta_x = 0.286527307369
tws0.beta_y = 0.286527307369
tws0.alpha_x = -0.838833736086
tws0.alpha_y = -0.838833736086


lat = MagneticLattice(i1.cell, i1.start_sim, method=MethodTM({"global": SecondTM}))


print(lat[i1.start_sim] is i1.start_sim)
tws_ref = twiss(lat, tws0)
print(tws_ref[-1])
plot_opt_func(lat, tws_ref)
plt.show()


for elem in lat.sequence:
    if elem.__class__ is Cavity:
        if "A1" in elem.id:
            elem.v = 0.018764
            elem.phi = 15.
        #if "AH1" in elem.id:
        #    elem.v = 0
lat.update_transfer_maps()


navi = Navigator(lat)


vars = [i1.qi_46_i1, i1.qi_47_i1, i1.qi_50_i1, i1.qi_52_i1]

beta_x  = 2.8317292131504344
beta_y  = 6.651738960640371
alpha_x = 0.2919751990869057
alpha_y = -1.9571969991015152
constr = {i1.stsub_62_i1: {'beta_x': beta_x, 'beta_y': beta_y,
                           "alpha_x": alpha_x, "alpha_y": alpha_y}}
#print(constr)


p_array = generate_parray(chirp=0.0, charge=5e-9, nparticles=20000, energy=0.0065, tws=tws0)
tw = get_envelope(p_array)
print(tw)
#for i, q in enumerate(vars):
#    q.k1 = res[i]
res = match(lat, constr, vars, tw, verbose=True, max_iter=1000, method='simplex')
print(res)
tws = twiss(lat, tw)
plot_opt_func(lat, tws, top_plot=["E"])
plt.show()



for i, q in enumerate(vars):
    q.k1 = res[i]

navi = Navigator(lat)
navi.unit_step = 0.05
sc = SpaceCharge()
sc.step = 1
sc.nmesh_xyz = [63, 63, 63]
# sc = LSC()

sc2 = SpaceCharge()
sc2.step = 5
sc2.nmesh_xyz = [63, 63, 63]

sc3 = SpaceCharge()
sc3.step = 10
sc3.nmesh_xyz = [63, 63, 63]

acc1_1_stop = i1.a1_1_stop
start_sim = i1.start_sim
acc1_wake_kick = i1.a1_sim_stop
#navi.add_physics_proc(sc, start_sim, acc1_1_stop)
#navi.add_physics_proc(sc2, acc1_1_stop, acc1_wake_kick)
#navi.add_physics_proc(sc3, i1.a1_sim_stop, i1.stsub_62_i1)

res = match_beam(lat, constr, vars, p_array, navi, verbose=True, max_iter=10, method='simplex')
print(res)
navi.reset_position()
for i, q in enumerate(vars):
    q.k1 = res[i]
lat.update_transfer_maps()
tws_track, _ = track(lat, p_array, navi)
print(get_envelope(p_array))
s = [tw.s for tw in tws_track]
bx = [tw.beta_x for tw in tws_track]
by = [tw.beta_y for tw in tws_track]

plt.plot(s, bx)
plt.plot(s, by)
plt.show()

fig, ax = plot_API(lat, legend=False)
ax.plot([tw.s for tw in tws_ref], [tw.beta_x for tw in tws_ref], "C1", lw=2, label=r"ref: $\beta_x$")
ax.plot([tw.s for tw in tws_ref], [tw.beta_y for tw in tws_ref], "C2", lw=2, label=r"ref: $\beta_y$")
ax.plot(s, bx, "C1--", lw=2, label=r"track: $\beta_x$")
ax.plot(s, by, "C2--", lw=2, label=r"track: $\beta_y$")
ax.legend()
ax.set_ylabel(r"$\beta_{x,y}$ [m]")
#plot_opt_func(lat, tws_track, top_plot=["E"])
plt.show()