"""
Linac Orbit Correction.

S.Tomin. 09.2019
"""

from ocelot import *
from ocelot.gui.accelerator import *
import dogleg_lattice as dl
from ocelot.cpbd.orbit_correction import *
from ocelot.cpbd.response_matrix import *
import seaborn as sns
import logging
logging.basicConfig(level=logging.INFO)

method = MethodTM()
method.global_method = SecondTM

# introduce misalignment
dl.qi_77_i1.dx = -100e-6
dl.qi_77_i1.dy = 100e-6

dl.qi_85_i1.dx = 100e-6
dl.qi_85_i1.dy = -100e-6


lat = MagneticLattice(dl.cell,  method=method)

tws = twiss(lat, tws0=dl.tws0)
plot_opt_func(lat, tws, top_plot=["Dy"])
plt.show()


orb = Orbit(lat)

method = LinacRmatrixRM(lattice=orb.lat, hcors=orb.hcors, vcors=orb.vcors, bpms=orb.bpms)
#drm_method = LinacDisperseSimRM

orb.response_matrix = ResponseMatrix(method=method)
# in that case the initial Twiss is needed only for initial energy
orb.response_matrix.calculate(tw_init=dl.tws0)

ax = sns.heatmap(orb.response_matrix.df, annot=True)
ax.set_title("Orbit response matrix")
plt.show()


s_bpm_b = np.array([p.s for p in orb.bpms])
x_bpm_b, y_bpm_b = method.read_virtual_orbit(p_init=Particle())
fig, ax = plot_API(lat)
ax.plot(s_bpm_b, x_bpm_b*1000., "ro-", label="X [mm]")
ax.plot(s_bpm_b, y_bpm_b*1000., "bo-", label="Y [mm]")
ax.legend()
plt.show()


orb.correction(beta=500)

x_bpm, y_bpm = method.read_virtual_orbit(p_init=Particle())

p_list = lattice_track(lat, method.particle0)
s = [p.s for p in p_list]
x = [p.x*1000. for p in p_list]
y = [p.y*1000. for p in p_list]

fig, ax = plot_API(lat)
ax.plot(s_bpm_b, x_bpm*1000., "ro" , label="X [mm]")
ax.plot(s, x, 'r')
ax.plot(s_bpm_b, y_bpm*1000., "bo", label="Y [mm]")
ax.plot(s, y, 'b')
plt.show()

rm = orb.response_matrix.extract(cor_list=['CIX.78.I1', 'CIY.80.I1'], bpm_list=['BPMA.85.I1', 'BPMA.87.I1'])
print(rm)

rm[0, 0] = 2
rm[2, 0] = 0.1
rm[1, 1] = -0.1
print(rm)
orb.response_matrix.inject(cor_list=['CIX.78.I1', 'CIY.80.I1'], bpm_list=['BPMA.85.I1', 'BPMA.87.I1'], inj_matrix=rm)
rm_check = orb.response_matrix.extract(cor_list=['CIX.78.I1', 'CIY.80.I1'], bpm_list=['BPMA.85.I1', 'BPMA.87.I1'])
print("is RMs equal: ", (np.equal(rm, rm_check)).all())
ax = sns.heatmap(orb.response_matrix.df, annot=True)
ax.set_title("Orbit response matrix")
plt.show()


