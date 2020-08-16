__author__ = 'Sergey Tomin'

from ocelot import *
from ocelot.gui.accelerator import *
from ocelot.cpbd.orbit_correction import *
from ocelot.cpbd.response_matrix import *
import seaborn as sns
import logging
logging.basicConfig(level=logging.INFO)

# **************************** LATTICE: START ***********************************************
D0 = Drift (l = 0., eid= "D0")
D1 = Drift (l = 1.49, eid= "D1")
D2 = Drift (l = 0.1035, eid= "D2")
D3 = Drift (l = 0.307, eid= "D3")
D4 = Drift (l = 0.33, eid= "D4")
D5 = Drift (l = 0.3515, eid= "D5")
D6 = Drift (l = 0.3145, eid= "D6")
D7 = Drift (l = 0.289, eid= "D7")
D8 = Drift (l = 0.399, eid= "D8")
D9 = Drift (l = 3.009, eid= "D9")

SF = Sextupole(l = 0.0001, k2 = 17673.786254063251*1, eid= "SF")
SD = Sextupole(l = 0.0001, k2 =-36169.817233025707*1, eid= "SD")

q1 = Quadrupole (l = 0.293, k1 = 2.62, eid= "Q1")
q1.dx = 0.001
q1.dy = 0.001
q2 = Quadrupole (l = 0.293, k1 = -3.1, eid= "Q2")
q3 = Quadrupole (l = 0.327, k1 = 2.8, eid= "Q3")
q4 = Quadrupole (l = 0.291, k1 = -3.7, eid= "Q4")
q5 = Quadrupole (l = 0.391, k1 = 4.0782, eid= "Q5")
q6 = Quadrupole (l = 0.291, k1 = -3.534859, eid= "D6")

q1s = Quadrupole (l = 0.293, k1 = 2.62, eid= "Q1")
q1s.dx = 0.001
q1s.dy = 0.001
q2s = Quadrupole (l = 0.293, k1 = -3.1, eid= "Q2")
q3s = Quadrupole (l = 0.327, k1 = 2.8, eid= "Q3")
q4s = Quadrupole (l = 0.291, k1 = -3.7, eid= "Q4")
q5s = Quadrupole (l = 0.391, k1 = 4.0782, eid= "Q5")
q6s = Quadrupole (l = 0.291, k1 = -3.534859, eid= "D6")

M1 = Monitor(eid="M1")
M2 = Monitor(eid="M2")
M3 = Monitor(eid="M3")
M4 = Monitor(eid="M4")
M5 = Monitor(eid="M5")
M6 = Monitor(eid="M6")
H1 = Hcor(eid="H1")
H2 = Hcor(eid="H2")
H3 = Hcor(eid="H3")
H4 = Hcor(eid="H4")
H5 = Hcor(eid="H5")
H6 = Hcor(eid="H6")
V1 = Vcor(eid="V1")
V2 = Vcor(eid="V2")
V3 = Vcor(eid="V3")
V4 = Vcor(eid="V4")
V5 = Vcor(eid="V5")
V6 = Vcor(eid="V6")

M1s = Monitor(eid="M1s")
M2s = Monitor(eid="M2s")
M3s = Monitor(eid="M3s")
M4s = Monitor(eid="M4s")
M5s = Monitor(eid="M5s")
M6s = Monitor(eid="M6s")
H1s = Hcor(eid="H1s")
H2s = Hcor(eid="H2s")
H3s = Hcor(eid="H3s")
H4s = Hcor(eid="H4s")
H5s = Hcor(eid="H5s")
H6s = Hcor(eid="H6s")
V1s = Vcor(eid="V1s")
V2s = Vcor(eid="V2s")
V3s = Vcor(eid="V3s")
V4s = Vcor(eid="V4s")
V5s = Vcor(eid="V5s")
V6s = Vcor(eid="V6s")

B1 = SBend(l = 0.23, angle = 0.23/19.626248, eid= "B1")
B2 = SBend(l = 1.227, angle = 1.227/4.906312, eid= "B2")

Q1 = [q1, M1, H1, V1]
Q2 = [q2, M2, H2, V2]
Q3 = [q3, M3, H3, V3]
Q4 = [q4, M4, H4, V4]
Q5 = [q5, M5, H5, V5]
Q6 = [q6, M6, H6, V6]
Q1s = [q1s, M1s, H1s, V1s]
Q2s = [q2s, M2s, H2s, V2s]
Q3s = [q3s, M3s, H3s, V3s]
Q4s = [q4s, M4s, H4s, V4s]
Q5s = [q5s, M5s, H5s, V5s]
Q6s = [q6s, M6s, H6s, V6s]

cell = ( D1,SF, D2,Q1,D3, Q2,D2,SD,D4,B1,B2,D5,Q3,D5,B2,B1,D6,Q4,D7,Q5,D8,Q6,D9,Q6s,D8,Q5s,D7,Q4s,D6,B1,B2,D5,Q3s,D5,B2,B1,D4,SD,D2,Q2s,D3,Q1s,D2,SF,D1)

# **************************** LATTICE: END ***********************************************


beam = Beam()
beam.E = 2.5
beam.sigma_E = 0.001
beam.I = 0.1

method = MethodTM()
#method.params[Sextupole] = "kick"
lat = MagneticLattice(cell, method=method)

tw0 = Twiss(beam)
tws=twiss(lat, tw0, nPoints=1000)

plot_opt_func(lat, tws, top_plot=["Dx"])
plt.show()


orb = NewOrbit(lat)

method = RingRM(lattice=orb.lat, hcors=orb.hcors, vcors=orb.vcors, bpms=orb.bpms)
orb.response_matrix = ResponseMatrix(method=method)
orb.response_matrix.calculate()
cor_names_ref = orb.response_matrix.cor_names
bpm_names_ref = orb.response_matrix.bpm_names
#print("ref = ", orb.response_matrix.matrix)
orb.response_matrix.dump("test.p")
orb.response_matrix.load("test.p")

ax = sns.heatmap(orb.response_matrix.df, annot=True)
ax.set_title("Orbit response matrix")
plt.show()

s_bpm_b = np.array([p.s for p in orb.bpms])
x_bpm_b, y_bpm_b = method.read_virtual_orbit()
fig, ax = plot_API(lat)
ax.plot(s_bpm_b, x_bpm_b*1000., "ro-")
ax.plot(s_bpm_b, y_bpm_b*1000., "bo-")
plt.show()


orb.correction(beta=500)

x_bpm, y_bpm = method.read_virtual_orbit()

p_list = lattice_track(lat, method.particle0)
s = [p.s for p in p_list]
x = [p.x*1000. for p in p_list]
y = [p.y*1000. for p in p_list]
fig, ax = plot_API(lat)
ax.plot(s_bpm_b, x_bpm*1000., "ro")
ax.plot(s, x, 'r')
ax.plot(s_bpm_b, y_bpm*1000., "bo")
ax.plot(s, y, 'b')
plt.show()