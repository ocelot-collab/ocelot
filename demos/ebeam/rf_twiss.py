from ocelot import *
from ocelot.gui import *
import numpy as np

phi_bc2 = 0.033646252962410

l0 = 0.5
l_bc2 = l0 *phi_bc2 /np.sin((phi_bc2))

ac_v = 0.02265625 # in GV

bb_393_b2 = Bend(l=l_bc2, angle=phi_bc2, e1=0.000000000, e2=phi_bc2, tilt=1.570796330, eid= 'bb_393_b2')
bb_402_b2 = Bend(l=l_bc2, angle=-phi_bc2, e1=-phi_bc2, e2=0.000000000, tilt=1.570796330, eid= 'bb_402_b2')
bb_404_b2 = Bend(l=l_bc2, angle=-phi_bc2, e1=0.000000000, e2=-phi_bc2, tilt=1.570796330,  eid= 'bb_404_b2')
bb_413_b2 = Bend(l=l_bc2, angle=phi_bc2, e1=phi_bc2, e2=0.000000000, tilt=1.570796330,  eid= 'bb_413_b2')

d10cm =   Drift(l=0.1, eid= 'd10cm')
cd850cm = Drift(l=8.5 / np.cos(phi_bc2), eid= 'cd850cm')
cd150cm = Drift(l=1.5, eid= 'cd150cm')
cd100cm = Drift(l=1, eid= 'cd100cm')
d34cm59 = Drift(l=0.3459, eid= 'd34cm59')
d13cm =   Drift(l=0.13, eid= 'd13cm')
d130cm =  Drift(l=1.3, eid= 'd130cm')

bc2  = (d10cm, bb_393_b2, cd850cm, bb_402_b2, cd150cm, bb_404_b2, cd850cm, bb_413_b2,  cd100cm)

qd_415_b2 = Quadrupole(l=0.2000000, k1=0.3, tilt=0.000000000, eid= 'qd_415_b2')
qd_417_b2 = Quadrupole(l=0.2000000, k1=-0.2, tilt=0.000000000, eid= 'qd_417_b2')
qd_418_b2 = Quadrupole(l=0.2000000, k1=-0.5, tilt=0.000000000, eid= 'qd_418_b2')
q_249_l2 =  Quadrupole(l=0.3000000, k1=0.25, tilt=0.000000000, eid= 'q_249_l2')
q_261_l2 =  Quadrupole(l=0.3000000, k1=-0.29711100, tilt=0.000000000, eid= 'q_261_l2')


c_a3 = Cavity(l=1.0377000, phi=0.0, v = ac_v, freq=1.300e+009, eid= 'c_a3')


l3  = (d13cm,qd_415_b2,d130cm, qd_417_b2, d130cm, qd_418_b2,d130cm,
        c_a3, d34cm59, c_a3, d34cm59, c_a3, d34cm59, c_a3, d34cm59,
        c_a3, d34cm59, c_a3, d34cm59, c_a3, d34cm59, c_a3, d13cm ,
        q_249_l2, d34cm59, c_a3, d34cm59, c_a3, d34cm59, c_a3, d34cm59, c_a3, d34cm59,
        c_a3, d34cm59, c_a3, d34cm59, c_a3, d34cm59, c_a3, d13cm, q_261_l2, d130cm)

bc2_l3  = (bc2,l3)



beam = Beam()
beam.E = 2.4
beam.beta_x = 41.1209
beam.beta_y = 86.3314
beam.alpha_x = 1.9630
beam.alpha_y = 4.0972


lat = MagneticLattice(bc2_l3)

tw0 = Twiss(beam)

tws=twiss(lat, tw0, nPoints = None)
#lat.update_transfer_maps()
#tws=twiss(lat, tw0, nPoints = None)
#lat.update_transfer_maps()
plot_opt_func(lat, tws, top_plot = ["E"])
plt.show()

f=plt.figure()
ax = f.add_subplot(211)
ax.set_xlim(0, lat.totalLen)

f.canvas.set_window_title('Betas [m]')
s = [p.s for p in tws]
p1, = plt.plot(s, [p.beta_x for p in tws], lw=2.0)
p2, = plt.plot(s, [p.beta_y for p in tws], lw=2.0)
plt.grid(True)
plt.legend([p1,p2], [r'$\beta_x$',r'$\beta_y$', r'$D_x$'])

ax2 = f.add_subplot(212)
plot_lattice(lat, ax2, alpha=0.5)

# add beam size (arbitrary scale)

s = np.array([p.s for p in tws])

scale = 5000

sig_x = scale * np.array([np.sqrt(p.beta_x*beam.emit_x) for p in tws]) # 0.03 is for plotting same scale
sig_y = scale * np.array([np.sqrt(p.beta_y*beam.emit_y) for p in tws])

x = scale * np.array([p.x for p in tws])
y = scale * np.array([p.y for p in tws])


plt.plot(s, x + sig_x, color='#0000AA', lw=2.0)
plt.plot(s, x-sig_x, color='#0000AA', lw=2.0)

plt.plot(s, sig_y, color='#00AA00', lw=2.0)
plt.plot(s, -sig_y, color='#00AA00', lw=2.0)

#f=plt.figure()
plt.plot(s, x, 'r--', lw=2.0)
#plt.plot(s, y, 'r--', lw=2.0)

plt.show()

