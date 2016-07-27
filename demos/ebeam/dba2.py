from ocelot import *
from ocelot.gui import*
from pylab import *

from copy import deepcopy
from scipy.optimize import *

D0 = Drift (l = 0., eid= "D0")
D1 = Drift (l = 1.49, eid= "D1")
D2 = Drift (l = 0.1035, eid= "D2")
D3 = Drift (l = 0.307, eid= "D3")
D4 = Drift (l = 0.33, eid= "D4")
D5 = Drift (l = 0.3515, eid= "D5")
D6 = Drift (l = 0.3145, eid= "D6")
D7 = Drift (l = 0.289, eid= "D7")
D8 = Drift (l = 0.399, eid= "D8")
D9 = Drift (l = 3.009/2., eid= "D9")

SF = Sextupole(l = 0.01, k2 = 176.73786254063251, eid= "SF")
SD = Sextupole(l = 0.01, k2 = -361.69817233025707, eid= "SD")


Q1 = Quadrupole (l = 0.293, k1 = 2.62, eid= "Q1")
Q2 = Quadrupole (l = 0.293, k1 = -3.1, eid= "Q2")
Q3 = Quadrupole (l = 0.327, k1 = 2.8, eid= "Q3")
Q4 = Quadrupole (l = 0.291, k1 = -3.7, eid= "Q4")
Q5 = Quadrupole (l = 0.391, k1 = 4.0782, eid= "Q5")
Q6 = Quadrupole (l = 0.291, k1 = -3.534859, eid= "D6")

B1 = SBend(l = 0.23, angle = 0.23/19.626248, eid= "B1")
B2 = SBend(l = 1.227, angle = 1.227/4.906312, eid= "B2")

#und = Undulator (nperiods=200,lperiod=0.07,Kx = 0.49, id = "und")
#und.field_map.units = "mm"
#und.ax = 0.05
#M1 = Monitor(id = "m1")
#H1 = Hcor(l = 0.0, angle = 0.00, id = "H1")
#V2 = Vcor(l = 0.0, angle = 0.00, id = "V2")


superperiod = ( D9,Q6,D8,Q5,D7,Q4,D6,B1,B2,D5,Q3,D5,B2,B1,D4,SD,D2,Q2,D3,Q1,D2,SF,D1,D1,SF, D2,Q1,D3, Q2,D2,SD,D4,B1,B2,D5,Q3,D5,B2,B1,D6,Q4,D7,Q5,D8,Q6,D9)



m1 = Monitor(eid="start")
m2 = Monitor(eid="end")

dba = (m1,superperiod,m2)

beam = Beam()
beam.E = 14.0
beam.sigma_E = 0.002
beam.emit_xn = 0.4e-6 
beam.emit_yn = 0.4e-6 
beam.gamma_rel = beam.E / (0.511e-3)
beam.emit_x = beam.emit_xn / beam.gamma_rel
beam.emit_y = beam.emit_yn / beam.gamma_rel
beam.beta_x = 4.03267210229
beam.beta_y = 0.573294833302
beam.alpha_x = 0.
beam.alpha_y = 0.



lat = MagneticLattice(dba)
tw0 = Twiss(beam)
tws=twiss(lat, Twiss(), nPoints = 1000)
print( tws[0].beta_x)
print( tws[0].beta_y)
print( tws[0].alpha_x)
print( tws[0].alpha_y)
constr = {'end':{'Dx':0.0, 'Dxp':0.0}, 'D1':{'Dx':0.55, 'Dxp':0.}}
#constr = {'end':{'Dx':0.0, 'Dxp':0.0}, 'start':{'beta_x':15.0, 'beta_y':30.0}}
#constr = {'end':{'Dx':0.0, 'Dxp':0.0}}


vars = [Q3]
#vars = [q1,q2,q5,q6, [tw0, 'beta_x'], [tw0, 'beta_y'], [tw0, 'alpha_x'], [tw0, 'alpha_y']]
#vars = [q1,q2,q3,q5,q6, [tw0, 'beta_'], [tw0, 'beta_y'], [tw0, 'alpha_x'], [tw0, 'alpha_y']]

match(lat, constr, vars, tw0)


#tws=twiss(lat, tw0, nPoints = 1000)
tws=twiss(lat, Twiss(), nPoints = 1000)
s = [p.s for p in tws]
f=plt.figure()
ax = f.add_subplot(211)
ax.set_xlim(0, lat.totalLen)

f.canvas.set_window_title('Betas [m]') 
p1, = plt.plot(s, [p.beta_x for p in tws], lw=2.0)
p2, = plt.plot(s, [p.beta_y for p in tws], lw=2.0)
plt.grid(True)

ax.twinx()
p3,=plt.plot(s, [p.Dx for p in tws], 'r',lw=2.0)

plt.legend([p1,p2,p3], [r'$\beta_x$',r'$\beta_y$', r'$D_x$'])

ax2 = f.add_subplot(212)
plot_lattice(lat, ax2, alpha=0.5)

# add beam size (arbitrary scale)



scale = 1000
#[np.sqrt(p.beta_x*beam.emit_x + (p.Dx*beam.sigma_E)**2) for p in tws]
#[np.sqrt(p.beta_y*beam.emit_y) for p in tws]
sig_x = scale * np.array([np.sqrt(p.beta_x*beam.emit_x + (p.Dx*beam.sigma_E)**2) for p in tws]) # 0.03 is for plotting same scale
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
