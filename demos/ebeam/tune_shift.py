__author__ = 'Sergey Tomin'

from ocelot.gui import *
from ocelot import *
import numpy as np
from time import time
from ocelot.cpbd.track import *
from ocelot.rad.radiation_py import und_field

D0 = Drift(l=0., eid= "D0")
D1 = Drift(l=1.49, eid= "D1")
D2 = Drift(l=0.1035, eid= "D2")
D3 = Drift(l=0.307, eid= "D3")
D4 = Drift(l=0.33, eid= "D4")
D5 = Drift(l=0.3515, eid= "D5")
D6 = Drift(l=0.3145, eid= "D6")
D7 = Drift(l=0.289, eid= "D7")
D8 = Drift(l=0.399, eid= "D8")
D9 = Drift(l=3.009, eid= "D9")

SF = Sextupole(l=0.0001, k2= 17673.786254063251*0, eid= "SF")
SD = Sextupole(l=0.0001, k2=-36169.817233025707*0, eid= "SD")

Q1 = Quadrupole(l=0.293, k1=2.62, eid= "Q1")
Q2 = Quadrupole(l=0.293, k1=-3.1, eid= "Q2")
Q3 = Quadrupole(l=0.327, k1=2.8, eid= "Q3")
Q4 = Quadrupole(l=0.291, k1=-3.7, eid= "Q4")
Q5 = Quadrupole(l=0.391, k1=4.0782, eid= "Q5")
Q6 = Quadrupole(l=0.291, k1=-3.534859, eid= "D6")

B1 = SBend(l = 0.23, angle = 0.23/19.626248, eid= "B1")
B2 = SBend(l = 1.227, angle = 1.227/4.906312, eid= "B2")


u = Undulator(lperiod=0.02, nperiods=50, Kx=2)
u.mag_field = lambda x, y, z: und_field(x, y, z, u.lperiod, u.Kx)
u.npoints = 500
cell_u = ( D1,SF, D2,Q1,D3, Q2,D2,SD,D4,B1,B2,D5,Q3,D5,B2,B1,D6,Q4,D7,Q5,D8,Q6,
           Drift(l=1.0045), u, Drift(l=1.0045),
           Q6,D8,Q5,D7,Q4,D6,B1,B2,D5,Q3,D5,B2,B1,D4,SD,D2,Q2,D3,Q1,D2,SF,D1)


cell = ( D1,SF, D2,Q1,D3, Q2,D2,SD,D4,B1,B2,D5,Q3,D5,B2,B1,D6,Q4,D7,Q5,D8,Q6,D9,Q6,D8,Q5,D7,Q4,D6,B1,B2,D5,Q3,D5,B2,B1,D4,SD,D2,Q2,D3,Q1,D2,SF,D1)

ring = 3*cell + cell_u + 2*cell

method = MethodTM()
method.params[Sextupole] = KickTM
#method.params[Undulator] = UndulatorTestTM
method.params[Undulator] = RungeKuttaTM
method.global_method = TransferMap
lat = MagneticLattice(ring, method=method)
beam = Beam()
beam.E = 0. #GeV
beam.I = 0.1 #A
tw0 = Twiss(beam)
tws = twiss(lat,tw0, nPoints=1000)
plot_opt_func(lat, tws)
plt.show()
mu_y_no_u = 1 - tws[-1].muy%(2*pi)/(2*pi)
print(mu_y_no_u)
nturns = 2048*1-1

start = time()

pxy_list = [Track_info(Particle(y=0.0001, E=2), 0.00, 0.0001)]

pxy_list = track_nturns( lat, nturns, pxy_list,  nsuperperiods = 1, save_track=True)
y = [p[2] for p in pxy_list[0].p_list]
py = [p[3] for p in pxy_list[0].p_list]

print("time exec = ", time() - start)
pxy_list = freq_analysis(pxy_list, lat, nturns, harm = True)
print("low resolution: mu_y = ", pxy_list[0].muy)


def dft(sample, freqs):
    n = len(sample)
    x_freq = freqs*n
    transf = np.zeros(len(freqs), dtype=np.complex)
    
    for i, ai in enumerate(sample):
        transf += ai*np.exp(-2*np.pi*i/n*1j*x_freq)
    return transf

x = np.linspace(0.3, 0.31, 2048+1)

f = np.abs(dft(y, x))
mu_y_h = x[np.argmax(f)]
print("high resolution: mu_y = ", mu_y_h)

print("simulation: dQy = ", mu_y_no_u - mu_y_h)

# analytical solution
beta = 0.57
E = 2e9
c = 3e8
L = u.lperiod*u.nperiods
Kx = u.Kx
lperiod = u.lperiod
Bx = Kx/93.4/lperiod
#print("Bx = ", Bx)
ro = E/(Bx*c)
#print("ro = ", ro)
delta_Q_y = L/(8*np.pi*ro**2)*(beta + L**2/(12*beta))
print("analytics: dQy = ", delta_Q_y)


"""
beta = 0.56
E = 2e9
c = 3e8
L = 1
Kx = 2
lperiod = 0.02
Bx = Kx/94.3/lperiod
Bx = 1.0709745594714459
print("Bx = ", Bx)
ro = E/(Bx*c)
print("ro = ", ro)
nu = L/(8*np.pi*ro**2)*(beta + L**2/(12*beta))
print("nu = ", nu)
print("lin = ", 0.304169921875 - 0.3034375)
print("1 mm: sim=", 0.304169921875 - 0.303413085937)
print("1 mm: rk = ", 0.304169921875 - 0.303432617187)
print("0.1 mm: sim=", 0.304169921875 - 0.30341796875)
print("0.1 mm: rk = ", 0.304169921875 - 0.3034375)
"""
