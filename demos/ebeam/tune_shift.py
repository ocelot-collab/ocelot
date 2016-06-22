__author__ = 'Sergey Tomin'

from ocelot.gui import *
from ocelot import *
import numpy as np
from time import time
from ocelot.cpbd.track import *
from ocelot.rad.radiation_py import und_field

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

SF = Sextupole(l = 0.0001, k2 = 17673.786254063251*0, eid= "SF")
SD = Sextupole(l = 0.0001, k2 =-36169.817233025707*0, eid= "SD")

Q1 = Quadrupole (l = 0.293, k1 = 2.62, eid= "Q1")
Q2 = Quadrupole (l = 0.293, k1 = -3.1, eid= "Q2")
Q3 = Quadrupole (l = 0.327, k1 = 2.8, eid= "Q3")
Q4 = Quadrupole (l = 0.291, k1 = -3.7, eid= "Q4")
Q5 = Quadrupole (l = 0.391, k1 = 4.0782, eid= "Q5")
Q6 = Quadrupole (l = 0.291, k1 = -3.534859, eid= "D6")

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
method.params[Undulator] = UndulatorTestTM
#method.params[Undulator] = RungeKuttaTM
method.global_method = TransferMap
lat = MagneticLattice(ring, method=method)
beam = Beam()
beam.E = 2. #GeV
beam.I = 0.1 #A
tw0 = Twiss(beam)
tws = twiss(lat,tw0, nPoints=1000)
plot_opt_func(lat, tws)
plt.show()
#compensate_chromaticity(lat, ksi_x_comp=0, ksi_y_comp=0,  nsuperperiod=1)

nturns = 2048*2
#nx = 100
#ny = 80

#x_array = np.linspace(-0.03, 0.03, nx)
#y_array = np.linspace(0.0001, 0.03, ny)
start = time()
#pxy_list = create_track_list(x_array, y_array, p_array=[0.])
pxy_list = [Track_info(Particle(y=0.001, E=2), 0.00, 0.0001)]

pxy_list = track_nturns( lat, nturns, pxy_list,  nsuperperiods = 1, save_track=True)
#print([p.turn for p in pxy_list])
y = [p[2] for p in pxy_list[0].p_list]
py = [p[3] for p in pxy_list[0].p_list]

print("time exec = ", time() - start)
pxy_list = freq_analysis(pxy_list, lat, nturns, harm = True)
print(pxy_list[0].muy)
sp = np.fft.fft(y)
freq = np.fft.fftfreq(len(y))

plt.plot( sp.real, sp.imag)
plt.show()

def dft(sample, freqs):
    transf = np.zeros(len(freqs), dtype=np.complex)
    n = len(sample)
    for i, ai in enumerate(sample):
        transf += ai*np.exp(-2*np.pi*i/n*1j*freqs)
    return transf

#transf = dft(y, np.arange(0, ))

#f = dft(y, np.linspace(0.27, 0.33, 500))
plt.plot(freq, sp.real, freq, sp.imag)
plt.show()
plt.plot(y, py, "r.")
plt.show()
#da = np.array([pxy.turn for pxy in pxy_list])
#show_da(da, x_array, y_array)

#da_contr = contour_da(pxy_list, nturns)

#da_mux = np.array([pxy.mux for pxy in pxy_list])
#da_muy = np.array([pxy.muy for pxy in pxy_list])
#show_mu(da_contr, da_mux, da_muy, x_array, y_array)