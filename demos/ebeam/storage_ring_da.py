__author__ = 'Sergey Tomin'

from ocelot.gui import *
from ocelot import *
import numpy as np
from time import time

Q1 = Quadrupole(l= 0.4, k1=-1.3, eid= "Q1")
Q2 = Quadrupole(l= 0.8, k1=1.4, eid= "Q2")
Q3 = Quadrupole(l= 0.4, k1=-1.7, eid= "Q3")
Q4 = Quadrupole(l= 0.5, k1=1.19250444829 , eid= "Q4")

B  = Bend(l=2.7, k1=-.06, angle=2*pi/16., e1=pi/16., e2=pi/16., eid= "B")

SF = Sextupole(l=0.01, k2 = 150, eid= "SF") #random value
SD = Sextupole(l=0.01, k2 =-150, eid= "SD") #random value

D1 = Drift(l=2., eid= "D1")
D2 = Drift(l=0.6, eid= "D2")
D3 = Drift(l=0.3, eid= "D3")
D4 = Drift(l=0.7, eid= "D4")
D5 = Drift(l=0.9, eid= "D5")
D6 = Drift(l=0.2, eid= "D6")


cell = (D1, Q1, D2, Q2, D3, Q3, D4, B, D5, SD, D5, SF, D6, Q4, D6, SF, D5, SD,D5, B, D4, Q3, D3, Q2, D2, Q1, D1)
ring = cell
method = MethodTM()
method.params[Sextupole] = KickTM
method.global_method = TransferMap
lat = MagneticLattice(ring, method=method)
#lat = MagneticLattice(ring, method=MethodTM({Sextupole: KickTM, "global": TransferMap}))

compensate_chromaticity(lat, ksi_x_comp=0, ksi_y_comp=0,  nsuperperiod=8)

nturns = 1000
nx = 100
ny = 80

x_array = np.linspace(-0.03, 0.03, nx)
y_array = np.linspace(0.0001, 0.03, ny)
start = time()
pxy_list = create_track_list(x_array, y_array, p_array=[0.])
pxy_list = track_nturns( lat, nturns, pxy_list,  nsuperperiods = 8, save_track=True)
#print([p.turn for p in pxy_list])

print("time exec = ", time() - start)
pxy_list = freq_analysis(pxy_list, lat, nturns, harm = True)

da = np.array([pxy.turn for pxy in pxy_list])
show_da(da, x_array, y_array)

da_contr = contour_da(pxy_list, nturns)

da_mux = np.array([pxy.mux for pxy in pxy_list])
da_muy = np.array([pxy.muy for pxy in pxy_list])
show_mu(da_contr, da_mux, da_muy, x_array, y_array)