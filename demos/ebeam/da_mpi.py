__author__ = 'Sergey Tomin'
import sys
ind = sys.path[0].find("ocelot")
sys.path.append(sys.path[0][:ind])

from ocelot import *
from time import time
import numpy as np
from mpi4py import MPI


mpi_comm = MPI.COMM_WORLD
size = mpi_comm.Get_size()
rank = mpi_comm.Get_rank()


Q1 = Quadrupole(l=0.4, k1=-1.3, eid= "Q1")
Q2 = Quadrupole(l=0.8, k1=1.4, eid= "Q2")
Q3 = Quadrupole(l=0.4, k1=-1.7, eid= "Q3")
Q4 = Quadrupole(l=0.5, k1=1.19250444829 , eid= "Q4")

B  = Bend(l=2.7, k1=-.06, angle=2*pi/16., e1=pi/16., e2=pi/16., eid= "B")

SF = Sextupole(l=0.01, k2 = 5.8914775395193555*100, eid= "SF") #random value
SD = Sextupole(l=0.01, k2 = -6.8036102026266558*100, eid= "SD") #random value

D1 = Drift(l=2., eid = "D1")
D2 = Drift(l=0.6, eid = "D2")
D3 = Drift(l=0.3, eid = "D3")
D4 = Drift(l=0.7, eid = "D4")
D5 = Drift(l=0.9, eid = "D5")
D6 = Drift(l=0.2, eid = "D6")


cell = (D1, Q1, D2, Q2, D3, Q3, D4, B, D5, SD, D5, SF, D6, Q4, D6, SF, D5, SD,D5, B, D4, Q3, D3, Q2, D2, Q1, D1)
ring = cell
method = MethodTM()
method.params[Sextupole] = KickTM
lat = MagneticLattice(ring, method=method)



#compensate_chromaticity(lat, ksi_x_comp = 0, ksi_y_comp = 0,  nsuperperiod = 8)

nturns = 1000
nx = 80
ny = 80

x_array = np.linspace(-0.03, 0.03, nx)
y_array = np.linspace(0.0001, 0.03, ny)
start = time()
pxy_list = create_track_list(x_array, y_array, p_array=[0.])
print("stop")
pxy_list = track_nturns_mpi( mpi_comm,lat, nturns, pxy_list,  nsuperperiods = 8, save_track=False)
if rank == 0:
    print( time() - start)
    da = np.array([ pxy.turn for pxy in pxy_list])
    np.savetxt("da.txt", (da))
    b = []
    for x, y in zip(x_array, y_array):
        a = [x, y]
        b.append(np.array([x,y]))
    np.savetxt("da_axis.txt", np.array(b))
    #from ocelot.gui.accelerator import *
    #show_da(da, x_array, y_array)
