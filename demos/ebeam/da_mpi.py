__author__ = 'Sergey Tomin'
import sys
ind = sys.path[0].find("ocelot")
sys.path.append(sys.path[0][:ind])
#from ocelot.cpbd.match import *
#from ocelot.gui.accelerator import *
from ocelot.cpbd.elements import *
from ocelot.cpbd.optics import *
from ocelot.cpbd.e_beam_params import *
from ocelot.cpbd.chromaticity import *
from ocelot.cpbd.track import *
from mpi4py import MPI


mpi_comm = MPI.COMM_WORLD
size = mpi_comm.Get_size()
rank = mpi_comm.Get_rank()


Q1 = Quadrupole(l= 0.4, k1=-1.3, id = "Q1")
Q2 = Quadrupole(l= 0.8, k1=1.4, id = "Q2")
Q3 = Quadrupole(l= 0.4, k1=-1.7, id = "Q3")
Q4 = Quadrupole(l= 0.5, k1=1.19250444829 , id = "Q4")

B  = Bend(l=2.7, k1=-.06, angle=2*pi/16., e1=pi/16., e2=pi/16., id = "B")

SF = Sextupole(l=0., ms = 5.8914775395193555, id = "SF") #random value
SD = Sextupole(l=0., ms =  -6.8036102026266558, id = "SD") #random value

D1 = Drift(l=2., id = "D1")
D2 = Drift(l=0.6, id = "D2")
D3 = Drift(l=0.3, id = "D3")
D4 = Drift(l=0.7, id = "D4")
D5 = Drift(l=0.9, id = "D5")
D6 = Drift(l=0.2, id = "D6")


cell = (D1, Q1, D2, Q2, D3, Q3, D4, B, D5, SD, D5, SF, D6, Q4, D6, SF, D5, SD,D5, B, D4, Q3, D3, Q2, D2, Q1, D1)
ring = cell

lat = MagneticLattice(ring)



#compensate_chromaticity(lat, ksi_x_comp = 0, ksi_y_comp = 0,  nsuperperiod = 8)

nturns = 1000
nx = 80
ny = 80

x_array = linspace(-0.03, 0.03, nx)
y_array = linspace(0.0001, 0.03, ny)
start = time()
pxy_list = create_track_list(x_array, y_array, p_array=[0.])
pxy_list = tracking_mpi( mpi_comm,lat, nturns, pxy_list,  nsuperperiods = 8, order = 3, save_track=False)
if rank == 0:
    print( time() - start)
    da = array(map(lambda pxy: pxy.turn, pxy_list))
    np.savetxt("da.txt", (da))
    b = []
    for x, y in zip(x_array, y_array):
        a = [x, y]
        b.append(array([x,y]))
    np.savetxt("da_axis.txt", array(b))
    #from ocelot.gui.accelerator import *
    #show_da(da, x_array, y_array)
