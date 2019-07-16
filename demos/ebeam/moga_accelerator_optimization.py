"""
MOGA optimization example

Simultaneously the beam emittance minimization and the dynamic aperture maximization of the storage ring
using MOGA NSGA2 algorithm

Launch this script as
python3 moga_accelerator_optimization.py                     - for MOGA single thread
mpiexec -n 10 python3 moga_accelerator_optimization.py       - for MOGA MPI with 10 threads
"""

from ocelot import *
from ocelot.cpbd.optics import *
from ocelot.cpbd.moga import *
from ocelot.cpbd.chromaticity import *


#  Lattice of Kurchatov Light Sourse "Siberia - 2"
# ********************************* START ********************************************
D1 = Drift (l = 1.49, eid = "D1")
D2 = Drift (l = 0.1035, eid = "D2")
D3 = Drift(l = 0.307, eid = "D3")
D4 = Drift(l = 0.33, eid = "D4")
D5 = Drift(l = 0.3515, eid = "D5")
D6 = Drift(l = 0.3145, eid = "D6")
D7 = Drift(l = 0.289, eid = "D7")
D8 = Drift(l = 0.399, eid = "D8")
D9_05 = Drift(l = 1.5045, eid = "D9_05")

SF = Sextupole(l = 0.001, k2 = 1.7673786254063251, eid = "SF")
SD = Sextupole(l = 0.001, k2 = -3.6169817233025707, eid = "SD")

Q1 = Quadrupole(l = 0.2931, k1 = 2.62, eid = "Q1")
Q2 = Quadrupole(l = 0.2931, k1 = -3.1, eid = "Q2")
Q3 = Quadrupole(l = 0.3268, k1 = 2.8, eid = "Q3")
Q4 = Quadrupole(l = 0.2910, k1 = -3.7, eid = "Q4")
Q5 = Quadrupole(l = 0.3912, k1 = 4.0782, eid = "Q5")
Q6 = Quadrupole(l = 0.2910, k1 = -3.534859, eid = "Q6")

B1 = SBend(l = 0.23, angle = 0.23/19.626248, eid  = "B1")
B2 = SBend(l = 1.227, angle = 1.227/4.906312, eid = "B2")

M1 = Monitor(eid = "m1")
M2 = Monitor(eid = "m2")

superperiod = (M1,D1,SF,D2,Q1,D3,Q2,D2,SD,D4,B1,B2,D5,Q3,D5,B2,B1,D6,Q4,D7,Q5,D8,Q6,D9_05,M2,D9_05,Q6,D8,Q5,D7,Q4,D6,B1,B2,D5,Q3,D5,B2,B1,D4,SD,D2,Q2,D3,Q1,D2,SF,D1,M1)
# ********************************* END ********************************************


# create lattice
method = MethodTM()
method.params[Sextupole] = KickTM
method.global_method = TransferMap

lattice = MagneticLattice(superperiod,  method=method)


# setting of infinit and zero values
INF = float("inf")
ZERO = 0.0

# setting of boundaries for all variable
bounds = ((0., 6.0), (-6.0, 0.), (0., 6.0), (-6.0, 0.), (0., 6.0), (-6.0, 0.))

# setting weights for fitness function
# example for one fitness function minimization: (-1.0,)        - comma is important
# example for two fitness functions minimization: (-1.0, -1.0)
weights = (-1.0, -1.0)

# setting of initial solutions (optional)
init_pop = []
init_pop.append((2.62, -3.1, 2.8, -3.7, 4.0782, -3.534859))
init_pop.append((2.62, -3.08, 3.6, -2.5, 3.09, -1.84))
init_pop.append((2.865, -3.16, 4.48, -3.74, 4.1, -3.39))

# setting of additional arguments for the fitness function
args = [[Q1, Q2, Q3, Q4, Q5, Q6], lattice]

# definition of fintess function
# in this example fitness function return emittance and dynamic aperture
def fit_func(x0, iter_data, args):
    """
    x0 contains variables
    iter_data contains information about current step (number_of_current_individual, number_of_current_iteration)
    args contains additional arguments

    fitness function have to return iterable data

    example for one fintess function optimization
    return (f1_result, )              - comma is important

    example for two fintess functions optimization
    return (f1_result, f2_result)
    """
    # update lattice for the new candidate solution
    lattice = args[1]
    
    vars = args[0]
    for i in range(len(args[0])):
        
        vars[i].k1 = x0[i]
        vars[i].transfer_map = lattice.method.create_tm(vars[i])
    
    beam = Beam()
    beam.E = 2.5

    # additional restriction for a periodic solution existence
    tw = Twiss()
    tw = periodic_twiss(tw, lattice_transfer_map(lattice, beam.E))
    if tw == None:
        return (INF, INF)

    # beam emittance calculation (fitness function parameter 1)
    eb = EbeamParams(lattice, Twiss(beam), nsuperperiod=6)
    
    # additional restriction for the beam emittance
    if eb.emittance < 0.0 or eb.emittance > 100.0e-9:
        return (INF, INF)

    # cromaticity compensation
    tws = twiss(lattice)
    ksi = natural_chromaticity(lattice, tws[0], nsuperperiod=6)
    compensate_chromaticity(lattice, ksi_x_comp=0, ksi_y_comp=0, nsuperperiod=6)

    # dynamic aperture calculation (fitness function parameter 2)
    nturns = 5
    nx = 100
    ny = 50

    x_array = np.linspace(-0.1, 0.1, nx)
    y_array = np.linspace(0.0001, 0.05, ny)
    
    pxy_list = create_track_list(x_array, y_array, p_array=[0.])
    pxy_list = track_nturns(lattice, nturns, pxy_list, nsuperperiods=6, save_track=True)

    da = np.array([pxy.turn for pxy in pxy_list])

    da_num = 1.0
    for i in da:
        if i >= nturns-1:
            da_num += 1.0

    # additional restriction for dynamic aperture
    if da_num < 100:
        return (eb.emittance, INF)
    
    return (eb.emittance, 1.0/da_num)

#
# MOGA optimization
#

# init MOGA
opt = Moga(bounds=bounds, weights=weights)

# set additionla parameters (optional)
opt.set_params(n_pop=10, n_gen=50)

# start optimization
result = opt.nsga2(fit_func=fit_func, fit_func_args=args, init_pop=init_pop)
