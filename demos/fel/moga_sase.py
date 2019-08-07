'''
Example for optimization parameters for FEL time dependent simulation

Launch this script as 
python3 moga_fit_new.pyplot        - for MOGA single thread

don't launch this script as MOGA MPI
because Genesis MPI is recomenned for FEL time dependent simulation
'''

import sys, os
#import matplotlib.pyplot as plt
#from ocelot.gui.accelerator import *
from ocelot.utils.xfel_utils import *
#from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.gui.genesis_plot import *
import time
#rom ocelot.common.globals import *  #import of constants like "h_eV_s" and "speed_of_light"
from ocelot.rad.undulator_params import *

from ocelot.cpbd.moga import *


# Setting of path parameters
param_dir = '/data/netapp/xfel/svitozar/projects/XFEL/parameters/'
exp_dir = '/data/netapp/xfel/yevgeniy/moga_test/sase_new_results/' #directory where the output will be stored
beam_fileName = param_dir + 'beams/beam_0.02nC_aftSASE1_gianluca.txt' #path to beam file

sys.path.append(param_dir + 'sase1/') 
from sase1 import * #import sase1 undulator parameters


# create experimental directory with iter_# and run_# subdirectories
if MPI_RANK == 0:
    if not os.path.exists(exp_dir): os.makedirs(exp_dir)
else:
    time.sleep(1.0)


# Setting of some parameters for undulator
beta_av = 25.0 # beta function for the ebeam in undulator
E_beam = 14.0   # Electon beam energy (e.g.14.4) [GeV]
E_photon = 8265.0   # FEL Resonance photon energy e.g. 8000.0 [eV]


# Range of statistical runs to execute, from 0 to n-1 - (0,n)
#run_ids = range(0,10) 


# Length of undulator 
#N1 = 11
# tapering:
#n0 = [0,5,60]
#a1 = [0.0, -0.0012]
#a2 = [0.0, -0.000]




#zsep=0 #z-separation between slices in radaition wavelengths (if 0 - then calculated based of beam peak current parameters)




#launcher = get_genesis_launcher()






# read and prepare beam file
beam =read_beam_file(beam_fileName)
# beam = cut_beam(beam,[-10e-6, 10e-6])
beam=set_beam_energy(beam, E_beam)
beam_pk=get_beam_peak(beam) #another object containing beam parameters at peak current position, get_beam_s() gets parameters at given position s


# Main undulator lattice preparation
lat = MagneticLattice(sase1_segment(n=20))
und_l=l_fodo # undulator cell length
und.Kx = Ephoton2K(E_photon, und.lperiod, E_beam)
up = UndulatorParameters(und,beam_pk.E)

rematch(beta_av, l_fodo, qdh, lat, extra_fodo, beam_pk, qf, qd) # jeez...
beam =transform_beam_twiss(beam,transform=[ [beam_pk.beta_x,beam_pk.alpha_x], [beam_pk.beta_y,beam_pk.alpha_y] ])
#plot_beam(beam,showfig=0,savefig=savefig)


#a0 = und.Kx
#taper_func = lambda n : f1(n, n0, a0, a1, a2 )

#lat1 = deepcopy(lat)
#lat1 = taper(lat1, taper_func)

# possibility to override the automatic zsep calculations (based on rho)
zsep = generate_input(up, beam_pk, itdp=True).zsep





###    
# Parameters setting for MOGA optimizarion
###


# Setting number of fitness function optimization
'''
Examples:
weights = (-1.0, -1.0)  # for two fitness functions minimization
weights = (-1.0,)       # for one fitness functions minimization - comma is important
'''
weights = (-1.0,)

# Setting of boundaries for all variable
'''
Example:
bounds = ((var_1_min, var_1_max), (var_2_min, var_2_max), (var_3_min, var_3_max))
'''
bounds = ((6,13), (und.Kx*0.999,und.Kx*1.001), (-0.005, 0.0), (1.0,1.5))
'''
here bounds mean:
1) n0
2) a0
3) a1
4) a2
'''

# Setting of additional arguments for the fitness function
args = [beam_pk,E_beam,up,lat,exp_dir,zsep,beam]

# Setting of initial solutions (optional)
init_pop = []
init_pop.append((np.random.uniform(6,13), np.random.uniform(und.Kx*0.999,und.Kx*1.001), np.random.uniform(-0.005, 0.0), np.random.uniform(1.0,1.5)))

# Fitness function definition
def fit_func(x0, iter_data, args):
    
    # x0 contains variables
    # iter_data contains information about current step (number_of_current_individual, number_of_current_iteration)
    # args contains additional arguments
    #
    # fitness function have to return iterable data
    #
    # example for one fintess function optimization
    # return (f1_result, )              - comma is important
    #
    # example for two fintess functions optimization
    # return (f1_result, f2_result)


    n0 = int(x0[0])
    a0 = x0[1]
    a1 = x0[2]
    a2 = x0[3]
    
    # Genesis simulation parameters
    npart = 8192 # number of macroparticles in Genesis
    dgrid= 6e-4
    ncar=101

    zstop = und_l*11
        
    beam_pk = args[0]  # e.g. 14.0
    E_ev   = args[1]   # e.g. 8000.0  
    up=args[2]
    lat = args[3]
    exp_dir = args[4] + 'iter_' + str(iter_data[1]) + '/'
    zsep = args[5]
    beam = args[6]
    
    run_id = iter_data[0]

    
    launcher = get_genesis_new_launcher()

    taper_func = lambda n : f1(n, n0, a0, a1, a2)
    lat1 = deepcopy(lat)
    #lat1 = taper(lat, taper_func)
    
    inp = generate_input(up, beam_pk, itdp=True)

    inp.stageid=1
    inp.runid = run_id
    inp.exp_dir = exp_dir
    inp.ipseed = 6123
        
    inp.lat=lat1
    inp.beam=beam
        
    inp.iphsty = 8 # Generate output in the main output file at each IPHSTYth integration step. To disable output set IPHSTY to zero. 
    inp.ishsty = 1 # Generate output in the main output file for each ISHSTYth slice. 
        
    inp.isravg = 0 #we assume it is compensated by proper tapering on top of ours
    inp.isrsig = 1 #use it wisely, due to causality problems (affects initial parameters?)
    inp.ncar = ncar
    inp.zstop = zstop
    inp.dgrid = dgrid
    inp.zsep=zsep
        
    inp.idmppar = 1
    inp.idmpfld = 1
        
    out = run_genesis(inp, launcher, read_level=2, assembly_ver='sys')
        
    return (1.0/out.p_int[0][-1],)



t0 = time.time()

###
# MOGA optimization
###

# init MOGA
opt = Moga(bounds=bounds, weights=weights)

# set additionla parameters (optional)
opt.set_params(n_gen=10, log_file=exp_dir+'moga_result.dat', plt_file=exp_dir+'moga_plot.dat')
# n_pop - numper of population
# n_gen - numper of iterations
# log_file - name with full path for log file
# plt_file - name with full path for file with data for Pareto frontier plotting of during optimiaztion 

# start optimization
results = opt.anneal(fit_func=fit_func, fit_func_args=args, init_pop=init_pop)

    
t1 = time.time()
print ('total simulation time = %.2f min' %((t1-t0)/60))
