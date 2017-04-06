from ocelot import *
from ocelot.gui.accelerator import *
from ocelot.adaptors import *
from ocelot.cpbd.moga import *

from scipy import stats
import numpy as np
from time import *
import copy

from ocelot.test.I2 import *


phi1 = 18.7268
V1 = 18.50662e-3 / np.cos(phi1 * np.pi / 180.0)

c_a1_1_1_i1.v = V1; c_a1_1_1_i1.phi = phi1
c_a1_1_2_i1.v = V1; c_a1_1_2_i1.phi = phi1
c_a1_1_3_i1.v = V1; c_a1_1_3_i1.phi = phi1
c_a1_1_4_i1.v = V1; c_a1_1_4_i1.phi = phi1
c_a1_1_5_i1.v = V1; c_a1_1_5_i1.phi = phi1
c_a1_1_6_i1.v = V1; c_a1_1_6_i1.phi = phi1
c_a1_1_7_i1.v = V1; c_a1_1_7_i1.phi = phi1
c_a1_1_8_i1.v = V1; c_a1_1_8_i1.phi = phi1


phi13 = 180.0
V13 = -20.2e-3 / 8.0 / np.cos(phi13 * np.pi / 180.0)

c3_ah1_1_1_i1.v = V13; c3_ah1_1_1_i1.phi = phi13
c3_ah1_1_2_i1.v = V13; c3_ah1_1_2_i1.phi = phi13
c3_ah1_1_3_i1.v = V13; c3_ah1_1_3_i1.phi = phi13
c3_ah1_1_4_i1.v = V13; c3_ah1_1_4_i1.phi = phi13
c3_ah1_1_5_i1.v = V13; c3_ah1_1_5_i1.phi = phi13
c3_ah1_1_6_i1.v = V13; c3_ah1_1_6_i1.phi = phi13
c3_ah1_1_7_i1.v = V13; c3_ah1_1_7_i1.phi = phi13
c3_ah1_1_8_i1.v = V13; c3_ah1_1_8_i1.phi = phi13


lattice = gun_5MeV + i1_150M

p_array_init = astraBeam2particleArray(filename='workshop/beam_6MeV.ast')

n = 100
p_array_init.particles = p_array_init.particles[::n]
p_array_init.q_array = p_array_init.q_array[::n]

s_start = 3.2
p_array_init.s = 0.0

method = MethodTM()
method.global_method = SecondTM

stop_element = i1_end_i1
lat = MagneticLattice(lattice, start=start_sim, stop=stop_element, method=method)

bounds = ((0.0, 30.0), (160.0, 200.0))

init_pop = []
init_pop.append([c_a1_1_1_i1.phi, c3_ah1_1_1_i1.phi])

v1 = [c_a1_1_1_i1, c_a1_1_2_i1, c_a1_1_3_i1, c_a1_1_4_i1, c_a1_1_5_i1, c_a1_1_6_i1, c_a1_1_7_i1, c_a1_1_8_i1]
v2 = [c3_ah1_1_1_i1, c3_ah1_1_2_i1, c3_ah1_1_3_i1, c3_ah1_1_4_i1, c3_ah1_1_5_i1, c3_ah1_1_6_i1, c3_ah1_1_7_i1, c3_ah1_1_8_i1]
args = [[v1, v2], lat, p_array_init]


def fit_func(x0, args):

    p_array = copy.deepcopy(args[2])
    navi = Navigator(args[1])
    navi.unit_step = 0.02

    vars = args[0][0]
    for i in range(len(vars)):
        vars[i].phi = x0[0]
        vars[i].v = 18.50662e-3 / np.cos(x0[0] * np.pi / 180.0)
        vars[i].transfer_map = args[1].method.create_tm(vars[i])

    vars = args[0][1]
    for i in range(len(vars)):
        vars[i].phi = x0[1]
        vars[i].v = 18.50662e-3 / np.cos(x0[1] * np.pi / 180.0)
        vars[i].transfer_map = args[1].method.create_tm(vars[i]) 


    tws_track, p_array = track(args[1], p_array, navi)

    sI1, I1 = get_current(p_array, charge=p_array.q_array[0], num_bins=25)
    tau = p_array.particles[4::6]
    dp = p_array.particles[5::6]
    slope, intercept, r_value, p_value, std_err = stats.linregress(-tau*1000, dp)

    f1 = 1.0 / np.max(I1)
    f2 = abs(slope)
    print(1./f1, f2)
    return f1, f2



opt = Moga(bounds=bounds)
opt.set_params(pop_num=100, ngen=10)
result = opt.nsga2(fit_func=fit_func, fit_func_args=args, init_pop=init_pop)
print(result)