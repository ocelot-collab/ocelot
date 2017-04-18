__author__ = 'Sergey Tomin'
from ocelot.adaptors.astra2ocelot import *
import numpy as np
from ocelot.cpbd.wake3D import *
from matplotlib.pyplot import *
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}

matplotlib.rc('font', **font)


def csrtrack2astra(infile, orient="H"):
    #%H z x y pz px py -> x y z px py pz
    #%V z y x pz py px -> x y -z px py -pz
    #if nargin<2, orient='H'; end;
    PD = np.loadtxt(infile)
    n = len(PD[:, 0])-1
    Q = sum(PD[1:n+1, 6])*1e9
    t0 = PD[0, 0]
    PD1 = np.zeros((n-1, 6))
    if orient == 'H':
        PD1[:, 0] = PD[1:n+1, 2-1]
        PD1[:, 1] = PD[1:n+1, 3-1]
        PD1[:, 2] = PD[1:n+1, 1-1]
        PD1[:, 3] = PD[1:n+1, 5-1]
        PD1[:, 4] = PD[1:n+1, 6-1]
        PD1[:, 5] = PD[1:n+1, 4-1]
    else:
        PD1[:, 0] =-PD[1:n+1, 3-1]
        PD1[:, 1] = PD[1:n+1, 2-1]
        PD1[:, 2] = PD[1:n+1, 1-1]
        PD1[:, 3] =-PD[1:n+1, 6-1]
        PD1[:, 4] = PD[1:n+1, 5-1]
        PD1[:, 5] = PD[1:n+1, 4-1]
    #%add a reference particle
    for i in range(6):
        PD1[1:n, i] = PD1[1:n, i] + PD1[0, i]
    return PD1, Q


ocelot_file = 'c:\\Users\\tomins\\Documents\\Dropbox\\XFEL\\Dohlus\\Igor\\ocelot.ast'

xtrack_file = 'c:\\Users\\tomins\\Documents\\Dropbox\\XFEL\\Dohlus\\Igor\\PP_XCode\\xcode.ast'

p_array_o = astraBeam2particleArray(ocelot_file)

p_array_x = astraBeam2particleArray(xtrack_file)

I_o = s2current(p_array_o.particles[4::6], p_array_o.q_array, 300, filter_order=5, mean_vel=speed_of_light)
I_x = s2current(p_array_x.particles[4::6], p_array_x.q_array, 300, filter_order=5, mean_vel=speed_of_light)
figure(1)
plot(-I_o[:, 0]*1000, I_o[:, 1], "r", lw=2)
plot(-I_x[:, 0]*1000, I_x[:, 1], "b", lw=2)
xlabel("s, mm")
ylabel("I, A")
legend(["ocelot", "xtrack"], loc=2)
grid(True)


figure(2)
plot(-p_array_o.particles[4::6]*1000, p_array_o.particles[5::6], "r.", lw=2)
plot(-p_array_x.particles[4::6]*1000, p_array_x.particles[5::6], "b.", lw=2)
xlabel("s, mm")
ylabel("dE/E")
grid(True)
legend(["ocelot", "xtrack"], loc=1)


show()
infile = 'c:\\Users\\tomins\\Documents\\Dropbox\\XFEL\\Dohlus\\_CSRtrack_comp\\test_bc0\\N5_BC3\\out\\out.fmt1'
pd, q = csrtrack2astra(infile, orient="H")

