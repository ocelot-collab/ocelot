"""
S Tomin
"""
import sys
path = sys.path[0]
indx = path.find("ocelot")
sys.path.append(path[:indx])


from ocelot import *
from ocelot.gui.accelerator import *
import time
from ocelot.optimizer.mint.xfel_interface import XFELMachineInterface
from ocelot.test.xfel_lat.I1B1 import *
#from ocelot.test.xfel_lat.B1D import *
mi = XFELMachineInterface()

def get_x(bpm):
    #x = mi.get_value("XFEL.DIAG/ORBIT/" + bpm + "/X.SA1")*0.001
    return (np.random.rand(1)[0]-0.5)*0.001

def get_y(bpm):
    #y = mi.get_value("XFEL.DIAG/ORBIT/" + bpm + "/Y.SA1") * 0.001
    return (np.random.rand(1)[0]-0.5)*0.001

def read_bpm(lat2):
    S = []
    X = []
    Y = []
    L = 0
    for elem in lat2.sequence:
        L += elem.l
        if elem.__class__ == Monitor:
            x = get_x(elem.id)
            y = get_y(elem.id)
            S.append(L)
            X.append(x)
            Y.append(y)
    return np.array(S), np.array(X), np.array(Y)


def get_energy():
    #E = mi.get_value('XFEL.DIAG/BEAM_ENERGY_MEASUREMENT/LH/ENERGY.SA1')
    return 130.
    #return mi.get_value('XFEL.DIAG/BEAM_ENERGY_MEASUREMENT/I1T/ENERGY.SA1')

# undu_49_i1.lperiod = 0.074
# undu_49_i1.nperiods = 10
# undu_49_i1.Kx = 1.935248
# undu_49_i1.l = 0.74
# Undulator(lperiod=0.074, nperiods=10, Kx=1.935248, Ky=0.0, eid='UNDU.49.I1')
lat = MagneticLattice(cell)

tws0 = Twiss()
tws0.E = 0.005000000
tws0.beta_x = 29.171000000
tws0.alpha_x = 10.955000000
tws0.beta_y = 29.171000
tws0.alpha_y = 10.955000
#L = 0
#for elem in lat.sequence:
#    L += elem.l
#    print(elem.id, L)
#    if elem.id == 'C3.AH1.1.8.I1':
#    	break

print("test")
#exit(0)
tws = twiss(lat, tws0)
plot_opt_func(lat, tws, top_plot=["Dx", "Dy"])
print("sdfsdf")
plt.show()
print("sdfsdf")
#bpmf_95_i1
#bpms_99_i1
# bpmf_103_i1
lat3 = MagneticLattice(cell, start=bpmf_95_i1, stop=bpms_99_i1)

R = lattice_transfer_map(lat3, energy=0.130)
print(R)


Ro11 =  -4.43466901e+00
Ro12 = 5.21261067e+00
    
Ri11 = -1.43615226e+00 
Ri12 = 2.61734725e+00
Ri16 = 1.82166844e-01

lat2 = MagneticLattice(cell, start=c3_ah1_1_8_i1, stop=c_a2_1_1_l1)


plt.ion()
fig, ax = plot_API(lat2)
ax.set_ylim(-1, 1)

timeout = 0.5
E = []
S, X, Y = read_bpm(lat2)

linex, = ax.plot(S, np.zeros_like(S), "ro-") 
liney, = ax.plot(S, np.zeros_like(S), "bo-")

Xref = np.array(X)
Yref = np.array(Y)
xa = get_x(bpmf_95_i1.id)
xc = get_x(bpmf_103_i1.id)
delta1 = (Ri11 * xa - Ri12 * (xc - Ro11*xa)/Ro12)/Ri16
E.append(get_energy())
time.sleep(timeout)
for i in range(40):
    xa = get_x(bpmf_95_i1.id)
    xc = get_x(bpmf_103_i1.id)
    delta = (Ri11 * xa - Ri12 * (xc - Ro11*xa)/Ro12)/Ri16
    print("delta=", delta)
    
    E2 = get_energy()
    E.append(E2)
    E0 = np.mean(E)
    S, X, Y = read_bpm(lat2)
    dx = np.zeros_like(Xref)
    E1 = E0*(delta + 1.)
    E2 = E0*(delta1 + 1.)
    #dy = (Yref/delta1 + Y/delta)/2.
    #dx = E0*(Xref - X)/(E1 - E2)*0.001
    dy = E0*(Yref - Y)/(E1 - E2)
    linex.set_ydata(dx)
    liney.set_ydata(dy)
    Xref = np.array(X)
    Yref = np.array(Y)
    delta1 = delta
    fig.canvas.draw()
    time.sleep(timeout)
#print(d)  
#XFEL.DIAG/BEAM_ENERGY_MEASUREMENT/LH/ENERGY.SA1
#energy = mi.get_value('XFEL.DIAG/BEAM_ENERGY_MEASUREMENT/I1T/ENERGY.SA1')

plt.show()
exit(0)
for i in range(100):
    energy = mi.get_value('XFEL.DIAG/BEAM_ENERGY_MEASUREMENT/LH/ENERGY.SA1')
    print(energy)
    time.sleep(0.1)
