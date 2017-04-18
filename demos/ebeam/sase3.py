'''
input deck for XFEL SASE3 beamline
'''
from ocelot import *
import numpy as np

und = Undulator(nperiods=73, lperiod=0.068, Kx=0.0, eid= "und"); voodoo = 1.5
#und = Undulator(nperiods=140, lperiod=0.036, Kx=0.0, id = "und"); voodoo = 0.8



d = Drift (l=1.0, eid= "d")

d1 = Drift (l=0.55, eid= "d1")
#d2 = Drift (l=0.45, id = "d2")
d2 = Drift (l=und.lperiod, eid= "d2")
#d3 = Drift (l=0.27, id = "d0.05nm3")
d3 = Drift (l=und.lperiod, eid= "d0.05nm3")

# phase shifter

#b1 = RBend (l=0.0575, angle=0.0, id = "b1")
#b2 = RBend (l=0.0575, angle=-0.0, id = "b2")
b1 = RBend (l=und.lperiod, angle=0.0, eid= "b1")
b2 = RBend (l=und.lperiod, angle=-0.0, eid= "b2")


#psu=(b1,b2,b2,b1)
psu= Drift (l=b1.l*2 + b2.l*2 + d3.l, eid= "d1")

#qf = Quadrupole (l=0.1, id = "qf")
#qd = Quadrupole (l=0.1, id = "qd")
qf = Quadrupole (l=und.lperiod*2, eid= "qf")
qd = Quadrupole (l=und.lperiod*2, eid= "qd")


#qfh = Quadrupole (l=0.05, id = "qfh")
#qdh = Quadrupole (l=0.05, id = "qdh")
qfh = Quadrupole (l=qf.l / 2., eid= "qfh")
qdh = Quadrupole (l=qd.l / 2., eid= "qdh")


cell_ps = (und, d2, qf, psu, und, d2, qd, psu)

def sase3_segment(n=0): return (und, d2, qd, psu) + n*cell_ps

sase3 = (und, d2, qd, psu) + 11*cell_ps
sase = sase3


# self-seeding setup
und2 = Undulator(nperiods=73, lperiod=0.068, Kx=0.0, eid= "und2")
sase3_ss_1 = (und2, d2, qd, psu, und2, d2, qf, psu, und2,d2, qd, psu, und2, d2, qf, psu, und2)

sase3_ss_2 = (d2, qd, psu, und, d2, qf, psu)

lc = d2.l + qd.l + psu.l + und.l +  d2.l + qf.l +  psu.l
lcm = psu.l +  und.l + d2.l


d1_c = Drift (l=0.1, eid= "d1_c")
b1_c = Hcor (l=0.2, angle=1.e-5, eid= "b1_c")
b2_c = Hcor (l=0.2, angle=-1.e-5, eid= "b2_c")

d2_c = Drift (l=(lcm - 2*d1_c.l - 2*b1_c.l - 2*b2_c.l)/3., eid= "d2_c")
d3_c = Drift (l=(lcm - 2*d1_c.l - 2*b1_c.l - 2*b2_c.l)/3., eid= "d3_c")

chic = (d2, qd, d1_c, b1_c, d2_c, b2_c, d3_c, b2_c, d2_c, b1_c, d1_c, qf, psu)
sase3_ss_2m = chic

sase3_ss_3 = (und, d2, qd, psu) + 3*cell_ps


sase3_ss = sase3_ss_1 + sase3_ss_2m + sase3_ss_3 + (d1_c, b1_c, d2_c, b2_c, d3_c, b2_c, d2_c, b1_c, d1_c, qf, psu) + sase3_ss_3 


# for matching
extra_fodo = (und, d2, qdh)
extra_fodo_2 = (qfh, psu, und, d2, qdh)
l_fodo = qf.l / 2 + (b1.l + b2.l + b2.l + b1.l + d3.l) + und.l + d2.l + qf.l / 2 

## uncomment this for simple estimates w/o focusing (SR)
#und = Undulator(nperiods=73*24, lperiod=0.068, Kx=0.0, id = "und") 
#sase3 = (und)

# example settings 28m beta, 1002.95987383 eV (1.23618298443e-09 m)
und.Kx = 9.52
qf.k1 = 0.72
qd.k1 = -0.72
b1.angle = 3.1e-05
b2.angle =-3.1e-05

# example settings 17.5 GeV
#und.Kx = 9.1398  # 1000eV
und.Kx = 7.4178  # 1500eV
#und.Kx = 6.3850  # 2000eV
#und.Kx = 5.1490  # 3000eV
#und.Kx = 3.5   # 14KeV

# example settings 14.0 GeV
#und.Kx = 4.03 # 3000eV
#und.Kx = 5.04 # 2000eV
und.Kx = 5.87 # 1500eV
#und.Kx = 9.06 # 650eV


beam = Beam()
beam.E = 17.5
beam.sigma_E = 0.002
beam.emit_xn = 0.4e-6 
beam.emit_yn = 0.4e-6 
beam.gamma_rel = beam.E / (0.511e-3)
beam.emit_x = beam.emit_xn / beam.gamma_rel
beam.emit_y = beam.emit_yn / beam.gamma_rel
beam.beta_x = 33.7
beam.beta_y = 23.218
beam.alpha_x = 1.219
beam.alpha_y = -0.842

beam.tpulse = 80    # electron bunch length in fs (rms)
beam.C = 1.0        # bunch charge (nC)
beam.I = 1.0e-9 * beam.C / ( np.sqrt(2*pi) * beam.tpulse * 1.e-15 ) 

#beam.emit = {0.02: [0.32e-6,0.32e-6], 0.1: [0.39e-6,0.39e-6], 0.25: [0.6e-6,0.6e-6], 0.5: [0.7e-6,0.7e-6], 1.0: [0.97e-6,0.97e-6]}
beam.emit = {0.02: [0.2e-6,0.18e-6], 0.1: [0.32e-6,0.27e-6], 0.25: [0.4e-6,0.36e-6], 0.5: [0.45e-6,0.42e-6], 1.0: [0.8e-6,0.84e-6]}

def f1(n, n0, a0, a1, a2):
    '''
    piecewise-quadratic tapering function
    '''
    for i in range(1,len(n0)):
        if n < n0[i]:
            return a0 + (n-n0[i-1])*a1[i-1] + (n-n0[i-1])**2 * a2[i-1]
        a0 += (n0[i]-n0[i-1])*a1[i-1] + (n0[i]-n0[i-1])**2 * a2[i-1]
    
    return 1.0


def f2(n, n0, a0, a1, a2):
    '''
    exponential tapering
    '''
    
    if n <= n0:
        return a0
    
    return a0 * (  1 + a1 * (n - n0)**a2 )



def get_taper_coeff(ebeam, ephoton):
    if ebeam == 14:
        if ephoton > 400 and ephoton < 1000:
            n0 = [0,6,25,35]
            a0 = 0.999
            a1 = [-0., -0.001,  -0.00 ]
            a2 = [0., -0.00012, -0.000 ]
            return n0, a0, a1, a2
        if ephoton >= 1000 and ephoton < 2000:
            n0 = [0,7, 25,35]
            a0 = 0.999
            a1 = [-0., -0.001,  -0.00 ]
            a2 = [0., -0.00012, -0.000 ]
            return n0, a0, a1, a2
        if ephoton >= 2000 and ephoton < 2999:
            n0 = [0,8, 25,35]
            #n0 = [0,10, 25,35] # 1nc
            a0 = 0.999
            a1 = [-0., -0.001,  -0.00 ]
            a2 = [0., -0.0001, -0.000 ]
            #a2 = [0., -0.00005, -0.000 ]
            return n0, a0, a1, a2
        if ephoton >= 2999:
            n0 = [0,10, 25,35]
            #n0 = [0,13, 25,35] # 1nc
            a0 = 0.999
            a1 = [-0., -0.001,  -0.00 ]
            a2 = [0., -0.0001, -0.000 ]
            return n0, a0, a1, a2

    if ebeam == 8.5:
        pass

