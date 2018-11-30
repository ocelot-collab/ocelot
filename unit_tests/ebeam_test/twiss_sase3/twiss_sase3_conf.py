"""Test parameters description"""

import pytest
import numpy as np

from ocelot import *

"""lattice elements description"""

und = Undulator(nperiods=73, lperiod=0.068, Kx=0.0, eid= "und")

d = Drift(l=1.0, eid="d")
d1 = Drift(l=0.55, eid="d1")
d2 = Drift(l=und.lperiod, eid="d2")
d3 = Drift(l=und.lperiod, eid="d0.05nm3")

b1 = RBend(l=und.lperiod, angle=0.0, eid="b1")
b2 = RBend(l=und.lperiod, angle=-0.0, eid="b2")

psu = Drift(l=2.0*b1.l + 2.0*b2.l + d3.l, eid="d1")

qf = Quadrupole(l=2.0*und.lperiod, eid="qf")
qd = Quadrupole(l=2.0*und.lperiod, eid="qd")

qfh = Quadrupole(l=qf.l/2.0, eid="qfh")
qdh = Quadrupole(l=qd.l/2.0, eid="qdh")

cell_ps = (und, d2, qf, psu, und, d2, qd, psu)


# self-seeding setup
und2 = Undulator(nperiods=73, lperiod=0.068, Kx=0.0, eid= "und2")
sase3_ss_1 = (und2, d2, qd, psu, und2, d2, qf, psu, und2, d2, qd, psu, und2, d2, qf, psu, und2)

sase3_ss_2 = (d2, qd, psu, und, d2, qf, psu)

lc = d2.l + qd.l + psu.l + und.l + d2.l + qf.l +  psu.l
lcm = psu.l +  und.l + d2.l


d1_c = Drift(l=0.1, eid="d1_c")
b1_c = Hcor(l=0.2, angle=1.e-5, eid="b1_c")
b2_c = Hcor(l=0.2, angle=-1.e-5, eid="b2_c")

d2_c = Drift(l=(lcm - 2*d1_c.l - 2*b1_c.l - 2*b2_c.l)/3.0, eid="d2_c")
d3_c = Drift(l=(lcm - 2*d1_c.l - 2*b1_c.l - 2*b2_c.l)/3.0, eid="d3_c")

chic = (d2, qd, d1_c, b1_c, d2_c, b2_c, d3_c, b2_c, d2_c, b1_c, d1_c, qf, psu)
sase3_ss_2m = chic
sase3_ss_3 = (und, d2, qd, psu) + 3*cell_ps
#sase3_ss = sase3_ss_1 + sase3_ss_2m + sase3_ss_3 + (d1_c, b1_c, d2_c, b2_c, d3_c, b2_c, d2_c, b1_c, d1_c, qf, psu) + sase3_ss_3 


# for matching
extra_fodo = (und, d2, qdh)
extra_fodo_2 = (qfh, psu, und, d2, qdh)
l_fodo = qf.l / 2 + (b1.l + b2.l + b2.l + b1.l + d3.l) + und.l + d2.l + qf.l / 2 


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


"""pytest fixtures descripteion"""

@pytest.fixture(scope='module')
def sase3_ss():
    return (sase3_ss_1 + sase3_ss_2m + sase3_ss_3 + (d1_c, b1_c, d2_c, b2_c, d3_c, b2_c, d2_c, b1_c, d1_c, qf, psu) + sase3_ss_3)

    
@pytest.fixture(scope='module')
def lattice(sase3_ss):
    return MagneticLattice(sase3_ss)


@pytest.fixture(scope='module')
def beam():

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
    beam.I = 1.0e-9 * beam.C / (np.sqrt(2.0*np.pi) * beam.tpulse * 1.0e-15)
    beam.emit = {0.02: [0.2e-6, 0.18e-6], 0.1: [0.32e-6, 0.27e-6], 0.25: [0.4e-6, 0.36e-6], 0.5: [0.45e-6, 0.42e-6], 1.0: [0.8e-6, 0.84e-6]}
    
    return beam
