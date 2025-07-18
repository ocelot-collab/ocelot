import pytest
from ocelot import *
from ocelot.gui import *
import matplotlib.pyplot as plt

tws0 = Twiss(beta_x=6.6, beta_y=16.4, emit_xn=0.5e-6, emit_yn=0.5e-6, E=1)

d = Drift(l=1)
qf = Quadrupole(l=0.5, k1=0.6)
qd = Quadrupole(l=0.25, k1=-0.6)
u = Undulator(lperiod=0.04, nperiods=50, Kx=4, Ky=0.)
m1 = Marker()
m2 = Marker()



@pytest.fixture(scope='module')
def lattice():
    return MagneticLattice((m1, qd, d, u, d, qf, d, u, d, qd, m2))



@pytest.fixture(scope='function')
def p_array():
    np.random.seed(10)
    return generate_parray(
                        sigma_tau=3e-6, sigma_p=1e-4, chirp=0.01,
                        charge=250e-12, nparticles=10000, tws=tws0, shape="gauss"
                        )