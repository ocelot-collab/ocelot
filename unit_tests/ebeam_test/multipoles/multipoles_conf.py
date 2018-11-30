"""Test parameters description"""

import pytest
import numpy as np

from ocelot import *

"""lattice elements description"""

D = Drift(l=1000/16/4, eid="D")
Qf = Multipole(kn=[0., 0.021/2.], eid="Qf")
Qd = Multipole(kn=[0., -0.02], eid="Qd")
B = Multipole(kn=2.0*np.pi/32.0)
Sf = Multipole(kn=(0., 0., 0.0), eid="Sf")
Sd = Multipole(kn=(0., 0., -0.0), eid="Sd")
F = Multipole(kn=[0., 0., 0., 0., 0.1])


"""pytest fixtures description"""

@pytest.fixture(scope='module')
def cell():
    return (Qf, Sf, D, F, B, D, Qd, Sd, D, B, D, Sf, Qf)


@pytest.fixture(scope='module')
def lattice(cell):
    return MagneticLattice(16*cell)
