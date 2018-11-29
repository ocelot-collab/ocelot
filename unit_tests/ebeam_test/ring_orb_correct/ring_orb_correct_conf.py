"""Test parameters description"""

import pytest

from ocelot import *

"""lattice elements description"""

D0 = Drift(l=0.0, eid="D0")
D1 = Drift(l=1.49, eid="D1")
D2 = Drift(l=0.1035, eid="D2")
D3 = Drift(l=0.307, eid="D3")
D4 = Drift(l=0.33, eid="D4")
D5 = Drift(l=0.3515, eid="D5")
D6 = Drift(l=0.3145, eid="D6")
D7 = Drift(l=0.289, eid="D7")
D8 = Drift(l=0.399, eid="D8")
D9 = Drift(l=3.009, eid="D9")

SF = Sextupole(l=0.0001, k2= 17673.786254063251, eid="SF")
SD = Sextupole(l=0.0001, k2=-36169.817233025707, eid="SD")

q1 = Quadrupole(l=0.293, k1=2.62, eid="Q1")
q1.dx = 0.001
q1.dy = 0.001
q2 = Quadrupole(l=0.293, k1=-3.1, eid="Q2")
q3 = Quadrupole(l=0.327, k1=2.8, eid="Q3")
q4 = Quadrupole(l=0.291, k1=-3.7, eid="Q4")
q5 = Quadrupole(l=0.391, k1=4.0782, eid="Q5")
q6 = Quadrupole(l=0.291, k1=-3.534859, eid="D6")

q1s = Quadrupole(l=0.293, k1=2.62, eid="Q1")
q1s.dx = 0.001
q1s.dy = 0.001
q2s = Quadrupole(l=0.293, k1=-3.1, eid="Q2")
q3s = Quadrupole(l=0.327, k1=2.8, eid="Q3")
q4s = Quadrupole(l=0.291, k1=-3.7, eid="Q4")
q5s = Quadrupole(l=0.391, k1=4.0782, eid="Q5")
q6s = Quadrupole(l=0.291, k1=-3.534859, eid="D6")

M1 = Monitor()
M2 = Monitor()
M3 = Monitor()
M4 = Monitor()
M5 = Monitor()
M6 = Monitor()
H1 = Hcor()
H2 = Hcor()
H3 = Hcor()
H4 = Hcor()
H5 = Hcor()
H6 = Hcor()
V1 = Vcor()
V2 = Vcor()
V3 = Vcor()
V4 = Vcor()
V5 = Vcor()
V6 = Vcor()

M1s = Monitor()
M2s = Monitor()
M3s = Monitor()
M4s = Monitor()
M5s = Monitor()
M6s = Monitor()
H1s = Hcor()
H2s = Hcor()
H3s = Hcor()
H4s = Hcor()
H5s = Hcor()
H6s = Hcor()
V1s = Vcor()
V2s = Vcor()
V3s = Vcor()
V4s = Vcor()
V5s = Vcor()
V6s = Vcor()

B1 = SBend(l=0.23, angle=0.23/19.626248, eid="B1")
B2 = SBend(l=1.227, angle=1.227/4.906312, eid="B2")

Q1 = [q1, M1, H1, V1]
Q2 = [q2, M2, H2, V2]
Q3 = [q3, M3, H3, V3]
Q4 = [q4, M4, H4, V4]
Q5 = [q5, M5, H5, V5]
Q6 = [q6, M6, H6, V6]
Q1s = [q1s, M1s, H1s, V1s]
Q2s = [q2s, M2s, H2s, V2s]
Q3s = [q3s, M3s, H3s, V3s]
Q4s = [q4s, M4s, H4s, V4s]
Q5s = [q5s, M5s, H5s, V5s]
Q6s = [q6s, M6s, H6s, V6s]


"""pytest fixtures description"""

@pytest.fixture(scope='module')
def cell():
    return (D1,SF,D2,Q1,D3,Q2,D2,SD,D4,B1,B2,D5,Q3,D5,B2,B1,D6,Q4,D7,Q5,D8,Q6,D9,Q6s,D8,Q5s,D7,Q4s,D6,B1,B2,D5,Q3s,D5,B2,B1,D4,SD,D2,Q2s,D3,Q1s,D2,SF,D1)


@pytest.fixture(scope='module')
def method():
    return MethodTM()


@pytest.fixture(scope='module')
def lattice(cell, method):
    return MagneticLattice(cell, method=method)
