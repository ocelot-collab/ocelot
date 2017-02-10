from ocelot import *

'''
Lattice of Kurchatov Light Sourse "Siberia - 2"
'''

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
