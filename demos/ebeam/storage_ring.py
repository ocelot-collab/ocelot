__author__ = 'Sergey Tomin'

from ocelot.gui import *
from ocelot import *

Q1 = Quadrupole(l=0.4, k1=-1.3, eid= "Q1")
Q2 = Quadrupole(l=0.8, k1=1.4, eid= "Q2")
Q3 = Quadrupole(l=0.4, k1=-1.7, eid= "Q3")
Q4 = Quadrupole(l=0.5, k1=1.19250444829 , eid= "Q4")

B  = Bend(l=2.7, k1=-.06, angle=2*pi/16., e1=pi/16., e2=pi/16., eid= "B")

SF = Sextupole(l=0.01, k2 = 150, eid= "SF") #random value
SD = Sextupole(l=0.01, k2 =-150, eid= "SD") #random value

D1 = Drift(l=2., eid= "D1")
D2 = Drift(l=0.6, eid= "D2")
D3 = Drift(l=0.3, eid= "D3")
D4 = Drift(l=0.7, eid= "D4")
D5 = Drift(l=0.9, eid= "D5")
D6 = Drift(l=0.2, eid= "D6")


cell = (D1, Q1, D2, Q2, D3, Q3, D4, B, D5, SD, D5, SF, D6, Q4, D6, SF, D5, SD,D5, B, D4, Q3, D3, Q2, D2, Q1, D1)
ring = cell

lat = MagneticLattice(ring)

beam = Beam()
beam.E = 2.5 #GeV
beam.I = 0.1 #A
tw0 = Twiss(beam)

tws = twiss(lat,tw0, nPoints=1000)
print("Qx = ", tws[-1].mux/2/pi, "  Qy = ", tws[-1].muy/2/pi)
match_tunes(lat, tws[-1], [Q1, Q2, Q3, Q4],  1.2, 0.91, ncells=1, print_proc=0)
tws = twiss(lat,tw0, nPoints=1000)
print("Qx = ", tws[-1].mux/2/pi, "  Qy = ", tws[-1].muy/2/pi)


plot_opt_func(lat, tws)
plt.show()


eb = EbeamParams(lat, beam, nsuperperiod=8)
print(eb)
