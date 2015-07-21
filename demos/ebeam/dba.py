__author__ = 'Sergey Tomin'

from ocelot.gui.accelerator import *
from ocelot.cpbd.match import *
#from matplotlib import pylab as plt
#plt.plot(range(10))
#plt.show()

Q1 = Quadrupole(l= 0.4, k1=-1.3, id = "Q1")
Q2 = Quadrupole(l= 0.8, k1=1.4, id = "Q2")
Q3 = Quadrupole(l= 0.4, k1=-1.7, id = "Q3")
Q4 = Quadrupole(l= 0.5, k1=1.3 , id = "Q4")

B  = Bend(l=2.7, k1=-.06, angle=2*pi/16., e1=pi/16., e2=pi/16., id = "B")

SF = Sextupole(l=0., ms = 1.5, id = "SF") #random value
SD = Sextupole(l=0., ms = -1.5, id = "SD")#random value

D1 = Drift(l=2., id = "D1")
D2 = Drift(l=0.6, id = "D2")
D3 = Drift(l=0.3, id = "D3")
D4 = Drift(l=0.7, id = "D4")
D5 = Drift(l=0.9, id = "D5")
D6 = Drift(l=0.2, id = "D6")


cell = (D1, Q1, D2, Q2, D3, Q3, D4, B, D5, SD, D5, SF, D6, Q4, D6, SF, D5, SD,D5, B, D4, Q3, D3, Q2, D2, Q1, D1)


lat = MagneticLattice(cell)

tw0 = Twiss()

tws=twiss(lat, tw0, nPoints=1000)
plot_opt_func(lat, tws)

plt.show()



constr = {D1:{'Dx':0.0, 'Dxp':0.0}, 'periodic':True}

vars = [Q4]

match(lat, constr, vars, tw0)

for element in lat.sequence:
    if element.id == "Q4":
        print( "new quadrupole strength: Q4.k1 = ", Q4.k1)

tws=twiss(lat, tw0, nPoints = 1000)
plot_opt_func(lat, tws)
plt.show()