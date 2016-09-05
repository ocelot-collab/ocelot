"""
This script was created by Sergey Tomin for Workshop: Designing future X-ray FELs. Source and license info is on GitHub.
August 2016.
"""

# import from Ocelot main modules and functions
from ocelot import *

# import from Ocelot graphical modules
from ocelot.gui.accelerator import *

# Creating simple cell
# you can play with element parameters
# defining of the drifts
D1 = Drift(l=2.)
D2 = Drift(l=0.6)
D3 = Drift(l=0.3)
D4 = Drift(l=0.7)
D5 = Drift(l=0.9)
D6 = Drift(l=0.2)

# defining of the quads
Q1 = Quadrupole(l= 0.4, k1=-1.3)
Q2 = Quadrupole(l= 0.8, k1=1.4)
Q3 = Quadrupole(l= 0.4, k1=-1.7)
Q4 = Quadrupole(l= 0.5, k1=1.3)

# defining of the bending magnet
B = Bend(l=2.7, k1=-.06, angle=2*pi/16., e1=pi/16., e2=pi/16.)

# defining of the sextupoles
SF = Sextupole(l=0.01, k2 = 1.5) #random value
SD = Sextupole(l=0.01, k2 = -1.5) #random value

# creating of the cell
cell = (D1, Q1, D2, Q2, D3, Q3, D4, B, D5, SD, D5, SF, D6, Q4,
        D6, SF, D5, SD,D5, B, D4, Q3, D3, Q2, D2, Q1, D1)


lat = MagneticLattice(cell)

# to see total lenth of the lattice
print("length of the cell: ", lat.totalLen, "m")

tws=twiss(lat, nPoints=1000)

# to see twiss paraments at the begining of the cell, uncomment next line
# print(tws[0])

# to see twiss paraments at the end of the cell, uncomment next line
# print(tws[-1])

# plot twiss paramentrs.
plot_opt_func(lat, tws, legend=False, font_size=10)


# you can play with quadrupole strength and try to make achromat
Q4.k1 = 1.5

# to make achromat uncomment next line
# Q4.k1 =  1.18543769836
# To use matching function, please see ocelot/demos/ebeam/dba.py

# updating trransfer maps after changing element parameters.
lat.update_transfer_maps()

# recalculate twiss parameters
tws=twiss(lat, nPoints=1000)

plot_opt_func(lat, tws, legend=False)
plt.show()