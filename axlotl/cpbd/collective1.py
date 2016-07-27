'''
tracking with collective effects
'''

from ocelot import *

# lattice definition

m1 = Marker()
m2 = Marker()
m3 = Marker()

lat = MagneticLattice( (m1,d1,b1,d1,m2,b1,d1,b1,d1,m3))

