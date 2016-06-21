'''
SASE
'''
from ocelot.gui.accelerator import plot_lattice
from ocelot.gui.misc import plot_beam
from ocelot.cpbd.elements import Element, Quadrupole, RBend, Drift, Undulator, Hcor, Vcor
from ocelot import MagneticLattice
from ocelot.cpbd.beam import Beam, ParticleArray
from ocelot.cpbd.optics import *
from ocelot.adaptors.genesis import *
from pylab import *

'''
x = np.linspace(0, 2*np.pi, 10)
y = np.sin(x)
xvals = np.linspace(0, 2*np.pi, 50)
yinterp = np.interp(xvals, x, y)

plt.plot(x, y, 'o')

plt.plot(xvals, yinterp, '-x')

plt.show()

sys.exit(0)
'''

beam = Beam()
beam.E = 14.0
beam.beta_x = 73.7
beam.beta_y = 23.218
beam.alpha_x = 1.219
beam.alpha_y = -0.842


beamf = 'data/test_beam.txt'
beam_new = transform_beam_file(beam_file = beamf,  
							transform = [ [beam.beta_x,beam.alpha_x], [beam.beta_y,beam.alpha_y] ], 
							energy_scale = beam.E / 17.5, n_interp = 1000)

		
print 'beam file size:', len(beam_new.z)

plot_beam(plt.figure(), beam_new)


plt.figure()
plt.plot(beam_new.z, beam_new.I, '-x')



plt.show()

