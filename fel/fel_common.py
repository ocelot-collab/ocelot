'''
SR calculations based on Monte Carlo photon generator
'''

import numpy as np
from numpy import *

import scipy.special as sf
import scipy.optimize as opt

import scipy.integrate as integrate

from numpy.polynomial.chebyshev import *

from ocelot.cpbd.elements import Element, Quadrupole, RBend, Drift, Undulator, MagneticLattice
from ocelot.cpbd.beam import Beam
from ocelot.common.screen import Screen
from ocelot.common.run import SRRunParameters
from ocelot.common.xio import XIO

from ocelot.adaptors import srwutil as srw
from ocelot.adaptors.genesis import *

import numpy as np

from ocelot.rad.sr import *
from ocelot.fel.fel import *

from ocelot.utils.launcher import *

from pylab import *
params = {'backend': 'ps', 'axes.labelsize': 18, 'text.fontsize': 16, 'legend.fontsize': 24, 'xtick.labelsize': 32,  'ytick.labelsize': 32, 'text.usetex': True}
rcParams.update(params)
rc('text', usetex=True) # required to have greek fonts on redhat


h = 4.135667516e-15
c = 299792458.0


