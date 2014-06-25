'''
common definition for electron optics
'''

from numpy import *
import numpy as np
import scipy.special as sf
import scipy.integrate as integrate
from numpy.polynomial.chebyshev import *

from cpbd.elements import *
from cpbd.beam import Beam, Particle, ParticleArray
from ocelot.cpbd.optics import *

from ocelot.common.screen import Screen
from ocelot.common.run import SRRunParameters
from ocelot.common.xio import XIO
from ocelot.common.globals import *
