'''
common definition for electron optics
'''

from numpy import *
import numpy as np
import scipy.special as sf
import scipy.integrate as integrate
from numpy.polynomial.chebyshev import *

from xframework.cpbd.elements import *
from xframework.cpbd.beam import Beam, Particle, ParticleArray
from xframework.cpbd.optics import *

from xframework.common.screen import Screen
from xframework.common.run import SRRunParameters
from xframework.common.xio import XIO
from xframework.common.globals import *
