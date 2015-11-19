'''
tuning using 4 corrector direct sase optimization
'''
import ocelot.utils.mint.mint as mint
import ocelot.utils.mint.swig.dcs as dcs
import os, sys

from pylab import *
from scipy.optimize import *
from time import sleep, time

from flash1_interface import FLASH1MachineInterface

mi = FLASH1MachineInterface()



while True:
	print "{0}\t{1}".format(time(), mi.get_sase(detector = 'gmd_fl1_slow'))
