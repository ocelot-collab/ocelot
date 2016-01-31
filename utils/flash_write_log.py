'''
tuning using 4 corrector direct sase optimization
'''
#import ocelot.utils.mint.mint as mint
#import ocelot.utils.mint.swig.dcs as dcs
import os, sys

from pylab import *
#from scipy.optimize import *
from time import sleep, time

from ocelot.utils.mint.flash1_interface_pydoocs import FLASH1MachineInterface

sys.path.append('../')

#from tune_common import *


h10 = 'H10SMATCH'
h12 = 'H12SMATCH'
v7 = 'V7SMATCH'
v14 = 'V14SMATCH'

hu1 = 'H3UND1'
hu2 = 'H3UND2'
hu3 = 'H3UND3'
hu4 = 'H3UND4'
hu5 = 'H3UND5'
hu6 = 'H3UND6'

qu13 = 'Q5UND1.3.5'
qu24 = 'Q5UND2.4'

#q12 = 'Q12SMATCH'
q13 = 'Q13SMATCH'
q14 = 'Q14SMATCH'
q15 = 'Q15SMATCH'


#names = ["H8TCOL", "V8TCOL"]
#names = [hu1, hu2, hu3, hu4, hu5, hu6]
#names = [qu13, qu24]
#names = [q13, q15]
#names = [h10, h12, v7, v14]
#names = ["H3DBC3", "V3DBC3", "H10ACC7", "V10ACC7"]
names = ["H8TCOL", "V8TCOL", hu1, hu2, hu3, hu4, hu5, hu6, qu13, qu24,h10, h12, v7, v14, "H3DBC3", "V3DBC3", "H10ACC7", "V10ACC7"]
#names =  ["H8TCOL", "V8TCOL", hu1, hu2, hu3, hu4, hu5, hu6, qu13, qu24,h10, h12, v7, v14]
#names = [  "V10ACC7"]
#values = [mi.get_galue(name) for name in names]


mi = FLASH1MachineInterface()

def print_string(values):
    print str(time())+"\t"+"\t".join([str(v) for v in values])


print "#time" + "\t" + "\t".join([str(v) for v in names]) + "\t" + "sase"
while True:
    sleep(1)
    #print time()
    values = [mi.get_value(name) for name in names]
    values.append( mi.get_sase(detector = 'gmd_fl1_slow'))

    print_string(values)


"""
def get_pyro():
    pyro_1_ch = 'TTF2.FEEDBACK/LONGITUDINAL/MONITOR5/MEAN_AVG'
    p1 = dcs.get_device_val(pyro_1_ch)
    #pyro_2_ch = 'TTF2.FEEDBACK/LONGITUDINAL/MONITOR8/MEAN_AVG'
    #p2 = dcs.get_device_val(pyro_2_ch)
    return p1
"""
