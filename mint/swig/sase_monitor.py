'''
tuning using 4 corrector direct sase optimization
'''
from time import time

import ocelot.mint.swig.dcs as dcs
from tune_common import *


def get_sase_old():
    #gmd_channel = 'TTF2.DAQ/PHFLUX/OUT33/VAL' 
    gmd_channel = 'TTF2.FEL/BKR.FLASH.STATE/BKR.FLASH.STATE/SLOW.INTENSITY' 

    n_avg = 100
    val = 0.0
    for i in xrange(n_avg):
	    val += dcs.get_device_val(gmd_channel) / n_avg


    if val <= 5.0: return val*0.00

    return val



def get_pyro():
    pyro_1_ch = 'TTF2.FEEDBACK/LONGITUDINAL/MONITOR5/MEAN_AVG'
    p1 = dcs.get_device_val(pyro_1_ch)
    #pyro_2_ch = 'TTF2.FEEDBACK/LONGITUDINAL/MONITOR8/MEAN_AVG'
    #p2 = dcs.get_device_val(pyro_2_ch)
    return p1

#open()

while True:
	print "{0}\t{1}\t{2}".format(time(), get_sase(detector = 'gmd_fl1_slow'), get_pyro())
