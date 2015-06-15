import ocelot.utils.mint.mint as mint
import ocelot.utils.mint.swig.dcs as dcs
import os, sys

from pylab import *


print 'corrector setting test...'

#mag_name = 'H3UND2'
mag_name = 'V14SMATCH'

mag_channel_rbv = 'TTF2.MAGNETS/STEERER/' + mag_name + '/PS.RBV'
mag_channel = 'TTF2.MAGNETS/STEERER/' + mag_name + '/PS'

val = dcs.get_device_val(mag_channel_rbv)
print 'before rbv:', val

res =  dcs.set_device_val(mag_channel, 1.255)

print 'setting result returned ', res

val = dcs.get_device_val(mag_channel_rbv)
print 'after. rbv:', val
