import ocelot.utils.mint.mint as mint
import ocelot.utils.mint.swig.dcs as dcs
import os, sys

from pylab import *


print 'bpm read:'

print 'bpm test ...'

bpm_names = ['2UND1', '2UND2', '2UND4', '2UND5', '2UND6']

blm_names = ['1R.UND3']

for bpm_name in bpm_names:
	bpm = dcs.BPM("TTF2.DIAG/BPM/" + bpm_name)
	dcs.get_bpm_val(bpm)
	print 'bpm read:', bpm.id, bpm.x, bpm.y, bpm.z_pos
	
for blm_name in blm_names:
	blm = dcs.Device("TTF2.DIAG/BLM/" + blm_name)
	dcs.get_device_info(blm)
	blm_channel = blm.id + '/CH00.TD'
	print 'blm info:', blm.id, bpm.z_pos
	h = np.array(dcs.get_device_td(blm_channel))
	print 'sum:', np.sum(h)

gmd_channel = 'TTF2.DAQ/PHFLUX/OUT33/VAL'
mag_channel = 'TTF2.MAGNETS/STEERER/H3UND1/PS.RBV'

blm_channel = 'TTF2.DIAG/BLM/1R.UND1/CH00.TD'

h = np.array(dcs.get_device_td(blm_channel))
print 'sum:', np.sum(h)
#plt.plot(h)
#plt.show()

val = dcs.get_device_val(gmd_channel)
print 'gmd:', val

val = dcs.get_device_val(mag_channel)
print 'steerer:', val
