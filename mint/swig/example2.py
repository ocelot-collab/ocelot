from pylab import *

import ocelot.mint.swig.dcs as dcs

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
	print 'max:', np.max(h)
	
	blm_alarm_ch = "TTF2.DIAG/BLM.ALARM/" + blm_name + '/THRFHI'
	val = dcs.get_device_val(blm_alarm_ch)
	print 'alarm:', val * 1.25e-3
	


gmd_channel = 'TTF2.DAQ/PHFLUX/OUT33/VAL'
mag_channel = 'TTF2.MAGNETS/STEERER/H3UND1/PS.RBV'
#blm_channel = 'TTF2.DIAG/BLM/1R.UND1/CH00.TD'
blm_channel = 'FLASH.DIAG/BLM/3.2FL2SASE3/SIGNAL.TD'
h = np.array(dcs.get_device_td(blm_channel))
print 'sum:', np.sum(h)
#plt.plot(h)
#plt.show()

val = dcs.get_device_val(gmd_channel)
print 'gmd:', val


gmd_bpm_x_channel = 'TTF2.FEL/GMDPOSMON/TUNNEL/IX.POS'
gmd_bpm_y_channel = 'TTF2.FEL/GMDPOSMON/TUNNEL/IY.POS'

val = dcs.get_device_val(gmd_bpm_x_channel)
print 'gmd tunnel bpm x:', val
val = dcs.get_device_val(gmd_bpm_y_channel)
print 'gmd tunnel bpm y:', val


gmd_bpm_x_channel = 'TTF2.FEL/GMDPOSMON/BDA/IX.POS'
gmd_bpm_y_channel = 'TTF2.FEL/GMDPOSMON/BDA/IY.POS'

val = dcs.get_device_val(gmd_bpm_x_channel)
print 'gmd bda bpm x:', val
val = dcs.get_device_val(gmd_bpm_y_channel)
print 'gmd bda bpm y:', val

val = dcs.get_device_val(mag_channel)
print 'steerer:', val

# probably incorrect
h = np.array(dcs.get_device_td('TTF2.FEL/BKR.FLASH.STATE/BKR.FLASH.STATE/ENERGY.CLIP.SPECT'))
print 'sum:', np.sum(h)
plt.plot(h)

spec = np.array(dcs.get_device_td('TTF2.EXP/PBD.PHOTONWL.ML/WAVE_LENGTH/VAL.TD'))
plt.figure()
plt.plot(spec)

mcp = np.array(dcs.get_device_td('FLASH.DIAG/MCP.ADC/FL2MCP/MCP1.TD'))
plt.figure()
plt.plot(mcp)

plt.show()
