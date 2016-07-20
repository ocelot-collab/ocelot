from pylab import *

import ocelot.mint.swig.dcs as dcs

print 'rf test...'

rf_channel = 'FLASH.RF/LLRF.DAQ/VS.ACC23/AMPL.TD'
rf_channel_ampl = 'FLASH.RF/LLRF.CONTROLLER/CTRL.ACC23/SP.AMPL'
rf_channel_ph = 'FLASH.RF/LLRF.CONTROLLER/CTRL.ACC23/SP.PHASE'
rf_ampl = dcs.get_device_val(rf_channel_ampl)
rf_ph = dcs.get_device_val(rf_channel_ph)
print 'rf:', rf_ampl, rf_ph


h = np.array(dcs.get_device_td(rf_channel))
print np.max( np.abs(h) ), len(h)

pyro_1_ch = 'TTF2.FEEDBACK/LONGITUDINAL/MONITOR5/MEAN_AVG'
print 'pyro1:', dcs.get_device_val(pyro_1_ch)
pyro_2_ch = 'TTF2.FEEDBACK/LONGITUDINAL/MONITOR8/MEAN_AVG'
print 'pyro2:', dcs.get_device_val(pyro_2_ch)


plt.plot(h)
plt.show()
