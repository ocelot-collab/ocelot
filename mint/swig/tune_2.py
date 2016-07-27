'''
n-stage tuning using action graph
'''
from pylab import *

import ocelot.mint.swig.dcs as dcs
from action_graph import *

blm_names = ['1L.UND1',
             '1R.UND1',
             '1L.UND2', 
             '1R.UND2', 
             '1L.UND3', 
             '1R.UND3', 
             '1L.UND4',
             '1R.UND4',
             '1L.UND5',
             '1R.UND5',
             '1L.UND6',
             '1R.UND6']

blms = []

for blm_name in blm_names:
	blm = dcs.Device("TTF2.DIAG/BLM/" + blm_name)
	dcs.get_device_info(blm)
	#print 'blm info:', blm.id, bpm.z_pos
	blms.append(blm)


def get_alarms():
    #print 'n alarms', len(blm_names)
    alarm_vals = np.zeros(len(blm_names))
    for i in xrange(len(blm_names)):
	blm_channel = blms[i].id + '/CH00.TD'
	h = np.array(dcs.get_device_td(blm_channel))
	#print 'max:', np.max( np.abs(h) )
        blm_alarm_ch  = blms[i].id.replace('BLM', 'BLM.ALARM') + '/THRFHI'
        alarm_val = dcs.get_device_val(blm_alarm_ch) * 1.25e-3 # alarm thr. in Volts
        #print 'alarm:', alarm_val
        alarm_vals[i] = np.max( np.abs(h) ) / alarm_val 

    return alarm_vals 


def variant1():
    print get_alarms()



if __name__ == "__main__":
    print 'starting optimization...'
    variant1()
