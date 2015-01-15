'''
tuning using 4 corrector direct sase optimization
'''
import ocelot.utils.mint.mint as mint
import ocelot.utils.mint.swig.dcs as dcs
import os, sys

from pylab import *
from scipy.optimize import *
from time import sleep


cor_names = ['V14SMATCH']

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


def get_sase():
    pass

def scan_one_cor(device_name, start_val = 0.0, step = 0.01, stop_alarm_level = 0.6):

    '''
    go up then down and then up again. Change direction when 
    stop_alarm_level of alarm level reached.
    '''
    mag_channel = 'TTF2.MAGNETS/STEERER/' + device_name + '/PS'
    mag_channel_rbv = 'TTF2.MAGNETS/STEERER/' + device_name + '/PS.RBV'

    val = start_val

    while True:

        val += step

        sleep(0.5)

        # set device
        rbv = dcs.get_device_val(mag_channel_rbv)
        #res =  dcs.set_device_val(mag_channel, val)

        alarm = np.max(get_alarms())
        print "\n\n rbv = {2} nval={0} alarm={1}\n\n".format(val, alarm, rbv)

        if alarm >= stop_alarm_level:
            step = - step


def match_orbit(correctors, bpm_target):


    bpm_names = []
    planes = []
    bpms = []
    targets = []

    for bpm in bpm_target.keys():
        plane  = bpm[-1:]
        bpm_name = bpm[:-2]
    
        planes.append(plane)
        bpm_names.append(bpm_name)

        targets.append(bpm_target[bpm])

	bpm = dcs.BPM("TTF2.DIAG/BPM/" + bpm_name)
	bpms.append(bpm)


    def error_func(x):
        
        pen_max = 100.0

        #print 'error_func: ', bpm_names, '->',  planes

        for i in xrange(len(x)):
            print 'x[{0}]={1}'.format(i, x[i])
            if abs(x[i]) > 3.5:
                return pen_max


        for i in xrange(len(correctors)):
            mag_channel = 'TTF2.MAGNETS/STEERER/' + correctors[i] + '/PS'
            print 'setting', mag_channel, '->',x[i]
            #res =  dcs.set_device_val(mag_channel, 1.255)

        sleep(0.0005)

        n_avg = 1000

        for i_bpm in xrange(len(bpms)):
            x = 0.0
            y = 0.0
            for i in xrange(n_avg):
                dcs.get_bpm_val(bpms[i_bpm])
                x += bpms[i_bpm].x / n_avg
                y += bpms[i_bpm].y / n_avg
            print 'bpm read:', bpms[i_bpm].id, x, y
            print 'target', planes[i_bpm], targets[i_bpm]


        alarm = np.max(get_alarms())

        pen = 0.0

        if alarm > 1.0:
            return pen_max
        if alarm > 0.7:
            return alarm * 10.0
        pen += alarm

        return pen
    
        
    for bpm in bpm_target.keys():
        plane  = bpm[-1:]
        bpm_name = bpm[:-2]
        print bpm_name, '->',  plane

    x = [0.0]*len(bpm_target.keys())
    res  = fmin(error_func,x,xtol=1e-5, maxiter=20, maxfun=20)


if __name__ == "__main__":
    print 'starting optimization...'
    #scan_one_cor(cor_names[0])
    match_orbit(['V14SMATCH'],{'2UND1.Y':0.0})
