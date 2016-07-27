'''
tuning using 4 corrector direct sase optimization
'''
from time import sleep

from pylab import *
from scipy.optimize import *

import ocelot.mint.swig.dcs as dcs
from rf import *



# for flash1
blm_names = ['14L.SMATCH',
	     '14R.SMATCH',
	     '1L.UND1',
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
             '1R.UND6',
	     '1SFUND1','1SFUND2','1SFUND3','1SFUND4',
             '1SFELC','3SFELC','4SFELC',
             '10SMATCH','3SDUMP']

# for sflash
blm_names = ['1SFUND1','1SFUND2','1SFUND3','1SFUND4',
             '1SFELC','3SFELC','4SFELC',
             '10SMATCH','3SDUMP']

# for flash2
blm_names = ['3.1FL2SASE3','3.2FL2SASE3',
	     '3.1FL2SASE4','3.2FL2SASE4',
	     '3.1FL2SASE5','3.2FL2SASE5',
	     '3.1FL2SASE6','3.2FL2SASE6',
	     '3.1FL2SASE7','3.2FL2SASE7',
	     '3.1FL2SASE8','3.2FL2SASE8',
	     '3.1FL2SASE9','3.2FL2SASE9',
	     '3.1FL2SASE10','3.2FL2SASE10',
	     '3.1FL2SASE11','3.2FL2SASE11',
	     '3.1FL2SASE12','3.2FL2SASE12',
	     '3.1FL2SASE13','3.2FL2SASE13',
	     '3.1FL2SASE14','3.2FL2SASE14']

blms = []

def init_blms(blms):
	for blm_name in blm_names:
		if blm_name.find('FL2') > 0:
			blm = dcs.Device("FLASH.DIAG/BLM/" + blm_name)
		else:
			blm = dcs.Device("TTF2.DIAG/BLM/" + blm_name)
		dcs.get_device_info(blm)
		print 'blm info:', blm.id, blm.z_pos
		blms.append(blm)


def get_alarms(blm_names):
    print 'n alarms', len(blm_names)
    alarm_vals = np.zeros(len(blm_names))
    for i in xrange(len(blm_names)):

	#blm_channel = 'TTF2.DIAG/BLM/'+blm_names[i]+'/CH00.TD' 
	#blm_channel = 'FLASH.DIAG/BLM/'+blm_names[i]+'/SIGNAL.TD' 
	#print 'reading blm channel', blm_channel
	#blm_channel = blms[i].id + '/CH00.TD'
	#h = np.array(dcs.get_device_td(blm_channel))
	#print 'max:', np.max( np.abs(h) )

	if blm_names[i].find('FL2') > 0:
		alarm_val = 2000
		blm_channel = 'FLASH.DIAG/BLM/'+blm_names[i]+'/SIGNAL.TD' 
	else:
		blm_channel = 'TTF2.DIAG/BLM/'+blm_names[i]+'/CH00.TD'
		blm_alarm_ch  = ('TTF2.DIAG/BLM/'+blm_names[i]).replace('BLM', 'BLM.ALARM') + '/THRFHI'
		blm_alarm_ch  = ('FLASH.DIAG/BLM/'+blm_names[i]).replace('BLM', 'BLM.ALARM') + '/THRFHI'
		print 'reading alarm channel', blm_alarm_ch
		alarm_val = dcs.get_device_val(blm_alarm_ch) * 1.25e-3 # alarm thr. in Volts
        print 'alarm:', alarm_val

	h = np.array(dcs.get_device_td(blm_channel))

        alarm_vals[i] = np.max( np.abs(h) ) / alarm_val 

    #print alarm_vals
    return alarm_vals 


def get_sase(detector='gmd_default'):

    if detector == 'mcp':
	    # incorrect
	    return dcs.get_device_val('TTF2.DIAG/MCP.HV/MCP.HV1/HV_CURRENT')
	    #return np.abs( np.mean(h) )

    if detector == 'flash2.mcp':
	    v = 0.0
	    n_av = 100
	    for i in range(0,n_av):
		    v += np.max( np.abs(dcs.get_device_td('FLASH.DIAG/MCP.ADC/FL2MCP/MCP1.TD')))
	    return v / n_av
	    #return np.abs( np.mean(h) )


    if detector == 'sflash.mcp':
	    val = 0.0
	    #h = dcs.get_device_td('TTF2.EXP/DIAGMCP.ADC/SFLASH2.0/CH02.TD')
	    #return np.max(h)
	    
	    for i in range(0,1000):
		    h = dcs.get_device_td('TTF2.EXP/DIAGMCP.ADC/SFLASH2.0/CH02.TD')
		    if val < np.max(h): val = np.max(h)
		    sleep(0.001)
		    #val += h[700] / 100.0
	    #return h[7000]
	    return val

    if detector == 'gmd_fl1_slow':
	     return dcs.get_device_val('TTF2.FEL/BKR.FLASH.STATE/BKR.FLASH.STATE/SLOW.INTENSITY' ) 
	    

    # default 'BKR' gmd
    h = np.array(dcs.get_device_td('TTF2.FEL/BKR.FLASH.STATE/BKR.FLASH.STATE/ENERGY.CLIP.SPECT'))
    val = np.mean(h)
    return val

    '''
    #gmd_channel = 'TTF2.DAQ/PHFLUX/OUT33/VAL' 
    gmd_channel = 'TTF2.FEL/BKR.FLASH.STATE/BKR.FLASH.STATE/SLOW.INTENSITY' 

    n_avg = 200
    val = 0.0
    for i in xrange(n_avg):
	    val += dcs.get_device_val(gmd_channel) / n_avg


    #if val <= 5.0: return val*0.00
    print 'sase:', val
    return val
    '''
    
def get_sase_pos():

    x1 = dcs.get_device_val('TTF2.FEL/GMDPOSMON/TUNNEL/IX.POS')
    y1 = dcs.get_device_val('TTF2.FEL/GMDPOSMON/TUNNEL/IY.POS')

    x2 = dcs.get_device_val('TTF2.FEL/GMDPOSMON/BDA/IX.POS')
    y2 = dcs.get_device_val('TTF2.FEL/GMDPOSMON/BDA/IY.POS')
    
    return [ (x1,y1), (x2,y2) ] 

def get_spectrum(f=None, detector='tunnel_default'):

    if detector == 'sflash':
	    spec = np.array(dcs.get_device_td('TTF2.EXP/CAM.SFELC/SFELCCD/SPECTRUM.X.TD'))
	    f = np.linspace(0, 1, len(spec))
	    return f, spec

    f_min = 13.0
    f_max = 14.0
    
    spec = np.array(dcs.get_device_td('TTF2.EXP/PBD.PHOTONWL.ML/WAVE_LENGTH/VAL.TD'))

    if f == None:
	    f = np.linspace(f_min, f_max, len(spec))

    return f, spec

def scan_one_cor(device_name, start_val = 0.0, step = 0.01, stop_alarm_level = 0.6):

    '''
    go up then down and then up again. Change direction when 
    stop_alarm_level of alarm level reached.
    '''
    mag_channel = 'TTF2.MAGNETS/STEERER/' + device_name + '/PS'
    mag_channel_rbv = 'TTF2.MAGNETS/STEERER/' + device_name + '/PS.RBV'

    val = start_val

    vals = []
    rbvs = []
    
    n_max = 5000
    i_set = 0
    while True:
        i_set += 1
        val += step

        # set device
        res =  dcs.set_device_val(mag_channel, val)
	sleep(2.0)
	rbv = dcs.get_device_val(mag_channel_rbv)

	rbvs.append(rbv)
	vals.append(val)

        alarm = np.max(get_alarms())
        print "\n\n rbv = {2} nval={0} alarm={1}\n\n".format(val, alarm, rbv)

        if alarm >= stop_alarm_level:
            step = - step
	if i_set > n_max:
           break

    return (vals, rbvs) 


def scan_one_qmover(device_name, step = 0.004, stop_alarm_level = 0.6):

    '''
    go up then down and then up again. Change direction when 
    stop_alarm_level of alarm level reached.
    '''
    mover_ch = 'TTF2.MAGNETS/QUAD.MOVER/' + device_name + '/FODOMETER.REFERENCE'
    mover_ch_rbv = 'TTF2.MAGNETS/QAUD.MOVER/' + device_name + '/FODOMETER'
    start_ch = 'TTF2.MAGNETS/QUAD.MOVER/' + device_name + '/CMD'

    val = dcs.get_device_val(mover_ch)

    vals = []
    rbvs = []
    
    n_max = 500
    i_set = 0
    while True:
        i_set += 1
        val += step

        # set device
        res =  dcs.set_device_val(mover_ch, val)
	dcs.set_device_val(start_ch, 7)

	sleep(1.0)

	#rbvs.append(rbv)
	vals.append(val)

        alarm = np.max(get_alarms())
        #print "\n\n rbv = {2} nval={0} alarm={1}\n\n".format(val, alarm, rbv)

        if alarm >= stop_alarm_level or val > 0.9 or val < -0.9:
            step = - step
	if i_set > n_max:
           break

    return vals


def scan_cor(device_name, minv, maxv, n, t_scan = 1):
    '''
    device scan (magnets) -- going up and down one time
    '''
    gmd = np.zeros(n)
    cor_val = np.zeros(n)


    #sleep(7.0)

    for i in range(len(gmd)):
        
        if i < n:
            cor_val[i] = minv + (maxv - minv) * i / float(n)
        else:
            cor_val[i] = maxv - (maxv - minv)*(i - n) / float(n)
        
        #gmd[i] = get_sase(cor_val[i], dt=t_scan)
	mag_channel = 'TTF2.MAGNETS/STEERER/' + device_name + '/PS'
	res =  dcs.set_device_val(mag_channel, cor_val[i])
	print 'setting...', cor_val[i]
	if i==0: sleep(10.0)
	sleep(0.5)
	gmd[i] = get_sase()

    return cor_val, gmd

def scan_rf(dev, minv, maxv, n, t_scan = 1):
    '''
    rf device (single phase/amplitude) scan -- going up and down one time
    '''
    gmd = np.zeros(n)
    cor_val = np.zeros(n)


    #sleep(7.0)

    for i in range(len(gmd)):
        
        if i < n:
            cor_val[i] = minv + (maxv - minv) * i / float(n)
        else:
            cor_val[i] = maxv - (maxv - minv)*(i - n) / float(n)
        
        #gmd[i] = get_sase(cor_val[i], dt=t_scan)

	rf_channel_ampl = 'FLASH.RF/LLRF.CONTROLLER/CTRL.'+dev+'/SP.AMPL'
	rf_channel_ph = 'FLASH.RF/LLRF.CONTROLLER/CTRL.'+dev+'/SP.PHASE'

	res =  dcs.set_device_val(rf_channel_ph, cor_val[i])
	print 'setting...', cor_val[i]
	if i==0: sleep(10.0)
	sleep(2.0)
	gmd[i] = get_sase()

    return cor_val, gmd


def scan_rf_knob1(dev, minv, maxv, n, t_scan = 1):
    '''
    scan rf knob using acc1 and acc39 (chirp, 2nd and 3rd order chirp)
    '''
    gmd = np.zeros(n)
    cor_val = np.zeros(n)


    #sleep(7.0)

    for i in range(len(gmd)):
        
        if i < n:
            cor_val[i] = minv + (maxv - minv) * i / float(n)
        else:
            cor_val[i] = maxv - (maxv - minv)*(i - n) / float(n)
        
        #gmd[i] = get_sase(cor_val[i], dt=t_scan)

	rf_channel_ampl = 'FLASH.RF/LLRF.CONTROLLER/CTRL.'+dev+'/SP.AMPL'
	rf_channel_ph = 'FLASH.RF/LLRF.CONTROLLER/CTRL.'+dev+'/SP.PHASE'

	res =  dcs.set_device_val(rf_channel_ph, cor_val[i])
	print 'setting...', cor_val[i]
	if i==0: sleep(10.0)
	sleep(2.0)
	gmd[i] = get_sase()

    return cor_val, gmd

def scan_rf_knob2(dev, minv, maxv, n, t_scan = 1):
    '''
    scan rf knob using acc2 (chirp)
    '''
    pass

def init_corrector_vals(correctors):
   vals = [0.0]*len(correctors)#np.zeros(len(correctors))
   for i in range(len(correctors)):
      mag_channel = 'TTF2.MAGNETS/STEERER/' + correctors[i] + '/PS'
      vals[i] = dcs.get_device_val(mag_channel)
   return vals

def parse_ref_orbit(file_name, bpm_mask):
   lines = open(file_name).read().split('\n')
   bpms = []
   for l in lines:
      tks = l.split()
      if len(tks)<2 : continue
      dev = tks[0].split('/')
      if dev[1] == 'ORBIT': 
         bpm = dev[2]
	 x, y, z = [float(i) for i in tks[2:5]]
         print bpm, x, y, z
	 if bpm in bpm_mask: bpms.append([bpm, x,y,z])

   return bpms

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
        
        pen_max = 10000.0

        #print 'error_func: ', bpm_names, '->',  planes

        for i in xrange(len(x)):
            print 'x[{0}]={1}'.format(i, x[i])
            if abs(x[i]) > 3.5:
                return pen_max


        for i in xrange(len(correctors)):
            mag_channel = 'TTF2.MAGNETS/STEERER/' + correctors[i] + '/PS'
            print 'setting', mag_channel, '->',x[i]
            #res =  dcs.set_device_val(mag_channel, x[i])

        sleep(0.0005)

	pen = 0.0

        n_avg = 20

        for i_bpm in xrange(len(bpms)):
            x = 0.0
            y = 0.0
            for i in xrange(n_avg):
                dcs.get_bpm_val(bpms[i_bpm])
                x += bpms[i_bpm].x / n_avg
                y += bpms[i_bpm].y / n_avg
            #print 'bpm read:', bpms[i_bpm].id, x, y
            #print 'target', planes[i_bpm], targets[i_bpm]
	    if planes[i_bpm] == 'Y': pen += (targets[i_bpm] - y)**2
	    if planes[i_bpm] == 'X': pen += (targets[i_bpm] - x)**2

        alarm = np.max(get_alarms())
	print 'alarm read...', alarm


        if alarm > 1.0:
            return pen_max
        if alarm > 0.7:
            return alarm * 10.0
        pen += alarm

	print 'pen=', pen

        return pen
    
        
    for bpm in bpm_target.keys():
        plane  = bpm[-1:]
        bpm_name = bpm[:-2]
        print bpm_name, '->',  plane

    x = init_corrector_vals(correctors)
    res  = fmin(error_func,x,xtol=1e-5, maxiter=20, maxfun=20)


def find_parameters_RF1(E,Es,Ess,Esss):
	'''
	determine RF amplitude and phase from 'knob' values. ACC1 + ACC39
	'''
	c = 299792458.0
	f = 1.3e9
	k=2*pi*f/c
	n=3 # harmonic cavity
	A0 = E; B0 = Es;  C0 = Ess; D0 = Esss;
	X1 = (C0+A0*k**2*n**2)/(k**2*(-1+n**2))
	Y1 = -((D0+B0*k**2*n**2)/(k**3*(-1+n**2)))
	X3 = -((C0+A0*k**2)/(k**2*(-1+n**2)))
	Y3 = (D0+B0*k**2)/(k**3*(-n+n**3))
	V1 = sqrt(X1**2 + Y1**2);
	f1 = arctan(Y1/X1) + pi/2*(1 - sign(X1)) 
	V3 = sqrt(X3**2 + Y3**2)
	f3 = arctan(Y3/X3) + pi/2*(1 - sign(X3))
	
	phi11=f1*180./pi; phi13=f3*180./pi-180.
	return (V1, phi11, V3, phi13)

def find_knob_val_RF1(V1,phi1,V39,phi39):
	'''
	determine RF knob values form RF phase/ampl parameters
	'''
	c = 299792458.0
	f = 1.3e9
	k=2*pi*f/c
	s = 0.0
	E = V1*cos(k*s + phi1 * pi / 180.) - V39*cos(3*k*s + phi39 * pi / 180.)
	Ep = -k*V1*sin(k*s + phi1* pi / 180.) + 3*k*V39*sin(3*k*s + phi39 * pi / 180.)
	Epp = -k**2*V1*cos(k*s + phi1 * pi / 180.) + 9*k**2*V39*cos(3*k*s + phi39 * pi / 180.)
	Eppp = k**3*V1*sin(k*s + phi1 * pi / 180.) - 27*k**3* V39*sin(3*k*s + phi39 * pi / 180.)

	return (E, Ep, Epp, Eppp)

def find_parameters_RF2(E,Es):
	'''
	determine RF amplitude and phase from 'knob' values. ACC2
	'''
	c = 299792458.0
	f = 1.3e9
	k=2*pi*f/c
	X1 = E
	Y1 = -Es/k
	V2 = sqrt(X1**2 + Y1**2)
	f2 = arctan(Y1/X1) + pi/2*(1 - sign(X1)) 
	
	phi2=f2*180./pi
	
	return (V2, phi2)

def find_knob_val_RF2(V2,phi2):
	'''
	determine RF knob values form RF phase/ampl parameters
	'''
	c = 299792458.0
	f = 1.3e9
	k=2*pi*f/c
	s = 0.0
	E = V2*cos(k*s + phi2 * pi / 180.)
	Ep = -k*V2*sin(k*s + phi2* pi / 180.)

	return (E, Ep)

def scan_rf_test():
	# scan acc1 and and acc39, keeping energy constant
	p1 = []
	p2 = []
	p3 = []
	p4 = []
	for e1 in np.linspace(-100, 100, 100):
		(phi11,V11,phi13,V13) = find_parameters_RF1(10.0,0.0,0.0,e1)
		p1.append(phi11)
		p2.append(V11)
		p3.append(phi13)
		p4.append(V13)
	plt.figure()
	plt.plot(p1)
	plt.plot(p2)
	plt.plot(p3)
	plt.plot(p4)
	plt.show()
