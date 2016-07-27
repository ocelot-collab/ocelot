'''
main tuning script
'''
from time import sleep

from pylab import *
from scipy.optimize import *
import scipy.optimize as opt

import ocelot.mint.swig.dcs as dcs
from tune_common import *

init_blms(blms)


def init_corrector_vals(correctors):
   vals = [0.0]*len(correctors)#np.zeros(len(correctors))
   for i in range(len(correctors)):
      mag_channel = 'TTF2.MAGNETS/STEERER/' + correctors[i] + '/PS'
      vals[i] = dcs.get_device_val(mag_channel)
   return vals

def init_gap_vals(unds):
   vals = [0.0]*len(unds)
   for i in range(len(unds)):
      #gap_ch_rbv = 'TTF2.EXP/SFLASH.UND.MOTOR/'+unds[i]+'/GAP.SET' #'FLASH.UTIL/FL2.UND.MOTOR/'+unds[i]+'/GAP.SET'
      gap_ch_rbv = 'FLASH.UTIL/FL2.UND.MOTOR/'+unds[i]+'/GAP'
      vals[i] = dcs.get_device_val(gap_ch_rbv)

   return vals

def init_ps_vals(unds):
   vals = [0.0]*len(unds)
   for i in range(len(unds)):
      #ps_ch_rbv = 'TTF2.MAGNETS/DIPOLE/'+unds[i]+'/PS' # 'TTF2.MAGNETS/STEERER/P3'+unds[i]+'/PS.RBV' 
      ps_ch_rbv = 'TTF2.MAGNETS/STEERER/P3'+unds[i]+'/PS.RBV' 
      vals[i] = dcs.get_device_val(ps_ch_rbv)

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
    res  = opt.fmin(error_func,x,xtol=1e-5, maxiter=20, maxfun=20)

def restore_orbit(ref_file):

    bpm_names = ['2UND1','4UND1','5UND1','2UND2','4UND2','5UND2','2UND3','4UND3','5UND3', '2UND4', '4UND4','5UND4',
                 '2UND5','4UND5','5UND5','2UND6','4UND6','5UND6']

    bpms = parse_ref_orbit('/home/iagapov/public/bpm_2014-12-15T14_30_47+01_149uj.csv', bpm_mask = bpm_names)
    
    z = [ i[3] for i in bpms]
    x = [ i[1] for i in bpms]
    y = [ i[2] for i in bpms]

    print bpms
    plt.figure()
    plt.grid(True)
    plt.plot(z,x,'rd-', lw=2)
    plt.plot(z,y,'gd-', lw=2)

    bpm_targx = {}
    for b in bpms:
      bpm_targx[b[0]+'.X'] = b[1] 

    bpm_targy = {} 
    for b in bpms:
      bpm_targy[b[0]+'.Y'] = b[2] 

    print bpm_targx
    print bpm_targy

    #match_orbit(['V14SMATCH','V7SMATCH'],bpm_targy)
    match_orbit(['H10SMATCH','H12SMATCH'],bpm_targy)

    plt.show()



def max_sase(correctors, opt_pointing=False):
    '''
    direct sase optimization with simplex, using correctors
    '''

    sase_ref = get_sase()

    if opt_pointing:
	    weight_gmd_bpm_1 = 10.0
	    weight_gmd_bpm_2 = 10.0
    else:
	    weight_gmd_bpm_1 = 0.0
	    weight_gmd_bpm_2 = 0.0

    def error_func(x):
        
        pen_max = 100.0

        #print 'error_func: ', bpm_names, '->',  planes

        for i in xrange(len(x)):
            print 'x[{0}]={1}'.format(i, x[i])
            if abs(x[i]) > 50.5:
                return pen_max


        for i in xrange(len(correctors)):
            mag_channel = 'TTF2.MAGNETS/STEERER/' + correctors[i] + '/PS'
            print 'setting', mag_channel, '->',x[i]
            res =  dcs.set_device_val(mag_channel, x[i])

        sleep(1.6)

	sase = get_sase(detector='gmd_fl1_slow')
        alarm = np.max(get_alarms())
	z1, z2 = get_sase_pos()



	print 'alarm:', alarm
	print 'sase:', sase
	print 'pointing', z1, z2, 'weights', weight_gmd_bpm_1, weight_gmd_bpm_2


        pen = 0.0

        if alarm > 1.0:
            return pen_max
        if alarm > 0.7:
            return alarm * 50.0
        pen += alarm

	pen -= sase	
	#pen += weight_gmd_bpm_1* 0*(z1[1]**2 + 0*z1[0]**2) + weight_gmd_bpm_2 * (z2[0]**2 + 0*z2[1]**2)
        #if sase < 0.9 * sase_ref: pen += 1000.0
        #pen += np.abs(z1[1])

        print 'penalty:', pen

        return pen
    
        
    x = init_corrector_vals(correctors)
    
    res  = fmin(error_func,x,xtol=1e-3, maxiter=150, maxfun=150)


def max_sase_new(unds, dev='ps'):
    '''
    direct sase optimization with simplex, using correctors
    '''

    sase_ref = get_sase(detector = 'sflash.mcp')


    def error_func_gap(x):
        
        pen_max = 100.0

        for i in xrange(len(x)):
            print 'x[{0}]={1}'.format(i, x[i])
            if x[i] < 10.0 or x[i] > 220.0:
                return pen_max


        #print 'error_func: ', bpm_names, '->',  planes

        for i in xrange(len(unds)):
              gap_ch = 'TTF2.EXP/SFLASH.UND.MOTOR/'+unds[i]+'/GAP.SET'
              gap_ch_rbv = 'TTF2.EXP/SFLASH.UND.MOTOR/'+unds[i]+'/GAP'
              gap_speed_ch = 'TTF2.EXP/SFLASH.UND.MOTOR/'+unds[i]+'/GAP.SPEED'
              gap_start_ch = 'TTF2.EXP/SFLASH.UND.MOTOR/'+unds[i]+'/CMD'
              print 'setting', gap_ch, '->',x[i]
              res =  dcs.set_device_val(gap_ch, x[i])
              #dcs.set_device_val(gap_start_ch, 1)

        sleep(1.0)

	sase = get_sase(detector = 'flash2.mcp')
        alarm = np.max(get_alarms())

	print 'alarm:', alarm
	print 'sase:', sase

        pen = 0.0

        if alarm > 1.0:
            return pen_max
        if alarm > 0.7:
            return alarm * 50.0
        pen += alarm

	pen -= sase	

        print 'penalty:', pen

        return pen
   
    def error_func_ps(x):
        
        pen_max = 100.0

        #print 'error_func: ', bpm_names, '->',  planes

        for i in xrange(len(unds)):
            #mag_channel = 'TTF2.MAGNETS/DIPOLE/'+unds[i]+'/PS'
            mag_channel = 'TTF2.MAGNETS/STEERER/P3'+unds[i]+'/PS' 
            print 'setting', mag_channel, '->',x[i]
            res =  dcs.set_device_val(mag_channel, x[i])

        sleep(1.0)

	sase = get_sase(detector = 'flash2.mcp')
        alarm = np.max(get_alarms())

	print 'alarm:', alarm
	print 'sase:', sase

        pen = 0.0

        if alarm > 1.0:
            return pen_max
        if alarm > 0.7:
            return alarm * 50.0
        pen += alarm
        
	pen -= sase	
        
        print 'penalty:', pen

        return pen   
     
    if dev == 'gap':
       print 'initial gaps:'
       x = init_gap_vals(unds)
       print x
       res  = fmin(error_func_gap,x,xtol=1e-3, maxiter=150, maxfun=150)
    if dev == 'ps':
       print 'initial ps:'
       x = init_ps_vals(unds)
       print x
       #res  = fmin(error_func_ps,x,xtol=1e-3, maxiter=150, maxfun=150)



def scan_gap(und, minv, maxv, n, t_scan = 1):
    '''
    gap scan -- going up and down one time
    '''
    gmd = np.zeros(n)
    gap_val = np.zeros(n)


    #sleep(7.0)

    for i in range(len(gmd)):
        
        if i < n:
            gap_val[i] = minv + (maxv - minv) * i / float(n)
        else:
            gap_val[i] = maxv - (maxv - minv)*(i - n) / float(n)
        
        gap_ch = 'TTF2.EXP/SFLASH.UND.MOTOR/'+und+'/GAP.SET'
        gap_ch_rbv = 'TTF2.EXP/SFLASH.UND.MOTOR/'+und+'/GAP'
        gap_speed_ch = 'TTF2.EXP/SFLASH.UND.MOTOR/'+und+'/GAP.SPEED'
        gap_start_ch = 'TTF2.EXP/SFLASH.UND.MOTOR/'+und+'/CMD'
	res =  dcs.set_device_val(gap_ch, gap_val[i])
        dcs.set_device_val(gap_start_ch, 1)

	print 'setting...', gap_val[i]
	if i==0: sleep(10.0)
	sleep(2.0)
	gmd[i] = get_sase(detector = 'sflash.mcp')

    return gap_val, gmd


def scan_ps_sflash(dev, minv, maxv, n, t_scan = 1):
    '''
    sflash phase shifter  scan -- going up and down one time
    '''
    gmd = np.zeros(n)
    cor_val = np.zeros(n)

    #sleep(7.0)

    for i in range(len(gmd)):
        
        if i < n:
            cor_val[i] = minv + (maxv - minv) * i / float(n)
        else:
            cor_val[i] = maxv - (maxv - minv)*(i - n) / float(n)

        ch_m = 'TTF2.MAGNETS/DIPOLE/'+dev+'M/PS'
        ch_c = 'TTF2.MAGNETS/DIPOLE/'+dev+'C/PS'

	#res =  dcs.set_device_val(ch1_m, cor_val[i])
        #res =  dcs.set_device_val(ch1_c, cor_val[i]*0.4)

	print 'setting...', cor_val[i]
	if i==0: sleep(10.0)
	sleep(1.0)
	gmd[i] = get_sase(detector = 'sflash.mcp')

    return cor_val, gmd


def scan_ps(dev, minv, maxv, n, t_scan = 1):
    '''
    phase shifter scan -- going up and down one time
    '''
    mcp = np.zeros(n)
    cor_val = np.zeros(n)

    #sleep(7.0)

    for i in range(len(mcp)):
        
        if i < n:
            cor_val[i] = minv + (maxv - minv) * i / float(n)
        else:
            cor_val[i] = maxv - (maxv - minv)*(i - n) / float(n)


        #ps_ch = 'TTF2.MAGNETS/STEERER/P3'+dev+'/PS' 
        ps_ch = 'TTF2.FEEDBACK/FL2.WAVELENGTHCONTROL/P' + dev + '/DELTA_PHI'

	res =  dcs.set_device_val(ps_ch, cor_val[i])

	print 'setting...', cor_val[i]
	if i==0: sleep(20.0)
	sleep(2.0)
	mcp[i] = get_sase(detector = 'flash2.mcp')

    return cor_val, mcp



def max_spec(rf_params, correctors):
    '''
    optimization of the spectrum (width/density???) with rf parameters
    use correctors to restore orbit after each rf step???
    '''
    rf_ch_amp_1 = 'FLASH.RF/LLRF.CONTROLLER/CTRL.ACC1/SP.AMPL'
    rf_ch_ph_1 = 'FLASH.RF/LLRF.CONTROLLER/CTRL.ACC1/SP.PHASE'

    rf_ch_amp_39 = 'FLASH.RF/LLRF.CONTROLLER/CTRL.ACC39/SP.AMPL'
    rf_ch_ph_39 = 'FLASH.RF/LLRF.CONTROLLER/CTRL.ACC39/SP.PHASE'

    rf_ch_amp_2 = 'FLASH.RF/LLRF.CONTROLLER/CTRL.ACC23/SP.AMPL'
    rf_ch_ph_2 = 'FLASH.RF/LLRF.CONTROLLER/CTRL.ACC23/SP.PHASE'

    init_amp_1 = dcs.get_device_val(rf_ch_amp_1)
    init_ph_1 = dcs.get_device_val(rf_ch_ph_1)
      
    init_amp_39 = dcs.get_device_val(rf_ch_amp_39)
    init_ph_39 = dcs.get_device_val(rf_ch_ph_39)
      
    init_amp_2 = dcs.get_device_val(rf_ch_amp_2)
    init_ph_2 = dcs.get_device_val(rf_ch_ph_2)
    

    def error_func(x):
        
        pen_max = 100.0

        #print 'error_func: ', bpm_names, '->',  planes

        for i in xrange(len(x)):

           print 'setting ->',x[i]

           
        v1, phi1, v39, phi39=find_parameters_RF1(c0,c1,x[0],c3)
        v2, phi2 = find_parameters_RF2(c4,c5)

        dcs.set_device_val(rf_ch_amp_1, v1)
        dcs.set_device_val(rf_ch_ph_1, phi1)
        dcs.set_device_val(rf_ch_amp_39, v39)
        dcs.set_device_val(rf_ch_ph_39, phi39)
            #dcs.set_device_val(rf_ch_amp_2, v2)
            #dcs.set_device_val(rf_ch_ph_2, phi2)

        print 'c:', c0, c1, c2, c3, c4, c5  


        sleep(1.0)

	sase = get_sase()
        alarm = np.max(get_alarms())
        f, spec = get_spectrum()

	print 'alarm:', alarm
	print 'sase:', sase
        print 'max spec:', np.max(spec) 

        pen = 0.0

        if alarm > 1.0:
            return pen_max
        if alarm > 0.7:
            return alarm * 50.0
        pen += alarm

        sase_weight = 1.0
        spec_weight = 1.0

	pen -= sase * sase_weight
        #pen -= np.max(spec) * spec_weight

        return pen
    
        
    print 'reading rf params', init_amp_1, init_ph_1, init_amp_39, init_ph_39, init_amp_2, init_ph_2
    c0,c1,c2,c3 = find_knob_val_RF1(init_amp_1, init_ph_1, init_amp_39, init_ph_39)
    c4, c5 = find_knob_val_RF2(init_amp_2, init_ph_2)
    print 'init knobs:', c0,c1,c2,c3, c4, c5

    x = [c2]

    res  = opt.fmin(error_func,x,xtol=1e-3, maxiter=100, maxfun=100)
    # max_sase(correctors)

def max_sase_2(cor_name, minv, maxv, niter=20):
    #rough_scan
    plt.figure()
    cor_val, gmd = scan_cor(cor_name, minv, maxv, niter)    
    plt.plot(cor_val, gmd, '.-')
    
    idx = (-gmd).argsort()[:5] # 5 largest values
    #TODO: this can be replaced by e.g. taking all data above threshold level
    
    print idx 
    i1, i2 = sorted(idx)[0], sorted(idx)[-1]
    print 'searching interval', cor_val[i1], cor_val[i2]
    
    cor_val, gmd = scan_cor(cor_name, cor_val[i1], cor_val[i2], niter) 
    plt.plot(cor_val, gmd, '.-')


    idx = (-gmd).argsort()[:5] # 5 largest values
    i1, i2 = sorted(idx)[0], sorted(idx)[-1]
    print 'searching interval', cor_val[i1], cor_val[i2]
    
    cor_val, gmd = scan_cor(cor_name, cor_val[i1], cor_val[i2], niter, t_scan = 20) 
    plt.plot(cor_val, gmd, '.-')

    cor_mean = mean(cor_val)
    print 'cor_val=', cor_mean
    
    scan_cor(cor_name, cor_mean, cor_mean + 1.e-3, 2)

    plt.show()



def seq_1():
    print 'starting optimization...'

    #vals, rbvs = scan_one_cor('Q13SMATCH',start_val = 42.5, step = 0.05 )
    #plt.plot(vals,rbvs, 'rp')
    #plt.show()

    #match_orbit(['V14SMATCH','V7SMATCH'],{'2UND1.Y':0.0, '4UND1.Y':0.0})
    #max_sase(['H10SMATCH','H12SMATCH','V7SMATCH','V14SMATCH'], opt_pointing = False)
    #max_sase(['H10SMATCH','H12SMATCH'], opt_pointing=False)
    #max_sase(['V14SMATCH','V7SMATCH'], opt_pointing=False)
    #max_sase(['H3UND1','H3UND2','H3UND3','H3UND4','H3UND5'], opt_pointing=True)
    #max_sase(['Q5UND1.3.5'])
    #max_sase(['Q5UND2.4'])
    #max_sase(['Q5UND1.3.5'])
    #max_sase(['Q5UND2.4'])
    #max_sase(['H3UND1','H3UND2'], opt_pointing=False)
    #max_sase(['H3UND3','H3UND4'], opt_pointing=False)
    #max_sase(['H10SMATCH','H12SMATCH'], opt_pointing=False)
    #max_sase(['V7SMATCH','V14SMATCH'], opt_pointing=True)
    #max_sase(['H8TCOL','V8TCOL'], opt_pointing=False)

    #max_sase(['Q13SMATCH'])
    #restore_orbit('/home/iagapov/public/bpm_2014-12-15T14_30_47+01_149uj.csv')
    
    '''
    device_name = 'Q5UND2.4'
    minv = -6.0
    maxv = -5.0

    corv, gmd = scan_cor(device_name, minv, maxv, 50)
    i1, = (-gmd).argsort()[:1]
    print i1, corv[i1]
    mag_channel = 'TTF2.MAGNETS/STEERER/' + device_name + '/PS'
    res =  dcs.set_device_val(mag_channel, minv)
    sleep(12.0)
    res =  dcs.set_device_val(mag_channel, corv[i1])

    #corv, gmd = scan_rf('GUN', 0.0, 3.0, 50)

    plt.plot(corv, gmd, 'rp')
    plt.show()
    '''

    '''
    max_sase_2('V7SMATCH',0.3, 0.5, niter=20)
    max_sase_2('V14SMATCH',0.6, 1.1, niter=20)
    '''

def seq_2():
   max_spec(None, None)

def seq_3():
   scan_one_qmover('Q5UND5.Y')

def seq_4():
   #scan sflash phase shifters

   ps = 'P1SFUND3'
   
   ps_m_ch = 'TTF2.MAGNETS/DIPOLE/'+ps+'M/PS' # main coil
   ps_m_ch_rbv = 'TTF2.MAGNETS/DIPOLE/'+ps+'M/PS.RBV' # maiin coil
   ps_c_ch = 'TTF2.MAGNETS/DIPOLE/'+ps+'C/PS' # correction coil
   ps_c_ch_rbv = 'TTF2.MAGNETS/DIPOLE/'+ps+'C/PS.RBV' # correction coil

   print 'ps main start ', dcs.get_device_val(ps_m_ch), '/', dcs.get_device_val(ps_m_ch_rbv)
   print 'ps cor start ', dcs.get_device_val(ps_c_ch), '/', dcs.get_device_val(ps_c_ch_rbv)
   res = dcs.set_device_val(ps_m_ch, 0.0)
   res = dcs.set_device_val(ps_c_ch, 0.0)
   #val = dcs.get_device_val(gap_ch_rbv)
   #print 'gap start ', val
   #dcs.set_device_val(gap_start_ch, 1)
   

def seq_5():
   #scan sflash undulator gaps

   und = '1SFUND2'
   
   gap_ch = 'TTF2.EXP/SFLASH.UND.MOTOR/'+und+'/GAP.SET'
   gap_ch_rbv = 'TTF2.EXP/SFLASH.UND.MOTOR/'+und+'/GAP'
   gap_speed_ch = 'TTF2.EXP/SFLASH.UND.MOTOR/'+und+'/GAP.SPEED'
   gap_start_ch = 'TTF2.EXP/SFLASH.UND.MOTOR/'+und+'/CMD'
   val = dcs.get_device_val(gap_ch_rbv)
   print 'gap start ', val
   print 'gap speed ', dcs.get_device_val(gap_speed_ch), 'mm/s'
   res = dcs.set_device_val(gap_ch, 210.0)
   val = dcs.get_device_val(gap_ch_rbv)
   print 'gap start ', val
   dcs.set_device_val(gap_start_ch, 1)


def seq_6():
   #flash2 phase shifters

   und = 'FL2SASE13'
   
   gap_ch = 'FLASH.UTIL/FL2.UND.MOTOR/'+und+'/GAP.SET'
   gap_ch_rbv = 'FLASH.UTIL/FL2.UND.MOTOR/'+und+'/GAP'
   gap_speed_ch = 'FLASH.UTIL/FL2.UND.MOTOR/'+und+'/GAP.SPEED'
   gap_start_ch = 'FLASH.UTIL/FL2.UND.MOTOR/'+und+'/CMD'
   print 'gap start ', dcs.get_device_val(gap_ch_rbv)
   print 'gap speed ', dcs.get_device_val(gap_speed_ch), 'mm/s'

   ps_ch = 'TTF2.MAGNETS/STEERER/P3'+und+'/PS' 
   ps_ch_rbv = 'TTF2.MAGNETS/STEERER/P3'+und+'/PS.RBV' 

   print 'ps start ', dcs.get_device_val(ps_ch), '/', dcs.get_device_val(ps_ch_rbv)

   res = dcs.set_device_val(ps_ch, 2.0)

   #res = dcs.set_device_val(gap_ch, 210.0)
   #val = dcs.get_device_val(gap_ch_rbv)
   #print 'gap start ', val
   #dcs.set_device_val(gap_start_ch, 1)

def seq_7():
   #chicane
   print 'chicane setting ', dcs.get_device_val('TTF2.FEEDBACK/CHICANECTRL/V2-6SFELC/ALPHA'), 'mrad'
   dcs.set_device_val('TTF2.FEEDBACK/CHICANECTRL/V2-6SFELC/ALPHA', 0.0)
   dcs.set_device_val('TTF2.FEEDBACK/CHICANECTRL/V2-6SFELC/WRITE', 1)


def seq_8():
   # sflash sase tuning
   print init_gap_vals(['1SFUND2', '1SFUND3'])
   print init_ps_vals(['P1SFUND2M','P1SFUND2C'])
   #scan_gap(und)
   #scan_ps(und)
   max_sase_new(['P1SFUND3C'], dev='ps')
   #max_sase_new(['1SFUND2', '1SFUND3'], dev='gap')

   #dev = '1SFUND4'
   #corv, gmd = scan_gap(dev, 11.0, 12.0, 50)
   
   dev = 'P1SFUND3'
   corv, gmd = scan_ps(dev, -1.0, 1.0, 50)

   i1, = (-gmd).argsort()[:1]
   print i1, corv[i1]

   #gap_ch = 'TTF2.EXP/SFLASH.UND.MOTOR/'+dev+'/GAP.SET'
   #gap_ch_rbv = 'TTF2.EXP/SFLASH.UND.MOTOR/'+dev+'/GAP'
   #gap_start_ch = 'TTF2.EXP/SFLASH.UND.MOTOR/'+dev+'/CMD'
   #print 'setting', gap_ch, '->',corv[i1]
   #res =  dcs.set_device_val(gap_ch, corv[i1])
   #dcs.set_device_val(gap_start_ch, 1)

   plt.plot(corv, gmd, 'rp')
   plt.show()
   
def seq_9():
    print 'starting optimization flash 2...'
    max_sase_new(['FL2SASE3','FL2SASE4'])

def seq_10():
   print 'phase shifter scan...'
   ps, mcp = scan_ps('FL2SASE5', 160,360 + 160, 36)
   plt.figure()
   plt.plot(ps, mcp)
   plt.show()


def seq_test():
	print 'testing...'
	print 'alarms', get_alarms()
	print 'sase gmd default:', get_sase()
	print 'sase mcp:', get_sase(detector='mcp')
	print 'sase bpm:', get_sase_pos()

	
	f, spec = get_spectrum()
	print 'spectrum:', spec
	plt.figure()
	plt.plot(f, spec)
	plt.show()
	

	V1, phi1, V3, phi3 = find_parameters_RF1(30.12, 10.0, -250.0, 5461.0)
	print 'rf knob test:', V1, phi1, V3, phi3
	print find_knob_val_RF1(V1, phi1, V3, phi3)

	V2, phi2 = find_parameters_RF2(25.0, 12.0)
	print 'rf knob test:', V2, phi2
	print find_knob_val_RF2(V2, phi2)


if __name__ == "__main__":
	seq_10()
 
