'''
tuning using 4 corrector direct sase optimization
'''
import ocelot.utils.mint.mint as mint
import ocelot.utils.mint.swig.dcs as dcs
import os, sys


from tune_common import *

from pylab import *
from scipy.optimize import *
import scipy.optimize as opt
from time import sleep

init_blms()    


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

    if opt_pointing:
	    weight_gmd_bpm_1 = 5.0
	    weight_gmd_bpm_2 = 5.0
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

        sleep(1.0)

	sase = get_sase()
        alarm = np.max(get_alarms())
	z1, z2 = get_sase_pos()


	print 'alarm:', alarm
	print 'sase:', sase
	print 'pointing', z1, z2, 'weights', weight_gmd_bpm_1, weight_gmd_bpm_2

	weight_gmd1 = 5.0
	weight_gmd2 = 5.0


        pen = 0.0

        if alarm > 1.0:
            return pen_max
        if alarm > 0.7:
            return alarm * 50.0
        pen += alarm

	pen -= sase
	
	pen += weight_gmd_bpm_1* (z1[1]**2 + z1[0]**2) + weight_gmd_bpm_2 * (z2[1]**2 + z2[0]**2)

        return pen
    
        
    x = init_corrector_vals(correctors)
    
    res  = fmin(error_func,x,xtol=1e-3, maxiter=100, maxfun=100)

def max_spec(rf_params, correctors):
    '''
    optimization of the spectrum (width/density???) with rf parameters
    use correctors to restore orbit after each rf step???
    '''
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

        sleep(1.0)

	sase = get_sase()
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

        return pen
    
        
    x = init_corrector_vals(correctors)
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
    max_sase(['H10SMATCH','H12SMATCH','V7SMATCH','V14SMATCH'], opt_pointing = True)
    #max_sase(['H10SMATCH','H12SMATCH'])
    #max_sase(['V7SMATCH','V14SMATCH'])
    #max_sase(['H3UND1','H3UND2','H3UND3','H3UND4','H3UND5'])
    #max_sase(['Q5UND2.4','Q5UND1.3.5'])
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
	seq_1()	
	#seq_test()
