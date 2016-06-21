'''
response function measurement
'''
from time import sleep
from pickle import dump

from pylab import *

import ocelot.mint.swig.dcs as dcs

steerers = []


class ResponseFunction:
   def __init__(self):
      pass

def get_sase():
    #gmd_channel = 'TTF2.DAQ/PHFLUX/OUT33/VAL' 
    gmd_channel = 'TTF2.FEL/BKR.FLASH.STATE/BKR.FLASH.STATE/SLOW.INTENSITY' 
    
    n_avg = 1000
    val = 0.0
    for i in xrange(n_avg):
	    val += dcs.get_device_val(gmd_channel) / n_avg
    return val


def measure_response(steerers, step = 0.01, n_steps = 3):
   '''
   response for blms, bpms, sase, 
   '''
   bpm_names = ['4UND3', '5UND3']


   bpm_names = ['1TCOL',
                '6TCOL',
                '8TCOL',
                '3ECOL',
                '5ECOL',
                '2ORS',
                '7ORS',
                '9ORS',
                '12ORS',
                '1SFUND2',
                '1SFUND3',
                '1SFUND4',
                '1SFELC',
                '1SMATCH',
                '6SMATCH',
                '13SMATCH',
                '14SMATCH',
                '2UND1',
                '4UND1',
                '5UND1',
                '2UND2',
                '4UND2',
                '5UND2',
                '2UND3',
                '4UND3',
                '5UND3',
                '2UND4',
                '4UND4',
                '5UND4',
                '2UND5',
                '4UND5',
                '5UND5',
                '2UND6',
                '4UND6',
                '5UND6']


   bpms = []
   for bpm_name in bpm_names:
	bpm = dcs.BPM("TTF2.DIAG/BPM/" + bpm_name)
	bpms.append(bpm)


   blms = []

   resp_funcs = []

   for s in steerers:

      rf = ResponseFunction()
      rf.device = s
      rf.set_vals = []
      rf.rb_vals = []

      rf.bpm_vals = {}
      for bpm in bpm_names:
         rf.bpm_vals[bpm] = []


      rf.blm_vals = {}
      rf.sase_vals = []
      
      mag_channel = 'TTF2.MAGNETS/STEERER/' + s + '/PS' 
      mag_channel_rbv = 'TTF2.MAGNETS/STEERER/' + s + '/PS.RBV' 
      
      rf_channel_ampl = 'FLASH.RF/LLRF.CONTROLLER/CTRL.ACC23/SP.AMPL'
      rf_channel_ph = 'FLASH.RF/LLRF.CONTROLLER/CTRL.ACC23/SP.PHASE'

      # TODO: rename mag_channel dev_channel
      #mag_channel = mag_channel_rbv = rf_channel_ph

      rbv = dcs.get_device_val(mag_channel_rbv)
      init_v = dcs.get_device_val(mag_channel)
      print 'reading ', mag_channel, init_v, rbv

      set_val = init_v

      for i in range(n_steps):
         set_val += step
	 # set value
         dcs.set_device_val(mag_channel, set_val)
         sleep(1.0)
	 rbv = dcs.get_device_val(mag_channel_rbv)
	 rf.set_vals.append(set_val)
	 rf.rb_vals.append(rbv)

	 for i_bpm in range(len(bpms)):
            #x, y = set_val**2, set_val**3	 
            dcs.get_bpm_val(bpms[i_bpm])
            rf.bpm_vals[bpm_names[i_bpm]].append((bpms[i_bpm].x,bpms[i_bpm].y)) 

         rf.sase_vals.append( get_sase() )

      for i in range(2*n_steps):
         set_val -= step
	 # set value
         dcs.set_device_val(mag_channel, set_val)
         sleep(1.0)
	 rbv = dcs.get_device_val(mag_channel_rbv)
	 rf.set_vals.append(set_val)
	 rf.rb_vals.append(rbv)
	 for i_bpm in range(len(bpms)):
            #x, y = set_val**2, set_val**3	 
            dcs.get_bpm_val(bpms[i_bpm])
            rf.bpm_vals[bpm_names[i_bpm]].append((bpms[i_bpm].x,bpms[i_bpm].y)) 
         rf.sase_vals.append( get_sase() )

      for i in range(n_steps):
         set_val += step
	 # set value
         dcs.set_device_val(mag_channel, set_val)
	 rbv = dcs.get_device_val(mag_channel_rbv)
         sleep(1.0)
	 rf.set_vals.append(set_val)
	 rf.rb_vals.append(rbv)
	 for i_bpm in range(len(bpms)):
            #x, y = set_val**2, set_val**3	 
            dcs.get_bpm_val(bpms[i_bpm])
            rf.bpm_vals[bpm_names[i_bpm]].append((bpms[i_bpm].x,bpms[i_bpm].y)) 
         rf.sase_vals.append( get_sase() )

      #dcs.set_device_val(mag_channel, init_v)

      resp_funcs.append(rf)
   return resp_funcs 
   
def measure_response_2(steerer1, steerer2, step1, step2, n_steps = 3):
   '''
   response for blms, bpms, sase, 
   '''
   bpm_names = ['4UND3', '5UND3']


   bpm_names = ['1TCOL',
                '6TCOL',
                '8TCOL',
                '3ECOL',
                '5ECOL',
                '2ORS',
                '7ORS',
                '9ORS',
                '12ORS',
                '1SFUND2',
                '1SFUND3',
                '1SFUND4',
                '1SFELC',
                '1SMATCH',
                '6SMATCH',
                '13SMATCH',
                '14SMATCH',
                '2UND1',
                '4UND1',
                '5UND1',
                '2UND2',
                '4UND2',
                '5UND2',
                '2UND3',
                '4UND3',
                '5UND3',
                '2UND4',
                '4UND4',
                '5UND4',
                '2UND5',
                '4UND5',
                '5UND5',
                '2UND6',
                '4UND6',
                '5UND6']


   bpms = []
   for bpm_name in bpm_names:
	bpm = dcs.BPM("TTF2.DIAG/BPM/" + bpm_name)
	bpms.append(bpm)


   blms = []

   resp_funcs = []


   rf = ResponseFunction()
   rf.device = steerer1 + ':' + steerer2
   rf.set_vals = []
   rf.rb_vals = []

   rf.bpm_vals = {}
   for bpm in bpm_names:
      rf.bpm_vals[bpm] = []


   rf.blm_vals = {}
   rf.sase_vals = []
      
   mag_channel_1 = 'TTF2.MAGNETS/STEERER/' + steerer1 + '/PS' 
   mag_channel_rbv_1 = 'TTF2.MAGNETS/STEERER/' + steerer1 + '/PS.RBV' 

   mag_channel_2 = 'TTF2.MAGNETS/STEERER/' + steerer2 + '/PS' 
   mag_channel_rbv_2 = 'TTF2.MAGNETS/STEERER/' + steerer2 + '/PS.RBV' 
      
   rbv_1 = dcs.get_device_val(mag_channel_rbv_1)
   init_v_1 = dcs.get_device_val(mag_channel_1)
   print 'reading ', mag_channel_1, init_v_1, rbv_1

   rbv_2 = dcs.get_device_val(mag_channel_rbv_2)
   init_v_2 = dcs.get_device_val(mag_channel_2)
   print 'reading ', mag_channel_2, init_v_2, rbv_2

   set_val_1 = init_v_1
   set_val_2 = init_v_2

   t_out = 10.0

   for i in range(n_steps):
      print '...', set_val_1, set_val_2, '...'
      set_val_1 += step1
      set_val_2 += step2
      dcs.set_device_val(mag_channel_1, set_val_1)
      dcs.set_device_val(mag_channel_2, set_val_2)
      sleep(t_out)
      rbv_1 = dcs.get_device_val(mag_channel_rbv_1)
      rbv_2 = dcs.get_device_val(mag_channel_rbv_2)
      rf.set_vals.append(set_val_1)
      rf.rb_vals.append(rbv_1)

      for i_bpm in range(len(bpms)):
         dcs.get_bpm_val(bpms[i_bpm])
         rf.bpm_vals[bpm_names[i_bpm]].append((bpms[i_bpm].x,bpms[i_bpm].y)) 

      rf.sase_vals.append( get_sase() )

   for i in range(2*n_steps):
      print '...', set_val_1, set_val_2, '...'
      set_val_1 -= step1
      set_val_2 -= step2
      dcs.set_device_val(mag_channel_1, set_val_1)
      dcs.set_device_val(mag_channel_2, set_val_2)
      sleep(t_out)
      rbv_1 = dcs.get_device_val(mag_channel_rbv_1)
      rbv_2 = dcs.get_device_val(mag_channel_rbv_2)
      rf.set_vals.append(set_val_1)
      rf.rb_vals.append(rbv_1)

      for i_bpm in range(len(bpms)):
         dcs.get_bpm_val(bpms[i_bpm])
         rf.bpm_vals[bpm_names[i_bpm]].append((bpms[i_bpm].x,bpms[i_bpm].y)) 

      rf.sase_vals.append( get_sase() )

   for i in range(n_steps):
      print '...', set_val_1, set_val_2, '...'
      set_val_1 += step1
      set_val_2 += step2
      dcs.set_device_val(mag_channel_1, set_val_1)
      dcs.set_device_val(mag_channel_2, set_val_2)
      sleep(t_out)
      rbv_1 = dcs.get_device_val(mag_channel_rbv_1)
      rbv_2 = dcs.get_device_val(mag_channel_rbv_2)
      rf.set_vals.append(set_val_1)
      rf.rb_vals.append(rbv_1)
      for i_bpm in range(len(bpms)):
         dcs.get_bpm_val(bpms[i_bpm])
         rf.bpm_vals[bpm_names[i_bpm]].append((bpms[i_bpm].x,bpms[i_bpm].y)) 
      rf.sase_vals.append( get_sase() )

      #dcs.set_device_val(mag_channel, init_v)

   resp_funcs.append(rf)
   return resp_funcs 


'''
steerers = ['Q13SMATCH']
resp_funcs = measure_response(steerers, step = 0.25, n_steps = 20)
'''

resp_funcs = measure_response_2('V7SMATCH','V14SMATCH', step1 = 0.005, step2=0.000, n_steps = 10)

print 'response fucntions'

try:
   f_name = sys.argv[1]
except:
   f_name = 'rf.dat'

dump(resp_funcs, open(f_name,'w'))

for rf in resp_funcs:
   print rf.device, rf.set_vals, rf.rb_vals, rf.bpm_vals['5UND3']
   plt.figure()
   print [s[0] for s in rf.bpm_vals['4UND3']]
   plt.plot(rf.set_vals, [s[0] for s in rf.bpm_vals['5UND3']], 'rp')
   plt.plot(rf.set_vals, [s[1] for s in rf.bpm_vals['5UND3']], 'bp')
   plt.figure()
   plt.plot(rf.set_vals, rf.sase_vals, 'gp')

plt.show()
