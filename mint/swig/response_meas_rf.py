'''
response function measurement
'''
from time import sleep
from pickle import dump

from pylab import *

import ocelot.mint.swig.dcs as dcs
from tune_common import *

init_blms()


class ResponseFunction:
   def __init__(self):
      pass


def measure_response(rf_params, init_v = None, step = 10.0, n_steps = 3, delay=0.5):
   '''
   response of sase, orbit, spectrum vs. rf knobs, 
   '''
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

   for s in rf_params:

      rf = ResponseFunction()
      rf.device = s
      rf.set_vals = []
      rf.rb_vals = []

      rf.bpm_vals = {}
      for bpm in bpm_names:
         rf.bpm_vals[bpm] = []


      rf.blm_vals = {}
      rf.sase_vals = []
      rf.spec_vals = []
            
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

      print 'reading rf params', init_amp_1, init_ph_1, init_amp_39, init_ph_39, init_amp_2, init_ph_2
      c0,c1,c2,c3 = find_knob_val_RF1(init_amp_1, init_ph_1, init_amp_39, init_ph_39)
      c4, c5 = find_knob_val_RF2(init_amp_2, init_ph_2)
      print 'init knobs:', c0,c1,c2,c3, c4, c5

      cur_step = step

      '''
      v1, phi1, v39, phi39=find_parameters_RF1(c0,-692.8,c2,c3)
      dcs.set_device_val(rf_ch_amp_1, v1)
      dcs.set_device_val(rf_ch_ph_1, phi1)
      dcs.set_device_val(rf_ch_amp_39, v39)
      dcs.set_device_val(rf_ch_ph_39, phi39)

      return None
      '''

      '''
      v2, phi2 = find_parameters_RF2(c4,-1588.0)
      dcs.set_device_val(rf_ch_amp_2, v2)
      dcs.set_device_val(rf_ch_ph_2, phi2)
      return None
      '''

      for i in range(n_steps*4):

         if i == n_steps:
            cur_step = - step
         if i == n_steps*3:
            cur_step = step

         if s == 'c0': c0 += cur_step
         if s == 'c1': c1 += cur_step
         if s == 'c2': c2 += cur_step
         if s == 'c3': c3 += cur_step
         if s == 'c4': c4 += cur_step
         if s == 'c5': c5 += cur_step

         v1, phi1, v39, phi39=find_parameters_RF1(c0,c1,c2,c3)
         v2, phi2 = find_parameters_RF2(c4,c5)

         if s in ['c0','c1','c2','c3']:
            dcs.set_device_val(rf_ch_amp_1, v1)
            dcs.set_device_val(rf_ch_ph_1, phi1)
            dcs.set_device_val(rf_ch_amp_39, v39)
            dcs.set_device_val(rf_ch_ph_39, phi39)
            print 'BC1'
         else:
            dcs.set_device_val(rf_ch_amp_2, v2)
            dcs.set_device_val(rf_ch_ph_2, phi2)
            print 'BC2'

         print 'knobs:', c0,c1,c2,c3, c4, c5
         print 'rf:', v1, phi1, v39, phi39, v2, phi2            
         sleep(delay)
	 rf.set_vals.append(locals()[s])

	 for i_bpm in range(len(bpms)):
            dcs.get_bpm_val(bpms[i_bpm])
            rf.bpm_vals[bpm_names[i_bpm]].append((bpms[i_bpm].x,bpms[i_bpm].y)) 

         rf.sase_vals.append( get_sase() )
         f, spec = get_spectrum()
         rf.f = f
         rf.spec_vals.append(spec)

      resp_funcs.append(rf)
   return resp_funcs 
   

resp_funcs = measure_response(['c1'], init_v=None, step = 0.05, n_steps = 30)
resp_funcs = measure_response(['c2'], init_v=None, step = 55.0, n_steps = 50)
resp_funcs = measure_response(['c3'], init_v=None, step = 15000.0, n_steps = 30)
#resp_funcs = measure_response(['c5'], init_v=None, step = 3.0, n_steps = 30)

print 'response fucntions'

try:
   f_name = sys.argv[1]
except:
   f_name = 'rf.dat'

dump(resp_funcs, open(f_name,'w'))

for rf in resp_funcs:
   #plt.figure()
   #print [s[0] for s in rf.bpm_vals['4UND3']]
   #plt.plot(rf.set_vals, [s[0] for s in rf.bpm_vals['5UND3']], 'rp')
   #plt.plot(rf.set_vals, [s[1] for s in rf.bpm_vals['5UND3']], 'bp')
   plt.figure()
   plt.plot(rf.set_vals, rf.sase_vals, 'gp')
   plt.figure()
   plt.plot(rf.set_vals, [max(s) for s in rf.spec_vals], 'rp')

plt.show()
