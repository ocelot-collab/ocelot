'''
response function measurement
'''

from pickle import load

from pylab import *

steerers = []


class ResponseFunction:
   def __init__(self):
      pass
try:
   f_name = sys.argv[1]
except:
   f_name = 'rf.dat'

resp_funcs = load( open(f_name,'r') )



print 'response functions'


for rf in resp_funcs:
   #print rf.device, rf.set_vals, rf.rb_vals, rf.bpm_vals['5UND3']
   plt.figure(rf.device + '.orb')
   #print [s[0] for s in rf.bpm_vals['2ORS']]
   p1, = plt.plot(rf.set_vals, [s[0] for s in rf.bpm_vals['2ORS']], 'rp')
   p2, = plt.plot(rf.set_vals, [s[1] for s in rf.bpm_vals['2ORS']], 'bp')
   legend([p1,p2], ['x','y'])
   plt.figure(rf.device + '.sase')
   plt.plot(rf.set_vals, rf.sase_vals, 'gp-')

plt.show()
