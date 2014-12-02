import os, sys
import dcs

from pylab import *

print 'testing info...'
print dcs.cvar.pi
print dcs.info()


print 'testing beam info...'
b = dcs.Beam(14)
print b.i
b.foo()
print b.i

print 'testing 1d functions...'
f = dcs.Func_1d(400)
print f[2]
f[2] = 4.1
print f[2]

print len(f)

#print 'np array', np.array(f)

dcs.test_func_1d(f, len(f))


plt.plot(f)

g = dcs.Func_1d(4)
g[0] = -1.1

print 'sum', f.sum()
print 'addition', f + g


print 'devce getter test:'

print dcs.get_device_val("bbb")

p = dcs.Parameters()
p.x = 37.6

'''
dcs.track(f, p)
print f

dcs.track(f, p)
print f
'''

plt.show()