from pylab import *

import ocelot.mint.swig.dcs as dcs

print 'testing info...'
print dcs.cvar.pi
print dcs.info()


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


h = dcs.test_func_1d_2(300)
plt.plot(h)

print "parameter getter test..."
p = dcs.get_parameters();
print p.energy

print 'device getter test:'

print dcs.get_device_val("bbb")


print 'device setter test:'

print 'bpm test ...'
bpm = dcs.BPM("bpm_id")

print bpm.id, bpm.x, bpm.y

dcs.get_bpm_val(bpm)

print bpm.id, bpm.x, bpm.y


print 'orbit test ...'

orb = dcs.get_orbit()
print orb[0].id, orb[0].x, orb[0].y
print orb[1].id, orb[1].x, orb[1].y

print orb.get("bpm1").id, orb.get("bpm1").x, orb.get("bpm1").y
print orb.get("bpm2").id, orb.get("bpm2").x, orb.get("bpm2").y
print orb.get("bpm48").id


plt.show()
