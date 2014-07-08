'''
example07 -- crystal
'''

from ocelot.optics.elements import *
from ocelot.optics.wave import *
from ocelot.optics.ray import Ray, trace as trace_ray
from ocelot.gui.optics import *

from numpy import *
import sys


m = 1.0
cm = 1.e-2
mm = 1.e-3
mum = 1.e-6



def init_geometry():
    global dt, t, wf, n

    g1 = Crystal(r=[0,0,0*cm], size=[13*cm,20*cm,10*cm], no=[0,0,-1], id="cr1")
    g2 = Crystal(r=[0,0,70*cm], size=[13*cm,20*cm,10*cm], no=[0,0,-1], id="cr2")
    g1.bw = 3.e-2
    g2.bw = 3.e-2

    geo = Geometry([g1,g2])

    return geo


geo = init_geometry()
scene = init_plots(['geometry:x'], geo)

r = Ray(r0=[0,0.0,-0.5], k=[0,0.0,1])
trace_ray(r, geo)
plot_rays(scene.ax[0], [r], proj='x')



class Signal(object):
    def __init__(self, n=100):
        self.t = np.linspace(-1,1, n)
        self.f = np.zeros_like(self.t, dtype=np.complex)
        self.n = n
        

w1 = Signal(n=1000)

v = 20 * 2 * pi
dt = 0.05
w1.f = np.exp(-1j*v*w1.t - (w1.t)**2 / (2*dt)**2) 

plt.figure()
plt.grid(True)

plt.plot(w1.t, np.real(w1.f))
plt.plot(w1.t, np.imag(w1.f))

w1.sp = np.fft.fft(w1.f)

plt.figure()
plt.grid(True)
w1.w = np.fft.fftfreq(w1.n, (w1.t[1] - w1.t[0]) )
plt.plot(w1.w, np.abs(w1.sp))


# filter
def transform_field(cr, wave):
    w0 = wave.w[np.argmax(np.abs(wave.sp))]
    bw = abs(cr.bw * w0)

    for i in xrange(wave.n):
        #print w1.w[i] , w0
        #print np.abs(w1.w[i] - w0), bw
        fact = np.exp(-(wave.w[i] - w0)**2/(2*bw)**2)
        #print fact
        #print w1.sp[i]
        wave.sp[i] = wave.sp[i] * fact
        #print w1.sp[i]


for i in xrange(len(r.r0)):
    print 'propagating thru', r.obj[i], r.s[i]
    
    if r.obj[i].__class__ == Drift:
        r.obj[i].l = r.s[i]
        print r.obj[i].__dict__
        w1.t += r.s[i]
     
    if r.obj[i].__class__ == Crystal:
        print 'applying crystal filter'
        
        transform_field(r.obj[i], w1)
            
        w1.t += r.s[i]
    

plt.plot(w1.w, np.abs(w1.sp))


w1.f = np.fft.ifft(w1.sp)

plt.figure()
plt.grid(True)


plt.plot(w1.t, np.real(w1.f))
plt.plot(w1.t, np.imag(w1.f))



plt.show()
