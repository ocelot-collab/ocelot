'''
example04 -- diffractive grating -- ray tracing
'''

from ocelot.optics.elements import *
from ocelot.optics.wave import *
from ocelot.optics.ray import Ray, trace as trace_ray
from ocelot.gui.optics import *

from numpy import *


m = 1.0
cm = 1.e-2
mm = 1.e-3
mum = 1.e-6



def init_geometry():
    global dt, t, wf, n

    g1 = Grating(r=[0,0,0*cm], size=[13*cm,20*cm,10*cm], no=[0.0,1,0], id="g1")

    geo = Geometry([g1])

    return geo


geo = init_geometry()
scene = init_plots(['geometry'], geo)

rays = []

sigma_y = 1.e-1 # rad

for i in xrange(10):
    r = Ray(r0=[0,0.5,-1], k=[0,-0.5,1])
    r.lamb = 1+i*0.1
    trace_ray(r, geo)
    print 'lamb:', r.lamb
    rays.append(r)

plot_rays(scene.ax[0], rays)

print rays[0].r0, rays[0].k, rays[0].s

plt.show()
