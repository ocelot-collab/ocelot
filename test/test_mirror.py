'''
example01 -- 4 mirrors -- ray tracing
UNDER DEVELOPMENT!
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

    m1 = Mirror(r=[0,0,-20*cm], size=[10*cm,10*cm,20*mm], no=[0.0,-1.0,-1], id="m1")

    geo = Geometry([m1])

    return geo


geo = init_geometry()
scene = init_plots(['geometry:y'], geo)

rays = []

r = Ray(r0=[0,0,-0.6], k=[0,0,1])
trace_ray(r, geo)
rays.append(r)

r = Ray(r0=[0,0,-0.6], k=[0,0.1,1])
trace_ray(r, geo)
rays.append(r)

r = Ray(r0=[0,0,-0.25], k=[0,-0.1,1])
trace_ray(r, geo)
rays.append(r)


plot_rays(scene.ax[0], rays, proj='y')

plt.show()
