'''
test: morror in zx plane
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

    m1 = Mirror(r=[0,0,-20*cm], size=[10*cm,10*cm,20*mm], yaw_ang=-pi/4., roll_ang=pi/2., id="m1")
    
    geo = Geometry([m1])

    return geo


geo = init_geometry()

print 'no', geo.find("m1").no

#exit(0)

scene = init_plots(['geometry:x'], geo)

rays = []

r = Ray(r0=[0,0,-0.6], k=[0,0,1])
trace_ray(r, geo)
rays.append(r)

r = Ray(r0=[0,0,-0.6], k=[0.1, 0., 1])
trace_ray(r, geo)
rays.append(r)

r = Ray(r0=[0,0,-0.6], k=[-0.1,0., 1])
trace_ray(r, geo)
rays.append(r)

plot_rays(scene.ax[0], rays, proj='x')

plt.show()
