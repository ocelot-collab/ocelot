'''
focusing in y
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

    m1 = EllipticMirror(r=[0,0,0*cm], size=[10*cm,50*cm,30*mm], a=[5*m, 1.*cm ], id="m1")
    a1 = Aperture(r=[0,0,5*m], size=[0.02*m,0.02*m,0.20*m], d=[10*mum, 10*mum], id="a1")
    geo = Geometry([m1, a1])

    return geo


geo = init_geometry()
scene = init_plots(['geometry:x','geometry:y'], geo)

rays = []


a = geo.find("m1").a

r = Ray(r0=[0.0,a[1],-sqrt(a[0]**2 - a[1]**2)], k=[0.0,-a[1]/a[0] * (1+ 0.0) ,1])
trace_ray(r, geo)
rays.append(r)

r = Ray(r0=[0.0,a[1],-sqrt(a[0]**2 - a[1]**2)], k=[0.0,-a[1]/a[0] * (1 - 0.05) ,1])
trace_ray(r, geo)
rays.append(r)

r = Ray(r0=[0.0,a[1],-sqrt(a[0]**2 - a[1]**2)], k=[0.0,-a[1]/a[0] * (1 + 0.05) ,1])
trace_ray(r, geo)
rays.append(r)


plot_rays(scene.ax[0], rays, proj='x', alpha=1.0)
plot_rays(scene.ax[1], rays, proj='y', alpha=1.0)


plt.show()
