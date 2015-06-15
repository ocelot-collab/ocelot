'''
example02 -- elliptical mirror (KB) -- ray tracing
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

    m1 = EllipticMirror(r=[0,0,1*m], size=[50*cm,50*cm,30*mm], a=[100*cm, 10.0*cm ], id="m1")
    m2 = EllipticMirror(r=[0,0,0*cm], size=[2*cm,50*cm,30*mm], a=[5*m, 1.*cm ], id="m2")
    a1 = Aperture(r=[0,0,6*m], size=[0.1*m,0.1*m,0.20*m], no=[0.0,0,1], d=[10*mum, 10*mum], id="a1")


    geo = Geometry([m2,a1])
    #geo = Geometry([m2])

    return geo


geo = init_geometry()
scene = init_plots(['geometry:x','geometry:y'], geo)

rays = []

a = geo.find("m2").a


for i in range(40):
    
    r = Ray(r0=[0.01,a[1],-sqrt(a[0]**2 - a[1]**2)], k=[0,-a[1]/a[0] * (1+ np.random.randn()*0.1) ,1])
    #test 1 -0.0359695
    #r = Ray(r0=[0,geo.find("m2").a[1],-sqrt(a[0]**2 - a[1]**2)], k=[0,-0.0359695 ,1])
    
    trace_ray(r, geo)
    rays.append(r)

'''
f = sqrt(a[0]**2 - a[1]**2)
r = Ray(r0=[0,a[1],-f], k=[0,-a[1]/a[0] ,1])
trace_ray(r, geo)
rays.append(r)
scene.ax[0].plot( [-f,f],[a[1],a[1]],color='red', lw=1)
'''

plot_rays(scene.ax[0], rays, proj='x')
plot_rays(scene.ax[1], rays, proj='y')


plt.show()
