'''
example01 -- 4 mirrors -- ray tracing
'''

from ocelot.optics.elements import *
from ocelot.optics.ray import Ray, trace as trace_ray
from ocelot.gui.optics import *

from numpy import *


m = 1.0
cm = 1.e-2
mm = 1.e-3
mum = 1.e-6



def init_geometry():
    global dt, t, wf, n

    m1 = Mirror(r=[0,0,-20*cm], size=[10*cm,10*cm,20*mm], pitch_ang=3*pi/4., id="m1")
    m2 = Mirror(r=[0,-50*cm,-20*cm], size=[10*cm,5*cm,20*mm], pitch_ang=-pi/4., id="m2")
    m3 = Mirror(r=[0,-50*cm,20*cm], size=[10*cm,10*cm,1*mm], pitch_ang=pi/4., id="m3")
    m4 = Mirror(r=[0,0,20*cm], size=[10*cm,8*cm,10*mm], pitch_ang=-3*pi/4., id="m4")

    geo = Geometry([m1,m2,m3,m4])

    return geo


geo = init_geometry()
scene = init_plots(['geometry:y'], geo)

rays = []

sigma_y = 0.1e-1 # rad

for i in range(100):
    r = Ray(r0=[0,0,-0.3], k=[0,np.random.randn()*sigma_y,1])
    trace_ray(r, geo)
    rays.append(r)

plot_rays(scene.ax[0], rays, proj='y')

plt.show()
