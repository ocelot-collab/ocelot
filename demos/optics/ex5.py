'''
diffraction effects on aperture -- wave
UNDER DEVELOPMENT!
'''

import matplotlib.pyplot as plt
import scipy.integrate as integrate


from ocelot.common.math_op import *

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

    a1 = Aperture(r=[0,0,-20*m], size=[1.9*mm,1.9*mm,0.20*mm], no=[0.0,0,1], d=[30*mum, 30*mum], id="a1")
    l1 = Lense(r=[0,0*cm,-19.9*m], D=0.5*mm, s1=0.1*mm, s2=0.1*mm, id="l1")
    l1.f = 2.5*m
    d1 = Detector(r=[0,0*cm,105.0*m], size=[0.2*mm, 0.2*mm, 1*mm], nx=111, ny=111, id="d1")
    d2 = Detector(r=[0,0*cm,130*m], size=[0.2*mm, 0.2*mm, 1*mm], nx=111, ny=111, id="d2")

    geo = Geometry([a1,l1,d1])

    return geo

        

def plot_field(of, title=None):
    
    x = np.linspace(-of.size_x, of.size_x, of.nx) / 1.e-3
    y = np.linspace(-of.size_y, of.size_y, of.ny) / 1.e-3
    #mu1, mu2, sig1, sig2, _ = fit_gauss_2d(x,y, np.abs(of.mesh.points))
    s2 = fwhm(x,np.abs(of.mesh.points[:,of.ny/2]))
    print('beam size:', s2, ' mm [fwhm]')

    
    fig = plt.figure(title)
    nx, ny = of.mesh.points.shape

    ax = fig.add_subplot(111)
    plt.plot(x,np.abs(of[:,ny/2]), lw=3)

    #ax = fig.add_subplot(212)
    #plt.plot(x,np.angle(of[:,ny/2]), lw=3)

    '''
    spec = fft.fft2(of[:,:])
    
    kx = np.fft.fftfreq(of.nx, d=2*of.size_x/of.nx)
    ky = np.fft.fftfreq(of.ny, d=2*of.size_y/of.ny)
    '''

    #ax = fig.add_subplot(325)
    #plt.plot(kx* of.lam, np.abs(spec[:,ny/2]), lw=3)
    #plt.imshow(np.abs(spec))
    
    '''
    propagate_fourier(of, dz=z)
    
    ax = fig.add_subplot(222)

    plt.plot(x,np.abs(of[:,ny/2]), lw=3)
    
    ax = fig.add_subplot(224)
    plt.plot(x,np.angle(of[:,ny/2]), lw=3)
    
    spec = fft.fft2(of[:,:])
    
    kx = np.fft.fftfreq(of.nx, d=2*of.size_x/of.nx)
    ky = np.fft.fftfreq(of.ny, d=2*of.size_y/of.ny)

    #ax = fig.add_subplot(326)
    #plt.plot(kx * of.lam, np.abs(spec[:,ny/2]), lw=5)

    
    s2 = fwhm(x,np.abs(of.mesh.points[:,of.ny/2]))
    print 'beam size (end):', s2, ' mm [fwhm]'
    '''

'''
aperture functions
'''

def f0(x,y, a=2):
    sig = 1.0
    r = sqrt(x**2 + y**2)
    return 1. / (sqrt(2*pi) * sig**2) * exp(-r**a / (2*sig**2))

def f02(x,y, a=2):
    sig = 1.0
    r = sqrt(x**2 + y**2)
    return 1. / (sqrt(2*pi) * sig) * exp(-r**a / (2*sig**2)) * exp(1j*r**2 ) 


def f1(x,y):
    nsig_cut = 10.0
    sig = 2.0e-5 * m
    
    if (x**2 + y**2) > (nsig_cut)**2 * sig**2: return 0.0
    
    return f0(x/sig,y/sig, a=2)
        

def f2(x,y):
    sig = 0.6e-4 * m
    return f02(x/sig,y/sig, a=1)
    return 0.0


def f3(x,y):
    sig = 1.e-4 * m
    
    phi = pi / 4
    sig1 = sig*1.05
    x1 = x*cos(phi) + y*sin(phi)
    y1 = -x*sin(phi) + y*cos(phi)
    
    if abs(x) < sig and abs(y) < sig and abs(x1) < sig1 and abs(y1) < sig1:
        return f0(x/sig,y/(2*sig))
        
    return 0.0


geo = init_geometry()
scene = init_plots(['geometry','detectors:d1','detectors:d2'], geo)
#scene = init_plots(['geometry'], geo)

rays = []

sigma_x = 30.e-2*mum
sigma_xp = 10.e-16 # rad
sigma_y = 30.e-2*mum # m
sigma_yp = 10.e-16 # rad



for i in range(1):
    r = Ray(r0=[np.random.randn()*sigma_x,np.random.randn()*sigma_y,-145*m], k=[np.random.randn()*sigma_xp,np.random.randn()*sigma_yp,1])
    print('tracing ray...')
    trace_ray(r, geo)
    rays.append(r)

plot_rays(scene.ax[0], rays)


#plot_detectors(scene, geo)

try:
    m1=geo.find("d1").matrix.transpose()
    scene.profile_im["d1"].imshow(m1, cmap='gist_heat',interpolation='none',extent=[0,1,0,1], vmin=0, vmax=10)

    m2=geo.find("d2").matrix.transpose()
    scene.profile_im["d2"].imshow(m2, cmap='gist_heat',interpolation='none',extent=[0,1,0,1], vmin=0, vmax=10)
except:
    pass


print('tracing nominal ray')
r = Ray(r0=[0,0,-145*m], k=[0,0,1])
trace_ray(r, geo)


print('initializing field...')
of = ParaxialFieldSlice(lam=4e-10*m, nx=251, ny=251, size_x=0.19*mm, size_y =0.19*mm)
of.init_field(f1)
x = np.linspace(-of.size_x, of.size_x, of.nx)
s2 = fwhm(x,np.abs(of.mesh.points[:,of.ny/2]))


print('w=', of.w , ' s^-1')
print('k=', of.k)
print('lam=', of.lam , ' m')
print('size=',s2 / 1.e-3, ' mm [fwhm]')
print('div~=',of.lam / s2 / 1.e-3, ' mrad [fwhm]')

plot_field(of, title="start")


for i in range(len(r.s)):
    print('propagating thru', r.obj[i], r.s[i] )
    propagate_fourier(of, obj=r.obj[i], dz=r.s[i])
    plot_field(of, title=r.obj[i].id + ".left")
    if r.obj[i] != None and r.obj[i].__class__ != OptDrift:
        propagate_fourier(of, obj=None, dz=r.s[i])
        plot_field(of, title= r.obj[i].id + ".right")

plt.show()

