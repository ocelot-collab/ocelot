'''
free space propagation -- wave 
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


        

def plot_field(of, title=None):
    
    x = np.linspace(-of.size_x, of.size_x, of.nx) / 1.e-3
    y = np.linspace(-of.size_y, of.size_y, of.ny) / 1.e-3
    #mu1, mu2, sig1, sig2, _ = fit_gauss_2d(x,y, np.abs(of.mesh.points))
    s2 = fwhm(x,np.abs(of.mesh.points[:,of.ny/2]))
    print ('beam size:', s2, ' mm [fwhm]')

    
    fig = plt.figure(title)
    nx, ny = of.mesh.points.shape

    ax = fig.add_subplot(211)
    plt.plot(x,np.abs(of[:,ny/2]), lw=3)

    ax = fig.add_subplot(212)
    plt.plot(x,np.angle(of[:,ny/2]), lw=3)

    
'''
aperture functions
'''

def f0(x,y, a=2):
    sig = 1.0
    r = sqrt(x**2 + y**2)
    return 1. / (sqrt(2*pi) * sig**2) * exp(-r**a / (2*sig**2))


def f1(x,y):
    nsig_cut = 10.0
    sig = 1.0e-5 * m
    
    if (x**2 + y**2) > (nsig_cut)**2 * sig**2: return 0.0
    
    return f0(x/sig,y/sig, a=2)
        


print('initializing field...')
of = ParaxialFieldSlice(lam=5e-9*m, nx=151, ny=151, size_x=0.15*mm, size_y =0.15*mm)
of.init_field(f1)
x = np.linspace(-of.size_x, of.size_x, of.nx)
s2 = fwhm(x,np.abs(of.mesh.points[:,of.ny/2]))


print('w=', of.w , ' s^-1')
print('k=', of.k)
print('lam=', of.lam , ' m')
print('size=',s2 / 1.e-3, ' mm [fwhm]')
print('div~=',of.lam / s2 / 1.e-3, ' mrad [fwhm]')

plot_field(of, title="start")

s = []
z = []

for i in range(10):
    dz = 1.0 * m
    propagate_fourier(of, obj=None, dz=dz)
    x = np.linspace(-of.size_x, of.size_x, of.nx) 
    y = np.linspace(-of.size_y, of.size_y, of.ny) 
    s2 = fwhm(x,np.abs(of.mesh.points[:,of.ny/2]))
    z.append(dz*i)
    s.append(s2)
    print('beam size:', s2 / 1.e-3, ' mm [fwhm]')
    
    if s2>of.size_x: 
        print ('warning: need rescaling', s2, of.size_x)
        rescale(of)

    #plot_field(of, title=str(i*dz) + 'm')

plot_field(of, title="end")

plt.figure()

plt.plot(z,s)

plt.show()

