'''
wave optics
'''

from numpy import sin, cos, pi, sqrt, log, exp, array, random, sign
from numpy.linalg import norm
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
import scipy.integrate as integrate
#import matplotlib.animation as animation

from ocelot.optics.elements import *
from ocelot.common.globals import *


class TransferFunction(object):
    '''
    data container for Fourier Optics transfer functions
    '''
    def __init__(self):
        self.k = None # wave vector - 2*pi/wavelength
        self.tr = None # complex value of transmission - modulus*exp(-i*phase)
        self.ref = None # .. of reflection
        self.xlamds = None # carrier wavelength
        self.mid_k = None # center of feature in spectrum
        self.dk = None # width of feature in spectrum
    
    def ev(self):
        return self.k* h_eV_s/2/pi * speed_of_light
    
    def __mul__(self, f):
        if f.__class__ == TransferFunction:
            f2 = TransferFunction()
            f2.k = f.k
            f2.ev = f.ev
            # TODO check data grid alignment
            
            f2.tr = self.tr * f.tr
            f2.ref = self.ref * f.ref
            
            return f2
        return None


class WaveFront:
    def __init__(self):
        pass

class Scene:
    def __init__(self):
        pass


def init():
    global scene
    scene.line_wf.set_data([], [])
    scene.time_text.set_text('')
    #scene.profile_im.set_data(np.ones([51,51])*10)
    res = []
    if 'geometry' in scene.views:
        res.append(scene.line_wf)
        res.append(scene.time_text)
    if 'detectors' in scene.views:
        res.append(scene.profile_im)
    return res


def normal(x1,x2):
    d = (x2[0] - x1[0]) / (x2[1] - x1[1])
    #print 'x1,x2,d=', x1, x2, d 
    n1, n2 = 1.0 / np.sqrt(1+d**2), -d / np.sqrt(1+d**2)
    return np.array([n1, n2])


def rotate_pi(v,n):
    vrot = -v + 2*n*np.dot(n,v)
    return vrot 



'''
data structures for optical field propagation
'''
class Mesh:
    def __init__(self, nx, ny, dtype=np.float):
        self.nx = nx
        self.ny = ny
        
        self.points = np.zeros([nx,ny], dtype=dtype)
                
        self.x = 0 
        self.y = 0 
         
        self.dx = 1.0 
        self.dy = 1.0 
        
    def __getitem__(self, idx):
        return self.points[idx]
    def __setitem__(self, idx, val):
        self.points[idx] = val
        
    def __str__(self):
        s = "Mesh " + str(self.nx) + 'x' + str(self.ny) + ' '
        s += 'xmin='+ str(self.x) + ' xmax=' + str(self.x + self.dx * (self.nx - 1) ) + ' '  
        s += 'ymin='+ str(self.y) + ' ymax=' + str(self.y + self.dy * (self.ny - 1) ) + ' '
        s+= '\n' + str(self.points)
        return s
    def idx(self, x):
        '''
        return mesh point idx of left-bottom corner of the cell a coordinate belongs to
        if coordinate outside mesh return -1
        '''
        ix = (x[0] - self.x) / self.dx
        if ix < 0 or ix > self.nx - 1:
            ix = -1

        iy = (x[1] - self.y) / self.dy
        if iy < 0 or iy > self.ny - 1:
            iy = -1

        return ix,iy
    
    def init(self, f = lambda x, y : 0):
        
        x=0; y=0
        
        for i1 in range(self.nx):
            for i2 in range(self.ny):
                x = self.x + self.dx * i1
                y = self.y + self.dy * i2
                self[i1,i2] = f(x,y)
                


class ParaxialFieldSlice():
    '''
    complex transverse electric field E_x + i E_y
    '''
    def __init__(self, lam=1.0, nx=31, ny=31, size_x=1.0, size_y = 1.0):
        '''
        lam -- wavelength (m)
        '''
        c = 1.0
        self.lam = lam     
        self.k = 2*pi / self.lam 
        self.w = self.k / c 
        
        self.nx = nx
        self.ny = ny

        self.size_x = size_x
        self.size_y = size_y
        
        self.x = np.zeros(nx)
        self.y = np.zeros(ny)
            
    def init_field(self, f = lambda x, y : 0):
        self.mesh = Mesh(nx=self.nx, ny=self.ny, dtype = np.complex)
        self.mesh.x = -self.size_x
        self.mesh.y = -self.size_y
        
        self.mesh.dx = 2*(-self.mesh.x) / ( self.mesh.nx -1)
        self.mesh.dy = 2*(-self.mesh.y) / ( self.mesh.ny -1)
        
        x=0; y=0
        
        for i1 in range(self.mesh.nx):
            for i2 in range(self.mesh.ny):
                x = self.mesh.x + self.mesh.dx * i1
                y = self.mesh.y + self.mesh.dy * i2
                self.mesh[i1,i2] = f(x,y)
                
                self.x[i1] = x
                self.y[i2] = y
                
                #print i1, i2, x, y, self.mesh[i1,i2] 

    def __getitem__(self, idx):
        return self.mesh[idx]
    def __setitem__(self, idx, val):
        self.mesh[idx] = val

def rescale(of, scale=2.0):
    of.size_x /= scale
    of.size_y /= scale
        
    of.x /= scale
    of.y /= scale

    of.mesh.x /= scale
    of.mesh.y /= scale

    of_old = np.copy(of.mesh.points)
    for i in range(of.nx):
        for j in range(of.ny):
            i_new = int(i*scale)
            j_new = int(j*scale)
            try:
                of[i,j] = of_old[ int(i*scale). int(j*scale)]
            except:
                of[i,j] = 0
             
            
            
def propagate_fourier(of, dz, obj=None, scale=1.0):
    '''
    wave propagator
    '''
    
    if obj == None or obj.__class__ == OptDrift:
        debug('wave propagator: drift')
        spec = fft.fft2(of[:,:])
        
        kx = np.fft.fftfreq(of.nx, d=2*of.size_x/of.nx)
        ky = np.fft.fftfreq(of.ny, d=2*of.size_y/of.ny)
        
        for i in range(of.nx):
            for j in range(of.ny):
                k = 2*pi / of.lam #of.w / c
                #print (kx[i]/k), (ky[j]/k)
                #phi = k * sqrt(1 - (kx[i]/k)**2 - (ky[j]/k)**2)
                phi = -pi * of.lam * ( (kx[i])**2 + (ky[j])**2 )
                #print phi*dz
                spec[i,j] *= exp(1j * phi*dz + 1j *k*dz)
            
        of[:,:] = fft.ifft2(spec)
        
    if obj.__class__ == Aperture:
        debug('wave propagator: aperture', obj.d)
        for i in range(of.nx):
            for j in range(of.ny):
                #print of.x[i], obj.d[0]
                if (of.x[i]/obj.d[0])**2 + (of.y[j]/obj.d[1])**2 >1:
                    of[i,j] = 0 

    if obj.__class__ == Lense:
        debug('wave propagator: lense, f=', obj.f, " [m]")
        for i in range(of.nx):
            for j in range(of.ny):
                phi = pi*( (of.x[i]/sqrt(of.lam*obj.f))**2 + (of.y[j]/sqrt(of.lam*obj.f))**2 )
                of[i,j] *= np.exp(-1j*phi) 
                #of[i,j] *= i
    

def propagate_fresnel(of, dz, scale=1.0):
    '''
    Propagate paraxial field slice in free space, Fresnel 
    '''
    k = 2*pi/ of.lam

    #of_old = np.copy(of.mesh.points)

    for i in range(of.nx):
        for j in range(of.ny):
            tmp = 0.0 + 0.0j
            print(i,j)
            for i1 in range(of.nx):
                for j1 in range(of.ny):
                    phi = 1j * k * ( (of.x[i1] - of.x[i])**2 + (of.y[j1] - of.y[j])**2 ) / (2.0 * dz)
                    #print phi
                    tmp = tmp +  of[i1,j1] * exp(phi) / (of.nx * of.ny)
            of[i,j] = tmp * exp(1j*k*dz) / (1j * of.lam * dz)
            print(of[i,j])
        
    
