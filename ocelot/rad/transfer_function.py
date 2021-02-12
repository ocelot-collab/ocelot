import logging
import time

import numpy as np
from copy import deepcopy
from math import factorial

from ocelot.optics.new_wave import *
from ocelot.common.globals import *
from ocelot import ocelog
from ocelot.common.ocelog import *
_logger = logging.getLogger(__name__) 

class Mask(Grid):
    """
    Mask element class
    The class represents a transfer function of an optical element. Note, several masks can correspond to one optical element.
    :param Grid: with following parameters  
        :params (dx, dy, dz): spatial step of the mesh, [m] 
        :param shape: number of point in each direction of 3D data mesh
        :param xlamds: carrier wavelength of the wavepacket, [m]
        :param used_aprox: the approximation used in solving the wave equation. OCELOT works in Slowly Varying Amplitude Approximation
    :param domain_z: longitudinal domain (t - time, f - frequency)   
    :param domain_x: transverse domain (s - space, k - inverse space)      
    :param domain_y: transverse domain (s - space, k - inverse space)       
    """
    def __init__(self, shape=(0, 0, 0)):
        Grid.__init__(self, shape=shape)
        self.domain_z = 't'   # longitudinal domain (t - time, f - frequency)
        self.domain_x = 's'   # transverse domain (s - space, k - inverse space)
        self.domain_y = 's'   # transverse domain (s - space, k - inverse space)
        
    def __mul__(self, other):#check this stuff when it will been needed
        """
        Multiplication operator for two masks 
        """
        m = deepcopy(self)
        if other.__class__ in [self] and self.mask is not None and other.mask is not None:
            m.mask = self.mask * other.mask
            return m
        else: ValueError("'other' must belong to Mask class") 

class PropMask(Mask):
    """
    Angular-spectrum propagation for fieldfile
    
    can handle wide spectrum
      (every slice in freq.domain is propagated
       according to its frequency)
    no kx**2+ky**2<<k0**2 limitation
    
    dfl is the RadiationField() object
    
    z0 is the propagation distance in [m]
    for 'kt' domain propagation
        no Fourier transform to frequency domain is done
        assumes no angular dispersion (true for plain FEL radiation)
        assumes narrow spectrum at center of xlamds (true for plain FEL radiation)
    'kf' propagation is a default option
    
    z>0 -> forward direction
    z<0 -> backward direction
    z=0 return original
    """   
    def __init__(self, z0, method='PropMask_kf'):
        Mask.__init__(self)
        self.z0 = z0 
        self.method = method
        self.mask = None

    def apply(self, dfl):
        """
        'apply' method for PropMask class 
        
        transform dfl to the domain of the propagation 'kf' or 'kt'
        get the transfer function calling 'get_mask' method, e.g. get H transfer function 
        multiply dfl and the mask, e.g. E = E_0 * H convolution in inverse space domain
        return the field to the original domain
        """
        
        _logger.info('propagating dfl file by %.2f meters' % (self.z0))
        
        if self.z0 == 0:
            _logger.debug(ind_str + 'z0=0, returning original')
            return dfl

        start = time.time()

        domains = dfl.domains()
        self.domain_x = 'k'
        self.domain_y = 'k'
        if self.method in ['PropMask', 'PropMask_kf']:
            dfl.to_domain('kf')
            self.domain_x = 'f'
        elif self.method == 'PropMask_kt':
            dfl.to_domain('kt')
            self.domain_x = 't'
        else:
            raise ValueError("propagation method should be 'PropMask', 'PropMask_kf' or 'PropMask_kt'")
                         
        if self.mask is None:
            self.get_mask(dfl) # get H transfer function 

        dfl.fld *= self.mask # E = E_0 * H convolution in inverse space domain
        
        dfl.to_domain(domains) # back to original domain  

        t_func = time.time() - start
        _logger.debug(ind_str + 'done in %.2f sec' % t_func)
        
        return dfl
        
    def get_mask(self, dfl):
        """
        'get_mask' method for PropMask class
        find the transfer function for propagation in 'kf' of 'ks' domains
        H = exp(iz0(k^2 - kx^2 - ky^2)^(1/2) - k) 
        """
        self.copy_grid(dfl)
       
        k_x, k_y = np.meshgrid(self.grid_kx(), self.grid_ky())
        
        if dfl.domain_z == 'f':
            k = self.grid_kz()
            self.mask = [np.exp(1j * self.z0 * (np.sqrt(k[i] ** 2 - k_x ** 2 - k_y ** 2) - k[i])) for i in range(dfl.Nz())] 
        elif dfl.domain_z == 't':
            k = 2 * np.pi / self.xlamds
            self.mask = [np.exp(1j * self.z0 * (np.sqrt(k ** 2 - k_x ** 2 - k_y ** 2) - k)) for i in range(dfl.Nz())]
        else: 
            raise ValueError('wrong field domain, domain must be ks or kf ')    
            
        return self.mask

class Prop_mMask(Mask):
    """
    Angular-spectrum propagation for fieldfile
    
    can handle wide spectrum
      (every slice in freq.domain is propagated
       according to its frequency)
    no kx**2+ky**2<<k0**2 limitation
    
    dfl is the RadiationField() object
    
    z0 is the propagation distance in [m]
    m is the output mesh size in terms of input mesh size (m = L_out/L_inp)
    for 'ks' domain propagation
        no Fourier transform to frequency domain is done
        assumes no angular dispersion (true for plain FEL radiation)
        assumes narrow spectrum at center of xlamds (true for plain FEL radiation)
    'kf' propagation is a default option
    
    z>0 -> forward direction
    z<0 -> backward direction
    z=0 return original
    """   
    def __init__(self, z0, mx, my, method):
        Mask.__init__(self)
        self.z0 = z0
        self.mx = mx
        self.my = my
        self.method = method
        self.mask = None

    def apply(self, dfl):
        """
        'apply' method for Prop_mMask class 
        
        transform dfl to the domain of the propagation 'kf' or 'kt'
        transforming the wavefront curvature transform by a phase factor z0/(1 - m)
        get the transfer function calling 'get_mask' method, e.g. get H transfer function 
        multiply dfl and the mask, e.g. E = E_0 * H convolution in inverse space domain
        transverce mesh step resizing 
        return the field to the original domain
        transforming the wavefront curvature transform by a phase factor -m*z0/(m - 1) --- inverse transformation
        """   
        _logger.info('propagating dfl file by %.2f meters' % (self.z0))

        if self.z0 == 0:
            _logger.debug(ind_str + 'z0=0, returning original')
            return dfl
        
        start = time.time()
    
        if self.mx != 1:
#            dfl.curve_wavefront(-self.z0 / (1 - self.mx), plane='x')
            dfl = QuadCurvMask(r = -self.z0 / (1 - self.mx), plane='x').apply(dfl)
        if self.my != 1:
#            dfl.curve_wavefront(-self.z0 / (1 - self.my), plane='y')
            dfl = QuadCurvMask(r = -self.z0 / (1 - self.my), plane='y').apply(dfl)

        domains = dfl.domains()      
        self.domain_x = 'k'
        self.domain_y = 'k'
        if self.method in ['PropMask', 'PropMask_kf']:
            dfl.to_domain('kf')
            self.domain_x = 'f'
        elif self.method == 'PropMask_kt':
            dfl.to_domain('kt')
            self.domain_x = 't'
        else:
            raise ValueError("propagation method should be 'PropMask', 'PropMask_kf' or 'PropMask_kt'")
        
        if self.mask is None:
            self.get_mask(dfl) # get H transfer function 
            
        dfl.fld *= self.mask #E = E_0 * H convolution in inverse space domain
               
        dfl.dx *= self.mx #transform grid to output mesh size
        dfl.dy *= self.my        
        self.dx *= self.mx #transform grid to output mesh size
        self.dy *= self.my   
        
        dfl.to_domain(domains) # back to the original domain
        
        if self.mx != 1:
#            dfl.curve_wavefront(-self.mx * self.z0 / (self.mx - 1), plane='x')
            dfl = QuadCurvMask(r = -self.mx * self.z0 / (self.mx - 1), plane='x').apply(dfl)
        if self.my != 1:
#            dfl.curve_wavefront(-self.my * self.z0 / (self.my - 1), plane='y')
            dfl = QuadCurvMask(r = -self.my * self.z0 / (self.my - 1), plane='y').apply(dfl)

        t_func = time.time() - start
        _logger.debug(ind_str + 'done in %.2f sec' % t_func)
        
        return dfl
    
    def get_mask(self, dfl):
        """
        'get_mask' method for Prop_mMask class
        find the transfer function for propagation in 'kf' of 'ks' domains
        H = exp(iz0(k^2 - kx^2 - ky^2)^(1/2) - k) 
        """
        self.copy_grid(dfl)
        
        k_x, k_y = np.meshgrid(self.grid_kx(), self.grid_ky())
        if dfl.domain_z == 'f':
            self.mask = np.ones(dfl.shape(), dtype=complex128)
            k = self.grid_kz()
            if self.mx != 0:
#                self.mask[i,:,:] *= [np.exp(1j * self.z0/self.mx * (np.sqrt(k[i] ** 2 - k_x ** 2) - k[i])) for i in range(dfl.Nz())] #Hx = exp(iz0/mx(k^2 - kx^2)^(1/2) - k)
                for i in range(self.Nz()):
                    self.mask[i,:,:] *= np.exp(1j * self.z0 / self.mx * (np.sqrt(k[i] ** 2 - k_x ** 2) - k[i]))
            if self.my != 0:
#                self.mask[i,:,:] *= [np.exp(1j * self.z0/self.my * (np.sqrt(k[i] ** 2 - k_y ** 2) - k[i])) for i in range(dfl.Nz())] #Hy = exp(iz0/my(k^2 - ky^2)^(1/2) - k)                   
                for i in range(self.Nz()):
                    self.mask[i,:,:] *= np.exp(1j * self.z0 / self.my * (np.sqrt(k[i] ** 2 - k_y ** 2) - k[i]))
        elif dfl.domain_z == 't':
            k = 2 * np.pi / self.xlamds
            Hx = [np.exp(1j * self.z0/self.mx * (np.sqrt(k ** 2 - k_x ** 2) - k)) for i in range(dfl.Nz())] 
            Hy = [np.exp(1j * self.z0/self.my * (np.sqrt(k ** 2 - k_y ** 2) - k)) for i in range(dfl.Nz())]          
            self.mask = Hx*Hy
        else: 
            raise ValueError("wrong field domain, domain must be 'ks' or 'kf'")    
            
        return self.mask    
    
class QuadCurvMask(Mask):
    """
    Quadratic curvature wavefront mask
    add quadratic x and y term in 's' space 
    
    :param r: radius of curveture
    :param plane: plane in which the curvature added
    :param mask: the transfer function itself
    """
    def __init__(self, r, plane):
        Mask.__init__(self)
        self.r = r 
        self.plane = plane 
        self.mask = None 

    def apply(self, dfl):
        """
        'apply' method for QuadCurvMask class      
        
        get the mask using 'get_mask' method
        convert the field to a spatial domain
        multiply the mask and the field
        convert the field to the original domain
        """
        domains = dfl.domain_z, dfl.domain_xy
        
        self.domain_x = 's'
        self.domain_y = 's'
        if self.domain_z is None:
            self.domain_z = dfl.domain_z
        
        _logger.debug('curving radiation wavefront by {}m in {} domain'.format(self.r, self.domain_z))
        
        if np.size(self.r) == 1:
            if self.r == 0:
                raise ValueError('radius of curvature cannot be zero')
            elif self.r == np.inf:
                _logger.debug(ind_str + 'radius of curvature is infinite, skipping')
                return
            else:
                pass
            
        if self.mask is None:
            self.get_mask(dfl) 
    
        dfl.to_domain(self.domain_z + 's')
        dfl.fld *= self.mask 
        dfl.to_domain(domains) 

        return dfl
    
    def get_mask(self, dfl):
        """
        'get_mask' method for QuadCurvMask class      
        
        TODO write a description 
        """        
        self.copy_grid(dfl)
       
        x, y = np.meshgrid(self.grid_x(), self.grid_y())
        if self.plane == 'xy' or self.plane == 'yx':
            arg2 = x ** 2 + y ** 2
        elif self.plane == 'x':
            arg2 = x ** 2
        elif self.plane == 'y':
            arg2 = y ** 2
        else:
            raise ValueError("'plane' should be in 'x', 'y'")

        if self.domain_z == 'f':
            k = 2 * np.pi /dfl.grid_z() 

            if np.size(self.r) == 1:
                self.mask = np.exp(-1j * k[:, np.newaxis, np.newaxis] / 2 * arg2[np.newaxis, :, :] / self.r) #H = exp(-i * k / 2 * (x^2 + y^2))
            elif np.size(self.r) == dfl.Nz():
                self.mask = np.exp(-1j * k[:, np.newaxis, np.newaxis] / 2 * arg2[np.newaxis, :, :] / self.r[:, np.newaxis, np.newaxis])
            else:
                raise ValueError('wrong dimensions of radius of curvature')    

        elif self.domain_z == 't':
            k = 2 * np.pi / dfl.xlamds

            if np.size(self.r) == 1:
                self.mask = np.exp(-1j * k / 2 * arg2 / self.r)
            elif np.size(self.r) == dfl.Nz():
                self.mask = np.exp(-1j * k / 2 * arg2[np.newaxis, :, :] / self.r[:, np.newaxis, np.newaxis])
            else:
                raise ValueError('wrong dimensions of radius of curvature')
        else:
            raise ValueError("domain_z should be in 'f', 't'")
            
        return self.mask

class LensMask(Mask):
    """
    A thin ideal Lens mask

    :param fx: focus length in x direction, [m]
    :param fy: focus length in y direction, [m]
    """
    def __init__(self, fx=np.inf, fy=np.inf):
        Mask.__init__(self)
        self.fx = fx
        self.fy = fy
        self.mask = None
              
    def apply(self, dfl):
        _logger.info('apply the lens mask')
        domains = dfl.domains()

        if self.mask is None:
            self.get_mask(dfl)
                
        dfl.to_domain('fs')
        dfl.fld *= self.mask
        dfl.to_domain(domains)
        
        return dfl
    
    def get_mask(self, dfl):

        _logger.info('get the lens mask')        
        H_fx = QuadCurvMask(r=self.fx, plane='x').get_mask(dfl)        
        H_fy = QuadCurvMask(r=self.fy, plane='y').get_mask(dfl)
        self.mask = H_fx * H_fy
      
        return self.mask
    
class ApertureRectMask(Mask):
    """
     Applies a rectangular appertire for a field
     
    :param lx: horizonatal size of the aperture, [m]
    :param ly: vertical size of the aperture, [m]
    :param cx: x coprdinate coordinate of the aperture center, [m] 
    :param cy: y coprdinate coordinate of the aperture center, [m] 
    :mask: the rectangular aperture transfer function
    """
    def __init__(self, lx=np.inf, ly=np.inf, cx=0, cy=0):
        Mask.__init__(self)
        self.lx = lx
        self.ly = ly
        self.cx = cx
        self.cy = cy
        self.mask = None

    def apply(self, dfl):
        
        domains = dfl.domain_z, dfl.domain_xy
        self.domain_x = 's'
        self.domain_y = 's'
        if self.domain_z is None:
            self.domain_z = dfl.domain_z
        
        dfl.to_domain(self.domain_z + 's')

        if self.mask is None:
            self.get_mask(dfl)
                
        mask_idx = np.where(self.mask == 0)

        dfl_energy_orig = dfl.E()
        dfl.fld[:, mask_idx[0], mask_idx[1]] = 0

        if dfl_energy_orig == 0:
            _logger.warn(ind_str + 'dfl_energy_orig = 0')
        elif dfl.E() == 0:
            _logger.warn(ind_str + 'done, %.2f%% energy lost' % (100))
        else:
            _logger.info(ind_str + 'done, %.2f%% energy lost' % ((dfl_energy_orig - dfl.E()) / dfl_energy_orig * 100))
        
        dfl.to_domain(domains) 
        return dfl

    def get_mask(self, dfl):
        """
        model a rectangular aperture to the radaition in 's' domain
        """
        _logger.info('applying square aperture to dfl')
        
        self.copy_grid(dfl)
            
        if np.size(self.lx) == 1:
            self.lx = [-self.lx / 2, self.lx / 2]
        if np.size(self.ly) == 1:
            self.ly = [-self.ly / 2, self.ly / 2]
        _logger.debug(ind_str + 'ap_x = {}'.format(self.lx))
        _logger.debug(ind_str + 'ap_y = {}'.format(self.ly ))
        print(self.lx, self.ly)
        
        idx_x = np.where((self.grid_x() >= self.lx[0]) & (self.grid_x() <= self.lx[1]))[0]
        idx_x1 = idx_x[0]# + 1 #magic number
        idx_x2 = idx_x[-1]  

        idx_y = np.where((self.grid_y() >= self.ly[0]) & (self.grid_y() <= self.ly[1]))[0]
        idx_y1 = idx_y[0]# + 1 #magic number
        idx_y2 = idx_y[-1]# + 1 #magic number

        _logger.debug(ind_str + 'idx_x = {}-{}'.format(idx_x1, idx_x2))
        _logger.debug(ind_str + 'idx_y = {}-{}'.format(idx_y1, idx_y2))

        self.mask = np.zeros_like(dfl.fld[0, :, :])
        self.mask[idx_y1:idx_y2, idx_x1:idx_x2] = 1
        
        return self.mask

class ApertureEllipsMask(Mask):
    """
     Applies a rectangular appertire for a field
     
    :param ax: ellipse x main axis, [m]
    :param ay: ellipse y main axis, [m]
    :param cx: ellipse x coordinate of center, [m] 
    :param cy: ellipse y coordinate of center, [m]  
    :mask: the elliptical aperture transfer function
    """
    def __init__(self, ax=np.inf, ay=np.inf, cx=0, cy=0):
        Mask.__init__(self)
        self.ax = ax
        self.ay = ay
        self.cx = cx
        self.cy = cy
        self.mask = None 
    
    def ellipse(self, dfl):    
        x, y = np.meshgrid(self.grid_x(), self.grid_y())
        xp =  (x - self.cx)*np.cos(pi) + (y - self.cy)*np.sin(pi)
        yp = -(x - self.cx)*np.sin(pi) + (y - self.cy)*np.cos(pi)
        return (2*xp/self.ax)**2 + (2*yp/self.ay)**2
    
    def apply(self, dfl):
        """
        apply elliptical aperture to the radaition in 's' domain
        """
        
        _logger.info('applying elliptical aperture to dfl')
        _logger.debug(ind_str + 'ap_x = {}'.format(self.ax) + 'cx = {}'.format(self.cx))
        _logger.debug(ind_str + 'ap_y = {}'.format(self.ay) + 'cy = {}'.format(self.cy))

        domains = dfl.domain_z, dfl.domain_xy
        self.domain_x = 's'
        self.domain_y = 's'
        if self.domain_z is None:
            self.domain_z = dfl.domain_z
        
        dfl.to_domain(self.domain_z + 's')
        if self.mask is None:
            self.get_mask(dfl)
        
        mask_idx = np.where(self.mask == 0)

        dfl_energy_orig = dfl.E()
        dfl.fld[:, mask_idx[0], mask_idx[1]] = 0
        
        if dfl_energy_orig == 0:
            _logger.warn(ind_str + 'dfl_energy_orig = 0')
        elif dfl.E() == 0:
            _logger.warn(ind_str + 'done, %.2f%% energy lost' % (100))
        else:
            _logger.info(ind_str + 'done, %.2f%% energy lost' % ((dfl_energy_orig - dfl.E()) / dfl_energy_orig * 100))
        
        dfl.to_domain(domains) 
        return dfl
    
    def get_mask(self, dfl):
        
        self.copy_grid(dfl)
        
        self.mask = np.zeros_like(dfl.fld[0, :, :])
        self.mask[self.ellipse(dfl) <= 1] = 1
        
        return self.mask

class PhaseDelayMask(Mask):
    """
    The function adds a phase shift to a fld object. The expression for the phase 
    -- coeff[0] + coeff[1]*(w - w0)/1! + coeff[2]*(w - w0)**2/2! + coeff[3]*(w - w0)**3/3!
    
    :param: coeff 
        coeff[0] =: measured in [rad]      --- phase
        coeff[1] =: measured in [fm s ^ 1] --- group delay
        coeff[2] =: measured in [fm s ^ 2] --- group delay dispersion (GDD)
        coeff[3] =: measured in [fm s ^ 3] --- third-order dispersion (TOD)
        ...
    :param E_ph0:  energy with respect to which the phase shift is calculated
    :mask: phase delay aperture transfer function
    """
    def __init__(self, coeff, E_ph0):
        Mask.__init__(self)
        self.coeff = coeff
        self.E_ph0 = E_ph0
        self.mask = None
        
    def apply(self, dfl):
        
        _logger.info('apply the frequency chirp')
        domains = dfl.domains()
        
        if self.mask is None: 
            self.get_mask(dfl)
        
        dfl.to_domain('f')
        dfl.fld *= self.mask
        dfl.to_domain(domains)
        
        _logger.info('done')
        return dfl

    def get_mask(self, dfl):

        _logger.info('get the frequency chirp mask')
        
        self.copy_grid(dfl)
        
        if self.E_ph0 == None:
            w0 = 2 * np.pi * speed_of_light / dfl.xlamds
    
        elif self.E_ph0 == 'center':
            _, k = np.mean([dfl.int_z(), dfl.grid_kz()], axis=(1))
            w0 = speed_of_light * k
    
        elif isinstance(self.E_ph0, str) is not True:
            w0 = 2 * np.pi * self.E_ph0 / h_eV_s
    
        else:
            raise ValueError("E_ph0 must be None or 'center' or some number")

        w = self.grid_kz() * speed_of_light
        delta_w = w - w0

        _logger.debug('calculating phase delay')
        _logger.debug(ind_str + 'coeffs for compression = {}'.format(self.coeff))

        coeff_norm = [ci / (1e15) ** i / factorial(i) for i, ci in enumerate(self.coeff)]
        coeff_norm = list(coeff_norm)[::-1]
        _logger.debug(ind_str + 'coeffs_norm = {}'.format(coeff_norm))
        delta_phi = np.polyval(coeff_norm, delta_w)

        _logger.debug(ind_str + 'delta_phi[0] = {}'.format(delta_phi[0]))
        _logger.debug(ind_str + 'delta_phi[-1] = {}'.format(delta_phi[-1]))

        self.mask = np.exp(-1j * delta_phi)[:, np.newaxis, np.newaxis]

        _logger.debug(ind_str + 'done')
        
        return self.mask
        
class MirrorMask(Mask):
    """
    Class for simulating HeightProfile of highly polished mirror surface
    
    TODO
    write documentation
    add logging
    """
    def __init__(self, height_profile=None, hrms=0, angle=2*np.pi/180, plane='x', lx=np.inf, ly=np.inf, eid=None, shape=(0, 0, 0)):
        Mask.__init__(self, shape=shape)
        self.height_profile = height_profile 
        self.hrms = hrms #must have a size equals 2 
        self.angle = angle
        self.mask = None
        self.eid = eid
        
        if plane is 'x':
            self.lx = lx/np.sin(self.angle)
            self.ly = ly
        elif plane is 'y':
            self.lx = lx
            self.ly = ly/np.sin(self.angle)
        else:
            raise ValueError(" 'plane' must be 'x' or 'y' ")

    def apply(self, dfl):
        
        _logger.info('apply HeightProfile errors')
        start = time.time()

        if self.mask is None: 
            self.get_mask(dfl)

        dfl.fld *= self.mask
        
        t_func = time.time() - start
        _logger.debug(ind_str + 'done in {}'.format(t_func))


    def get_mask(self, dfl):
        _logger.info('getting HeightProfile errors')
        
        if self.mask is None:
            heightErr_x = HeightErrorMask_1D(height_profile=self.height_profile, hrms=self.hrms, axis='x', angle=self.angle)
            heightErr_x.get_mask(dfl) 
            heightErr_y = HeightErrorMask_1D(height_profile=self.height_profile, hrms=self.hrms, axis='y', angle=self.angle)
            heightErr_y.get_mask(dfl) 
            RectApp = ApertureRectMask()
            RectApp.lx = self.lx
            RectApp.ly = self.ly           
            RectApp.get_mask(dfl)
            self.mask = heightErr_x.mask * heightErr_y.mask * RectApp.mask
        return self.mask
        
        
class HeightErrorMask_1D(Mask):
    """
    Mask for simulating HeightProfile of highly polished mirror surface

    :param hrms: [meters] height errors root mean square
    :param length: [meters] length of the surface
    :param points_number: number of points (pixels) at the surface
    :param wavevector_cutoff: [1/meters] point on k axis for cut off small wavevectors (large wave lengths) in the PSD
                                    (with default value 0 effects on nothing)
    :param psd: [meters^3] 1d array; power spectral density of surface (if not specified, will be generated)
            (if specified, must have shape = (points_number // 2 + 1, ), otherwise it will be cut to appropriate shape)
    :param seed: seed for np.random.seed() to allow reproducibility
    :return: HeightProfile object
    """

    def __init__(self, height_profile=None, hrms=0, axis='x', angle=np.pi * 2 / 180, seed=None):
        Mask.__init__(self)
        self.mask = None
        self.height_profile = height_profile 
        self.hrms = hrms
        self.axis = axis
        self.angle = angle
        self.seed = seed
        
    def apply(self, dfl):
        
        _logger.info('apply HeightProfile errors')
        start = time.time()

        if self.mask is None: 
            self.get_mask(dfl)

        dfl.fld *= self.mask
        
        t_func = time.time() - start
        _logger.debug(ind_str + 'done in {}'.format(t_func))

        return dfl
        
    def get_mask(self, dfl):
        
        self.copy_grid(dfl)
        
        dict_axes = {'z': 0, 'y': 1, 'x': 2}
        dl = {0: dfl.dz, 1: dfl.dy, 2: dfl.dx}
        
        if isinstance(self.axis, str):
            axis = dict_axes[self.axis]
 
        n = dfl.fld.shape[axis]
        eff_l = n * dl[axis] / np.sin(self.angle)
        if self.height_profile is None:
            if self.hrms is None:
                _logger.error('hrms and height_profile not specified')
                raise ValueError('hrms and height_profile not specified')
            self.height_profile = HeightProfile()
            self.height_profile = self.height_profile.generate_1d_profile(self.hrms, L=eff_l, N=n, seed=self.seed)
        
        elif self.height_profile.L != eff_l or self.height_profile.N != n:
            if self.height_profile.L < eff_l:
                _logger.warning(
                    'provided height profile length {} is smaller than projected footprint {}'.format(height_profile.L,
                                                                                                      eff_l))  # warn and pad with zeroes
    
            # interpolation of height_profile to appropriate sizes
            _logger.info(ind_str + 'given height_profile was interpolated to length and '.format(dfl.shape()))
            s = np.linspace(-eff_l / 2, eff_l / 2, n)
            h = np.interp(s, height_profile.s, height_profile.h, right=height_profile.h[0], left=height_profile.h[-1])
            self.height_profile = HeightProfile()
            self.height_profile.points_number = dfl.fld.shape[axis]
            self.height_profile.L = eff_l
            self.height_profile.s = s
            self.height_profile.h = h
    
        phase_delay = 2 * 2 * np.pi * np.sin(self.angle) * self.height_profile.h / dfl.xlamds
        
        if self.axis == 'x':
            self.mask = np.exp(1j * phase_delay)[np.newaxis, np.newaxis, :]
        elif self.axis == 'y':
            self.mask = np.exp(1j * phase_delay)[np.newaxis, :, np.newaxis]
        elif self.axis == 'z':
            self.mask = np.exp(1j * phase_delay)[:, np.newaxis, np.newaxis]             
        return self.mask
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    