import logging
import time

import numpy as np
from copy import deepcopy
from math import factorial

from ocelot.common.globals import *
from ocelot import ocelog
from ocelot.common.ocelog import *
_logger = logging.getLogger(__name__) 

class Grid:
    
    def __init__(self, shape=(0, 0, 0)):
        self.dx = []
        self.dy = []
        self.dz = []
        self.shape = shape 
        self.xlamds = None #picky point!!!
    def Lz(self):  
        '''
        full longitudinal mesh size
        '''
        return self.dz * self.Nz()

    def Ly(self):  
        '''
        full transverse vertical mesh size
        '''
        return self.dy * self.Ny()

    def Lx(self):  
        '''
        full transverse horizontal mesh size
        '''
        return self.dx * self.Nx()

    def Nz(self):
        '''
        number of points in z
        '''
        return self.shape[0]

    def Ny(self):
        '''
        number of points in y
        '''
        return self.shape[1]

    def Nx(self):
        '''
        number of points in x
        '''
        return self.shape[2]
    
    def x(self):
        return np.linspace(-self.Lx() / 2, self.Lx() / 2, self.Nx())
    
    def kx(self):
        k = 2 * np.pi / self.dx
        return np.linspace(-k / 2, k / 2, self.Nx())
    
    def y(self):
        return np.linspace(-self.Ly() / 2, self.Ly() / 2, self.Ny())
    
    def ky(self):
        k = 2 * np.pi / self.dy
        return np.linspace(-k / 2, k / 2, self.Ny())

    def z(self):
        return np.linspace(0, self.Lz(), self.Nz())

    def kz(self):
        dk = 2 * np.pi / self.Lz()
        k = 2 * np.pi / self.xlamds #picky point!!!
        return np.linspace(k - dk / 2 * self.Nz(), k + dk / 2 * self.Nz(), self.Nz())
    
class Mask(Grid):
    
    def __init__(self, shape=(0, 0, 0)):
        Grid.__init__(self, shape=shape)

    def apply(self, dfl):
        return dfl

    def get_mask(self, dfl):
        pass

    def __mul__(self, other):
        m = deepcopy(self)
        if other.__class__ in [self] and self.mask is not None and other.mask is not None:
            m.mask = self.mask * other.mask
            return m

class RectMask(Mask):
    def __init__(self, shape=(0, 0, 0)):
        Mask.__init__(self, shape=shape)
        self.lx = np.inf
        self.ly = np.inf
        self.cx = 0
        self.cy = 0
        self.mask = None

    def apply(self, dfl):
        if self.mask is None:
            self.get_mask(dfl)
        mask_idx = np.where(self.mask == 0)

        # dfl_out = deepcopy(dfl)
        dfl_energy_orig = dfl.E()
        dfl.fld[:, mask_idx[0], mask_idx[1]] = 0

        if dfl_energy_orig == 0:
            _logger.warn(ind_str + 'dfl_energy_orig = 0')
        elif dfl.E() == 0:
            _logger.warn(ind_str + 'done, %.2f%% energy lost' % (100))
        else:
            _logger.info(ind_str + 'done, %.2f%% energy lost' % ((dfl_energy_orig - dfl.E()) / dfl_energy_orig * 100))
        # tmp_fld = dfl.fld[:,idx_x1:idx_x2,idx_y1:idx_y2]
        return dfl

    def get_mask(self, dfl):
        """
        model rectangular aperture to the radaition in either domain
        """
        _logger.info('applying square aperture to dfl')
        
        self.dx = dfl.dx
        self.dy = dfl.dy
        self.dz = dfl.dz
        self.shape = dfl.shape()
            
        if np.size(self.lx) == 1:
            self.lx = [-self.lx / 2, self.lx / 2]
        if np.size(self.ly) == 1:
            self.ly = [-self.ly / 2, self.ly / 2]
        _logger.debug(ind_str + 'ap_x = {}'.format(self.lx))
        _logger.debug(ind_str + 'ap_y = {}'.format(self.ly ))

        idx_x = np.where((self.x() >= self.lx[0]) & (self.x() <= self.lx[1]))[0]
        idx_x1 = idx_x[0]
        idx_x2 = idx_x[-1]

        idx_y = np.where((self.y() >= self.ly [0]) & (self.y() <= self.ly [1]))[0]
        idx_y1 = idx_y[0]
        idx_y2 = idx_y[-1]

        _logger.debug(ind_str + 'idx_x = {}-{}'.format(idx_x1, idx_x2))
        _logger.debug(ind_str + 'idx_y = {}-{}'.format(idx_y1, idx_y2))

        self.mask = np.zeros_like(dfl.fld[0, :, :])
        self.mask[idx_y1:idx_y2, idx_x1:idx_x2] = 1
        return self.mask

    def __mul__(self, other):
        m = RectMask()
        if other.__class__ in [RectMask, ] and self.mask is not None and other.mask is not None:
            m.mask = self.mask * other.mask
            return m

class EllipsMask(Mask):   
    def __init__(self, shape):
        Mask.__init__(self, shape=shape)
        self.ax = np.inf
        self.ay = np.inf
        self.cx = 0
        self.cy = 0
        self.mask = None 
    
    def ellipse(self, dfl):    
        x, y = np.meshgrid(self.x(), self.y())
        xp =  (x - self.cx)*np.cos(pi) + (y - self.cy)*np.sin(pi)
        yp = -(x - self.cx)*np.sin(pi) + (y - self.cy)*np.cos(pi)
        return (2*xp/self.ax)**2 + (2*yp/self.ay)**2
    
    def apply(self, dfl):
        """
        apply elliptical aperture to the radaition in either domain
        """
        
        _logger.info('applying elliptical aperture to dfl')
        _logger.debug(ind_str + 'ap_x = {}'.format(self.ax) + 'cx = {}'.format(self.cx))
        _logger.debug(ind_str + 'ap_y = {}'.format(self.ay) + 'cy = {}'.format(self.cy))

        
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
            
        return dfl
    
    def get_mask(self, dfl):
        
        self.dx = dfl.dx
        self.dy = dfl.dy
        self.dz = dfl.dz
        self.shape = dfl.shape()
        
        self.mask = np.zeros_like(dfl.fld[0, :, :])
        self.mask[self.ellipse(dfl) <= 1] = 1
        return self.mask
        
        
    
class MaxwPropagator_m(Mask):
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
    def __init__(self, z0, mx, my):
        Mask.__init__(self)
        self.z0 = z0
        self.mx = mx
        self.my = my
        self.mask = None
        self.domains = 'kf'

    def apply(self, dfl):
        
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

        if self.domains == 'kf' or self.domains == 'kf' or self.domains == 'k':
            dfl.to_domain(self.domains) #the field is transformed in inverce space domain and, optionaly, in 'f' or 't' domains
        else:
            raise ValueError("domains value should be 'kf' or 'kt' or 'k'")
        
        if self.mask is None:
            self.get_mask(dfl) # get H transfer function 
            
        dfl.fld *= self.mask #E = E_0 * H convolution in inverse space domain
        
        dfl.dx *= self.mx #transform grid to output mesh size
        dfl.dy *= self.my        
        
        dfl.to_domain(domains) # back to original domain
        
        if self.mx != 1:
#            dfl.curve_wavefront(-self.mx / (self.mx - 1), plane='x')
            dfl = QuadCurvMask(r = -self.mx * self.z0 / (self.mx - 1), plane='x').apply(dfl)
        if self.my != 1:
#            dfl.curve_wavefront(-self.my / (self.my - 1), plane='y')
            dfl = QuadCurvMask(r = -self.my * self.z0 / (self.my - 1), plane='y').apply(dfl)

        t_func = time.time() - start
        _logger.debug(ind_str + 'done in %.2f sec' % t_func)
        
        return dfl
        
    def get_mask(self, dfl):
     
        k_x, k_y = np.meshgrid(self.kx(), self.ky())
        k = self.kz()

        if self.domains == 'kf':
            k = self.kz()
            Hx = [np.exp(1j * self.z0/self.mx * (np.sqrt(k[i] ** 2 - k_x ** 2) - k[i])) for i in range(dfl.Nz())][0] #Hx = exp(iz0/mx(k^2 - kx^2)^(1/2) - k)
            Hy = [np.exp(1j * self.z0/self.my * (np.sqrt(k[i] ** 2 - k_y ** 2) - k[i])) for i in range(dfl.Nz())][0] #Hy = exp(iz0/my(k^2 - ky^2)^(1/2) - k)                  
            self.mask = Hx*Hy
        elif self.domains == 'ks':
            k = 2 * np.pi / dfl.xlamds
            Hx = [np.exp(1j * self.z0/self.mx * (np.sqrt(k ** 2 - k_x ** 2) - k)) for i in range(dfl.Nz())][0] 
            Hy = [np.exp(1j * self.z0/self.my * (np.sqrt(k ** 2 - k_y ** 2) - k)) for i in range(dfl.Nz())][0]          
            self.mask = Hx*Hy
        else: 
            raise ValueError('wrong field domain, domain must be ks or kf ')    
            
        return self.mask
    
class MaxwPropagator(Mask):
    """
    Angular-spectrum propagation for fieldfile
    
    can handle wide spectrum
      (every slice in freq.domain is propagated
       according to its frequency)
    no kx**2+ky**2<<k0**2 limitation
    
    dfl is the RadiationField() object
    
    z0 is the propagation distance in [m]
    for 'ks' domain propagation
        no Fourier transform to frequency domain is done
        assumes no angular dispersion (true for plain FEL radiation)
        assumes narrow spectrum at center of xlamds (true for plain FEL radiation)
    'kf' propagation is a default option
    
    z>0 -> forward direction
    z<0 -> backward direction
    z=0 return original
    """   
    def __init__(self, z0):
        Mask.__init__(self)
        self.z0 = z0 
        self.mask = None
        self.domains = 'kf'

    def apply(self, dfl):
        _logger.info('propagating dfl file by %.2f meters' % (self.z0))
        
        if self.z0 == 0:
            _logger.debug(ind_str + 'z0=0, returning original')
            return dfl

        start = time.time()

        domains = dfl.domains()

        if self.domains == 'kf' or self.domains == 'kf' or self.domains == 'k':
            dfl.to_domain(self.domains) #the field is transformed in inverce space domain and, optionaly, in 'f' or 't' domains
        else:
            raise ValueError("domains value should be 'kf' or 'kt' or 'k'")
            
        if self.mask is None:
            self.get_mask(dfl) # get H transfer function 

        dfl.fld *= self.mask # E = E_0 * H convolution in inverse space domain
        
        dfl.to_domain(domains) # back to original domain !!! mast be modified for an optimization procedure?

        t_func = time.time() - start
        _logger.debug(ind_str + 'done in %.2f sec' % t_func)
        
        return dfl
        
    def get_mask(self, dfl):
     
        k_x, k_y = np.meshgrid(self.kx(), self.ky())
        
        if dfl.domain_z is 'f':
            k = self.kz()
            self.mask = [np.exp(1j * self.z0 * (np.sqrt(k[i] ** 2 - k_x ** 2 - k_y ** 2) - k[i])) for i in range(dfl.Nz())] #H = exp(iz0(k^2 - kx^2 - ky^2)^(1/2) - k)
        else:
            k = 2 * np.pi / dfl.xlamds
            self.mask = [np.exp(1j * self.z0 * (np.sqrt(k ** 2 - k_x ** 2 - k_y ** 2) - k)) for i in range(dfl.Nz())]  #H = exp(iz0(k^2 - kx^2 - ky^2)^(1/2) - k)

        return self.mask


class DriftMask(Mask):

    def __init__(self, z0, mx, my):
        Mask.__init__(self)
        self.z0 = z0
        self.mx = mx #
        self.my = my #
        self.mask = None
        self.type = 'MaxwPropagator'  #type of propagation. also may be 
        # 'Fraunhofer_Propagator'
        # 'Fresnel_Propagator' . . .
    def apply(self, dfl): 
        
        if self.mask is None:         
            if self.type == 'MaxwPropagator' and self.mx == 1 and self.my == 1:
                dfl = MaxwPropagator(z0 = self.z0).apply(dfl)
            elif self.type == 'MaxwPropagator' and (self.mx != 1 or self.my != 1):
                dfl = MaxwPropagator_m(z0 = self.z0, mx = self.mx, my = self.my).apply(dfl)
            elif self.type == 'Fraunhofer_Propagator':
                pass
            elif self.type == 'Fresnel_Propagator':
                pass
            else:
                raise ValueError("check a propagation type")
            
        return dfl
        
    def get_mask(self, dfl):
        if self.type == 'MaxwPropagator' and self.mx == 1 and self.my == 1:
            self.mask = MaxwPropagator(z0 = self.z0).get(dfl)
        elif self.type == 'MaxwPropagator' and self.mx != 1 and self.my != 1:
            self.mask = MaxwPropagator_m(z0 = self.z0, mx = self.mx, my = self.my).get(dfl)
        elif self.type == 'Fraunhofer_Propagator':
            pass
        elif self.type == 'Fresnel_Propagator':
            pass
        else:
            raise ValueError("check a propagation type")
     
        
        return self.mask

class QuadCurvMask(Mask):

    def __init__(self, r, plane):
        Mask.__init__(self)
        self.r = r #is the radius of curvature
        self.plane = plane #is the plane in which the curvature is applied
        self.mask = None #is the transfer function itself
        self.domains = 's' #is the domain in which the transfer function is calculated
        self.domain_z = None #is the domain in which wavefront curvature is introduced 

    def apply(self, dfl):
        
        domains = dfl.domain_z, dfl.domain_xy

        if self.domain_z is None:
            self.domain_z = dfl.domain_z
        
        _logger.debug('curving radiation wavefront by {}m in {} domain'.format(self.r, self.domain_z))
        
        if np.size(self.r) == 1:
            if self.r == 0:
                raise ValueError('radius of curvature should not be zero')
            elif self.r == np.inf:
                _logger.debug(ind_str + 'radius of curvature is infinite, skipping')
                return
            else:
                pass
            
        if self.mask is None:
            self.get_mask(dfl) 
    
        dfl.to_domain(self.domain_z + 's') #the field must be in 's' domain
        dfl.fld *= self.mask 
        dfl.to_domain(domains) #return to the original domain

        return dfl
    
    def get_mask(self, dfl):
        
        x, y = np.meshgrid(self.x(), self.y())
        if self.plane == 'xy' or self.plane == 'yx':
            arg2 = x ** 2 + y ** 2
        elif self.plane == 'x':
            arg2 = x ** 2
        elif self.plane == 'y':
            arg2 = y ** 2
        else:
            raise ValueError('"plane" should be in ["x", "y"]')

        if self.domain_z == 'f':
            k = 2 * np.pi /dfl.scale_z() # <- change on somethenin' more reliable

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
            ValueError('domain_z should be in ["f", "t", None]')
            
        return self.mask
    

class LensMask(Mask):

    def __init__(self, fx, fy):
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
        H_fx = QuadCurvMask(self.fx, plane='x')
        H_fy = QuadCurvMask(self.fy, plane='y')
        H_fx.domain_z = 'f'
        H_fy.domain_z = 'f'
        
        H_fx.get_mask(dfl)
        H_fy.get_mask(dfl)
        
        self.mask = H_fy.mask*H_fx.mask
         
        return self.mask
    
        
class PhaseDelayMask(Mask):
    """
    The function adds a phase shift to a fld object. The expression for the phase see in the calc_phase_delay function
    dfl   --- is a fld object
    coeff --- coefficients in phase (see in the calc_phase_delay function)
    E_ph0 --- energy with respect to which the phase shift is calculated
    """
    def __init__(self, coeff):
        Mask.__init__(self)
        self.coeff = coeff
        self.E_ph0 = None
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
        """
        expression for the phase -- coeff[0] + coeff[1]*(w - w0)/1! + coeff[2]*(w - w0)**2/2! + coeff[3]*(w - w0)**3/3!
        coeff is a list with
        coeff[0] =: measured in [rad]      --- phase
        coeff[1] =: measured in [fm s ^ 1] --- group delay
        coeff[2] =: measured in [fm s ^ 2] --- group delay dispersion (GDD)
        coeff[3] =: measured in [fm s ^ 3] --- third-order dispersion (TOD)
        ...
        """
        _logger.info('get the frequency chirp mask')
        
        self.dx = dfl.dx
        self.dy = dfl.dy
        self.dz = dfl.dz
        self.shape = dfl.shape()
        self.xlamds = dfl.xlamds
        
        if self.E_ph0 == None:
            w0 = 2 * np.pi * speed_of_light / dfl.xlamds
    
        elif self.E_ph0 == 'center':
            _, lamds = np.mean([dfl.int_z(), dfl.scale_z()], axis=(1))
            w0 = 2 * np.pi * speed_of_light / lamds
    
        elif isinstance(self.E_ph0, str) is not True:
            w0 = 2 * np.pi * self.E_ph0 / h_eV_s
    
        else:
            raise ValueError("E_ph0 must be None or 'center' or some value")

        w = self.kz() * speed_of_light
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
        #####################################
'''
class HightErrorMask(Mask):
    """
    Class for generating HeightProfile of highly polished mirror surface

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

    def __init__(self):
        Mask.__init__(self)
        self.mask = None
        self.hrms = None
        self.axis = 'x'
        self.seed = None
        
    def apply(self, dfl):
        
        _logger.info('apply HeightProfile errors')
        domains = dfl.domains()
        start = time.time()

        if self.mask is None: 
            self.get_mask(dfl)
        
        if self.axis == 'x':
            dfl.fld *= self.mask[np.newaxis, np.newaxis, :]
        elif self.axis == 'y':
            dfl.fld *= self.mask[np.newaxis, :, np.newaxis]
        elif self.axis == 'z':
            dfl.fld *= self.mask[:, np.newaxis, np.newaxis]     

         t_func = time.time() - start
        _logger.debug(ind_str + 'done in {}'.format(t_func))

        return dfl
        
    def get_mask(self, dfl):
 
        dict_axes = {'z': 0, 'y': 1, 'x': 2}
        dlength = {0: dfl.dz, 1: dfl.dy, 2: dfl.dx}
        
        if isinstance(self.axis, str):
            axis = dict_axes[self.axis]
 
        points_number = dfl.fld.shape[axis]
        footprint_len = points_number * dlength[axis] / np.sin(angle)
        
        if height_profile is None:
            if hrms is None:
                _logger.error('hrms and height_profile not specified')
                raise ValueError('hrms and height_profile not specified')
            height_profile = HeightProfile()
            height_profile.generate_1d_profile(hrms, length=footprint_len, points_number=points_number, seed=seed)
        
        elif height_profile.length != footprint_len or height_profile.points_number != points_number:
            pass
    
        phase_delay = 2 * 2 * np.pi * np.sin(angle) * height_profile.h / dfl.xlamds

        return np.exp(1j * phase_delay)
'''
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    