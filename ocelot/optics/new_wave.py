"""
wave optics
"""

from numpy import random
from numpy.linalg import norm
import numpy as np
from math import factorial
from numpy import inf, complex128, complex64
import scipy
import numpy.fft as fft
from copy import deepcopy
import time
import os

# from ocelot.optics.elements import *
from ocelot.common.globals import *
from ocelot.common.math_op import find_nearest_idx, fwhm, std_moment, bin_scale, bin_array, mut_coh_func
from ocelot.common.py_func import filename_from_path
# from ocelot.optics.utils import calc_ph_sp_dens
# from ocelot.adaptors.genesis import *
# import ocelot.adaptors.genesis as genesis_ad
# GenesisOutput = genesis_ad.GenesisOutput
from ocelot.common.ocelog import *
_logger = logging.getLogger(__name__)

import multiprocessing
nthread = multiprocessing.cpu_count()

try:
    import pyfftw
    fftw_avail = True
except ImportError:
    print("wave.py: module PYFFTW is not installed. Install it if you want speed up dfl wavefront calculations")
    fftw_avail = False

__author__ = "Svitozar Serkez, Andrei Trebushinin, Mykola Veremchuk"

### just must to be here for generating dfl :_)

def generate_gaussian_dfl(xlamds=1e-9, shape=(51, 51, 100), dgrid=(1e-3, 1e-3, 50e-6), power_rms=(0.1e-3, 0.1e-3, 5e-6),
                          power_center=(0, 0, None), power_angle=(0, 0), power_waistpos=(0, 0), wavelength=None,
                          zsep=None, freq_chirp=0, en_pulse=None, power=1e6, **kwargs):
    """
    generates RadiationField object
    narrow-bandwidth, paraxial approximations

    xlamds [m] - central wavelength
    shape (x,y,z) - shape of field matrix (reversed) to dfl.fld
    dgrid (x,y,z) [m] - size of field matrix
    power_rms (x,y,z) [m] - rms size of the radiation distribution (gaussian)
    power_center (x,y,z) [m] - position of the radiation distribution
    power_angle (x,y) [rad] - angle of further radiation propagation
    power_waistpos (x,y) [m] downstrean location of the waist of the beam
    wavelength [m] - central frequency of the radiation, if different from xlamds
    zsep (integer) - distance between slices in z as zsep*xlamds
    freq_chirp dw/dt=[1/fs**2] - requency chirp of the beam around power_center[2]
    en_pulse, power = total energy or max power of the pulse, use only one
    """

    start = time.time()

    if dgrid[2] is not None and zsep is not None:
        if shape[2] == None:
            shape = (shape[0], shape[1], int(dgrid[2] / xlamds / zsep))
        else:
            _logger.error(ind_str + 'dgrid[2] or zsep should be None, since either determines longiduninal grid size')

    _logger.info('generating radiation field of shape (nz,ny,nx): ' + str(shape))
    if 'energy' in kwargs:
        _logger.warn(ind_str + 'rename energy to en_pulse, soon arg energy will be deprecated')
        en_pulse = kwargs.pop('energy', 1)

    dfl = RadiationField((shape[2], shape[1], shape[0]))

    k = 2 * np.pi / xlamds

    dfl.xlamds = xlamds
    dfl.domain_z = 't'
    dfl.domain_xy = 's'
    dfl.dx = dgrid[0] / dfl.Nx()
    dfl.dy = dgrid[1] / dfl.Ny()

    if dgrid[2] is not None:
        dz = dgrid[2] / dfl.Nz()
        zsep = int(dz / xlamds)
        if zsep == 0:
            _logger.warning(
                ind_str + 'dgrid[2]/dfl.Nz() = dz = {}, which is smaller than xlamds = {}. zsep set to 1'.format(dz,
                                                                                                                 xlamds))
            zsep = 1
        dfl.dz = xlamds * zsep
    elif zsep is not None:
        dfl.dz = xlamds * zsep
    else:
        _logger.error('dgrid[2] or zsep should be not None, since they determine longiduninal grid size')

    rms_x, rms_y, rms_z = power_rms  # intensity rms [m]
    _logger.debug(ind_str + 'rms sizes = [{}, {}, {}]m (x,y,z)'.format(rms_x, rms_y, rms_z))
    xp, yp = power_angle
    x0, y0, z0 = power_center
    zx, zy = power_waistpos

    if z0 == None:
        z0 = dfl.Lz() / 2

    xl = np.linspace(-dfl.Lx() / 2, dfl.Lx() / 2, dfl.Nx())
    yl = np.linspace(-dfl.Ly() / 2, dfl.Ly() / 2, dfl.Ny())
    zl = np.linspace(0, dfl.Lz(), dfl.Nz())
    z, y, x = np.meshgrid(zl, yl, xl, indexing='ij')

    qx = 1j * np.pi * (2 * rms_x) ** 2 / xlamds + zx
    qy = 1j * np.pi * (2 * rms_y) ** 2 / xlamds + zy
    qz = 1j * np.pi * (2 * rms_z) ** 2 / xlamds

    if wavelength.__class__ in [list, tuple, np.ndarray] and len(wavelength) == 2:
        domega = 2 * np.pi * speed_of_light * (1 / wavelength[0] - 1 / wavelength[1])
        dt = (z[-1, 0, 0] - z[0, 0, 0]) / speed_of_light
        freq_chirp = domega / dt / 1e30 / zsep
        # freq_chirp = (wavelength[1] - wavelength[0]) / (z[-1,0,0] - z[0,0,0])
        _logger.debug(ind_str + 'difference wavelengths {} {}'.format(wavelength[0], wavelength[1]))
        _logger.debug(ind_str + 'difference z {} {}'.format(z[-1, 0, 0], z[0, 0, 0]))
        _logger.debug(ind_str + 'd omega {}'.format(domega))
        _logger.debug(ind_str + 'd t     {}'.format(dt))
        _logger.debug(ind_str + 'calculated chirp {}'.format(freq_chirp))
        wavelength = np.mean([wavelength[0], wavelength[1]])

    if wavelength == None and xp == 0 and yp == 0:
        phase_chirp_lin = 0
    elif wavelength == None:
        phase_chirp_lin = x * np.sin(xp) + y * np.sin(yp)
    else:
        phase_chirp_lin = (z - z0) / dfl.dz * (dfl.xlamds - wavelength) / wavelength * xlamds * zsep + x * np.sin(
            xp) + y * np.sin(yp)

    if freq_chirp == 0:
        phase_chirp_quad = 0
    else:
        # print(dfl.scale_z() / speed_of_light * 1e15)
        # phase_chirp_quad = freq_chirp *((z-z0)/dfl.dz*zsep)**2 * xlamds / 2# / pi**2
        phase_chirp_quad = freq_chirp / (speed_of_light * 1e-15) ** 2 * (zl - z0) ** 2 * dfl.xlamds  # / pi**2
        # print(phase_chirp_quad.shape)

    # if qz == 0 or qz == None:
    #     dfl.fld = np.exp(-1j * k * ( (x-x0)**2/2/qx + (y-y0)**2/2/qy - phase_chirp_lin + phase_chirp_quad ) )
    # else:
    arg = np.zeros_like(z).astype('complex128')
    if qx != 0:
        arg += (x - x0) ** 2 / 2 / qx
    if qy != 0:
        arg += (y - y0) ** 2 / 2 / qy
    if abs(qz) == 0:
        idx = abs(zl - z0).argmin()
        zz = -1j * np.ones_like(arg)
        zz[idx, :, :] = 0
        arg += zz
    else:
        arg += (z - z0) ** 2 / 2 / qz
        # print(zz[:,25,25])

    if np.size(phase_chirp_lin) > 1:
        arg -= phase_chirp_lin
    if np.size(phase_chirp_quad) > 1:
        arg += phase_chirp_quad[:, np.newaxis, np.newaxis]
    dfl.fld = np.exp(-1j * k * arg)  # - (grid[0]-z0)**2/qz
    # dfl.fld = np.exp(-1j * k * ( (x-x0)**2/2/qx + (y-y0)**2/2/qy + (z-z0)**2/2/qz - phase_chirp_lin + phase_chirp_quad) ) #  - (grid[0]-z0)**2/qz

    if en_pulse != None and power == None:
        dfl.fld *= np.sqrt(en_pulse / dfl.E())
    elif en_pulse == None and power != None:
        dfl.fld *= np.sqrt(power / np.amax(dfl.int_z()))
    else:
        _logger.error('Either en_pulse or power should be defined')
        raise ValueError('Either en_pulse or power should be defined')

    dfl.filePath = ''

    t_func = time.time() - start
    _logger.debug(ind_str + 'done in %.2f sec' % (t_func))

    return dfl


class Grid:
    def __init__(self, shape=(0,0,0)):
        self.dx = []
        self.dy = []
        self.dz = []
        self.shapes = shape
        
        self.xlamds = None
        self.used_aprox = 'SVEA'
    
    def copy_grid(self, other, version=2):
    
        if version == 1:
            self.dx = other.dx
            self.dy = other.dy
            self.dz = other.dz
            self.shape = other.shape
            
            self.xlamds = other.xlamds
            self.used_aprox = other.used_aprox
            
        elif version == 2: #copy the same attributes of Mask and RadiationField objects
            attr_list = np.intersect1d(dir(self),dir(other))
            for attr in attr_list:
                if attr.startswith('__') or callable(getattr(self, attr)):
                    continue
                setattr(self, attr, getattr(other, attr))
        else:
            raise ValueError        
    
    def shape(self):
        return self.shapes
    
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
        return self.shapes[0]

    def Ny(self):
        '''
        number of points in y
        '''
        return self.shapes[1]

    def Nx(self):
        '''
        number of points in x
        '''
        return self.shapes[2]
    
    def grid_x(self):
        return np.linspace(-self.Lx() / 2, self.Lx() / 2, self.Nx())
    
    def grid_kx(self):
        k = 2 * np.pi / self.dx
        return np.linspace(-k / 2, k / 2, self.Nx())
    
    def grid_y(self):
        return np.linspace(-self.Ly() / 2, self.Ly() / 2, self.Ny())
    
    def grid_ky(self):
        k = 2 * np.pi / self.dy
        return np.linspace(-k / 2, k / 2, self.Ny())

    def grid_t(self):
        return np.linspace(0, self.Lz()/speed_of_light, self.Nz())
    
    def grid_w(self):
        dw = 2 * np.pi * speed_of_light / self.Lz() # self.Lz() must be replaced with self.Tz()
        
        if self.xlamds is None:        
            return np.linspace(-dw / 2 * self.Nz(), dw / 2 * self.Nz(), self.Nz())
        
        elif self.used_aprox == 'SVEA' and self.xlamds is not None:
            w0 = 2 * np.pi * speed_of_light / self.xlamds
            return np.linspace(w0 - dw / 2 * self.Nz(), w0 + dw / 2 * self.Nz(), self.Nz())
        
        else:
            raise ValueError

    def grid_z(self):
        return self.grid_t() * speed_of_light

    def grid_kz(self):
        return self.grid_w() / speed_of_light
            
class RadiationField(Grid):
    """
    3d or 2d coherent radiation distribution, *.fld variable is the same as Genesis dfl structure
    """

    def __init__(self, shape=(0,0,0)):
        super().__init__(shape=shape)
        self.fld = np.zeros(shape, dtype=complex128)  # (z,y,x)
        self.xlamds = None    # carrier wavelength [m]
        self.domain_z = 't'   # longitudinal domain (t - time, f - frequency)
        self.domain_x = 's'   # transverse domain (s - space, k - inverse space)
        self.domain_y = 's'   # transverse domain (s - space, k - inverse space)
        self.domain_xy = 's'  # transverse domain (s - space, k - inverse space)
        self.filePath = ''
                      
    def fileName(self):
        return filename_from_path(self.filePath)

    def copy_param(self, dfl1, version=1):
        if version == 1:
            self.dx = dfl1.dx
            self.dy = dfl1.dy
            self.dz = dfl1.dz
            self.xlamds = dfl1.xlamds
            self.domain_z = dfl1.domain_z
            self.domain_x = dfl1.domain_x
            self.domain_y = dfl1.domain_y

            self.domain_xy = dfl1.domain_xy
            self.filePath = dfl1.filePath
        elif version == 2: #does it link the address of these two objects only? : _) then exactly what we need for grid copying
            attr_list = dir(dfl1)
            for attr in attr_list:
                if attr.startswith('__') or callable(getattr(self, attr)):
                    continue
                if attr == 'fld':
                    continue
                setattr(self, attr, getattr(dfl1, attr))

    def __getitem__(self, i):
        return self.fld[i]

    def __setitem__(self, i, fld):
        self.fld[i] = fld

    def shape(self):
        '''
        returns the shape of fld attribute
        '''
        return self.fld.shape

    def domains(self):
        '''
        returns domains of the radiation field
        '''
        return self.domain_z, self.domain_xy

    def intensity(self):
        '''
        3d intensity, abs(fld)**2
        '''
        return self.fld.real ** 2 + self.fld.imag ** 2 # calculates faster

    def int_z(self):
        '''
        intensity projection on z
        power [W] or spectral density [arb.units]
        '''
        return np.sum(self.intensity(), axis=(1, 2))

    def ang_z_onaxis(self):
        '''
        on-axis phase
        '''
        xn = int((self.Nx() + 1) / 2)
        yn = int((self.Ny() + 1) / 2)
        fld = self[:, yn, xn]
        return np.angle(fld)

    def int_y(self):
        '''
        intensity projection on y
        '''
        return np.sum(self.intensity(), axis=(0, 2))

    def int_x(self):
        '''
        intensity projection on x
        '''
        return np.sum(self.intensity(), axis=(0, 1))

    def int_xy(self):
        # return np.swapaxes(np.sum(self.intensity(), axis=0), 1, 0)
        return np.sum(self.intensity(), axis=0)

    def int_zx(self):
        return np.sum(self.intensity(), axis=1)

    def int_zy(self):
        return np.sum(self.intensity(), axis=2)

    def E(self):
        '''
        energy in the pulse [J]
        '''
        if self.Nz() > 1:
            return np.sum(self.intensity()) * self.Lz() / self.Nz() / speed_of_light
        else:
            return np.sum(self.intensity())

#   old scales for versions compatibility
#   propper scales in meters or 2 pi / meters
    def scale_kx(self):  # scale in meters or meters**-1
        _logger.warning('"scale_kx" will be deprecated, use "grid_x and grid_kx" instead')
        if 's' in [self.domain_xy, self.domain_x]:    # space domain
            return self.grid_x()
        elif 'k' in [self.domain_xy, self.domain_x]:  # inverse space domain
            return self.grid_kx()
        else:
            raise AttributeError('Wrong domain_xy attribute')

    def scale_ky(self):  # scale in meters or meters**-1
        _logger.warning('"scale_ky" will be deprecated, use "grid_y and grid_ky" instead')
        if 's' in [self.domain_xy, self.domain_y]:    # space domain
            return self.grid_y()
        elif 'k' in [self.domain_xy, self.domain_y]:  # inverse space domain
            return self.grid_ky()
        else:
            raise AttributeError('Wrong domain_xy attribute')
            
    def scale_kz(self):  # scale in meters or meters**-1
        _logger.warning('"scale_kz" will be deprecated, use "grid_z and grid_kz" instead')        
        if self.domain_z == 't':  # time domain
            return self.grid_z()
        elif self.domain_z == 'f':  # frequency domain
            return self.grid_kz()
        else:
            raise AttributeError('Wrong domain_z attribute')

    def scale_x(self):  # scale in meters or radians
        _logger.warning('"scale_x" will be deprecated, use "grid_x and grid_kx" instead')        
        if self.domain_xy == 's':  # space domain
            return self.scale_kx()
        elif self.domain_xy == 'k':  # inverse space domain
            return self.scale_kx() * self.xlamds / 2 / np.pi
        else:
            raise AttributeError('Wrong domain_xy attribute')

    def scale_y(self):  # scale in meters or radians
        if self.domain_xy == 's':  # space domain
            return self.scale_ky()
        elif self.domain_xy == 'k':  # inverse space domain
            return self.scale_ky() * self.xlamds / 2 / np.pi
        else:
            raise AttributeError('Wrong domain_xy attribute')

    def scale_z(self):  # scale in meters
        if self.domain_z == 't':  # time domain
            return self.scale_kz()
        elif self.domain_z == 'f':  # frequency domain
            return 2 * pi / self.scale_kz()
        else:
            raise AttributeError('Wrong domain_z attribute')

    def ph_sp_dens(self):
        if self.domain_z == 't':
            dfl = deepcopy(self)
            dfl.fft_z()
        else:
            dfl = self
        pulse_energy = dfl.E()
        spec0 = dfl.int_z()
        freq_ev = h_eV_s * speed_of_light / dfl.scale_z()
        freq_ev_mean = np.sum(freq_ev * spec0) / np.sum(spec0)
        n_photons = pulse_energy / q_e / freq_ev_mean
        spec = calc_ph_sp_dens(spec0, freq_ev, n_photons)
        return freq_ev, spec

    def to_domain(self, domains='ts', **kwargs):
        """
        tranfers radiation to specified domains
        *domains is a string with one or two letters:
            ("t" or "f") and ("s" or "k")
        where
            't' (time); 'f' (frequency); 's' (space); 'k' (inverse space);
        e.g.
            't'; 'f'; 's'; 'k'; 'ts'; 'fs'; 'tk'; 'fk'
        order does not matter

        **kwargs are passed down to self.fft_z and self.fft_xy
        """
        _logger.debug('transforming radiation field to {} domain'.format(str(domains)))
        dfldomain_check(domains)

        for domain in domains:
            domain_o_z, domain_o_xy = self.domain_z, self.domain_xy
            if domain in ['t', 'f'] and domain is not domain_o_z:
                self.fft_z(**kwargs)
            if domain in ['s', 'k'] and domain is not domain_o_xy:
                self.fft_xy(**kwargs)

    def fft_z(self, method='mp', nthread=multiprocessing.cpu_count(),
              **kwargs):  # move to another domain ( time<->frequency )
        _logger.debug('calculating dfl fft_z from ' + self.domain_z + ' domain with ' + method)
        start = time.time()
        orig_domain = self.domain_z
        
        if nthread < 2:
            method = 'np'
        
        if orig_domain == 't':
            if method == 'mp' and fftw_avail:
                fft_exec = pyfftw.builders.fft(self.fld, axis=0, overwrite_input=True, planner_effort='FFTW_ESTIMATE',
                                               threads=nthread, auto_align_input=False, auto_contiguous=False,
                                               avoid_copy=True)
                self.fld = fft_exec()
            else:
                self.fld = np.fft.fft(self.fld, axis=0)
            # else:
            #     raise ValueError('fft method should be "np" or "mp"')
            self.fld = np.fft.ifftshift(self.fld, 0)
            self.fld /= np.sqrt(self.Nz())
            self.domain_z = 'f'
        elif orig_domain == 'f':
            self.fld = np.fft.fftshift(self.fld, 0)
            if method == 'mp' and fftw_avail:
                fft_exec = pyfftw.builders.ifft(self.fld, axis=0, overwrite_input=True, planner_effort='FFTW_ESTIMATE',
                                                threads=nthread, auto_align_input=False, auto_contiguous=False,
                                                avoid_copy=True)
                self.fld = fft_exec()
            else:
                self.fld = np.fft.ifft(self.fld, axis=0)
                
                # else:
                # raise ValueError("fft method should be 'np' or 'mp'")
            self.fld *= np.sqrt(self.Nz())
            self.domain_z = 't'
        else:
            raise ValueError("domain_z value should be 't' or 'f'")
        
        t_func = time.time() - start
        if t_func < 60:
            _logger.debug(ind_str + 'done in %.2f sec' % (t_func))
        else:
            _logger.debug(ind_str + 'done in %.2f min' % (t_func / 60))

    def fft_xy(self, method='mp', nthread=multiprocessing.cpu_count(),
               **kwargs):  # move to another domain ( spce<->inverse_space )
        _logger.debug('calculating fft_xy from ' + self.domain_xy + ' domain with ' + method)
        start = time.time()
        domain_orig = self.domain_xy

        if nthread < 2:
            method = 'np'
        
        if domain_orig == 's':
            if method == 'mp' and fftw_avail:
                fft_exec = pyfftw.builders.fft2(self.fld, axes=(1, 2), overwrite_input=False,
                                                planner_effort='FFTW_ESTIMATE', threads=nthread, auto_align_input=False,
                                                auto_contiguous=False, avoid_copy=True)
                self.fld = fft_exec()
            else:
                self.fld = np.fft.fft2(self.fld, axes=(1, 2))
                # else:
                # raise ValueError("fft method should be 'np' or 'mp'")
            self.fld = np.fft.fftshift(self.fld, axes=(1, 2))
            self.fld /= np.sqrt(self.Nx() * self.Ny())
            self.domain_xy = 'k'
        elif domain_orig == 'k':
            self.fld = np.fft.ifftshift(self.fld, axes=(1, 2))
            if method == 'mp' and fftw_avail:
                fft_exec = pyfftw.builders.ifft2(self.fld, axes=(1, 2), overwrite_input=False,
                                                 planner_effort='FFTW_ESTIMATE', threads=nthread,
                                                 auto_align_input=False, auto_contiguous=False, avoid_copy=True)
                self.fld = fft_exec()
            else:
                self.fld = np.fft.ifft2(self.fld, axes=(1, 2))
            # else:
            #     raise ValueError("fft method should be 'np' or 'mp'")
            self.fld *= np.sqrt(self.Nx() * self.Ny())
            self.domain_xy = 's'
        
        else:
            raise ValueError("domain_xy value should be 's' or 'k'")
        
        t_func = time.time() - start
        if t_func < 60:
            _logger.debug(ind_str + 'done in %.2f sec' % (t_func))
        else:
            _logger.debug(ind_str + 'done in %.2f min' % (t_func / 60))
    

    def mut_coh_func(self, norm=1, jit=1):
        '''
        calculates mutual coherence function
        consider downsampling the field first
        '''
        if jit:
            J = np.zeros([self.Ny(), self.Nx(), self.Ny(), self.Nx()]).astype(np.complex128)
            mut_coh_func(J, self.fld, norm=norm)
        else:
            I = self.int_xy() / self.Nz()
            J = np.mean(
                self.fld[:, :, :, np.newaxis, np.newaxis].conjugate() * self.fld[:, np.newaxis, np.newaxis, :, :],
                axis=0)
            if norm:
                J /= (I[:, :, np.newaxis, np.newaxis] * I[np.newaxis, np.newaxis, :, :])
        return J
    
    def coh(self, jit=0):
        '''
        calculates degree of transverse coherence
        consider downsampling the field first
        '''
        I = self.int_xy() / self.Nz()
        J = self.mut_coh_func(norm=0, jit=jit)
        coh = np.sum(abs(J) ** 2) / np.sum(I) ** 2
        return coh
        
    def tilt(self, angle=0, plane='x', return_orig_domains=True):
        '''
        deflects the radaition in given direction by given angle
        by introducing transverse phase chirp
        '''
        _logger.info('tilting radiation by {:.4e} rad in {} plane'.format(angle, plane))
        _logger.warn(ind_str + 'in beta')
        angle_warn = ind_str + 'deflection angle exceeds inverse space mesh range'
        
        k = 2 * pi / self.xlamds
        domains = self.domains()
        
        self.to_domain('s')
        if plane == 'y':
            if np.abs(angle) > self.xlamds / self.dy / 2:
                _logger.warning(angle_warn)
            dphi =  angle * k * self.scale_y()
            self.fld = self.fld * np.exp(1j * dphi)[np.newaxis, :, np.newaxis]
        elif plane == 'x':
            if np.abs(angle) > self.xlamds / self.dx / 2:
                _logger.warning(angle_warn)
            dphi =  angle * k * self.scale_x()
            self.fld = self.fld * np.exp(1j * dphi)[np.newaxis, np.newaxis, :]
        else:
            raise ValueError('plane should be "x" or "y"')
            
        if return_orig_domains:
            self.to_domain(domains)
    
            
    def disperse(self, disp=0, E_ph0=None, plane='x', return_orig_domains=True):
        '''
        introducing angular dispersion in given plane by deflecting the radaition by given angle depending on its frequency
        disp is the dispertion coefficient [rad/eV]
        E_ph0 is the photon energy in [eV] direction of which would not be changed (principal ray)
        '''
        _logger.info('introducing dispersion of {:.4e} [rad/eV] in {} plane'.format(disp, plane))
        _logger.warn(ind_str + 'in beta')
        angle_warn = ind_str + 'deflection angle exceeds inverse space mesh range'
        if E_ph0 == None:
            E_ph0 = 2 *np.pi / self.xlamds * speed_of_light * hr_eV_s
        
        dk = 2 * pi / self.Lz()
        k = 2 * pi / self.xlamds        
        phen = np.linspace(k - dk / 2 * self.Nz(), k + dk / 2 * self.Nz(), self.Nz()) * speed_of_light * hr_eV_s
        angle = disp * (phen - E_ph0)
        
        if np.amax([np.abs(np.min(angle)), np.abs(np.max(angle))]) > self.xlamds / self.dy / 2:
            _logger.warning(angle_warn)
        
        domains = self.domains()
        self.to_domain('sf')
        if plane =='y':
            dphi =  angle[:,np.newaxis] * k * self.scale_y()[np.newaxis, :]
            self.fld = self.fld * np.exp(1j *dphi)[:, :, np.newaxis]
        elif plane == 'x':
            dphi =  angle[:,np.newaxis] * k * self.scale_x()[np.newaxis, :]
            self.fld = self.fld * np.exp(1j *dphi)[:, np.newaxis, :]
        
        if return_orig_domains:
            self.to_domain(domains)
            
def imitate_sase_dfl(xlamds, rho=2e-4, seed=None, **kwargs):
    """
    imitation of SASE radiation in 3D

    xlamds - wavelength of the substracted fast-varying component
    rho - half of the expected FEL bandwidth
    **kwargs identical to generate_dfl()

    returns RadiationField object
    """

    _logger.info('imitating SASE radiation')
    if kwargs.get('wavelength', None) is not None:
        E0 = h_eV_s * speed_of_light / kwargs.pop('wavelength')
        _logger.debug(ind_str + 'using wavelength')
    else:
        E0 = h_eV_s * speed_of_light / xlamds
        _logger.debug(ind_str + 'using xlamds')
    dE = E0 * 2 * rho
    _logger.debug(ind_str + 'E0 = {}'.format(E0))
    _logger.debug(ind_str + 'dE = {}'.format(dE))
    dfl = generate_gaussian_dfl(xlamds, **kwargs)

    _logger.debug(ind_str + 'dfl.shape = {}'.format(dfl.shape()))
    td_scale = dfl.scale_z()
    _logger.debug(ind_str + 'time domain range = [{},  {}]m'.format(td_scale[0], td_scale[-1]))

    dk = 2 * np.pi / dfl.Lz()
    k = 2 * np.pi / dfl.xlamds
    fd_scale_ev = h_eV_s * speed_of_light * (
        np.linspace(k - dk / 2 * dfl.Nz(), k + dk / 2 * dfl.Nz(), dfl.Nz())) / 2 / np.pi
    fd_env = np.exp(-(fd_scale_ev - E0) ** 2 / 2 / (dE) ** 2)
    _logger.debug(ind_str + 'frequency domain range = [{},  {}]eV'.format(fd_scale_ev[0], fd_scale_ev[-1]))
    
    for key in imitate_1d_sase_like.__code__.co_varnames:
        kwargs.pop(key, None)
        
    _, td_envelope, _, _ = imitate_1d_sase_like(td_scale=td_scale, td_env=np.ones_like(td_scale), fd_scale=fd_scale_ev,
                                                fd_env=fd_env, td_phase=None, fd_phase=None, phen0=None, en_pulse=1,
                                                fit_scale='td', n_events=1, seed=seed, **kwargs)

    dfl.fld *= td_envelope[:, :, np.newaxis]

    return dfl

def calc_ph_sp_dens(spec, freq_ev, n_photons, spec_squared=1):
    """
    calculates number of photons per electronvolt
    """
    # _logger.debug('spec.shape = {}'.format(spec.shape))
    if spec.ndim == 1:
        axis = 0
    else:
        if spec.shape[0] == freq_ev.shape[0]:
            spec = spec.T
        axis = 1
        #     axis=0
        # elif spec.shape[1] == freq_ev.shape[0]:
        #     axis=1
        # else:
        #     raise ValueError('operands could not be broadcast together with shapes ', spec.shape, ' and ', freq_ev.shape)
    # _logger.debug('spec.shape = {}'.format(spec.shape))

    if spec_squared:
        spec_sum = np.trapz(spec, x=freq_ev, axis=axis)
    else:
        spec_sum = np.trapz(abs(spec) ** 2, x=freq_ev, axis=axis)

    if np.size(spec_sum) == 1:
        if spec_sum == 0:  # avoid division by zero
            spec_sum = np.inf
    else:
        spec_sum[spec_sum == 0] = np.inf  # avoid division by zero

    if spec_squared:
        norm_factor = n_photons / spec_sum
    else:
        norm_factor = np.sqrt(n_photons / spec_sum)

    if spec.ndim == 2:
        norm_factor = norm_factor[:, np.newaxis]
    # _logger.debug('spec.shape = {}'.format(spec.shape))
    # _logger.debug('norm_factor.shape = {}'.format(norm_factor.shape))
    spec = spec * norm_factor
    if axis == 1:
        spec = spec.T
    # _logger.debug('spec.shape = {}'.format(spec.shape))
    return spec

def imitate_1d_sase_like(td_scale, td_env, fd_scale, fd_env, td_phase=None, fd_phase=None, phen0=None, en_pulse=None,
                         fit_scale='td', n_events=1, **kwargs):
    """
    Models FEL pulse(s) based on Gaussian statistics
    td_scale - scale of the pulse on time domain [m]
    td_env - expected pulse envelope in time domain [W]
    fd_scale - scale of the pulse in frequency domain [eV]
    fd_env - expected pulse envelope in frequency domain [a.u.]
    td_phase - additional phase chirp to be added in time domain
    fd_phase - additional phase chirp to be added in frequency domain
    phen0 - sampling wavelength expressed in photon energy [eV]
    en_pulse - expected average energy of the pulses [J]
    fit_scale - defines the scale in which outputs should be returned:
        'td' - time domain scale td_scale is used for the outputs, frequency domain phase and envelope will be re-interpolated
        'fd' - frequency domain scale fd_scale is used for the outputs, time domain phase and envelope will be re-interpolated
    n_events - number of spectra to be generated

    returns tuple of 4 arguments: (ph_en, fd, s, td)
    fd_scale - colunm of photon energies in eV
    fd - matrix of radiation in frequency domain with shape, normalized such that np.sum(abs(fd)**2) is photon spectral density, i.e: np.sum(abs(fd)**2)*fd_scale = N_photons
    td - matrix of radiation in time domain, normalized such that abs(td)**2 = radiation_power in [w]
    """

    _logger.info('generating 1d radiation field imitating SASE')
    
    seed = kwargs.get('seed', None)
    if seed is not None:
        np.random.seed(seed)

    if fit_scale == 'td':

        n_points = len(td_scale)
        s = td_scale
        Ds = (td_scale[-1] - td_scale[0])
        ds = Ds / n_points

        td = np.random.randn(n_points, n_events) + 1j * np.random.randn(n_points, n_events)
        td *= np.sqrt(td_env[:, np.newaxis])
        fd = np.fft.ifftshift(np.fft.fft(np.fft.fftshift(td, axes=0), axis=0), axes=0)
        # fd = np.fft.ifft(td, axis=0)
        # fd = np.fft.fftshift(fd, axes=0)

        if phen0 is not None:
            e_0 = phen0
        else:
            e_0 = np.mean(fd_scale)

        # internal interpolated values
        fd_scale_i = h_eV_s * np.fft.fftfreq(n_points, d=(
                ds / speed_of_light)) + e_0  # internal freq.domain scale based on td_scale
        fd_scale_i = np.fft.fftshift(fd_scale_i, axes=0)
        fd_env_i = np.interp(fd_scale_i, fd_scale, fd_env, right=0, left=0)

        if fd_phase is None:
            fd_phase_i = np.zeros_like(fd_env_i)
        else:
            fd_phase_i = np.interp(fd_scale_i, fd_scale, fd_phase, right=0, left=0)

        fd *= np.sqrt(fd_env_i[:, np.newaxis]) * np.exp(1j * fd_phase_i[:, np.newaxis])

        # td = np.fft.ifftshift(fd, axes=0)
        # td = np.fft.fft(td, axis=0)
        td = np.fft.ifftshift(np.fft.ifft(np.fft.fftshift(fd, axes=0), axis=0), axes=0)

        td_scale_i = td_scale

    elif fit_scale == 'fd':

        n_points = len(fd_scale)
        Df = abs(fd_scale[-1] - fd_scale[0]) / h_eV_s
        df = Df / n_points

        fd = np.random.randn(n_points, n_events) + 1j * np.random.randn(n_points, n_events)
        fd *= np.sqrt(fd_env[:, np.newaxis])
        td = np.fft.ifftshift(np.fft.ifft(np.fft.fftshift(fd, axes=0), axis=0), axes=0)

        td_scale_i = np.fft.fftfreq(n_points, d=df) * speed_of_light
        td_scale_i = np.fft.fftshift(td_scale_i, axes=0)
        td_scale_i -= np.amin(td_scale_i)
        td_env_i = np.interp(td_scale_i, td_scale, td_env, right=0, left=0)

        if td_phase is None:
            td_phase_i = np.zeros_like(td_env_i)
        else:
            td_phase_i = np.interp(td_scale_i, td_scale, td_phase, right=0, left=0)

        td *= np.sqrt(td_env_i[:, np.newaxis]) * np.exp(1j * td_phase_i[:, np.newaxis])

        fd = np.fft.ifftshift(np.fft.fft(np.fft.fftshift(td, axes=0), axis=0), axes=0)

        fd_scale_i = fd_scale

    else:
        raise ValueError('fit_scale should be either "td" of "fd"')

    # normalization for pulse energy
    if en_pulse == None:
        _logger.debug(ind_str + 'no en_pulse provided, calculating from integral of td_env')
        en_pulse = np.trapz(td_env, td_scale / speed_of_light)

    pulse_energies = np.trapz(abs(td) ** 2, td_scale_i / speed_of_light, axis=0)
    scale_coeff = en_pulse / np.mean(pulse_energies)
    td *= np.sqrt(scale_coeff)

    # normalization for photon spectral density
    spec = np.mean(np.abs(fd) ** 2, axis=1)
    spec_center = np.sum(spec * fd_scale_i) / np.sum(spec)

    n_photons = pulse_energies * scale_coeff / q_e / spec_center
    fd = calc_ph_sp_dens(fd, fd_scale_i, n_photons, spec_squared=0)
    td_scale, fd_scale = td_scale_i, fd_scale_i

    np.random.seed()

    return (td_scale, td, fd_scale, fd)


def imitate_1d_sase(spec_center=500, spec_res=0.01, spec_width=2.5, spec_range=(None, None), pulse_length=6,
                    en_pulse=1e-3, flattop=0, n_events=1, spec_extend=5, **kwargs):
    """
    Models FEL pulse(s) based on Gaussian statistics
    spec_center - central photon energy in eV
    spec_res - spectral resolution in eV
    spec_width - width of spectrum in eV (fwhm of E**2)
    spec_range = (E1, E2) - energy range of the spectrum. If not defined, spec_range = (spec_center - spec_width * spec_extend, spec_center + spec_width * spec_extend)
    pulse_length - longitudinal size of the pulse in um (fwhm of E**2)
    en_pulse - expected average energy of the pulses in Joules
    flattop - if true, flat-top pulse in time domain is generated with length 'pulse_length' in um
    n_events - number of spectra to be generated

    return tuple of 4 arguments: (s, td, ph_en, fd)
    ph_en - colunm of photon energies in eV with size (spec_range[2]-spec_range[1])/spec_res
    fd - matrix of radiation in frequency domain with shape ((spec_range[2]-spec_range[1])/spec_res, n_events), normalized such that np.sum(abs(fd)**2) is photon spectral density, i.e: np.sum(abs(fd)**2)*spec_res = N_photons
    s - colunm of longitudinal positions along the pulse in yime domain in um
    td - matrix of radiation in time domain with shape ((spec_range[2]-spec_range[1])/spec_res, n_events), normalized such that abs(td)**2 = radiation_power
    """

    if spec_range == (None, None):
        spec_range = (spec_center - spec_width * spec_extend, spec_center + spec_width * spec_extend)
    elif spec_center == None:
        spec_center = (spec_range[1] + spec_range[0]) / 2

    pulse_length_sigm = pulse_length / (2 * np.sqrt(2 * np.log(2)))
    spec_width_sigm = spec_width / (2 * np.sqrt(2 * np.log(2)))

    fd_scale = np.arange(spec_range[0], spec_range[1], spec_res)
    n_points = len(fd_scale)
    _logger.debug(ind_str + 'N_points * N_events = %i * %i' % (n_points, n_events))

    fd_env = np.exp(-(fd_scale - spec_center) ** 2 / 2 / spec_width_sigm ** 2)
    td_scale = np.linspace(0, 2 * np.pi / (fd_scale[1] - fd_scale[0]) * hr_eV_s * speed_of_light, n_points)

    if flattop:
        td_env = np.zeros_like(td_scale)
        il = find_nearest_idx(td_scale, np.mean(td_scale) - pulse_length * 1e-6 / 2)
        ir = find_nearest_idx(td_scale, np.mean(td_scale) + pulse_length * 1e-6 / 2)
        td_env[il:ir] = 1
    else:
        s0 = np.mean(td_scale)
        td_env = np.exp(-(td_scale - s0) ** 2 / 2 / (pulse_length_sigm * 1e-6) ** 2)

    result = imitate_1d_sase_like(td_scale, td_env, fd_scale, fd_env, phen0=spec_center, en_pulse=en_pulse,
                                  fit_scale='fd', n_events=n_events, **kwargs)

    return result

            
def dfldomain_check(domains, both_req=False):
    err = ValueError(
        'domains should be a string with one or two letters from ("t" or "f") and ("s" or "k"), not {}'.format(
            str(domains)))

    # if type(domains) is not str:
    #     raise err
    if len(domains) < 1 or len(domains) > 2:
        raise err
    if len(domains) < 2 and both_req == True:
        raise ValueError('please provide both domains, e.g. "ts" "fs" "tk" "fk"')

    domains_avail = ['t', 'f', 's', 'k']
    for letter in domains:
        if letter not in domains_avail:
            raise err

    if len(domains) == 2:
        D = [['t', 'f'], ['s', 'k']]
        for d in D:
            if domains[0] in d and domains[1] in d:
                raise err

        """
        tranfers radiation to specified domains
        *domains is a string with one or two letters: 
            ("t" or "f") and ("s" or "k")
        where 
            't' (time); 'f' (frequency); 's' (space); 'k' (inverse space); 
        e.g.
            't'; 'f'; 's'; 'k'; 'ts'; 'fs'; 'tk'; 'fk'
        order does not matter
        
        **kwargs are passed down to self.fft_z and self.fft_xy
        """            
#%%    
### script itself ###

#optics elements check
#dfl = RadiationField()
#E_pohoton = 1239.8#200 #central photon energy [eV]
#kwargs={'xlamds':(h_eV_s * speed_of_light / E_pohoton), #[m] - central wavelength
#        'rho':1.0e-4, 
#        'shape':(101,101,11),             #(x,y,z) shape of field matrix (reversed) to dfl.fld
#        'dgrid':(400e-5,400e-5,35e-6), #(x,y,z) [m] - size of field matrix
#        'power_rms':(25e-5,25e-5,4e-6),#(x,y,z) [m] - rms size of the radiation distribution (gaussian)
#        'power_center':(0,0,None),     #(x,y,z) [m] - position of the radiation distribution
#        'power_angle':(0,0),           #(x,y) [rad] - angle of further radiation propagation
#        'power_waistpos':(0,0),        #(Z_x,Z_y) [m] downstrean location of the waist of the beam
#        'wavelength':None,             #central frequency of the radiation, if different from xlamds
#        'zsep':None,                   #distance between slices in z as zsep*xlamds
#        'freq_chirp':0,                #dw/dt=[1/fs**2] - requency chirp of the beam around power_center[2]
#        'en_pulse':None,                #total energy or max power of the pulse, use only one
#        'power':1e6,
#        }

#dfl = generate_gaussian_dfl(**kwargs);  #Gaussian beam defenition











