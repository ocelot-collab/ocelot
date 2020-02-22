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
    def __init__(self, shape):
        self.dx = []
        self.dy = []
        self.dz = []
        self.shape = shape
        
        self.xlamds = None
        self.used_aprox = 'SVEA'
        
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
    
    def grid_f(self):
        df = 2 * np.pi * speed_of_light / self.Lz() # self.Lz() must be replaced with self.Tz()
        
        if self.xlamds is None:        
            return np.linspace(-df / 2 * self.Nz(), df / 2 * self.Nz(), self.Nz())
        
        elif self.used_aprox == 'SVEA' and self.xlamds is not None:
            f0 = 2 * np.pi * speed_of_light / self.xlamds
            return np.linspace(f0 - df / 2 * self.Nz(), f0 + df / 2 * self.Nz(), self.Nz())
        
        else:
            raise ValueError

    def grid_z(self):
        return self.scale_t() * speed_of_light

    def grid_kz(self):
        return self.grid_f() / speed_of_light
            
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
            
            
#%%    
### script itself ###

#optics elements check
dfl = RadiationField()
E_pohoton = 1239.8#200 #central photon energy [eV]
kwargs={'xlamds':(h_eV_s * speed_of_light / E_pohoton), #[m] - central wavelength
        'rho':1.0e-4, 
        'shape':(101,101,11),             #(x,y,z) shape of field matrix (reversed) to dfl.fld
        'dgrid':(400e-5,400e-5,35e-6), #(x,y,z) [m] - size of field matrix
        'power_rms':(25e-5,25e-5,4e-6),#(x,y,z) [m] - rms size of the radiation distribution (gaussian)
        'power_center':(0,0,None),     #(x,y,z) [m] - position of the radiation distribution
        'power_angle':(0,0),           #(x,y) [rad] - angle of further radiation propagation
        'power_waistpos':(0,0),        #(Z_x,Z_y) [m] downstrean location of the waist of the beam
        'wavelength':None,             #central frequency of the radiation, if different from xlamds
        'zsep':None,                   #distance between slices in z as zsep*xlamds
        'freq_chirp':0,                #dw/dt=[1/fs**2] - requency chirp of the beam around power_center[2]
        'en_pulse':None,                #total energy or max power of the pulse, use only one
        'power':1e6,
        }

dfl = generate_gaussian_dfl(**kwargs);  #Gaussian beam defenition











