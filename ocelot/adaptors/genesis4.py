import time
import os
import socket
import copy
import h5py
import numpy as np
from ocelot.optics.wave import calc_ph_sp_dens, RadiationField
from ocelot.common.globals import *

class Genesis4Input:
    '''
    Genesis input files storage object
    '''
    def __init__(self):
        pass
    

class GenesisParticlesDump:
    '''
    Genesis particle *.dpa files storage object
    Each particle record in z starts with the energy of all particles 
    followed by the output of the particle phases, 
    positions in x and y and the momenta in x and y. 
    The momenta are normalized to mc
    '''

    def __init__(self):
        self.e = []
        self.ph = []
        self.x = []
        self.y = []
        self.px = []
        self.py = []

        # self.fileName = ''
        self.filePath = ''

    def fileName(self):
        return filename_from_path(self.filePath)


class Genesis4Output:
    '''
    Genesis input files storage object
    '''
    def __init__(self):
        self.h5 = None #hdf5 pointer
        
    def close(self):
        self.h5.close()
       
    @property
    def filePath(self):
        return self.h5.filename
    
    def fileName(self):
        return os.path.basename(self.h5.filename)

    @property
    def nZ(self):
        return self.z.size
    
    @property
    def nSlices(self):
        return self.h5['Beam/current'].size
    
    @property
    def lambdaref(self):
        return self.h5['Global/lambdaref'][0]
    
    @property
    def phenref(self):
        return h_eV_s * speed_of_light / self.lambdaref
    
    @property
    def I(self):
        return self.h5['Beam/current'][0]
    
    @property
    def beam_charge(self):
        return np.trapz(self.I, self.t)
    
    @property
    def rad_power(self):
        return self.h5['Field/power']
    
    @property
    def rad_energy(self):
        return np.trapz(self.rad_power, self.t)
    
    @property
    def n_photons(self):
        return self.rad_energy / q_e / self.phenref
        
    @property
    def t(self):
        if not self.tdp:
            return None
        else:
            return self.s / speed_of_light
        
    def rad_field(self, zi=None, loc='near'):
        if loc == 'far':
            intens = self.h5['Field/intensity-farfield']
            phase = self.h5['Field/phase-farfield']
        elif loc == 'near':
            intens = self.h5['Field/intensity-nearfield']
            phase = self.h5['Field/phase-nearfield']
        else:
            raise ValueError('loc should be either "far" or "near"')
            
        if zi is not None:
            intens = intens[zi,:]
            phase = phase[zi,:]
            
        #not scaled properly!!!!
        field = np.sqrt(intens[:]) * np.exp(1j * phase[:])
        return field
    
    def calc_spec(self, zi=None, loc='near', npad=1, estimate_ph_sp_dens=1):
        
        field = self.rad_field(zi=zi, loc=loc)
        axis = field.ndim - 1
        
        spec = np.abs(np.fft.fft(field, axis=axis))**2
        spec = np.fft.fftshift(spec, axes=axis)
        
        scale_ev = h_eV_s * speed_of_light * (np.fft.fftfreq(self.nSlices, d=self.s[1]-self.s[0]) + 1 / self.lambdaref)
        scale_ev = np.fft.fftshift(scale_ev)
        
        if estimate_ph_sp_dens:
            tt=np.trapz(spec, scale_ev, axis=axis)
            if axis==1:
                tt[tt==0] = np.inf
                spec *= (self.n_photons / tt)[:, np.newaxis]
            else:
                if tt==0:
                    tt = np.inf
                spec *= (self.n_photons[zi] / tt)
        
        return scale_ev, spec

        

def get_genesis4_launcher(launcher_program='genesis4', launcher_argument=''):
    '''
    Returns MpiLauncher() object for given program
    '''
    host = socket.gethostname()
    
    launcher = MpiLauncher()
    launcher.program = launcher_program
    launcher.argument = launcher_argument
    # launcher.program = '/data/netapp/xfel/products/genesis/genesis'
    # launcher.argument = ' < tmp.cmd | tee log'
    
    return launcher

def read_gout4(file_path):
    out = Genesis4Output()
    out.h5 = h5py.File(file_path, 'r')
    
    out.z = out.h5['Lattice/zplot'][:]
    out.zlat = out.h5['Lattice/z'][:]
    
    
    if 'time' in out.h5['Global'] and out.h5['Global/time'][0] == 1:
        out.tdp = True
    else:
        out.tdp = False
    
    if out.tdp:
        if 's0' in out.h5['Global']:
            s0 = out.h5['Global/s0']
        else:
            s0 = 0
        
        sn = out.h5['Beam/current'].size
        out.s = np.linspace(s0, out.h5['Global/slen'][()]+s0, sn)
    
    return out

def read_dfl4(file_path):
    
    fld = h5py.File(file_path, 'r')
    
    nslice = fld.get('slicecount')[0]
    lambdaref = fld.get('wavelength')[0]
    sepslice = fld.get('slicespacing')[0]
    gridsize = fld.get('gridsize')[0]
    N = int(np.sqrt(fld.get('slice000001/field-real').size))
    
    field_real = []
    field_imag = []
    for dset in fld:
        if dset.startswith('slice0'):
            field_real.append(fld[dset]['field-real'][:].reshape(N,N))
            field_imag.append(fld[dset]['field-imag'][:].reshape(N,N))
    
    dfl = RadiationField()
    dfl.fld = np.array(field_real) + 1j * np.array(field_imag)
    dfl.dx = gridsize
    dfl.dy = gridsize
    dfl.dz = sepslice
    dfl.xlamds = lambdaref
    dfl.domain_z = 't'  # longitudinal domain (t - time, f - frequency)
    dfl.domain_xy = 's'  # transverse domain (s - space, k - inverse space)
    dfl.filePath = fld.filename
    
    return dfl




#def read_dpa4(file_path):
#par = h5py.File(file_path, 'r')
#
#nslice = par.get('slicecount')[0]
#lslice = par.get('slicelength')[0]
#sepslice = par.get('slicespacing')[0]
#npart = par.get('slice000001/gamma').size
#zsep = int(sepslice / lslice)
#
##fill_gaps=1
#
#x = []
#y = []
#px=[]
#py=[]
#g=[]
#ph = []
#s0 = 0
#s = []
#I = out.I
#
#for dset in par:
#        if dset.startswith('slice0'):
#            I.append(par[dset]['current'][:])
#
##if fill_gaps:
##    for dset in par:
##        if dset.startswith('slice0'):
##            ts = par[dset]['theta'][:]
##            t.append(ts.repeat(zsep))
##            x.append(par[dset]['x'][:].repeat(zsep))
##            px.append(par[dset]['px'][:].repeat(zsep))
##            y.append(par[dset]['y'][:].repeat(zsep))
##            py.append(par[dset]['py'][:].repeat(zsep))
##            g.append(par[dset]['gamma'][:].repeat(zsep))
##            for sl in range(zsep):
##                s.append(s0 + ts / 2 / np.pi * lslice)            
##                s0 =+ lslice
##else:
#for dset in par:
#    if dset.startswith('slice0'):
#        ph0 = par[dset]['theta'][:]
#        s.append(s0 + ts / 2 / np.pi * lslice)
#        x.append(par[dset]['x'][:])
#        px.append(par[dset]['px'][:])
#        y.append(par[dset]['y'][:])
#        py.append(par[dset]['py'][:])
#        g.append(par[dset]['gamma'][:])
#        ph.append(ph0)
#        s0 += sepslice
#
#nbins=4
#npartpb=int(npart/nbins)
#dpa.x = np.reshape(x, (nslice, nbins, npartpb), order='F')
#dpa.px = np.reshape(px, (nslice, nbins, npartpb), order='F')
#dpa.y = np.reshape(y, (nslice, nbins, npartpb), order='F')
#dpa.py = np.reshape(py, (nslice, nbins, npartpb), order='F')
#dpa.ph = np.reshape(py, (nslice, nbins, npartpb), order='F')
#dpa.e = np.reshape(g, (nslice, nbins, npartpb), order='F')
#
#dpa.filePath = par.filename