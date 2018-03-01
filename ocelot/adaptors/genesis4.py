import time
import os
import socket
import copy
import h5py
import numpy as np
from ocelot.optics.wave import calc_ph_sp_dens, RadiationField
from ocelot.common.globals import *
from ocelot.adaptors.genesis import GenesisElectronDist #tmp

class Genesis4Input:
    '''
    Genesis input files storage object
    '''
    def __init__(self):
        pass
    

class Genesis4ParticlesDump:
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
    
    with h5py.File(file_path, 'r') as h5:
        
        nslice = h5.get('slicecount')[0]
        lambdaref = h5.get('wavelength')[0]
        sepslice = h5.get('slicespacing')[0]
        gridsize = h5.get('gridsize')[0]
        N = int(np.sqrt(h5.get('slice000001/field-real').size))
        
        field_real = []
        field_imag = []
        for dset in h5:
            if dset.startswith('slice0'):
                field_real.append(h5[dset]['field-real'][:].reshape(N,N))
                field_imag.append(h5[dset]['field-imag'][:].reshape(N,N))
        
        dfl = RadiationField()
        dfl.fld = np.array(field_real) + 1j * np.array(field_imag)
        dfl.dx = gridsize
        dfl.dy = gridsize
        dfl.dz = sepslice
        dfl.xlamds = lambdaref
        dfl.domain_z = 't'  # longitudinal domain (t - time, f - frequency)
        dfl.domain_xy = 's'  # transverse domain (s - space, k - inverse space)
#        dfl.h5 = h5
        dfl.filePath = h5.filename
    
    return dfl




def read_dpa4(file_path):
    h5 = h5py.File(file_path, 'r')
    
    nslice = int(h5.get('slicecount')[0])
    lslice = h5.get('slicelength')[0]
    sepslice = h5.get('slicespacing')[0]
    npart = int(h5.get('slice000001/gamma').size)
    nbins = int(h5.get('beamletsize')[0])
    zsep = int(sepslice / lslice)
    
    fill_gaps=0
    
    x = []
    y = []
    px=[]
    py=[]
    g=[]
    ph = []
    s0 = []
    s = []
    I = []
    
    for dset in h5:
            if dset.startswith('slice0'):
                I.append(h5[dset]['current'][:])
    
    #if fill_gaps:
    #    for dset in h5:
    #        if dset.startswith('slice0'):
    #            ph0 = h5[dset]['theta'][:]
    #            ph.append(ph0.repeat(zsep))
    #            x.append(h5[dset]['x'][:].repeat(zsep))
    #            px.append(h5[dset]['px'][:].repeat(zsep))
    #            y.append(h5[dset]['y'][:].repeat(zsep))
    #            py.append(h5[dset]['py'][:].repeat(zsep))
    #            g.append(h5[dset]['gamma'][:].repeat(zsep))
    #            for sl in range(zsep):
    #                s.append(s0 + ph0 / 2 / np.pi * lslice)            
    #                s0 =+ lslice
    else:
        for dset in h5:
            if dset.startswith('slice0'):
                ph0 = h5[dset]['theta'][:]
    #            s.append(s0 + ph0 / 2 / np.pi * lslice)
                x.append(h5[dset]['x'][:])
                px.append(h5[dset]['px'][:])
                y.append(h5[dset]['y'][:])
                py.append(h5[dset]['py'][:])
                g.append(h5[dset]['gamma'][:])
                ph.append(ph0)
                s0 += sepslice
            
    
    npartpb=int(npart/nbins)
    
    dpa = Genesis4ParticlesDump()
    
    dpa.x = np.reshape(x, (nslice, nbins, npartpb), order='F')
    dpa.px = np.reshape(px, (nslice, nbins, npartpb), order='F')
    dpa.y = np.reshape(y, (nslice, nbins, npartpb), order='F')
    dpa.py = np.reshape(py, (nslice, nbins, npartpb), order='F')
    dpa.ph = np.reshape(ph, (nslice, nbins, npartpb), order='F')
    dpa.g = np.reshape(g, (nslice, nbins, npartpb), order='F')
    dpa.I = np.array(I).flatten()
    #dpa.s = np.array(s).flatten()
    dpa.nslice = nslice
    dpa.lslice = lslice
    dpa.npart = npart
    dpa.nbins = nbins
    dpa.zsep = zsep
    dpa.filePath = h5.filename
    
    return dpa





def dpa42edist(dpa, n_part=None, fill_gaps=1, debug=1):
    fill_gaps = 0
    '''
    Convert dpa to edist objects
    reads GenesisParticlesDump() object
    returns GenesisElectronDist() object
    num_part - desired approximate number of particles in edist
    smear - whether to shuffle macroparticles smearing microbunching
    '''
    import random
    start_time = time.time()
    #if debug > 0:
    #    print ('    transforming particle to distribution file')
    
    #assert out('itdp') == True, '! steadystate Genesis simulation, dpa2dist() not implemented yet!'
    
    npart = dpa.npart
    # nslice=int(out('nslice'))
    nslice = dpa.nslice
    nbins = dpa.nbins
    xlamds = dpa.lslice
    zsep = dpa.zsep
    gen_I = dpa.I
    
#    npart_bin = int(npart / nbins)
        
    
    if fill_gaps:
        s0 = np.linspace(0, nslice * zsep * xlamds, nslice)
        s = np.linspace(0, nslice * zsep * xlamds, nslice * zsep)
        I = np.interp(s, s0, dpa.I)
        dt = (s[1] - s[0]) / speed_of_light
    else:
        s = np.linspace(0, nslice * zsep * xlamds, nslice)
        I = dpa.I
        dt = (s[1] - s[0]) / speed_of_light
        
    C = np.sum(I) * dt
        
    n_part_max = np.sum(I / I.max() * dpa.npart)
#    print(n_part)
    
    if n_part is not None:
        if n_part > n_part_max:
            n_part = n_part_max
    else:
        n_part = int(np.floor(n_part_max))
#    print(n_part)
    
    
    #n_part_bin = (I / np.sum(I) * n_part / nbins).astype(int)
    #n_part_bin[n_part_bin > npart_bin] = npart_bin
    #print(n_part_bin.max())
    n_part_slice = (I / np.sum(I) * n_part).astype(int)
    n_part_slice[n_part_slice > npart] = npart
#    print(n_part_slice.max())
    
    
    #pick_i = random.sample(range(n_part), n_part_slice[i])

    g = np.reshape(dpa.g, (nslice, npart))
    x = np.reshape(dpa.x, (nslice, npart))
    y = np.reshape(dpa.y, (nslice, npart))
    px = np.reshape(dpa.px, (nslice, npart)) / g
    py = np.reshape(dpa.py, (nslice, npart)) / g
    ph = np.reshape(dpa.ph, (nslice, npart))
    t1 = ph * xlamds / speed_of_light
    t0 = np.arange(nslice)[:,np.newaxis] * xlamds * zsep / speed_of_light
    t = t1+t0
    
    
    
    edist = GenesisElectronDist()
    #g1 = np.array([])
    #x1 = np.array([])
    #y1 = np.array([])
    
    
    for i in np.arange(I.size):
    #    for ii in np.arange(nbins):
    #    pick_i = random.sample(range(npart_bin), n_part_bin[i])
        pick_i = random.sample(range(npart), n_part_slice[i])
    
    #    t.append(dpa.t, t[i, pick_i])
    #    g = np.append(g, dpa.e[i, ii, pick_i])
        
        edist.g = np.append(edist.g, g[i, pick_i])
        edist.xp = np.append(edist.xp, px[i, pick_i])
        edist.yp = np.append(edist.yp, py[i, pick_i])
        edist.x = np.append(edist.x, x[i, pick_i])  
        edist.y = np.append(edist.y, y[i, pick_i])
        edist.t = np.append(edist.t, t[i, pick_i])
        
    #    if fill_gaps and zsep>1:
    #        for ii in range(zsep):
                
    
    edist.part_charge = C / n_part
    
    return edist