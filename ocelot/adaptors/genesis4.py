import time
import os
import socket
import copy
    
try:
    import h5py
    h5py_avail = True
except ImportError:
    print("wave.py: module h5py is not installed. Install it if you want to use genesis4 adaptor")
    h5py_avail = False

import numpy as np
from ocelot import ParticleArray
from ocelot.optics.wave import calc_ph_sp_dens, RadiationField
from ocelot.common.globals import *
from ocelot.adaptors.genesis import GenesisElectronDist #tmp
from ocelot.common.logging import *
from ocelot.utils.launcher import *

_logger = logging.getLogger('ocelot.gen4') 


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

def read_gout4(filePath):
    
    _logger.info('reading gen4 .out file')
    _logger.debug(ind_str + 'reading from ' + filePath)
    
    out = Genesis4Output()
    try:
        out.h5 = h5py.File(filePath, 'r')
    except Exception:
        _logger.error(ind_str + 'no such file ' + filePath)
        raise
    
    out.z = out.h5['Lattice/zplot'][:]
    out.zlat = out.h5['Lattice/z'][:]
    
    
    if 'time' in out.h5['Global'] and out.h5['Global/time'][0] == 1:
        out.tdp = True
        _logger.debug(ind_str + 'tdp=True')
    else:
        out.tdp = False
        _logger.debug(ind_str + 'tdp=False')
    
    if out.tdp:
        
        if 's0' in out.h5['Global']:
            s0 = out.h5['Global/s0']
        else:
            s0 = 0
        
        sn = out.h5['Beam/current'].size
        out.s = np.linspace(s0, out.h5['Global/slen'][()]+s0, sn)
    
    _logger.debug(ind_str + 'done')
    
    return out

def read_dfl4(filePath):
    
    _logger.info('reading gen4 .dfl file')
    _logger.debug(ind_str + 'reading from ' + filePath)
    
    with h5py.File(filePath, 'r') as h5:
        
        nslice = h5.get('slicecount')[0]
        lambdaref = h5.get('wavelength')[0]
        sepslice = h5.get('slicespacing')[0]
        gridsize = h5.get('gridsize')[0]
        N = int(np.sqrt(h5.get('slice000001/field-real').size))
        _logger.warn(ind_str + 'tbd')
        
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
    
    _logger.debug(ind_str + 'done')
    
    return dfl




def read_dpa4(filePath):

    _logger.info('reading gen4 .dpa file')
    _logger.debug(ind_str + 'reading from ' + filePath)
    
    h5 = h5py.File(filePath, 'r')
    
    nslice = int(h5.get('slicecount')[0])
    lslice = h5.get('slicelength')[0]
    sepslice = h5.get('slicespacing')[0]
    npart = int(h5.get('slice000001/gamma').size)
    nbins = int(h5.get('beamletsize')[0])
    zsep = int(sepslice / lslice)
    _logger.warn(ind_str + 'tbd')
    
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
    
    h5.close()
    
    _logger.debug(ind_str + 'done')
    
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
    
    _logger.info('converting gen4 .dpa to edist')
    
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

def read_dpa42parray(filePath, N_part=None, fill_gaps=True):
    
    _logger.info('reading gen4 .dpa file into parray')
    _logger.debug(ind_str + 'reading from ' + filePath)
    
    import random
    N_part = None
    
    #N_part = 100000
    fill_gaps=True
    _logger.debug('fill_gaps = ' + str(fill_gaps))
    
    h5 = h5py.File(filePath, 'r')
    
    nslice = int(h5.get('slicecount')[0])
    lslice = h5.get('slicelength')[0]
    sepslice = h5.get('slicespacing')[0]
    npart = int(h5.get('slice000001/gamma').size)
    nbins = int(h5.get('beamletsize')[0])
    zsep = int(sepslice / lslice)
    
    _logger.debug('nslice = ' + str(nslice))
    _logger.debug('lslice = ' + str(lslice) + 'm')
    _logger.debug('sepslice = ' + str(sepslice)+'m')
    _logger.debug('zsep = ' + str(zsep))
    _logger.debug('npart = ' + str(npart))
    _logger.debug('nbins = ' + str(nbins))
    
    I = []
    for dset in h5:
            if dset.startswith('slice0'):
                I.append(h5[dset]['current'][0])
    I = np.array(I)
    
    dt = zsep * lslice / speed_of_light
    _logger.debug('dt = ' + str(dt) + 'sec')
        
    N_part_max = np.sum(I / I.max() * npart) # total maximum reasonable number of macroparticles of the same charge that can be extracted
    
    if N_part is not None:
        if N_part > N_part_max:
            N_part = int(np.floor(N_part_max))
    else:
        N_part = int(np.floor(N_part_max))
    _logger.debug('Number of particles max= ' + str(N_part))
    
    n_part_slice = (I / np.sum(I) * N_part).astype(int) #array of number of particles per new bin
    n_part_slice[n_part_slice > npart] = npart
    
    N_part_act = np.sum(n_part_slice) #actual number of particles
    _logger.debug('Number of particles actual= ' + str(N_part_act))
    
    dt = zsep * lslice / speed_of_light
    C = np.sum(I) * dt #total charge
    c = C / N_part_act # particle charge
    
    pick_i = [random.sample(range(npart), n_part_slice[i]) for i in range(nslice)]
    
    x = []
    y = []
    px = []
    py = []
    g = []
    ph = []
    s = []
    
    for dset in h5:
        if dset.startswith('slice0'):
            i = int(dset.strip('/slice'))-1
            if len(pick_i[i]) > 0:
                ph0 = h5[dset]['theta'].value[pick_i[i]]
                s.append(i*zsep*lslice + ph0 / 2 / np.pi * lslice)
                x.append(np.array(h5[dset]['x'].value[pick_i[i]]))
                px.append(np.array(h5[dset]['px'].value[pick_i[i]]))
                y.append(np.array(h5[dset]['y'].value[pick_i[i]]))
                py.append(np.array(h5[dset]['py'].value[pick_i[i]]))
                g.append(np.array(h5[dset]['gamma'].value[pick_i[i]]))
    
    p_array = ParticleArray()
    p_array.rparticles = np.empty((6, N_part_act))
    
    g = np.concatenate(g).ravel()
    g0 = np.mean(g) # average gamma
    p_array.E = g0 * m_e_GeV # average energy in GeV
    p0 = sqrt(g0**2-1) * m_e_eV / speed_of_light
    
    p_array.rparticles[0] = np.concatenate(x).ravel() # position in x in meters
    p_array.rparticles[1] = np.concatenate(px).ravel() / g0  # divergence in x
    p_array.rparticles[2] = np.concatenate(y).ravel() # position in x in meters
    p_array.rparticles[3] = np.concatenate(py).ravel() / g0  # divergence in x
    p_array.rparticles[4] = -np.concatenate(s).ravel()
    p_array.rparticles[5] = (g - g0) * m_e_eV / p0 / speed_of_light
    
    if fill_gaps:
        p_array.rparticles[4] -= np.random.randint(0, zsep, N_part_act) * lslice
    
    p_array.q_array = np.ones(N_part_act) * c
    
    h5.close()
    return p_array


def write_gen4_lat(lat, filePath, line_name='LINE', l=np.inf):
    from ocelot.cpbd.elements import Undulator, Drift, Quadrupole, UnknownElement
    _logger.info('writing genesis4 lattice')
    _logger.debug(ind_str + 'writing to ' + filePath)
    
    lat_str = []
    beamline = []
    ll=0
    
    lat_str.append('# generated with Ocelot\n')
    
    for element in lat.sequence:
        
        if ll >= l:
            break
        
        element_num = str(len(beamline) + 1).zfill(3)
        
        if hasattr(element,'l'):
            ll += element.l
        
        if isinstance(element, Undulator):
            element_name = element_num + 'UND'
            s = '{:}: UNDULATOR = {{lambdau = {:}, nwig = {:}, aw = {:.6f}}};'.format(element_name, element.lperiod, element.nperiods, element.Kx/sqrt(2))
#            print(s)
    
        elif isinstance(element, Drift):
            element_name = element_num + 'DR'
            s = '{:}: DRIFT = {{l={:}}};'.format(element_name, element.l)
    
        elif isinstance(element, Quadrupole):
            if element.k1>=0:
                element_name = element_num + 'QF'
            else:
                element_name =  element_num + 'QD'
            s = '{:}: QUADRUPOLE = {{l = {:}, k1 = {:.6f} }};'.format(element_name, element.l, element.k1)
    
        else:
            _logger.debug('Unknown element with length '+ str(element.l))
            continue
        
        beamline.append(element_name)
        lat_str.append(s)
        
    lat_str.append('')
    lat_str.append('{:}: LINE = {{{:}}};'.format(line_name, ','.join(beamline)))
    lat_str.append('\n# end of file\n')
                   
    with open(filePath, 'w') as f:
        f.write("\n".join(lat_str))

def write_edist_hdf5(edist, filepath):
    f = h5py.File(filepath, 'w')
    f.create_dataset('p', data=edist.g)
    f.create_dataset('t', data=edist.t)
    f.create_dataset('x', data=edist.x)
    f.create_dataset('y', data=edist.y)
    f.create_dataset('xp', data=edist.xp)
    f.create_dataset('yp', data=edist.yp)
    f.create_dataset('charge', data = edist.charge())
    f.close()
    
def read_edist_hdf5(filepath, charge=None):
    
    edist = GenesisElectronDist()
    with h5py.File(filepath, 'r') as h5:
        
        edist.g =  h5.get('p')[:]
        edist.t =  h5.get('t')[:]
        edist.x =  h5.get('x')[:]
        edist.y =  h5.get('y')[:]
        edist.xp =  h5.get('xp')[:]
        edist.yp =  h5.get('yp')[:]
        
        if charge is not None:
            charge = h5.get('charge')[:]
    
    edist.part_charge = charge / edist.g.size
    return edist
