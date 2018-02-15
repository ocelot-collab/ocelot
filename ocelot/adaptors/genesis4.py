import time
import os
import socket
import copy
import h5py

from ocelot.common.globals import *

class Genesis4Input:
    '''
    Genesis input files storage object
    '''
    def __init__(self):
        pass
    


class Genesis4Output:
    '''
    Genesis input files storage object
    '''
    def __init__(self):
        self.h5 = None #hdf5 pointer
        
    def close(self):
        self.h5.close()
        
    @property
    def filename(self):
        return self.h5.filename
    
    @property
    def lambdaref(self):
        return self.h5['Global/lambdaref'][0]
    
    @property
    def phenref(self):
        return h_eV_s * speed_of_light / self.lambdaref
    
    @property
    def I(self):
        return self.h5['Beam/current'][:]
    
    @property
    def beam_charge(self):
        return np.trapz(self.I, self.t)[0]
    
    @property
    def rad_power(self):
        return self.h5['Field/power'][:]
    
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

def read_gout4(filename):
    out = Genesis4Output()
    out.h5 = h5py.File(file_path, 'r')
    
    out.z = out.h5['Lattice/zplot'][()]
    
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
        out.s = np.linspace(s0, out.h5['Global/slen'][()]+s0, sn)[np.newaxis,:]
    
    return out


