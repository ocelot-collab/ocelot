from ocelot.common.xio import *
import numpy as np
import numpy

class Xdb:
    def __init__(self, index_file=None, mode='r'):
        self.index_file = index_file
        self.file = h5py.File(self.index_file, mode)
    
    def create_index(self):

        #self.file = h5py.File(self.index_file,'w')
        
        g1 = self.file.create_group('Undulators')
        g1.attrs['Description'] = 'information on undulator section layout'
        
        g2 = self.file.create_group('Beams')
        g2.attrs['Description'] = 'information on beam parameters'
        
        g3 = self.file.create_group('SR')
        g3.attrs['Description'] = 'spontaneous SR calculations'
        
        g4 = self.file.create_group('FEL')
        g4.attrs['Description'] = 'SASE FEL calculations'
        
        #self.file.close()
        
    def add_beam(self, name, params):
        '''
        beam parameters
        '''
        g = self.file['Beams']
        gb = g.create_group(name)
        
        '''
        for p in params.keys():
            if p == 'input_file': 
                gb.attrs[p] = open(params[p]).read()
        '''
        
    def add_undulator(self, name, params):
        '''
        physical undulator line
        '''
        #self.file = h5py.File(self.index_file,'r+')
        g = self.file['Undulators']
        gu = g.create_group(name)
        for p in params.keys():
            if p == 'input_file': 
                gu.attrs[p] = open(params[p]).read()
            
        #self.file.close()
        
    def add_undulator_config(self, und_name, config_name, config):
        '''
        particular setup, i.e. settings, minor hardware modifications, etc.
        '''
        #self.file = h5py.File(self.index_file,'r+')
        g = self.file['Undulators']
        gu = g[und_name]
        gu2 = gu.create_group(config_name)
        for p in config.keys():
            if p == 'input_file': 
                gu2.attrs[p] = open(config[p]).read()

        #self.file.close()

    def read_undulator_config(self, path):
        g = self.file['Undulators/' + path]
        return g.attrs['input_file']
        #path = path.split('/')


    def submit_1d(self, g, params, name):
        print 'creating dataset ', name
        data = params[name]
        print data.shape, data.dtype

        if name in g: del g[name]
        #dset = g.create_dataset(name, (len(data),), dtype='f')
        dset = g.create_dataset(name, data.shape, dtype=data.dtype)
            
        dset[...] = data
        
        
    def add_fel_calculation(self, path, params, root='FEL'):
        g = self.file[root]
        
        if path in g: del g[path]
        gu = g.create_group(path)
        
        
        # TODO: add validation        
        for p in params.keys(): 
            data = params[p]
            print data.__class__
            if data.__class__ in (str, float, complex, int):
                gu.attrs[p] = data
            if data.__class__ == numpy.ndarray:
                self.submit_1d(gu, params, p)
        
    def list_available_fel_parameters(self, root='FEL'):
        g = self.file[root]
        pars = {}
        for p in g.keys():
            pars[p] = g[p].keys()
        return pars

        
        # parameters genesis_input, power_exit, sig_power_exit

    def add_sr_calculation(self, genesis_output, genesis_input):
        pass

    def list_available_sr_parameters(self, root='SR'):
        g = self.file[root]
        pars = {}
        for p in g.keys():
            pars[p] = g[p].keys()
        return pars
        
    
def list_available_sr_parameters(index_file):
    pass

    
