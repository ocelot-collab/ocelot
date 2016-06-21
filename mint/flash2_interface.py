'''
machine interface
includes online optimizer, response measurement and other stuff
'''

import swig.dcs as dcs
from pylab import *

'''
Any machine interface class should implement following methods  
'''
class MachineInterface:
    def __init__(self):
        pass
    def get_alarms(self):
        '''
        return an array of values for all active bpms
        values range between (0,1) 1 corresponds to alarm level
        '''
        pass
    def get_sase(self, detector_name='default'):
        '''
        return SASE pulse energy
        units and averaging can be arbitrary        '''

        pass
    def get_value(self, device_name):
        pass
    def set_value(self, device_name, val):
        pass
    
    
class FLASH1MachineInterface():
    def __init__(self):
        
        self.debug = False
        
        self.blm_names = ['14L.SMATCH','14R.SMATCH',
                          '1L.UND1', '1R.UND1',
                          '1L.UND2', '1R.UND2', 
                          '1L.UND3', '1R.UND3', 
                          '1L.UND4', '1R.UND4',
                          '1L.UND5', '1R.UND5',
                          '1L.UND6', '1R.UND6',
                          '1SFUND1', '1SFUND2', 
                          '1SFUND3','1SFUND4',
                          '1SFELC','3SFELC','4SFELC',
                          '10SMATCH','3SDUMP']

    def get_alarms(self):
        alarm_vals = np.zeros(len(self.blm_names))
        for i in xrange(len(self.blm_names)):
            blm_channel = 'TTF2.DIAG/BLM/'+self.blm_names[i]+'/CH00.TD'
            blm_alarm_ch  = ('TTF2.DIAG/BLM/'+self.blm_names[i]).replace('BLM', 'BLM.ALARM') + '/THRFHI'
            blm_alarm_ch  = ('FLASH.DIAG/BLM/'+self.blm_names[i]).replace('BLM', 'BLM.ALARM') + '/THRFHI'
            if (self.debug): print ('reading alarm channel', blm_alarm_ch)
            alarm_val = dcs.get_device_val(blm_alarm_ch) * 1.25e-3 # alarm thr. in Volts
            if (self.debug): print ('alarm:', alarm_val)
    
            h = np.array(dcs.get_device_td(blm_channel))
    
            alarm_vals[i] = np.max( np.abs(h) ) / alarm_val 
            
        return alarm_vals
    
    def get_sase_pos(self):

        x1 = dcs.get_device_val('TTF2.FEL/GMDPOSMON/TUNNEL/IX.POS')
        y1 = dcs.get_device_val('TTF2.FEL/GMDPOSMON/TUNNEL/IY.POS')

        x2 = dcs.get_device_val('TTF2.FEL/GMDPOSMON/BDA/IX.POS')
        y2 = dcs.get_device_val('TTF2.FEL/GMDPOSMON/BDA/IY.POS')
    
        return [ (x1,y1), (x2,y2) ] 

    def get_spectrum(self, f=None, detector='tunnel_default'):

        f_min = 13.0 # spectrum window (nm). TODO: replace with readout
        f_max = 14.0
        
        spec = np.array(dcs.get_device_td('TTF2.EXP/PBD.PHOTONWL.ML/WAVE_LENGTH/VAL.TD'))
    
        if f == None:
            f = np.linspace(f_min, f_max, len(spec))
    
        return f, spec
 
    def get_value(self, device_name):
        pass
    
    def set_value(self, device_name, val):
        pass
    
    
class FLASH2MachineInterface():
    def __init__(self):
        self.blm_names = ['3.1FL2SASE3','3.2FL2SASE3',
                          '3.1FL2SASE4','3.2FL2SASE4',
                          '3.1FL2SASE5','3.2FL2SASE5',
                          '3.1FL2SASE6','3.2FL2SASE6',
                          '3.1FL2SASE7','3.2FL2SASE7',
                          '3.1FL2SASE8','3.2FL2SASE8',
                          '3.1FL2SASE9','3.2FL2SASE9',
                          '3.1FL2SASE10','3.2FL2SASE10',
                          '3.1FL2SASE11','3.2FL2SASE11',
                          '3.1FL2SASE12','3.2FL2SASE12',
                          '3.1FL2SASE13','3.2FL2SASE13',
                          '3.1FL2SASE14','3.2FL2SASE14']
 
    def get_alarms(self):
        alarm_vals = np.zeros(len(self.blm_names))
        for i in range(len(self.blm_names)):
            alarm_val = 2000
            blm_channel = 'FLASH.DIAG/BLM/'+self.blm_names[i]+'/SIGNAL.TD'     
            h = np.array(dcs.get_device_td(blm_channel))
            alarm_vals[i] = np.max( np.abs(h) ) / alarm_val 
        
        return alarm_vals 



'''
flight simulator implementation of the machine interface
'''
class FlightSimulator:
    def __init__(self, lattice, beam):
        self.lattice = lattice
        self.beam = beam
    def get_alarms(self):
        return 0.0
    def get_sase(self, detector_name='default'):
        return 0.0
