'''
machine interface
LCLS
'''

from pylab import *
import epics


class LCLSMachineInterface:
    
    def __init__(self):
        pass
  
    def get_alarms(self):
        return [0.0,]
    
    def get_sase(self, detector_name='default'):
        '''
        at lcls the repetition is  120Hz and the readout buf size is 2400
        last 120 entries correspond to pulse energies over past 1 second  
        '''
        data = epics.caget("GDET:FEE1:241:ENRCHSTBR")[-120:]
        print ("gmd averaging over n points: ", len(data) )
        return np.mean( data ) 

    
    def get_value(self, device_name):
        return epics.caget(device_name + ':BCTRL')
    
    def set_value(self, device_name, val):
        return epics.caput(device_name + ':BCTRL', val)
    
    def init_corrector_vals(self, correctors):
        vals = [0.0]*len(correctors)#np.zeros(len(correctors))
        for i in range(len(correctors)):
            mag_channel = correctors[i] + ':BCTRL'
            vals[i] = epics.caget(mag_channel)
        return vals

    
class LCLSDeviceProperties:
    def __init__(self):
        pass
    def get_limits(self, device):
        if device.startswith("XCOR:"): return [-0.006, 0.006]
        if device.startswith("YCOR:"): return [-0.006, 0.006]
        
        return [-100, 100]