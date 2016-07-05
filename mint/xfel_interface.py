'''
machine interface
includes online optimizer, response measurement and other stuff
'''
try:
    # in server "doocsdev12" set environment
    #  $ export PYTHONPATH=/home/ttflinac/user/python-2.7/Debian/
    import pydoocs
except:
    print ('error importing doocs library')
"""
import sys
sys.path.append("/home/ttflinac/user/tomins/repository/pyqtgraph-0.9.10/build/lib")
sys.path.append("/home/ttflinac/user/tomins/repository")
"""
import re
from pylab import *
import time
import pickle
from threading import Thread, Lock
from ocelot.utils.db import *


class XFELMachineInterface():
    def __init__(self):
        
        self.debug = False

        self.blm_names = ['BLM.23.I1', 'BLM.25L.I1',
                          'BLM.25R.I1', 'BLM.48.I1',
                          'BLM.49.1.I1', 'BLM.49.2.I1',
                          'BLM.51.1.I1', 'BLM.51.2.I1',
                          'BLM.55.I1', 'BLM.56.I1',
                          'BLM.60.I1', 'BLM.63.I1',
                          'BLM.65U.I1','BLM.65D.I1',
                          "BLM.65L.I1", 'BLM.65R.I1',
                          #'BLM.66.I1'
                          ]
        self.mutex = Lock()


    def init_corrector_vals(self, correctors):
        vals = np.zeros(len(correctors))
        for i in range(len(correctors)):
            #mag_channel = 'TTF2.MAGNETS/STEERER/' + correctors[i] + '/PS'
            #self.mutex.acquire()
            vals[i] = self.get_value(correctors[i])#pydoocs.read(mag_channel)["data"]

            #self.mutex.release()
        return vals

    def get_cav_ampl(self, cav):
        return pydoocs.read("FLASH.RF/LLRF.CONTROLLER/PVS." + cav + "/AMPL.SAMPLE")["data"]

    def get_cav_phase(self, cav):
        return pydoocs.read("FLASH.RF/LLRF.CONTROLLER/PVS." + cav + "/PHASE.SAMPLE")["data"]

    def get_cavity_info(self, cavs):
        ampls = [0.0]*len(cavs)#np.zeros(len(correctors))
        phases = [0.0]*len(cavs)#np.zeros(len(correctors))
        for i in range(len(cavs)):
            #ampl_channel = 'FLASH.RF/LLRF.CONTROLLER/CTRL.' + cavs[i] + '/SP.AMPL'
            #phase_channel = 'FLASH.RF/LLRF.CONTROLLER/CTRL.' + cavs[i] + '/SP.PHASE'
            ampls[i] = self.get_cav_ampl(cavs[i])
            phases[i] = self.get_cav_phase(cavs[i])
        return ampls, phases

    def get_gun_energy(self):
        gun_energy = pydoocs.read("FLASH.RF/LLRF.ENERGYGAIN.ML/GUN/ENERGYGAIN.FLASH1")['data']
        gun_energy = gun_energy*0.001 # MeV -> GeV
        return gun_energy

    def get_bpms_xy(self, bpms):
        X = [0.0]*len(bpms)#np.zeros(len(correctors))
        Y = [0.0]*len(bpms)
        for i in range(len(bpms)):
            ch = 'XFEL.DIAG/ORBIT/' + bpms[i]
            #print(ch)
            X[i] = pydoocs.read(ch + "/X.SA1")['data']*0.001 # mm -> m
            Y[i] = pydoocs.read(ch + "/Y.SA1")['data']*0.001 # mm -> m
        return X, Y

    #def get_quads_current(self, quads):
    #    vals = np.zeros(len(quads))
    #    for i in range(len(quads)):
    #        mag_channel = 'TTF2.MAGNETS/QUAD/' + quads[i]# + '/PS'
    #        vals[i] = pydoocs.read(mag_channel + "/PS")['data']
    #    return vals
    #
    #def get_bends_current(self, bends):
    #    vals = [0.0]*len(bends)#np.zeros(len(correctors))
    #    for i in range(len(bends)):
    #        mag_channel = 'TTF2.MAGNETS/DIPOLE/' + bends[i]# + '/PS'
    #        vals[i] = pydoocs.read(mag_channel + "/PS")['data']
    #    return vals
    #
    #def get_sext_current(self, sext):
    #    vals = [0.0]*len(sext)#np.zeros(len(correctors))
    #    for i in range(len(sext)):
    #        mag_channel = "TTF2.MAGNETS/SEXT/" + sext[i]
    #        vals[i] = pydoocs.read(mag_channel + "/PS")['data']
    #    return vals


    def get_alarms(self):
        alarm_vals = np.zeros(len(self.blm_names))
        for i in range(len(self.blm_names)):
            blm_channel = 'XFEL.DIAG/BLM/'+self.blm_names[i]+'/SIGNAL.SA1'
            blm_alarm_ch  = 'XFEL.DIAG/BLM/'+self.blm_names[i] + '/SINGLE_ALARM_THRESHOLD'
            if (self.debug): print('reading alarm channel', blm_alarm_ch)
            #print('reading alarm channel', blm_alarm_ch)
            alarm_val = pydoocs.read(blm_alarm_ch)['data']
            if (self.debug): print ('alarm:', alarm_val)
            #print("blm channel: ", blm_channel)
            val = pydoocs.read(blm_channel)['data']
            #h = np.array([x[1] for x in sample])
            #alarm_vals[i] = np.max( np.abs(h) ) / alarm_val
            alarm_vals[i] = val/alarm_val
        return alarm_vals


    def get_charge(self):
        charge = pydoocs.read('TTF2.FEEDBACK/LONGITUDINAL/MONITOR2/MEAN_AVG')
        #print("charge = ", charge["data"], " nQ")
        return charge["data"]


    #def get_final_energy(self):
    #    final_energy = pydoocs.read("TTF2.FEEDBACK/LONGITUDINAL/MONITOR11/MEAN_AVG")
    #    #print("final_energy = ", final_energy["data"], "MeV")
    #    return final_energy["data"]

    #def get_sol_value(self):
    #    sol = pydoocs.read("TTF2.MAGNETS/SOL/1GUN/PS")
    #    #print("sol = ", sol["data"], "A")
    #    return sol["data"]

    def get_nbunches(self):
        nbunches = pydoocs.read("FLASH.DIAG/TIMER/FLASHCPUTIME1.0/NUMBER_BUNCHES.1")
        #print("nbunches = ", nbunches["data"])
        return nbunches["data"]

    def get_value_ps(self, device_name):
        ch = 'TTF2.MAGNETS/STEERER/' + device_name + '/PS'
        self.mutex.acquire()
        val = pydoocs.read(ch)['data']
        self.mutex.release()
        return val

    def get_value(self, device_name):
        #ch = 'XFEL.MAGNETS/MAGNET.ML/' + device_name + '/CURRENT.RBV'
        ch = 'XFEL.MAGNETS/MAGNET.ML/' + device_name + '/KICK_MRAD.SP'
        #print("getting value = ", ch)

        #self.mutex.acquire()
        val = pydoocs.read(ch)
        #print(val)
        #['data']
        #self.mutex.release()
        return val["data"]
    
    def set_value(self, device_name, val):
        #ch = 'XFEL.MAGNETS/MAGNET.ML/' + device_name + '/CURRENT.SP'
        ch = 'XFEL.MAGNETS/MAGNET.ML/' + device_name + '/KICK_MRAD.SP'
        print (ch, val)
        #self.mutex.acquire()
        pydoocs.write(ch, str(val))
        #self.mutex.release()
        return 0
 
class XFELDeviceProperties:
    def __init__(self):
        self.stop_exec = False
        self.save_machine = False
        self.patterns = {}
        self.limits = {}
        """
        self.patterns['launch_steerer'] = re.compile('[HV][0-9]+SMATCH')
        self.limits['launch_steerer'] = [-4,4]
        
        self.patterns['intra_steerer'] = re.compile('H3UND[0-9]')
        self.limits['intra_steerer'] = [-5.0,-2.0]
        
        self.patterns['QF'] = re.compile('Q5UND1.3.5')
        self.limits['QF'] = [1,7]
        """

    def set_limits(self, dev_name, limits):
        self.patterns[dev_name] = re.compile(dev_name)
        #print(self.patterns[dev_name])
        self.limits[dev_name] = limits
        #print("inside dp set = ", self.patterns[dev_name], self.limits)

    def get_limits(self, device):
        #print(self.limits)
        for k in self.patterns.keys():
            #print('testing', k)
            if self.patterns[k].match(device) != None:
                #print("inside dp get = ", device, self.limits[k])
                return self.limits[k]
        return [-2, 2]

    def get_polarity(self, quads):
        vals = [0.0]*len(quads)#np.zeros(len(correctors))
        for i in range(len(quads)):
            mag_channel = 'TTF2.MAGNETS/QUAD/' + quads[i]# + '/PS'
            vals[i] = pydoocs.read(mag_channel + "/PS.Polarity")['data']
        return vals

    def get_type_magnet(self, quads):
        vals = [0.0]*len(quads)#np.zeros(len(correctors))
        for i in range(len(quads)):
            mag_channel = 'TTF2.MAGNETS/QUAD/' + quads[i]# + '/PS'
            vals[i] = pydoocs.get(mag_channel + "/DEVTYPE")['data']
        return vals


'''
test interface
'''
class TestInterface:
    def __init__(self):
        pass
    def get_alarms(self):
        return np.random.rand(4)#0.0, 0.0, 0.0, 0.0]

    def get_sase(self, detector=None):
        return 0.0

    def init_corrector_vals(self, correctors):
        vals = [0.0]*len(correctors)
        return vals

    def get_cor_value(self, devname):
        return np.random.rand(1)[0]

    def get_value(self, device_name):
        print("get value", device_name)
        return np.random.rand(1)[0]

    def set_value(self, device_name, val):
        return 0.0

    def get_quads_current(self, device_names):
        return np.random.rand(len(device_names))

    def get_bpms_xy(self, bpms):
        X = np.zeros(len(bpms))
        Y = np.zeros(len(bpms))
        return X, Y

    def get_sase_pos(self):
        return [(0,0), (0, 0)]

    def get_gun_energy(self):
        return 0.

    def get_cavity_info(self, devnames):
        X = np.zeros(len(devnames))
        Y = np.zeros(len(devnames))
        return X, Y

    def get_charge(self):
        return 0.

    def get_wavelangth(self):
        return 0.

    def get_bc2_pyros(self):
        return (0, 0)

    def get_bc3_pyros(self):
        return (0, 0)

    def get_final_energy(self):
        return 0

    def get_sol_value(self):
        time.sleep(0.1)
        return 0

    def get_nbunches(self):
        return 0
    def get_value_ps(self, device_name):
        return 0.
