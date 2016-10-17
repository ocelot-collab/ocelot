"""
S.Tomin
Object function and devices
"""

import numpy as np
from ocelot.optimizer.mint.lcls_interface import *
import time


mi = TestLCLSMachineInterface()
dp = TestLCLSDeviceProperties()

#mi = LCLSMachineInterface()
#dp = LCLSDeviceProperties()

class Device(object):
    def __init__(self, eid=None):
        self.eid = eid
        self.id = eid

    def set_value(self, x):
        pass

    def get_value(self):
        return

    def check_limits(self, value):
        return False

    def state(self):
        return


class SLACDevice(Device):
   def __init__(self, eid=None):
       super(SLACDevice, self).__init__(eid=eid)
       self.test_value = 0.
       self.data = []
       self.time = []
       self.nsets = 0

   def set_value(self, x):
       self.data.append(x)
       self.time.append(time.time())
       mi.set_value(self.eid, x)

   def get_value(self):
       return mi.get_value(self.eid)

   def check_limits(self, value):
       return False

   def state(self):
       return

# for testing
class TestDevice(Device):
    def __init__(self, eid=None):
        super(TestDevice, self).__init__(eid=eid)
        self.test_value = 0.
        self.data = []
        self.time = []
        self.nsets = 0
        pass

    def get_value(self):
        return self.test_value

    def set_value(self, value):
        self.data.append(value)
        self.nsets += 1
        self.time.append(self.nsets)
        self.test_value = value

    def check_limits(self, value):
        limits = self.get_limits()

        if value < limits[0] or value > limits[1]:
            print ('limits exceeded', value, limits[0] , value, limits[1])
            return True
        return False

    def get_limits(self):
        return [-100, 100]


class Target(object):
    def __init__(self,  eid=None):
        """
        :param mi: Machine interface
        :param dp: Device property
        :param eid: ID
        """
        self.eid = eid
        self.id = eid

    def get_value(self):
        return 0

    def get_penalty(self):
        pen = - self.get_value()
        return pen

class SLACTarget(Target):
    def __init__(self, mi=None, dp=None, eid=None):
        super(SLACTarget, self).__init__(eid=eid)
        """
        :param mi: Machine interface
        :param dp: Device property
        :param eid: ID
        """
        self.secs_to_ave = 2
        self.debug = False
        self.kill = False
        self.pen_max = 100
        self.niter = 0
        self.y = []
        self.x = []

    def get_penalty(self):
        sase = self.get_value()
        alarm = self.get_alarm()
        if self.debug: print('alarm:', alarm)
        if self.debug: print('sase:', sase)
        pen = 0.0
        if alarm > 1.0:
            return self.pen_max
        if alarm > 0.7:
            return alarm * 50.0
        pen += alarm
        pen -= sase
        if self.debug: print('penalty:', pen)
        self.niter += 1
        print("niter = ", self.niter)
        self.y.append(pen)
        self.x.append(self.niter)
        return pen

    def get_value(self, seconds=None):
        """
        Returns data for the ojective function (sase) from the selected detector PV.

        At lcls the repetition is  120Hz and the readout buf size is 2800.
        The last 120 entries correspond to pulse energies over past 1 second.

        Args:
                seconds (float): Variable input on how many seconds to average data

        Returns:
                Float of SASE or other detecor measurment
        """
        datain = mi.get_value(self.eid)
        try: #try to average over and array input
            if seconds == None: #if a resquested seconds is passed
                dataout = np.mean(datain[-(self.secs_to_ave*120):])
                sigma   = np.std( datain[-(self.secs_to_ave*120):])
            else:
                dataout = np.mean(datain[-(seconds*120):])
                sigma   = np.std( datain[-(seconds*120):])
        except: #if average fails use the scaler input
            print ("Detector is not a waveform PV, using scalar value")
            dataout = datain
            sigma   = -1
        return dataout

    def get_stat_params(self):
        #get the current mean and std of the chosen detector
        obj_func = mi.get_value(self.eid)
        try:
            std = np.std(  obj_func[(2599-5*120):-1])
            ave = np.mean( obj_func[(2599-5*120):-1])
        except:
            print "Detector is not a waveform, Using scalar for hyperparameter calc"
            ave = obj_func
            # Hard code in the std when obj func is a scalar
            # Not a great way to deal with this, should probably be fixed
            std = 0.1

        print ("DETECTOR AVE", ave)
        print ("DETECTOR STD", std)

        return ave, std

    def get_alarm(self):
        return 0

    def get_energy(self):
        return mi.get_energy()


class TestTarget(Target):
    def __init__(self, mi=None, dp=None, eid=None):
        super(TestTarget, self).__init__(eid=eid)
        """
        :param mi: Machine interface
        :param dp: Device property
        :param eid: ID
        """
        self.mi = mi
        self.dp = dp
        self.debug = False
        self.kill = False
        self.pen_max = 100
        self.niter = 0
        self.y = []
        self.x = []

    def get_penalty(self):
        sase = self.get_value()
        alarm = self.get_alarm()
        if self.debug: print('alarm:', alarm)
        if self.debug: print('sase:', sase)
        pen = 0.0
        if alarm > 1.0:
            return self.pen_max
        if alarm > 0.7:
            return alarm * 50.0
        pen += alarm
        pen -= sase
        if self.debug: print('penalty:', pen)
        self.niter += 1
        print("niter = ", self.niter)
        self.y.append(pen)
        self.x.append(self.niter)
        return pen

    def get_value(self):
        values = np.array([dev.get_value() for dev in self.devices])
        return np.sum(np.exp(-np.power((values - np.ones_like(values)), 2) / 5.))

    def get_spectrum(self):
        return [0, 0]

    def get_stat_params(self):
        #spetrum = self.get_spectrum()
        #ave = np.mean(spetrum[(2599 - 5 * 120):-1])
        #std = np.std(spetrum[(2599 - 5 * 120):-1])
        ave = self.get_value()
        std = 0.1
        return ave, std

    def get_alarm(self):
        return 0

    def get_energy(self):
        return 3