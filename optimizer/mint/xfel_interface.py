"""
XFEL machine interface
S.Tomin, 2017
"""
try:
    # in server "doocsdev12" set environment
    #  $ export PYTHONPATH=/home/ttflinac/user/python-2.7/Debian/
    import pydoocs
except:
    print ('error importing doocs library')

#import re
import numpy as np
from threading import Thread, Lock
from ocelot.optimizer.mint.opt_objects import Device

class AlarmDevice(Device):
    """
    Devices for getting information about Machine status
    """
    def __init__(self, eid=None):
        super(AlarmDevice, self).__init__(eid=eid)


class XFELMachineInterface():
    """
    Machine Interface for European XFEL
    """
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

    def get_value(self, channel):
        """
        Getter function for XFEL.

        :param channel: (str) String of the devices name used in doocs
        :return: Data from pydoocs.read(), variable data type depending on channel
        """
        val = pydoocs.read(channel)
        return val["data"]

    def set_value(self, channel, val):
        """
        Method to set value to a channel

        :param channel: (str) String of the devices name used in doocs
        :param val: value
        :return: None
        """
        #print("SETTING")
        pydoocs.write(channel, val)
        return

    def get_bpms_xy(self, bpms):
        """
        Method for getting bmps data

        :param bpms: list of string. BPMs names
        :return: X, Y - two arrays in [m]
        """
        X = [0.0]*len(bpms)
        Y = [0.0]*len(bpms)
        for i in range(len(bpms)):
            ch = 'XFEL.DIAG/ORBIT/' + bpms[i]
            X[i] = pydoocs.read(ch + "/X.SA1")['data']*0.001 # mm -> m
            Y[i] = pydoocs.read(ch + "/Y.SA1")['data']*0.001 # mm -> m
        return X, Y


    def get_alarms(self):
        """
        Method for getting BLMs level. BLMs are predefined

        :return: arrays of the BLMs level
        """
        alarm_vals = np.zeros(len(self.blm_names))
        for i in range(len(self.blm_names)):
            blm_channel = 'XFEL.DIAG/BLM/'+self.blm_names[i]+'/SIGNAL.SA1'
            blm_alarm_ch  = 'XFEL.DIAG/BLM/'+self.blm_names[i] + '/SINGLE_ALARM_THRESHOLD'
            if (self.debug): print('reading alarm channel', blm_alarm_ch)
            alarm_val = pydoocs.read(blm_alarm_ch)['data']
            if (self.debug): print ('alarm:', alarm_val)
            val = pydoocs.read(blm_channel)['data']
            alarm_vals[i] = val/alarm_val
        return alarm_vals

    def get_charge(self):
        return self.get_value("XFEL.DIAG/CHARGE.ML/TORA.25.I1/CHARGE.SA1")


class XFELDeviceProperties:
    """
    Class for getting device properties
    """

    def __init__(self, ui=None):
        """

        :param ui:  ui=None, user interface
        """
        self.limits = {}
        self.ui = ui

    def get_limits(self, dev_name):
        """
        Method for getting limits of the devices from the GUI

        :param dev_name: str. Device name as in the GUI
        :return: [min, max]
        """
        if self.ui != None:
            lims = self.ui.get_limits(dev_name)
            if lims == None:
                return [0., 0.]
        else:
            return [-100., 100.]
        return lims


# test interface

class TestMachineInterface:
    """
    Machine interface for testing
    """
    def __init__(self):
        self.data = 1.
        pass
    def get_alarms(self):
        return np.random.rand(4)#0.0, 0.0, 0.0, 0.0]

    def get_value(self, device_name):
        """
        Testing getter function for XFEL.

        :param channel: (str) String of the devices name used in doocs
        :return: Data from pydoocs.read(), variable data type depending on channel
        """
        #if "QUAD" in device_name:
        #    return 0
        return np.random.rand(1)[0]-0.5 #self.data

    def set_value(self, device_name, val):
        """
        Testing Method to set value to a channel

        :param channel: (str) String of the devices name used in doocs
        :param val: value
        :return: None
        """
        #print("set:", device_name,  "-->", val)
        self.data += np.sqrt(val**2)
        return 0.0

    def get_bpms_xy(self, bpms):
        """
        Testing method for getting bmps data

        :param bpms: list of string. BPMs names
        :return: X, Y - two arrays in [m]
        """
        X = np.zeros(len(bpms))
        Y = np.zeros(len(bpms))
        return X, Y

    def get_charge(self):
        return 0

class TestDeviceProperties:
    def __init__(self, ui=None):
        """
        For testing.
        Class for getting device properties

        :param ui:  ui=None, user interface
        """
        self.patterns = {}
        self.limits = {}
        self.ui = ui

    def get_limits(self, dev_name):
        """
        Testing method for getting limits of the devices from the GUI

        :param dev_name: str. Device name as in the GUI
        :return: [min, max]
        """
        if self.ui != None:
            lims = self.ui.get_limits(dev_name)
            if lims == None:
                return [0., 0.]
        else:
            return [-100., 100.]
        return lims

