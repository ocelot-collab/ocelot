"""
Objective Function

S.Tomin, 2017
"""

from ocelot.optimizer.mint.opt_objects import Target
import numpy as np
import time


class XFELTarget(Target):
    """
    Objective function

    :param mi: Machine interface
    :param dp: Device property
    :param pen_max: 100, maximum penalty
    :param niter: 0, calls number get_penalty()
    :param penalties: [], appending penalty
    :param times: [], appending the time evolution of get_penalty()
    :param nreadings: 1, number of objective function readings
    :param interval: 0 (secunds), interval between readings
    """
    def __init__(self, mi=None, dp=None, eid="x57**2 + y57**2 + x59**2 + y59"):
        super(XFELTarget, self).__init__(eid=eid)

        self.mi = mi
        self.dp = dp
        self.debug = False
        self.kill = False
        self.pen_max = 100
        self.clean()
        self.nreadings = 1
        self.interval = 0.0

    def get_alarm(self):
        """
        Method to get alarm level (e.g. BLM value).

        alarm level must be normalized: 0 is min, 1 is max

        :return: alarm level
        """
        return 0

    def get_value(self):
        """
        Method to get signal of target function (e.g. SASE signal).

        :return: value
        """
        x57 = self.mi.get_value("XFEL.DIAG/ORBIT/BPMA.57.I1/X.SA1")
        y57 = self.mi.get_value("XFEL.DIAG/ORBIT/BPMA.57.I1/Y.SA1")
        x59 = self.mi.get_value("XFEL.DIAG/ORBIT/BPMA.59.I1/X.SA1")
        y59 = self.mi.get_value("XFEL.DIAG/ORBIT/BPMA.59.I1/Y.SA1")
        return -np.sqrt(x57 ** 2 + y57 ** 2 + x59 ** 2 + y59 ** 2)

        # values = np.array([dev.get_value() for dev in self.devices])
        # return 2*np.sum(np.exp(-np.power((values - np.ones_like(values)), 2) / 5.))
        # value = self.mi.get_value(self.eid)


    def get_value_test(self):
        """
        For testing

        :return:
        """
        values = np.array([dev.get_value() for dev in self.devices])
        value = 2*np.sum(np.exp(-np.power((values - np.ones_like(values)), 2) / 5.))
        value = value * (1. + (np.random.rand(1)[0] - 0.5) * 0.001)
        return value 


    def get_penalty(self):
        """
        Method to calculate the penalty on the basis of the value and alarm level.

        penalty = -get_value() + alarm()


        :return: penalty
        """
        sase = 0.
        for i in range(self.nreadings):
            sase += self.get_value()
            time.sleep(self.interval)
        sase = sase/self.nreadings
        print("SASE", sase)
        alarm = self.get_alarm()
        if self.debug: print('alarm:', alarm)
        if self.debug: print('sase:', sase)
        pen = 0.0
        if alarm > 1.0:
            return self.pen_max
        if alarm > 0.7:
            return alarm * self.pen_max / 2.
        pen += alarm
        pen -= sase
        if self.debug: print('penalty:', pen)
        self.niter += 1
        # print("niter = ", self.niter)
        self.penalties.append(pen)
        self.times.append(time.time())
        self.values.append(sase)
        self.alarms.append(alarm)
        return pen

    def get_spectrum(self):
        return [0, 0]

    def get_stat_params(self):
        # spetrum = self.get_spectrum()
        # ave = np.mean(spetrum[(2599 - 5 * 120):-1])
        # std = np.std(spetrum[(2599 - 5 * 120):-1])
        ave = self.get_value()
        std = 0.1
        return ave, std

    def get_energy(self):
        return 3

    def clean(self):
        self.niter = 0
        self.penalties = []
        self.times = []
        self.alarms = []
        self.values = []
