from ocelot.optimizer.mint.opt_objects import Target
import numpy as np
import time


class XFELTarget(Target):
    def __init__(self, mi=None, dp=None, eid=None):
        super(XFELTarget, self).__init__(eid=eid)
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
        self.penalties = []
        self.times = []
        self.alarms = []
        self.values = []

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
        self.penalties.append(pen)
        self.times.append(time.time())
        self.values.append(sase)
        self.alarms.append(alarm)
        return pen

    def get_value(self):
        """
        changeable
        :return:
        """
        values = np.array([dev.get_value() for dev in self.devices])
        return 2*np.sum(np.exp(-np.power((values - np.ones_like(values)), 2) / 5.))
        #value = self.mi.get_value(self.eid)

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
        """
        changeable
        :return:
        """
        return 0

    def get_energy(self):
        return 3
