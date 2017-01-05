from ocelot.optimizer.mint.opt_objects import *


class TestTarget_new(Target):
    def __init__(self, mi=None, dp=None, eid=None):
        super(TestTarget_new, self).__init__(eid=eid)
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
        self.x.append(time.time())
        return pen

    def get_value(self):
        values = np.array([dev.get_value() for dev in self.devices])
        #print("hi4")
        return 2*np.sum(np.exp(-np.power((values - np.ones_like(values)), 2) / 5.))

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
