"""
Main Ocelot optimization file
Contains the setup for using the scipy.optimize package run simplex and other algorothms
Modifiedi for use at LCLS from Ilya's version
"""

from time import sleep

import numpy as np
from ocelot.optimizer.mint.opt_objects import *
from scipy import optimize


class Logger(object):
    def __init__(self, log_file):
        self.log_file = log_file
        self.f = open(self.log_file, 'a')

    def log_start(self, dev_ids, method, x_init, target_ref):
        self.f.write('\n*** optimization step ***\n')
        self.f.write(str(dev_ids) + '\n')
        self.f.write(method + '\n')
        self.f.write('x_init =' + str(x_init) + '\n')
        self.f.write('target_ref =' + str(target_ref) + '\n')

    def log(self, data):
        self.f.write(data)

    def log_fin(self, target_new):
        self.f.write('target_new=' + str(target_new) + '\n')
        self.f.close()



class Minimizer(object):
    def __init__(self):
        self.max_iter = 100

    def minimize(self, error_func, x):
        pass


class Simplex(Minimizer):
    def __init__(self):
        super(Simplex, self).__init__()
        self.xtol = 1e-5

    def minimize(self,  error_func, x):
        optimize.fmin(error_func, x, maxiter=self.max_iter, maxfun=self.max_iter, xtol=self.xtol)


class Optimizer:
    def __init__(self, normalize=False):
        self.debug   = False
        self.minimizer = Simplex()
        self.logging = False
        self.kill    = False #intructed by tmc to terminate thread of this class
        self.log_file = "log.txt"
        self.logger = Logger(self.log_file)
        self.devices = []
        self.target = None
        self.timeout = 0
        self.x_data = []
        self.y_data = []

    def eval(self, seq, logging=False, log_file=None):
        """
        Run the sequence of tuning events
        """
        for s in seq:
            s.apply()

    def exceed_limits(self, x):
        for i in range(len(x)):
            print('{0} x[{1}]={2}'.format(self.devices[i].id, i, x[i]))
            if self.devices[i].check_limits(x[i]):
                return True
        return False

    def set_values(self, x):
        for i in range(len(self.devices)):
            print('setting', self.devices[i].id, '->', x[i])
            self.devices[i].set_value(x[i])

    def error_func(self, x):

        self.minimizer.kill = self.kill
        if self.kill:
            print('Killed from external process')
            # NEW CODE - to kill if run from outside thread
            return
        # check limits
        if self.exceed_limits(x):
            return self.target.pen_max
        # set values
        self.set_values(x)

        print('sleeping ' + str(self.timeout))
        sleep(self.timeout)
        print ('done sleeping')
        pen = self.target.get_penalty()
        print('penalty:', pen)
        if self.debug: print('penalty:', pen)
        self.x_data.append(x)
        self.y_data.append(pen)
        return pen

    def max_target_func(self, target, devices, params = {}):
        """
        Direct target function optimization with simplex/GP, using Devices as a multiknob
        """
        self.target = target
        self.devices = devices
        # testing
        self.target.devices = self.devices

        dev_ids = [dev.eid for dev in self.devices]
        if self.debug: print('starting multiknob optimization, devices = ', dev_ids)

        target_ref = self.target.get_penalty()

        x = [dev.get_value() for dev in self.devices] # self.mi.init_corrector_vals(devices)
        x_init = x

        if self.logging:
            self.logger.log_start(dev_ids, method=self.minimizer.__class__.__name__, x_init=x_init, target_ref=target_ref)

        self.minimizer.minimize(self.error_func, x)
        # set best solution
        x = self.x_data[np.argmin(self.y_data)]
        if self.exceed_limits(x):
            return self.target.pen_max
        self.set_values(x)

        target_new = self.target.get_penalty()

        print ('step ended changing sase from/to', target_ref, target_new)

        if self.logging:
            self.logger.log_fin(target_new=target_new)



class Action:
    def __init__(self, func, args = None, id = None):
        self.func = func
        self.args = args

        self.id = id

    def apply(self):
        print ('applying...', self.id)
        self.func(*self.args)
    #def to_JSON(self):
        #print "hoo"
    #def __repr__(self):
        #return json.dumps(self.__dict__)




def test_simplex():
    """
    test simplex method
    :return:
    """
    d1 = TestDevice(eid="d1")
    d2 = TestDevice(eid="d2")
    d3 = TestDevice(eid="d3")
    target = TestTarget()

    opt = Optimizer()
    opt.timeout = 0
    minimizer = Simplex()
    minimizer.max_iter = 300
    opt.minimizer = minimizer
    #opt.debug = True

    seq = [Action(func=opt.max_target_func, args=[ target, [d1, d2, d3]])]
    opt.eval(seq)



from itertools import chain
import scipy
from ocelot.optimizer.GP.OnlineGP import OGP
from ocelot.optimizer.GP.bayes_optimization import BayesOpt, HyperParams
import pandas as pd

def test_GP():
    """
    test GP method
    :return:
    """

    d1 = TestDevice(eid="d1")
    d2 = TestDevice(eid="d2")
    d3 = TestDevice(eid="d3")

    devices = [d1, d2]
    target = TestTarget()

    opt = Optimizer()
    opt.timeout = 0

    opt_smx = Optimizer()
    opt_smx.timeout = 0
    minimizer = Simplex()
    minimizer.max_iter = 3
    opt_smx.minimizer = minimizer
    #opt.debug = True

    seq = [Action(func=opt_smx.max_target_func, args=[target, devices])]
    opt_smx.eval(seq)
    s_data = np.append(np.vstack(opt_smx.x_data), np.transpose(-np.array([opt_smx.y_data])), axis=1)
    print(s_data)

    # -------------- GP config setup -------------- #
    #GP parameters
    numBV = 30
    xi = 0.01
    #no input bounds on GP selection for now

    pvs = [dev.eid for dev in devices]
    hyp_params = HyperParams(pvs=pvs, filename="../parameters/hyperparameters.npy")
    ave = np.mean(-np.array(opt_smx.y_data))
    std = np.std(-np.array(opt_smx.y_data))
    noise = hyp_params.calcNoiseHP(ave, std=0.)
    coeff = hyp_params.calcAmpCoeffHP(ave, std=0.)
    len_sc_hyps = []
    for dev in devices:
        ave = 10
        std = 3
        len_sc_hyps.append(hyp_params.calcLengthScaleHP(ave, std))
    print(len_sc_hyps )

    bnds = None
    #hyps = hyp_params.loadHyperParams(energy=3, detector_stat_params=target.get_stat_params())
    hyps1 = (np.array([len_sc_hyps]), coeff, noise) #(np.array([hyps]), coeff, noise)

    #init model
    dim = len(pvs)

    model = OGP(dim, hyps1, maxBV=numBV, weighted=False)

    minimizer = BayesOpt(model, target_func=target, xi=0.01, acq_func='EI', bounds=bnds, prior_data=pd.DataFrame(s_data))
    minimizer.devices = devices
    minimizer.max_iter = 300
    opt.minimizer = minimizer

    seq = [Action(func=opt.max_target_func, args=[ target, devices])]
    opt.eval(seq)


if __name__ == "__main__":
    test_simplex()
    #test_GP()

