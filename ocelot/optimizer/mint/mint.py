"""
Main Ocelot optimization file
Contains the setup for using the scipy.optimize package run simplex and other algorothms
Modified for use at LCLS from Ilya's version

The file was modified and were introduced new Objects and methods.
S. Tomin, 2017

"""
from __future__ import print_function, absolute_import
from time import sleep
from scipy.optimize import OptimizeResult
import scipy
import numpy as np
from ocelot.optimizer.mint.opt_objects import *
from scipy import optimize
from ocelot.optimizer.GP.bayes_optimization import *
from ocelot.optimizer.GP.OnlineGP import OGP
import pandas as pd
from threading import Thread
import sklearn
sklearn_version = sklearn.__version__
if sklearn_version >= "0.18":
    from ocelot.optimizer.GP import gaussian_process as gp_sklearn

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
        self.maximize = False

    def minimize(self, error_func, x):
        pass


class Simplex(Minimizer):
    def __init__(self):
        super(Simplex, self).__init__()
        self.xtol = 1e-5
        self.dev_steps = None

    def minimize(self,  error_func, x):
        #print("start seed", np.count_nonzero(self.dev_steps))
        if self.dev_steps == None or len(self.dev_steps) != len(x):
            print("initial simplex is None")
            isim = None
        elif np.count_nonzero(self.dev_steps) != len(x):
            print("There is zero step. Initial simplex is None")
            isim = None
        else:
            #step = np.ones(len(x))*0.05
            isim = np.zeros((len(x) + 1, len(x)))
            isim[0, :] = x
            for i in range(len(x)):
                vertex = np.zeros(len(x))
                vertex[i] = self.dev_steps[i]
                isim[i + 1, :] = x + vertex
            print("ISIM = ", isim)
        #res = optimize.minimize(error_func, x, method='Nelder-Mead',  tol=self.xtol,
        #                        options = {'disp': False, 'initial_simplex': [0.05, 0.05], 'maxiter': self.max_iter})
        if scipy.__version__ < "0.18":
            res = optimize.fmin(error_func, x, maxiter=self.max_iter, maxfun=self.max_iter, xtol=self.xtol)
        else:
            res = optimize.fmin(error_func, x, maxiter=self.max_iter, maxfun=self.max_iter, xtol=self.xtol, initial_simplex=isim)

        #print("finish seed")
        return res


class GaussProcess(Minimizer):
    def __init__(self):
        super(GaussProcess, self).__init__()
        self.seed_iter = 5
        self.seed_timeout = 0.1
        #BayesOpt.__init__(self, model=model, target_func=target_func, acq_func=acq_func, xi=xi, alt_param=alt_param,
        #                  m=m, bounds=bounds, iter_bound=iter_bound, prior_data=prior_data)
        self.target = None
        self.devices = []
        self.energy = 3
        #GP parameters
        self.numBV = 30
        self.xi = 0.01
        self.bounds = None
        self.acq_func = 'EI'
        self.alt_param = -1
        self.m = 200
        self.iter_bound = False
        self.hyper_file = "../parameters/hyperparameters.npy"
        self.max_iter = 50
        self.norm_coef = 0.1

    def seed_simplex(self):
        opt_smx = Optimizer()
        opt_smx.normalization = True
        opt_smx.norm_coef = self.norm_coef
        opt_smx.timeout = self.seed_timeout
        minimizer = Simplex()
        minimizer.max_iter = self.seed_iter
        opt_smx.minimizer = minimizer
        # opt.debug = True
        seq = [Action(func=opt_smx.max_target_func, args=[self.target, self.devices])]
        opt_smx.eval(seq)

        seed_data = np.append(np.vstack(opt_smx.opt_ctrl.dev_sets), np.transpose(-np.array([opt_smx.opt_ctrl.penalty])), axis=1)
        self.prior_data = pd.DataFrame(seed_data)
        self.seed_y_data = opt_smx.opt_ctrl.penalty


    def preprocess(self):
        # -------------- GP config setup -------------- #

        #no input bounds on GP selection for now
        #print(self.devices)
        pvs = [dev.eid for dev in self.devices]
        #print(pvs)
        hyp_params = HyperParams(pvs=pvs, filename=self.hyper_file)
        ave = np.mean(-np.array(self.seed_y_data))
        std = np.std(-np.array(self.seed_y_data))
        noise = hyp_params.calcNoiseHP(ave, std=0.)
        coeff = hyp_params.calcAmpCoeffHP(ave, std=0.)

        len_sc_hyps = []
        for dev in self.devices:
            ave = 10
            std = 3
            len_sc_hyps.append(hyp_params.calcLengthScaleHP(ave, std))

        #hyps = hyp_params.loadHyperParams(energy=self.energy, detector_stat_params=target.get_stat_params())
        hyps1 = (np.array([len_sc_hyps]), coeff, noise)
        #(np.array([hyp_params.calcLengthScaleHP(ave, std)]), coeff=hyp_params.calcAmpCoeffHP(ave, std=0.), noise=hyp_params.calcAmpCoeffHP(ave, std=0.))
        #print("hyps1", hyps1)
        #init model
        dim = len(pvs)

        self.model = OGP(dim, hyps1, maxBV=self.numBV, weighted=False)
        self.scanner = BayesOpt(model=self.model, target_func=self.target, acq_func=self.acq_func, xi=self.xi,
                           alt_param=self.alt_param, m=self.m, bounds=self.bounds, iter_bound=self.iter_bound,
                                prior_data=self.prior_data)

        self.scanner.max_iter = self.max_iter

    def minimize(self,  error_func, x):
        #self.target_func = error_func

        self.seed_simplex()
        self.preprocess()
        x = [dev.get_value() for dev in self.devices]
        print("start GP")
        self.scanner.minimize(error_func, x)
        print("finish GP")
        return


class GaussProcessSKLearn(Minimizer):
    def __init__(self):
        super(GaussProcessSKLearn, self).__init__()
        self.seed_iter = 5
        self.seed_timeout = 0.1

        self.target = None
        self.devices = []

        self.x_obs = []
        self.y_obs = []
        #GP parameters

        self.max_iter = 50
        self.norm_coef = 0.1
        self.kill = False
        self.opt_ctrl = None

    def seed_simplex(self):
        opt_smx = Optimizer()
        opt_smx.normalization = True
        opt_smx.norm_coef = self.norm_coef
        opt_smx.timeout = self.seed_timeout
        opt_smx.opt_ctrl = self.opt_ctrl
        minimizer = Simplex()
        minimizer.max_iter = self.seed_iter
        opt_smx.minimizer = minimizer
        # opt.debug = True
        seq = [Action(func=opt_smx.max_target_func, args=[self.target, self.devices])]
        opt_smx.eval(seq)
        print(opt_smx.opt_ctrl.dev_sets)
        self.x_obs = np.vstack(opt_smx.opt_ctrl.dev_sets)
        self.y_obs = np.array(opt_smx.opt_ctrl.penalty)
        self.y_sigma_obs = np.zeros(len(self.y_obs))

    def load_seed(self, x_sets, penalty, sigma_pen=None):
        self.x_obs = np.vstack(x_sets)
        self.y_obs = np.array(penalty)
        if sigma_pen == None:
            self.y_sigma_obs = np.zeros(len(self.y_obs))
        else:
            self.y_sigma_obs = sigma_pen

    def preprocess(self):

        self.scanner = gp_sklearn.GP()
        self.scanner.opt_ctrl = self.opt_ctrl
        devs_std = []
        devs_search_area = []
        for dev in self.devices:
            lims = dev.get_limits()
            devs_std.append((lims[-1] - lims[0])/3.)
            x_vec = np.atleast_2d(np.linspace(lims[0], lims[-1], num=50)).T
            devs_search_area.append(x_vec)

        self.scanner.x_search = np.hstack(devs_search_area)
        self.scanner.x_obs = self.x_obs
        self.scanner.y_obs = self.y_obs
        self.scanner.y_sigma_obs = self.y_sigma_obs

        self.scanner.ck_const_value = (0.5*np.mean(self.scanner.y_obs))**2 + 0.1
        #self.scanner.ck_const_value_bounds = (self.scanner.ck_const_value,self.scanner.ck_const_value)
        self.scanner.rbf_length_scale = np.array(devs_std)/2. + 0.01
        #self.scanner.rbf_length_scale_bounds = (self.scanner.rbf_length_scale, self.scanner.rbf_length_scale)
        self.scanner.max_iter = self.max_iter

    def minimize(self,  error_func, x):
        #self.target_func = error_func

        self.seed_simplex()
        if self.opt_ctrl.kill:
            return
        self.preprocess()
        x = [dev.get_value() for dev in self.devices]
        print("start GP")
        self.scanner.minimize(error_func, x)
        print("finish GP")
        return


class CustomMinimizer(Minimizer):
    def __init__(self):
        super(CustomMinimizer, self).__init__()
        self.dev_steps = [0.05]

    def minimize(self,  error_func, x):
        def custmin(fun, x0, args=(), maxfev=None, stepsize=[0.1],
                    maxiter=self.max_iter, callback=None, **options):

            print("inside ", stepsize)

            if np.size(stepsize) != np.size(x0):
                stepsize = np.ones(np.size(x0))*stepsize[0]
            print("inside ", stepsize)
            bestx = x0
            besty = fun(x0)
            print("BEST", bestx, besty)
            funcalls = 1
            niter = 0
            improved = True
            stop = False

            while improved and not stop and niter < maxiter:
                improved = False
                niter += 1
                for dim in range(np.size(x0)):
                    for s in [bestx[dim] - stepsize[dim], bestx[dim] + stepsize[dim]]:
                        print("custom", niter, dim, s)
                        testx = np.copy(bestx)
                        testx[dim] = s
                        testy = fun(testx, *args)
                        funcalls += 1
                        if testy < besty:
                            besty = testy
                            bestx = testx
                            improved = True
                    if callback is not None:
                        callback(bestx)
                    if maxfev is not None and funcalls >= maxfev:
                        stop = True
                        break

            return OptimizeResult(fun=besty, x=bestx, nit=niter,
                                  nfev=funcalls, success=(niter > 1))
        res = optimize.minimize(error_func, x, method=custmin, options=dict(stepsize=self.dev_steps))
        return res


class MachineStatus:
    def __init__(self):
        pass

    def is_ok(self):
        return True


class OptControl:
    def __init__(self):
        self.penalty = []
        self.dev_sets = []
        self.devices = []
        self.nsteps = 0
        self.m_status = MachineStatus()
        self.pause = False
        self.kill = False
        self.is_ok = True
        self.timeout = 0.1
        self.alarm_timeout = 0

    def wait(self):
        while 1:
            if self.m_status.is_ok():
                self.is_ok = True
                time.sleep(self.alarm_timeout)
                return 1
            time.sleep(self.timeout)
            self.is_ok = False
            print(".",)

    def stop(self):
        self.kill = True

    def start(self):
        self.kill = False

    def back_nsteps(self, n):
        if n <= self.nsteps:
            n = -1 - n
        else:
            print("OptControl: back_nsteps n > nsteps. return last step")
            n = -1
        return self.dev_sets[-n]


    def save_step(self, pen, x):
        self.penalty.append(pen)
        self.dev_sets.append(x)
        self.nsteps = len(self.penalty)

    def best_step(self):
        #if len(self.penalty)== 0:
        #    print("No ")
        x = self.dev_sets[np.argmin(self.penalty)]
        return x


class Optimizer(Thread):
    def __init__(self, normalize=False):
        super(Optimizer, self).__init__()
        self.debug = False
        self.minimizer = Simplex()
        self.logging = False
        #self.kill = False #intructed by tmc to terminate thread of this class
        self.log_file = "log.txt"
        self.logger = Logger(self.log_file)
        self.devices = []
        self.target = None
        self.timeout = 0
        self.opt_ctrl = OptControl()
        self.seq = []
        self.set_best_solution = True
        self.normalization = False
        self.norm_coef = 0.05

    def eval(self, seq=None, logging=False, log_file=None):
        """
        Run the sequence of tuning events
        """
        if seq != None:
            self.seq = seq
        for s in self.seq:
            s.apply()

    def exceed_limits(self, x):
        for i in range(len(x)):
            if self.devices[i].check_limits(x[i]):
                return True
        return False

    def set_values(self, x):
        for i in range(len(self.devices)):
            print('setting', self.devices[i].id, '->', x[i])
            self.devices[i].set_value(x[i])

    def calc_scales(self):
        """
        calculate scales for normalized simplex

        :return: np.array() - device_delta_limits * norm_coef
        """
        self.norm_scales = np.zeros(np.size(self.devices))
        for i, dev in enumerate(self.devices):
            lims = dev.get_limits()
            delta = lims[-1] - lims[0]
            self.norm_scales[i] = delta*self.norm_coef

        return self.norm_scales

    def error_func(self, x):
        if self.normalization:
            delta_x = x
            delta_x_scaled = delta_x/0.00025*self.norm_scales
            x = self.x_init + delta_x_scaled
            print("delta_x = ", delta_x, "delta_x_scaled = ", delta_x_scaled)

        if self.opt_ctrl.kill:
            #self.minimizer.kill = self.opt_ctrl.kill
            print('Killed from external process')
            # NEW CODE - to kill if run from outside thread
            return self.target.pen_max

        self.opt_ctrl.wait()

        # check limits
        if self.exceed_limits(x):
            return self.target.pen_max
        # set values
        self.set_values(x)

        print('sleeping ' + str(self.timeout))
        sleep(self.timeout)
        #print ('done sleeping')

        pen = self.target.get_penalty()
        print('penalty:', pen)
        if self.debug: print('penalty:', pen)

        self.opt_ctrl.save_step(pen, x)
        return pen

    def max_target_func(self, target, devices, params = {}):
        """
        Direct target function optimization with simplex/GP, using Devices as a multiknob
        """
        [dev.clean() for dev in devices]
        target.clean()
        self.target = target
        #print(self.target)
        self.devices = devices
        # testing
        self.minimizer.devices = devices
        self.minimizer.target = target
        self.minimizer.opt_ctrl = self.opt_ctrl
        self.target.devices = self.devices
        dev_ids = [dev.eid for dev in self.devices]
        if self.debug: print('starting multiknob optimization, devices = ', dev_ids)

        target_ref = self.target.get_penalty()

        x = [dev.get_value() for dev in self.devices]
        x_init = x

        if self.logging:
            self.logger.log_start(dev_ids, method=self.minimizer.__class__.__name__, x_init=x_init, target_ref=target_ref)

        if self.normalization:
            self.x_init = x_init
            x = np.zeros_like(x)
            self.calc_scales()

        res = self.minimizer.minimize(self.error_func, x)
        print("result", res)

        # set best solution
        if self.set_best_solution:
            print("SET the best solution", x)
            x = self.opt_ctrl.best_step()
            if self.exceed_limits(x):
                return self.target.pen_max
            self.set_values(x)

        target_new = self.target.get_penalty()

        print ('step ended changing sase from/to', target_ref, target_new)

        if self.logging:
            self.logger.log_fin(target_new=target_new)


    def run(self):
        self.opt_ctrl.start()
        self.eval(self.seq)
        print("FINISHED")
        #self.opt_ctrl.stop()
        return 0


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

    def get_limits():
        return [-100, 100]

    d1.get_limits = get_limits
    d2.get_limits = get_limits
    d3.get_limits = get_limits

    devices = [d1, d2, d3]
    target = Target_test()

    # init Optimizer
    opt = Optimizer()
    opt.timeout = 0

    # init Minimizer
    minimizer = Simplex()
    minimizer.max_iter = 300

    opt.minimizer = minimizer

    seq = [Action(func=opt.max_target_func, args=[target, devices])]
    opt.eval(seq)


def test_gauss_process():
    """
    test simplex method
    :return:
    """
    d1 = TestDevice(eid="d1")
    d2 = TestDevice(eid="d2")
    d3 = TestDevice(eid="d3")

    def get_limits():
        return [-100, 100]

    d1.get_limits = get_limits
    d2.get_limits = get_limits
    d3.get_limits = get_limits

    devices = [d1, d2, d3]
    target = TestTarget()

    # init Optimizer
    opt = Optimizer()
    opt.timeout = 0

    # init Minimizer
    minimizer = GaussProcess()
    minimizer.seed_iter = 3
    minimizer.max_iter = 300

    opt.minimizer = minimizer

    seq = [Action(func=opt.max_target_func, args=[ target, devices])]
    opt.eval(seq)


#from itertools import chain
#import scipy
#from ocelot.optimizer.GP.OnlineGP import OGP
#from ocelot.optimizer.GP.bayes_optimization import BayesOpt, HyperParams


def test_GP():
    """
    test GP method
    :return:
    """

    def get_limits():
        return [-100, 100]
    d1 = TestDevice(eid="d1")
    d1.get_limits = get_limits
    d2 = TestDevice(eid="d2")
    d2.get_limits = get_limits
    d3 = TestDevice(eid="d3")
    d3.get_limits = get_limits

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
    s_data = np.append(np.vstack(opt_smx.opt_ctrl.dev_sets), np.transpose(-np.array([opt_smx.opt_ctrl.penalty])), axis=1)
    print(s_data)

    # -------------- GP config setup -------------- #
    #GP parameters
    numBV = 30
    xi = 0.01
    #no input bounds on GP selection for now

    pvs = [dev.eid for dev in devices]
    hyp_params = HyperParams(pvs=pvs, filename="../parameters/hyperparameters.npy")
    ave = np.mean(-np.array(opt_smx.opt_ctrl.penalty))
    std = np.std(-np.array(opt_smx.opt_ctrl.penalty))
    noise = hyp_params.calcNoiseHP(ave, std=0.)
    coeff = hyp_params.calcAmpCoeffHP(ave, std=0.)
    len_sc_hyps = []
    for dev in devices:
        ave = 10
        std = 3
        len_sc_hyps.append(hyp_params.calcLengthScaleHP(ave, std))
    print("y_data", opt_smx.opt_ctrl.penalty)
    print("pd.DataFrame(s_data)", pd.DataFrame(s_data))
    print("len_sc_hyps", len_sc_hyps )

    bnds = None
    #hyps = hyp_params.loadHyperParams(energy=3, detector_stat_params=target.get_stat_params())
    hyps1 = (np.array([len_sc_hyps]), coeff, noise) #(np.array([hyps]), coeff, noise)
    print("hyps1", hyps1)
    #exit(0)
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
    #test_gauss_process()
    #test_GP()

