'''
machine interface
includes online optimizer, response measurement and other stuff
'''

from pylab import *
from scipy.optimize import *
import scipy.optimize as opt
from time import sleep, time
import json

"""
# Any machine interface class should implement following methods

class MachineInterface:
    def __init__(self):
        pass
    def get_alarms(self):
        '''
        return an array of values for all active bpms
        values range between (0,1) 1 corresponds to alarm level
        '''
        pass
    def get_sase(self, detector='default'):
        '''
        return SASE pulse energy
        units and averaging can be arbitrary        '''

        pass
    def get_value(self, device_name):
        pass
    def set_value(self, device_name, val):
        pass
    def init_corrector_vals(self, correctors):
        pass
    

# placeholder for magnet field ranges etc.

class DeviceProperties:
    def __init__(self):
        self.limits = [-1000, 1000]
"""

class Optimizer:
    def __init__(self, mi, dp, save_opt_params=None):
        self.debug = False
        self.mi = mi
        self.dp = dp
        self.timeout = 1.0
        self.logging = False
        self.log_file = "test.log"
        self.detector = 'gmd_default'
        self.wasSaved = False
        self.isRunning = False
        self.niter = 0
        self.maxiter = None
        self.sop = save_opt_params
        self.penalty = 0.


    def eval(self, seq, logging = False, log_file = None):
        self.isRunning = True
        self.wasSaved = False

        for s in seq:
            if self.sop != None: self.save_action(s.args, flag="start")
            self.wasSaved = False

            s.apply()
            if not self.wasSaved and self.sop != None:
                self.save_action(s.args, flag="stop")
        self.isRunning = False

    def save_action(self, args, flag):
        print("SAVE MACHINE",  flag)
        if self.sop is not None:
            self.wasSaved = self.sop.save(args, time=time(), niter=self.niter, flag=flag)
        self.niter = 0

    def run(self, sequence):
        if self.isRunning:
            print("Optimization is still running")
        else:
            self.eval(sequence)


    def max_sase(self, devices, method = 'simplex', params = {}, opt_pointing = False):
        '''
        direct sase optimization with simplex, using correctors as a multiknob
        '''

        if self.debug: print('starting multiknob optimization, correctors = ', devices)

        def error_func(x, x_init, tols):
            if self.debug: print("X_relative = ", x, x_init, tols)
            x = x_init + (x-1)*tols/(10.*0.05) # relative to absolute

            if self.debug: print("X_absolute = ", x, x_init)
            if self.debug: print("isRunning:", self.isRunning)
            self.niter += 1

            print("number of func eval: ", self.niter)

            if not self.isRunning:
                print("save machine parameters and kill optimizer")
                if self.sop != None: self.save_action([devices, method, params], flag="force stop")
                sleep(1) # in order to give time to save parameters before next evaluation of error_func

            pen_max = 100.0
    
            for i in range(len(x)):
                if self.debug: print('{0} x[{1}]={2}'.format(devices[i], i, x[i]))
                limits = self.dp.get_limits(devices[i])
                if self.debug: print('limits=[{0}, {1}]'.format(limits[0], limits[1]))
                if x[i] < limits[0] or x[i] > limits[1]:
                    print('limits exceeded')
                    return pen_max

            for i in range(len(devices)):
                print ('setting', devices[i], '->',x[i])
                self.mi.set_value(devices[i], x[i])

            sleep(self.timeout)

            sase = self.mi.get_sase(detector=self.detector)
            alarm = np.max(self.mi.get_alarms())

            if self.debug: print ('alarm:', alarm)
            if self.debug: print ('sase:', sase)

            pen = 0.0
            if alarm > 1.0:
                return pen_max
            if alarm > 0.7:
                return alarm * 50.0

            pen += alarm
            pen -= sase

            if self.debug: print ('penalty:', pen)
            self.penalty = pen
            return pen

        sase_ref = self.mi.get_sase(detector=self.detector)
    
        x = np.array([self.mi.get_value(dev) for dev in devices])
        tols = np.zeros(len(devices))
        for i, dev in enumerate(devices):
            limits = self.dp.get_limits(dev)
            tols[i] = (limits[1] - limits[0])/2.
        
        x_init = x
        
        if self.logging: 
            f = open(self.log_file,'a')
            f.write('\n*** optimization step ***\n')
            f.write(str(devices) + '\n')
            f.write(method + '\n')
            f.write('x0=' + str(x_init) + '\n')
            f.write('sase0=' + str(sase_ref) + '\n')

        try:
            max_iter = params['maxiter']
        except KeyError:
            max_iter = 10 * len(x)

        if max_iter == None:
            max_iter = 10 * len(x)

        try:
            xtol = params['xtol']
        except KeyError:
            xtol = 1.e-3

        self.maxiter = max_iter
        if method == 'simplex':
            print('using simplex optimizer, params:', params)
            # opt.fmin(error_func,x,xtol=xtol, maxiter=max_iter)
            opt.fmin(error_func, np.ones(len(x)), args=(x_init, tols), xtol=xtol, maxfun=max_iter)
            # opt.fmin(error_func, x, args=(x_init,tols), xtol=xtol, maxfun=max_iter)
        
        if method == 'cg':
            print ('using CG optimizer, params:', params )

            try:
                epsilon = params['epsilon']
            except KeyError:
                epsilon = 0.1

            opt.fmin_cg(error_func, x/x_init, args =(x_init,), gtol=xtol, epsilon = epsilon, maxiter=max_iter)

        if method == 'powell': 
            print ('using powell optimizer, params:', params)
            opt.fmin_powell(error_func,x/x_init, args=(x_init, ), xtol=xtol, maxiter=max_iter)


        if method == 'fancy_stuff_from': 
            print ('using fancy optimizer, params:', params)
            pass

        sase_new = self.mi.get_sase(detector=self.detector)
        #self.save_machine_set()

        #print ('step ended changing sase from/to', sase_ref, sase_new)
        #if sase_new <= sase_ref:
        #    for i in range(len(devices)):
        #        print ('reverting', devices[i], '->',x_init[i])
        #        self.mi.set_value(devices[i], x_init[i])

        if self.logging:
             f.write('sase_new=' + str(sase_new) + '\n')
             f.close()

    def max_sase_bump(self, bump, alpha, method = 'simplex', params = {}, opt_pointing = False):
        '''
        direct sase optimization with simplex, using correctors as a multiknob
        '''
        if self.debug: print ('starting multiknob optimization, correctors = ', bump)

        if opt_pointing:
            weight_gmd_bpm_1 = 10.0
            weight_gmd_bpm_2 = 10.0
        else:
            weight_gmd_bpm_1 = 0.0
            weight_gmd_bpm_2 = 0.0

        def error_func(x):

            print (self.dp)

            pen_max = 100.0

            if abs(x) >1:
                return pen_max

            dI = bump["dI"]
            currents = bump["currents"]
            #print 'error_func: ', bpm_names, '->',  planes
            correctors_ = bump["correctors"]
            for i in range(len(correctors_)):
                print ("alpha = ", x)
                print ('{0} x[{1}]={2}'.format(correctors_[i], i, currents[i] + dI[i]*x))
                limits = self.dp.get_limits(correctors_[i])
                print  ('limits=[{0}, {1}]'.format(limits[0], limits[1]))
                if currents[i] + dI[i]*x < limits[0] or currents[i] + dI[i]*x > limits[1]:
                    print ('limits exceeded')
                    return pen_max


            for i in range(len(correctors_)):
                print ('setting', correctors_[i], '->', currents[i] + dI[i]*x)
                self.mi.set_value(correctors_[i], currents[i] + dI[i]*x)

            sleep(self.timeout)

            sase = self.mi.get_sase(detector=self.detector)
            alarm = np.max(self.mi.get_alarms())
            #z1, z2 = get_sase_pos()

            if self.debug: print ('alarm:', alarm)
            if self.debug: print ('sase:', sase)
            #print 'pointing', z1, z2, 'weights', weight_gmd_bpm_1, weight_gmd_bpm_2

            pen = 0.0

            if alarm > 1.0:
                return pen_max
            if alarm > 0.7:
                return alarm * 50.0
            pen += alarm

            pen -= sase

            if self.debug: print ('penalty:', pen)

            return pen


        sase_ref = self.mi.get_sase(detector=self.detector)

        x = alpha
        x_init = x

        if self.logging:
            f = open(self.log_file,'a')
            f.write('\n*** optimization step ***\n')
            f.write(str(bump["correctors"]) + '\n')
            f.write(method + '\n')
            f.write('x0=' + str(x_init) + '\n')
            f.write('sase0=' + str(sase_ref) + '\n')


        if method == 'cg':
            print ('using CG optimizer, params:', params)

            try:
                max_iter = params['maxiter']
            except KeyError:
                max_iter = 10 * len(x)

            try:
                epsilon = params['epsilon']
            except KeyError:
                epsilon = 0.1

            try:
                gtol = params['gtol']
            except KeyError:
                gtol = 1.e-3

            opt.fmin_cg(error_func,x,gtol=gtol, epsilon = epsilon, maxiter=max_iter)

        if method == 'simplex':
            print ('using simplex optimizer, params:', params)

            try:
                max_iter = params['maxiter']
            except KeyError:
                max_iter = 10 * len(bump["correctors"])

            try:
                xtol = params['xtol']
            except KeyError:
                xtol = 1.e-3

            opt.fmin(error_func,x,xtol=xtol, maxiter=max_iter)

        if method == 'powell':
            print ('using powell optimizer, params:', params)

            try:
                max_iter = params['maxiter']
            except KeyError:
                max_iter = 10 * len(x)

            try:
                xtol = params['xtol']
            except KeyError:
                xtol = 1.e-3

            opt.fmin_powell(error_func,x,xtol=xtol, maxiter=max_iter)


        if method == 'fancy_stuff_from':
            print ('using fancy optimizer, params:', params)
            pass

        sase_new = self.mi.get_sase(detector=self.detector)

        print ('step ended changing sase from/to', sase_ref, sase_new)
        if sase_new <= sase_ref:
            for i in range(len(bump["correctors"])):
                print ('reverting', bump["correctors"][i], '->',bump["currents"][i])
                self.mi.set_value(bump["correctors"][i], bump["currents"][i])

        if self.logging:
             f.write('sase_new=' + str(sase_new) + '\n')
             f.close()


    def min_orbit(self, orbit, method = 'simplex', params = {}, opt_pointing = False):
        '''
        direct sase optimization with simplex, using correctors as a multiknob
        '''
        correctors = orbit["correctors"]
        bpms = orbit["bpms"]

        if self.debug: print ('starting multiknob optimization, correctors = ', correctors)

        def get_rms(bpms):
            #print
            #print "bpms= ",  bpms.keys()
            n = len(bpms.keys())
            x = 0.
            y = 0.
            for bpm in bpms.keys():
                X, Y = self.mi.get_bpms_xy([bpm])
                x += (X[0] - bpms[bpm]["x"])**2
                y += (Y[0] - bpms[bpm]["y"])**2
            print ("rms = ", np.sqrt(x/n)*1000. + np.sqrt(y/n)*1000.)
            return np.sqrt(x/n)*1000. + np.sqrt(y/n)*1000.

        def error_func(x):

            #print self.dp

            pen_max = 100.0

            #print 'error_func: ', bpm_names, '->',  planes

            for i in range(len(x)):
                print ('{0} x[{1}]={2}'.format(correctors[i], i, x[i]))
                limits = self.dp.get_limits(correctors[i])
                print  ('limits=[{0}, {1}]'.format(limits[0], limits[1]))
                if x[i] < limits[0] or x[i] > limits[1]:
                    print ('limits exceeded')
                    return pen_max


            for i in range(len(correctors)):
                print ('setting', correctors[i], '->',x[i])
                self.mi.set_value(correctors[i], x[i])

            sleep(self.timeout)


            orb_min = get_rms(bpms)
            alarm = np.max(self.mi.get_alarms())
            #z1, z2 = get_sase_pos()

            if self.debug: print ('alarm:', alarm)
            if self.debug: print ('sase:', orb_min)
            #print 'pointing', z1, z2, 'weights', weight_gmd_bpm_1, weight_gmd_bpm_2

            pen = 0.0

            if alarm > 1.0:
                return pen_max
            if alarm > 0.7:
                return alarm * 50.0
            pen += alarm

            pen += orb_min

            if self.debug: print ('penalty:', pen)

            return pen


        orb_ref = get_rms(bpms)

        x = self.mi.init_corrector_vals(correctors)
        print ("initial X", x)
        x_init = x

        if self.logging:
            f = open(self.log_file,'a')
            f.write('\n*** optimization step ***\n')
            f.write(str(correctors) + '\n')
            f.write(method + '\n')
            f.write('x0=' + str(x_init) + '\n')
            f.write('sase0=' + str(orb_ref) + '\n')




        if method == 'simplex':
            print ('using simplex optimizer, params:', params)

            try:
                max_iter = params['maxiter']
            except KeyError:
                max_iter = 10 * len(x)

            try:
                xtol = params['xtol']
            except KeyError:
                xtol = 1.e-3

            opt.fmin(error_func,x,xtol=xtol, maxiter=max_iter)



        orb_new = get_rms(bpms)

        print ('step ended changing sase from/to', orb_ref, orb_new)

        if self.logging:
             f.write('sase_new=' + str(orb_new) + '\n')
             f.close()

class Action:
    def __init__(self, func, args = None, id = None):
        self.func = func
        self.args = args
        self.id = id
    def apply(self):
        print ('applying...', self.id)
        n_iter = self.func(*self.args)
        return n_iter
    def to_JSON(self):
        print ("hoo")
    def __repr__(self):
        return json.dumps(self.__dict__)


'''
test interface
'''
"""
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
        sleep(0.1)
        return 0

    def get_nbunches(self):
        return 0
    def get_value_ps(self, device_name):
        return 0.

"""
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
