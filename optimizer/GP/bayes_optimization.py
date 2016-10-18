"""
Contains the Bayes optimization class.
Initialization parameters:
    model: an object with methods 'predict', 'fit', and 'update'
    interface: an object which supplies the state of the system and
        allows for changing the system's x-value.
        Should have methods '(x,y) = intfc.getState()' and 'intfc.setX(x_new)'.
        Note that this interface system is rough, and used for testing and
            as a placeholder for the machine interface.
    acq_func: specifies how the optimizer should choose its next point.
        'EI': uses expected improvement. The interface should supply y-values.
        'testEI': uses EI over a finite set of points. This set must be
            provided as alt_param, and the interface need not supply
            meaningful y-values.
    xi: exploration parameter suggested in some Bayesian opt. literature
    alt_param: currently only used when acq_func=='testEI'
    m: the maximum size of model; can be ignored unless passing an untrained
        SPGP or other model which doesn't already know its own size
    bounds: a tuple of (min,max) tuples specifying search bounds for each
        input dimension. Generally leads to better performance.
        Has a different interpretation when iter_bounds is True.
    iter_bounds: if True, bounds the distance that can be moved in a single
        iteration in terms of the length scale in each dimension. Uses the
        bounds variable as a multiple of the length scales, so bounds==2
        with iter_bounds==True limits movement per iteration to two length
        scales in each dimension. Generally a good idea for safety, etc.
    prior_data: input data to train the model on initially. For convenience,
        since the model can be trained externally as well.
        Assumed to be a pandas DataFrame of shape (n, dim+1) where the last
            column contains y-values.

Methods:
    acquire(): Returns the point that maximizes the acquisition function.
        For 'testEI', returns the index of the point instead.
        For normal acquisition, currently uses the bounded L-BFGS optimizer.
            Haven't tested alternatives much.
    best_seen(): Uses the model to make predictions at every observed point,
        returning the best-performing (x,y) pair. This is more robust to noise
        than returning the best observation, but could be replaced by other,
        faster methods.
    OptIter(): The main method for Bayesian optimization. Maximizes the
        acquisition function, then uses the interface to test this point and
        update the model.
"""

import operator as op
import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize
import time
from copy import deepcopy

class BayesOpt:
    def __init__(self, model, target_func, acq_func='EI', xi=0.0, alt_param=-1, m=200, bounds=None, iter_bound=False, prior_data=None):
        self.model = model
        self.m = m
        #self.delay = 0
        self.bounds = bounds
        self.iter_bound = iter_bound
        #self.interface = interface
        self.target_func = target_func
        #self.devices = devices
        self.acq_func = (acq_func, xi, alt_param)
        self.max_iter = 100
        self.check = None
        self.alpha = 1.

        #if(acq_func=='testEI'):
        #    (x_init, y_init) = np.array(alt_param.iloc[0, :-1],ndmin=2), alt_param.iloc[0, -1]
        #else:
        #    x_init = np.array([dev.get_value for dev in self.devices], ndmin=2)
        #    y_init = np.array([[target_func.get_penalty()]])
        #    #(x_init, y_init) = interface.getState()
        #
        #self.X_obs = np.array(x_init)
        #self.Y_obs = [y_init]
        #self.current_x = x_init

        # initialize model on prior data
        if(prior_data is not None):
            p_X = prior_data.iloc[:, :-1]
            p_Y = prior_data.iloc[:, -1]
            num = len(prior_data.index)
            self.model.fit(p_X, p_Y, min(m, num))


    def terminate(self, devices):
        """
        Sets the position back to the location that seems best in hindsight.
        It's a good idea to run this at the end of the optimization, since
        Bayesian optimization tries to explore and might not always end in
        a good place.
        """
        print("TERMINATE", self.x_best)
        if(self.acq_func[0] == 'EI'):
            # set position back to something reasonable
            for i, dev in enumerate(devices):
                dev.set_value(self.x_best[i])
            #error_func(self.x_best)
        if(self.acq_func[0] == 'UCB'):
            # UCB doesn't keep track of x_best, so find it
            (x_best, y_best) = self.best_seen()
            for i, dev in enumerate(devices):
                dev.set_value(x_best[i])


    def minimize(self, error_func, x):
        # weighting for exploration vs exploitation in the GP at the end of scan, alpha array goes from 1 to zero
        #alpha = [1.0 for i in range(40)]+[np.sqrt(50-i)/3.0 for i in range(41,51)]
        inverse_sign = -1
        self.current_x = np.array(np.array(x).flatten(), ndmin=2)
        self.X_obs = np.array(self.current_x)
        self.Y_obs = [np.array([[inverse_sign*error_func(x)]])]
        # iterate though the GP method
        for i in range(self.max_iter):
            # get next point to try using acquisition function
            x_next = self.acquire(self.alpha)

            #check for problems with the beam
            if self.check != None: self.check.errorCheck()

            y_new = error_func(x_next.flatten())
            if self.kill:
                #disable so user does not start another scan while the data is being saved
                break
            y_new = np.array([[inverse_sign *y_new]])

            #advance the optimizer to the next iteration
            #self.opt.OptIter(alpha=alpha[i])
            #self.OptIter() # no alpha



            # change position of interface and get resulting y-value


            x_new = deepcopy(x_next)
            #(x_new, y_new) = self.interface.getState()
            self.current_x = x_new

            # add new entry to observed data
            self.X_obs = np.concatenate((self.X_obs, x_new), axis=0)
            self.Y_obs.append(y_new)

            # update the model (may want to add noise if using testEI)
            self.model.update(x_new, y_new)# + .5*np.random.randn())


    def best_seen(self):
        """
        Checks the observed points to see which is predicted to be best.
        Probably safer than just returning the maximum observed, since the
        model has noise. It takes longer this way, though; you could
        instead take the model's prediction at the x-value that has
        done best if this needs to be faster.
        """
        (mu, var) = self.model.predict(self.X_obs)

        (ind_best, mu_best) = max(enumerate(mu), key=op.itemgetter(1))
        return (self.X_obs[ind_best], mu_best)

    def acquire(self, alpha=None):
        """
        Computes the next point for the optimizer to try by maximizing
        the acquisition function. If movement per iteration is bounded,
        starts search at current position.
        """
        if(self.acq_func[0] == 'EI'):
            (x_best, y_best) = self.best_seen()
            self.x_best = x_best
            x_start = x_best

            if(self.iter_bound):
                x_start = self.current_x
                if(self.bounds is None):
                    self.bounds = 1.0
                lengths = 1/np.sqrt(np.exp(self.model.covar_params[0]))
                iter_bounds = [(x_start[:,i] - self.bounds*lengths[:,i],x_start[:,i] + self.bounds*lengths[:,i]) for i in x_start.shape[1]]
            else:
                iter_bounds = self.bounds

            # maximize the EI (by minimizing negative EI)
            try:
                res = minimize(negExpImprove, x_start, args=(self.model, y_best, self.acq_func[1], alpha),
                                bounds=iter_bounds, method='L-BFGS-B', options={'maxfun':100})
            except:
                raise
            # return resulting x value as a (1 x dim) vector
            return np.array(res.x,ndmin=2)

        if(self.acq_func[0] == 'UCB'):
            mult = 0
            #curr_x = self.interface.getState()[0]
            curr_x = self.current_x
            res = minimize(negUCB, curr_x, args=(self.model, mult), bounds=self.bounds, method='L-BFGS-B')

            return np.array(res.x,ndmin=2)

        elif(self.acq_func[0] == 'testEI'):
            # collect all possible x values
            options = np.array(self.acq_func[2].iloc[:, :-1])
            (x_best, y_best) = self.best_seen()

            # find the option with best EI
            best_option_score = (-1,1e12)
            for i in range(options.shape[0]):
                result = negExpImprove(options[i],self.model,y_best,self.acq_func[1])
                if(result < best_option_score[1]):
                    best_option_score = (i, result)

            # return the index of the best option
            return best_option_score[0]
        else:
            print('Unknown acquisition function.')
            return 0


class HyperParams:
    def __init__(self, pvs, filename):
        self.pvs = pvs
        self.filename = filename
        pass

    def extract_hypdata(self, energy):
        key = str(energy)
        #energy = str(round(float(self.mi.get_energy())))
        #if len(energy) is 3: key = energy[0:1]
        #if len(energy) is 4: key = energy[0:2]
        #print ("Loading raw data for", key, "GeV from", self.filename)
        #print()
        f = np.load(str(self.filename))
        filedata = f[0][key]
        return filedata

    def loadHyperParams(self, energy, detector_stat_params):
        """
        Method to load in the hyperparameters from a .npy file.

        Sorts data, ordering parameters with this objects pv list.
        Formats data into tuple format that the GP model object can accept.
        ( [device_1, ..., device_N ], coefficent, noise)


        Args:
                filename (str): String for the file directory.
                energy:

        Returns:
                List of hyperparameters, ordered using the UI's "self.pvs" list.
        """

        #Load in a npy file containing hyperparameters binned for every 1 GeV of beam energy
        extention = self.filename[-4:]
        if extention == ".npy":
            filedata = self.extract_hypdata(energy)

        #sort list to match the UIs PV list order
        #if they are loaded in the wrong order, the optimzer will get the wrong params for a device
        keys = []
        hyps = []
        match_count = 0
        for pv in self.pvs:
            names = filedata.keys()
            if pv in names:
                keys.append(pv)
                ave = float(filedata[pv][0])
                std = float(filedata[pv][1])
                hyp = self.calcLengthScaleHP(ave, std)
                hyps.append(hyp)
                print ("calculate hyper params", pv, ave, std, hyp)
                match_count+=1

        if match_count != len(self.pvs):
            # TODO: what is it?
            # self.parent.scanFinished()
            raise Exception("Number of PVs in list does not match PVs found in hyperparameter file")


        ave, std = detector_stat_params
        print ("DETECTOR AVE", ave)
        print ("DETECTOR STD", std)

        coeff = self.calcAmpCoeffHP(ave, std)
        noise = self.calcNoiseHP(ave, std)

        dout = ( np.array([hyps]), coeff, noise )
        #prints for debug
        print()
        print ("Calculated Hyperparameters ( [device_1, ..., device_N ], amplitude coefficent, noise coefficent)")
        print()
        for i in range(len(hyps)):
            print(self.pvs[i], hyps[i])
        print ("AMP COEFF   = ", coeff)
        print ("NOISE COEFF = ", noise)
        print()
        return dout

    def calcLengthScaleHP(self, ave, std, c = 1.0, pv = None):
        """
        Method to calculate the GP length scale hyperparameters using history data

        Formula for hyperparameters are from Mitch and some papers he read on the GP.

        Args:
                ave (float): Mean of the device, binned around current machine energy
                std (float): Standard deviation of the device
                c   (float): Scaling factor to change the output to be larger or smaller, determined empirically
                pv  (str): PV input string to scale hyps depending on pv, not currently used

        Returns:
                Float of the calculated length scale hyperparameter
        """
        #for future use
        if pv is not None:
            #[pv,val]
            pass
        #+- 1 std around the mean
        hi  = ave+std
        lo  = ave-std
        hyp = -2*np.log( ( ( c*(hi-lo) ) / 4.0 ) + 0.01 )
        return hyp

    def calcAmpCoeffHP(self, ave, std, c = 0.5):
        """
        Method to calculate the GP amplitude hyperparameter

        Formula for hyperparameters are from Mitch and some papers he read on the GP.

        First we tried using the standard deviation to calc this but we found it needed to scale with mean instead


        Args:
                ave (float): Mean of of the objective function (GDET or something else)
                std (float): Standard deviation of the objective function
                c (float): Scaling factor to change the output to be larger or smaller, determined empirically

        Returns:
                Float of the calculated amplitude hyperparameter
        """
        #We would c = 0.5 to work well, could get changed at some point
        hyp2 = np.log( ( ((c*ave)**2) + 0.1 ) )
        return hyp2

    def calcNoiseHP(self, ave, std, c = 1.0):
        """
        Method to calculate the GP noise hyperparameter

        Formula for hyperparameters are from Mitch and some papers he read on the GP.

        Args:
                ave (float): Mean of of the objective function (GDET or something else)
                std (float): Standard deviation of the objective function
                c (float): Scaling factor to change the output to be larger or smaller, determined empirically

        Returns:
                Float of the calculated noise hyperparameter
        """
        hyp = np.log((c*std / 4.0) + 0.01)
        return hyp



def negExpImprove(x_new, model, y_best, xi, alpha=1.0):
    """
    The common acquisition function, expected improvement. Returns the
    negative for the minimizer (so that EI is maximized). Alpha attempts
    to control the ratio of exploration to exploitation, but seems to not
    work well in practice. The terminate() method is a better choice.
    """
    (y_new, var) = model.predict(np.array(x_new, ndmin=2))
    diff = y_new - y_best - xi
    if(var == 0):
        return 0
    else:
        Z = diff / np.sqrt(var)

    EI = diff * norm.cdf(Z) + np.sqrt(var) * norm.pdf(Z)
    #print(x_new, EI)
    return alpha * (-EI) + (1 - alpha) * (-y_new)


def negUCB(x_new, model, mult):
    """
    The upper confidence bound acquisition function. Currently only partially
    implemented. The mult parameter specifies how wide the confidence bound
    should be, and there currently is no way to compute this parameter. This
    acquisition function shouldn't be used until there is a proper mult.
    """
    (y_new, var) = model.predict(np.array(x_new,ndmin=2))

    UCB = y_new + mult * np.sqrt(var)
    return -UCB

def negProbImprove(x_new, model, y_best, xi):
    """
    The probability of improvement acquisition function. Untested.
    Performs worse than EI according to the literature.
    """
    (y_new, var) = model.predict(np.array(x_new,ndmin=2))
    diff = y_new - y_best - xi
    if(var == 0):
        return 0
    else:
        Z = diff / np.sqrt(var)

    return -norm.cdf(Z)
