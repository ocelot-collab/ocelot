import numpy as np
from matplotlib import pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
from scipy import optimize
from copy import deepcopy
import time
import matplotlib.cm as cm

np.random.seed(1)

class GP:
    def __init__(self):
        self.x_search = np.array([])
        self.x_obs = np.array([])
        self.y_obs = np.array([])
        self.y_sigma_obs = np.array([])
        self.y_pred = np.array([])
        self.sigma_y = 0

        # RBF kernel
        self.rbf_length_scale = 1
        self.rbf_length_scale_bounds = (0.01, 100)
        # ConstantKernel
        self.ck_const_value = 1.0
        self.ck_const_value_bounds = (1e-05, 100000.0)

        self.n_restarts_optimizer = 10

        self.max_iter = 40
        self.pen_max = 100

        self.ytol = 0.001
        self.xtol = 0.001

        self.opt_ctrl = None

    def append_new_data(self, x_new, y_obs, sigma_y_obs):
        # add new data
        #if np.size(x_cur) == 1:
        #    x_cur = [x_cur]
        #self.x_obs = np.concatenate((self.x_obs, x_new), axis=0)
        self.x_obs = np.append(self.x_obs, [x_new], axis=0)
        self.y_obs = np.append(self.y_obs, y_obs)
        self.y_sigma_obs = np.append(self.y_sigma_obs, sigma_y_obs)
        #print(self.y_sigma_init)

        #self.fit()

    def fit(self):
        """
        RBF(length_scale=1.0, length_scale_bounds=(1e-05, 100000.0))
        k(x_i, x_j) = exp(-1 / 2 d(x_i / length_scale, x_j / length_scale)^2)

        :return:
        """
        # Instanciate a Gaussian Process model
        kernel = ConstantKernel(self.ck_const_value, self.ck_const_value_bounds)\
                 * RBF(self.rbf_length_scale, self.rbf_length_scale_bounds)
        # Instanciate a Gaussian Process model
        if self.sigma_y != 0:
            self.alpha = (self.y_sigma_obs / self.y_obs) ** 2
        else:
            self.alpha = 1e-10
        self.gp = GaussianProcessRegressor(kernel=kernel, alpha=self.alpha,
                                           n_restarts_optimizer=self.n_restarts_optimizer)

        # Fit to data using Maximum Likelihood Estimation of the parameters
        self.gp.fit(self.x_obs, self.y_obs)

    #def find_max(self):
    #    # Make the prediction on the meshed x-axis (ask for MSE as well)
    #    y_pred, sigma = self.gp.predict(self.x_search, return_std=True)
    #    # print("y_pred", y_pred)
    #    self.y_pred = y_pred
    #    self.sigma = sigma
    #    return self.x_search[np.argmax(y_pred)]

    def acquire_simplex(self):
        # Make the prediction on the meshed x-axis (ask for MSE as well)
        def func(x):
            for i, xi in enumerate(x):
                if self.x_search[0][i] > xi or xi > self.x_search[-1][i]:
                    print("exceed limits ")
                    return self.pen_max

            y_pred, sigma = self.gp.predict(np.atleast_2d(x), return_std=True)
            self.sigma = sigma
            return y_pred

        y_pred, sigma = self.gp.predict(self.x_obs, return_std=True)
        x = self.x_obs[np.argmin(y_pred)]
        res = optimize.fmin(func, x)
        return res #self.x_search[np.argmin(y_pred)]


    def acquire(self):
        # Make the prediction on the meshed x-axis (ask for MSE as well)

        y_pred, sigma = self.gp.predict(self.x_search, return_std=True)
        x = self.x_search[np.argmin(y_pred)]
        #res = optimize.fmin(func, x)
        return x #self.x_search[np.argmin(y_pred)]

    def minimize(self, error_func, x):
        # weighting for exploration vs exploitation in the GP at the end of scan, alpha array goes from 1 to zero
        self.fit()
        for i in range(self.max_iter):
            # get next point to try using acquisition function
            if self.opt_ctrl != None and self.opt_ctrl.kill == True:
                print('GP: Killed from external process')
                break
            start = time.time()
            x_next = self.acquire()
            print("acquire ", start - time.time(), " sec")

            y_new = error_func(x_next.flatten())

            self.append_new_data(x_next, y_new, sigma_y_obs=self.sigma_y)

            # update the model (may want to add noise if using testEI)
            self.fit()
            #print("niter = ", i, self.x_obs[-3] , self.x_obs[-1], self.x_obs[-3] - self.x_obs[-1])
            if i>3 and np.linalg.norm((self.x_obs[-3] - self.x_obs[-1])) <= self.xtol:
                #print("ytol = ", np.linalg.norm((self.x_obs[-3] - self.x_obs[-1])))
                break
        return self.x_obs[-1]


def f(x):
    """The function to predict."""
    sigma = 2
    y = -np.exp(-np.sum((np.ones(np.size(x))*9 - x)**2)/(2*sigma**2))

    return y


def seed_simplex(error_func, x):

    res = optimize.fmin(error_func, x)

    return res


def test():
    xv, yv = np.meshgrid(np.arange(0, 10, 0.1), np.arange(0, 10, 0.1))
    n = np.shape(xv)

    xv = xv.reshape(n[0]*n[1], 1)
    yv = yv.reshape(n[0]*n[1], 1)
    X = np.hstack((xv, yv))
    y = np.array([f(x) for x in X])
    fig, ax = plt.subplots()
    im = ax.imshow(y.reshape(n), extent=[0, 10, 10, 0],  interpolation='bilinear', cmap=cm.RdYlGn)
    fig.colorbar(im)
    plt.show()



    #  First the noiseless case
    X1 = np.atleast_2d([0., 0.5, 0.5, 3., 4]).T
    X2 = np.atleast_2d([0., 0., 1., 2., 3]).T
    X = np.hstack([X1, X2])
    #print(X)

    # Observations
    y = [f(x) for x in X]

    sigma_y = 0.001*np.ones_like(y)
    noise = np.random.normal(0, sigma_y)
    print(noise)
    y += noise

    gp = GP()
    gp.ck_const_value=5
    gp.rbf_length_scale = np.sqrt(2)*4
    gp.rbf_length_scale_bounds = (gp.rbf_length_scale, gp.rbf_length_scale)
    gp.x_obs = X
    gp.y_obs = y
    gp.y_sigma_obs = sigma_y
    gp.sigma_y = 0.001

    xv, yv = np.meshgrid(np.arange(0, 10, 0.1), np.arange(0, 10, 0.1))
    n = np.shape(xv)
    xv = xv.reshape(n[0]*n[1], 1)
    yv = yv.reshape(n[0]*n[1], 1)
    gp.x_search = np.hstack((xv, yv))

    print(gp.x_search, X, X[-1])
    x = gp.minimize(error_func=f, x=X[-1])
    print(x)
    #for i in range(15):
    #    fig, ax = plt.subplots()
    #    gp.fit()
    #
    #    x_cur = gp.acquire()
    #    y_cur = f(x_cur)
    #
    #    sigma_y_cur = sigma_y[0]
    #    print("next", x_cur, y_cur, sigma_y[0])
    #    gp.append_new_data(x_cur, y_cur, sigma_y_cur)
    #
    #    y_pred, sigma = gp.gp.predict(gp.x_search, return_std=True)
    #    ax.plot( gp.x_obs[:,0], gp.x_obs[:,1], 'k.-', markersize=10)
    #    im = ax.imshow(y_pred.reshape(n), extent=[0, 10, 10, 0], cmap=cm.RdYlGn,  interpolation='bilinear')
    #    fig.colorbar(im)
    #    plt.show()

def one_d():

    X = np.arange(0, 10, 0.1)
    y = np.array([f(x) for x in X])
    fig, ax = plt.subplots()
    ax.plot(X, y)
    plt.show()

    gp = GP()
    #  First the noiseless case
    X = np.atleast_2d([7., 4., 5., 1., 0]).T


    # Observations
    y = np.array([f(x) for x in X])

    sigma_y = 0.01 * np.abs(y)
    noise = np.random.normal(0, sigma_y)
    #print(noise)
    y += noise
    gp.ck_const_value = 5
    gp.rbf_length_scale = np.sqrt(2) * 4
    # gp.rbf_length_scale_bounds = (gp.rbf_length_scale, gp.rbf_length_scale)
    gp.x_obs = X
    gp.y_obs = y
    gp.y_sigma_obs = sigma_y

    gp.x_search = np.atleast_2d(np.arange(0, 10, 0.1)).T

    #print(gp.x_search)

    for i in range(15):
        fig, ax = plt.subplots()
        gp.fit()
        #gp.find_max()
        x_cur = gp.acquire()
        y_cur = f(x_cur)

        plt.plot(gp.x_search, np.array([f(x) for x in gp.x_search]), 'r:', label=u'$f(x) = x\,\sin(x)$')
        plt.plot(gp.x_obs, gp.y_obs, 'r.', markersize=10, label=u'Observations')
        #print(gp.y_pred)

        y_pred, sigma = gp.gp.predict(gp.x_search, return_std=True)
        plt.plot(gp.x_search, y_pred, "b:")
        plt.fill(np.concatenate([gp.x_search, gp.x_search[::-1]]),
            np.concatenate([y_pred - 1.9600 * sigma,
                           (y_pred + 1.9600 * sigma)[::-1]]),
            alpha=.5, fc='b', ec='None', label='95% confidence interval')


        sigma_y_cur = sigma_y[0]
        print("next", x_cur, y_cur, sigma_y[0])
        gp.append_new_data(x_cur, y_cur, sigma_y_cur)

        plt.show()

def two_d():
    xv, yv = np.meshgrid(np.arange(0, 10, 0.1), np.arange(0, 10, 0.1))
    n = np.shape(xv)

    xv = xv.reshape(n[0]*n[1], 1)
    yv = yv.reshape(n[0]*n[1], 1)
    X = np.hstack((xv, yv))
    y = np.array([f(x) for x in X])
    fig, ax = plt.subplots()
    im = ax.imshow(y.reshape(n), extent=[0, 10, 10, 0],  interpolation='bilinear', cmap=cm.RdYlGn)
    fig.colorbar(im)
    plt.show()


    gp = GP()
    #  First the noiseless case
    X1 = np.atleast_2d([0., 0.5, 0.5, 3., 4]).T
    X2 = np.atleast_2d([0., 0., 1., 2., 3]).T
    X = np.hstack([X1, X2])
    #print(X)

    # Observations
    y = [f(x) for x in X]

    sigma_y = 0.001*np.ones(len(y))
    noise = np.random.normal(0, sigma_y)
    print(noise, sigma_y)
    y += noise
    gp.ck_const_value=5
    gp.rbf_length_scale = np.sqrt(2)*4
    gp.rbf_length_scale_bounds = (gp.rbf_length_scale, gp.rbf_length_scale)
    gp.x_obs = X
    gp.y_obs = y
    gp.y_sigma_obs = sigma_y
    gp.sigma_y = 0.001
    xv, yv = np.meshgrid(np.arange(0, 10, 0.1), np.arange(0, 10, 0.1))
    n = np.shape(xv)

    xv = xv.reshape(n[0]*n[1], 1)
    yv = yv.reshape(n[0]*n[1], 1)
    gp.x_search = np.hstack((xv, yv))

    print(gp.x_search)

    for i in range(15):
        fig, ax = plt.subplots()
        gp.fit()

        x_cur = gp.acquire()
        y_cur = f(x_cur)

        #sigma_y_cur = sigma_y[0]
        print("next", x_cur, y_cur, sigma_y[0])
        gp.append_new_data(x_cur, y_cur, gp.sigma_y)

        y_pred, sigma = gp.gp.predict(gp.x_search, return_std=True)
        ax.plot( gp.x_obs[:,0], gp.x_obs[:,1], 'k.-', markersize=10)
        im = ax.imshow(y_pred.reshape(n), extent=[0, 10, 10, 0], cmap=cm.RdYlGn,  interpolation='bilinear')
        fig.colorbar(im)
        plt.show()


if __name__ == "__main__":
    one_d()
    #two_d()
    #test()