# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 19:46:52 2015

@author: Mitch
"""

import numpy as np
import numpy.linalg as alg
from minimize import minimize
from GP_utils import *

def SPGP_train(X,Y,num_pseudo_inputs,num_starts=1):
    """
    Trains a sparse Gaussian process on the input data.
    X -- DataFrame with training data (n x dim)  
    Y -- Labels for training data (n x 1)
    num_pseudo_inputs -- number of points used to fill sparse model
    num_starts -- number of attempts at minimization. Increases runtime linearly.
    
    Returns:
    xb -- pseudo-inputs as ndarray (m x dim)
    hyperparams -- tuple containing GP parameters
    
    Translated to python from Edward Snelson's matlab code by Mitchell McIntire.
    """
    
    (n, dim) = X.shape
    m = np.min([num_pseudo_inputs,n])
    
    # center data
    mu_y = np.mean(Y)
    y0 = Y - mu_y
    
    min_lik = np.inf
    for i in range(num_starts):    
        # randomly choose initial points
        #   should randomly sample, but hacking this in for the ACR since
        #   the pandas version is older
        #xb_init = np.array(X.sample(m))
        xb_init = np.array(X.iloc[:m,:])
    
        # initialize hyperparameters
        hyp_ARD = np.array([-2*np.log((X.max() - X.min() + 0.1) / 2)])
        hyp_coeff = np.array([[np.log(Y.var() + 0.1)]])
        hyp_noise = np.array([[np.log(Y.var() / 4 + 0.01)]])
        hyperparams = pack_hyps(xb_init, hyp_ARD, hyp_coeff, hyp_noise)
        
        # minimize neg. log likelihood
        # min_result = minimize(SPGP_likelihood, hyperparams, args=(y0,np.array(X),m), method='BFGS', jac=True)
        #iter_res = np.reshape(min_result.x, (1,(m+1)*dim + 2))
        #lik = SPGP_likelihood(iter_res,y0,np.array(X),m,compute_deriv=False)
        #st = time.time()
        (iter_res, lik, i) = minimize(hyperparams, SPGP_likelihood, args=(y0,np.array(X), m), maxnumfuneval=200)
        #print(time.time() - st)
        if(lik[0] < min_lik):
            min_lik = lik[0]
            opt_res = iter_res
    
    # extract minimizing hyperparameters
    (xb, hyp_ARD, hyp_coeff, hyp_noise) = unpack_hyps(opt_res, m, dim)
    
    hyperparams = (hyp_ARD, hyp_coeff, hyp_noise)
    
    return xb, hyperparams #, mu_y
    
    
def SPGP_predict(X, y, xb, xt, hyperparams):
    (N, dim) = X.shape
    X = np.array(X)
    xt = np.array(xt)
    xb = np.array(xb)
    y = np.reshape(np.array(y), (N,1))
    m = xb.shape[0]
    
    (hyp_ARD, hyp_coeff, hyp_noise) = hyperparams
    sigma = np.exp(hyp_noise)
    coeff = np.exp(hyp_coeff)
    
    K = RBF_kernel(xb, xb, hyp_ARD, hyp_coeff, is_self=True)
    L = alg.cholesky(K)
    K = RBF_kernel(xb, X, hyp_ARD, hyp_coeff)
    V = alg.solve(L, K)
    
    ep = 1 + np.reshape(coeff - np.sum(V * V, axis=0), (1,N)) / sigma
    ep_sqrt = np.sqrt(ep)
    V = V / ep_sqrt
    y = y / ep_sqrt.transpose()
    
    Lm = alg.cholesky(sigma * np.eye(m) + np.dot(V, V.transpose()))
    bet = alg.solve(Lm, np.dot(V, y))
    
    K = RBF_kernel(xb, xt, hyp_ARD, hyp_coeff)
    lst = alg.solve(L, K)
    lmst = alg.solve(Lm, lst)
    
    mu = np.dot(bet.transpose(), lmst).transpose()
    lst_cols = np.sum(lst * lst, axis=0).transpose()
    lmst_cols = np.sum(lmst * lmst, axis=0).transpose()
    s2 = coeff - lst_cols + sigma * lmst_cols
    
    return mu, s2
    
    
