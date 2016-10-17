# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 15:46:42 2015

@author: Mitch
"""

import numpy as np
import numpy.linalg as alg

# These translate hyperparameters to vector form for minimization

def pack_hyps(xb, hyp_ARD, hyp_coeff, hyp_noise):
    (m,dim) = xb.shape
    xb_row = np.reshape(xb,(1,m*dim))
    return np.concatenate((xb_row,hyp_ARD,hyp_coeff,hyp_noise), axis=1)
    
def unpack_hyps(vec, m, dim):
    vec = np.reshape(vec, (1,(m+1)*dim + 2))
    xb = np.reshape(vec[0][:m*dim],(m,dim))
    hyp_ARD = np.reshape(vec[0][m*dim:-2],(1,dim))
    hyp_coeff = vec[0][-2]
    hyp_noise = vec[0][-1]
    
    return (xb, hyp_ARD, hyp_coeff, hyp_noise)
    
#################################################
    
# Computes kernel on input points
    
def RBF_kernel(x1, x2, hyp_ARD, hyp_coeff, is_self=False):    
    (n1, dim) = x1.shape
    n2 = x2.shape[0]
    
    b = np.exp(hyp_ARD)
    coeff = np.exp(hyp_coeff)
    
    # use ARD to scale
    b_sqrt = np.sqrt(b)
    x1 = x1 * b_sqrt
    x2 = x2 * b_sqrt
    
    try:
        x1_sum_sq = np.reshape(np.sum(x1 * x1, axis=1), (n1,1))
        x2_sum_sq = np.reshape(np.sum(x2 * x2, axis=1), (1,n2))
    except:
        raise
    K = -2 * np.dot(x1, x2.transpose())
    K = K + x1_sum_sq + x2_sum_sq
    K = coeff * np.exp(-.5 * K)
    
    if(is_self):
        jitter = 1e-6
        K = K + jitter * np.eye(n1)
    
    return K
    
#########################################
    
def r2(y, ytrue):
    n = y.shape[0]
    y = np.reshape(y,(n,1))
    ytrue = np.reshape(ytrue,(n,1))
    return 1 - np.sum((y - ytrue)**2) / np.sum((ytrue - ytrue.mean())**2)
    
#############################################
    
# Inverts based on cholesky factorization inverse
    
def chol_invert(factor):
    inv_fact = alg.inv(factor)
    return np.dot(inv_fact.transpose(), inv_fact)
    
#############################################
    
# finds index of closest point in list to input point
    
def closestPoint(x, X):
    x = np.array(x)
    points = np.array(X)
    deltas = points - x
    dist2 = np.sum(deltas * deltas, axis=1)
    return np.argmin(dist2)
    
#############################################
    
# computes pairwise distance between vector elements
    
def pair_dist(x1, x2):
    n1 = max(x1.shape)
    n2 = max(x2.shape)
    return np.reshape(np.subtract.outer(x1,x2),(n1,n2))
    
##############################################
    
# Computes negative log likelihood and derivatives for minimization
    
def SPGP_likelihood(params, y, x, m, compute_deriv=True):
    (N, dim) = x.shape
    y = np.reshape(np.array(y),(N,1))
    (xb, hyp_ARD, hyp_coeff, hyp_noise) = unpack_hyps(np.array([params]), m, dim)

    jitter = 1e-6
    b = np.exp(hyp_ARD)[0]
    b_sqrt = np.sqrt(b)
    coeff = np.exp(hyp_coeff)
    sigma = np.exp(hyp_noise)
    
    # compute Q matrix, covariance between pseudo-inputs
    Q = RBF_kernel(xb, xb, hyp_ARD, hyp_coeff, is_self=True)
    
    # compute K matrix, cov between pseudo-inputs and data
    K = RBF_kernel(xb, x, hyp_ARD, hyp_coeff)
    
    L = alg.cholesky(Q)
    
    V = alg.solve(L, K)
    ep = 1 + np.reshape(coeff - np.sum(V * V, axis=0), (1,N)) / sigma
    ep_sqrt = np.sqrt(ep)
    K = K / ep_sqrt
    V = V / ep_sqrt
    y = y / ep_sqrt.transpose()
    
    Lm = alg.cholesky(sigma * np.eye(m) + np.dot(V, V.transpose()))
        
    invLmV = alg.solve(Lm, V)
    bet = np.dot(invLmV, y)
    
    # compute negative log likelihood
    partial = np.log(Lm.diagonal()).sum()
    partial += (N - m) * hyp_noise / 2
    partial += float(np.dot(y.transpose(), y) - np.dot(bet.transpose(), bet)) / (2 * sigma)
    partial += np.log(ep).sum() / 2
    partial += N * np.log(2 * np.pi) / 2
    
    lik = partial
    
    if(not compute_deriv):
        return lik
    
    # now compute derivates for minimization
    
    # ugly precomputations
    Lt = np.dot(L, Lm)
    B1 = alg.solve(Lt.transpose(), invLmV)
    b1 = alg.solve(Lt.transpose(), bet)
    invLV = alg.solve(L.transpose(), V)
    invQ = chol_invert(L)
    invA = chol_invert(Lt)
    mu = np.dot(alg.solve(Lm.transpose(),bet).transpose(), V).transpose()
    sumVsq = np.reshape(np.sum(V * V, axis=0), (N,1))
    bigsum = 0.5 + y * np.dot(bet.transpose(), invLmV).transpose() / sigma
    bigsum -= np.reshape(np.sum(invLmV * invLmV, axis=0),(N,1)) / 2
    bigsum -= (y*y + mu*mu) / (2 * sigma)
    TT = np.dot(invLV, invLV.transpose() * bigsum)
    
    xb_deriv = np.zeros((m,dim))
    ARD_deriv = np.zeros((1,dim))
    for i in range(dim):
        dnnQ = pair_dist(xb[:,i], xb[:,i]) * Q
        dNnK = pair_dist(-xb[:,i],-x[:,i]) * K
        
        epdot = -2 * dNnK * invLV / sigma
        epPmod = -np.reshape(np.sum(epdot, axis=0), (N,1))
        
        dxb = -b1 * (np.dot(dNnK, (y - mu)) / sigma + np.dot(dnnQ, b1))
        dxb += np.reshape(np.sum((invQ - sigma * invA) * dnnQ,axis=1),(m,1))
        dxb += np.dot(epdot, bigsum)
        dxb -= 2 * np.reshape(np.sum(dnnQ * TT, axis=1),(m,1)) / sigma
        
        dARD = ((y - mu).transpose() * np.dot(b1.transpose(),dNnK)) / sigma
        dARD = dARD + (epPmod * bigsum).transpose()
        dARD = np.dot(dARD, np.reshape(x[:,i],(N,1)))
        
        dNnK = dNnK * B1
        dxb = dxb + np.reshape(np.sum(dNnK,axis=1),(m,1))
        dARD = dARD - np.dot(np.reshape(np.sum(dNnK,axis=0),(1,N)), np.reshape(x[:,i],(N,1)))
        
        dxb = b_sqrt[i] * dxb
        
        dARD = dARD / b_sqrt[i]
        dARD = dARD + np.dot(dxb.transpose(),np.reshape(xb[:,i],(m,1))) / b[i]
        dARD = b_sqrt[i] * dARD / 2
        
        xb_deriv[:,i] = dxb[:,0]
        ARD_deriv[0][i] = dARD
        
    ep = ep.transpose()
    epc = (coeff / ep - sumVsq - jitter * np.reshape(np.sum(invLV * invLV,axis=0),(N,1))) / sigma
    
    coeff_deriv = (m + jitter*np.trace(invQ - sigma * invA))
    coeff_deriv -= sigma * np.sum(invA * Q.transpose()) / 2
    coeff_deriv -= np.dot(mu.transpose(),y - mu) / sigma
    coeff_deriv += np.dot(np.dot(b1.transpose(), Q - jitter * np.eye(m)),b1) / 2
    coeff_deriv += np.dot(epc.transpose(), bigsum)
    
    noise_deriv = np.reshape(np.sum(bigsum / ep),(1,1))
    
    deriv = pack_hyps(xb_deriv, ARD_deriv, coeff_deriv, noise_deriv)
       
    return lik, deriv[0] / alg.norm(deriv)