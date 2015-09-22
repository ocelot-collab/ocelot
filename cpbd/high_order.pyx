__author__ = 'Sergey Tomin'
import numpy as np
from ocelot.common.globals import *
def sym_map(z, X, h, k1, k2, energy=0.):

    if h != 0. or k1 != 0. or k2 != 0.:
        step = 0.005
    else:
        step = z
    if step > z:
        step = z
    if step == 0.:
        return X
    n = int(z/step) + 1
    gamma = energy/m_e_GeV
    g_inv = 0.
    beta = 1.
    if gamma != 0.:
        g_inv = 1./gamma
        beta = np.sqrt(1. - g_inv*g_inv)
    z_array = np.linspace(0., z, num=n)
    step = z_array[1] - z_array[0]

    x =     X[0::6]
    px =    X[1::6]
    y =     X[2::6]
    py =    X[3::6]
    sigma = X[4::6]
    ps =    X[5::6]

    ps_beta = ps/beta
    #vec = [x, px, y, py, sigma, ps]
    c1 = h*h + k1
    c2 = h*k1 + k2/2.
    c3 = h*k1 + k2
    c4 = g_inv*g_inv/(beta*(1. + beta))
    for z in xrange(len(z_array) - 1):
        #vec = verlet(vec, step, h, k1, k2, beta=beta, g_inv=g_inv)
        px2_py2 = px*px + py*py
        x = (x + step*px*(1. - ps_beta))/(1. - step*h*px)
        y = y + step*py*(1. + h*x - ps_beta)
        sigma = sigma + step*(-h*x/beta - px2_py2/(2.*beta) - c4)

        px = px + step*(h*ps_beta + (-h*px2_py2 + c3*y*y)/2. - (c1 + c2*x)*x)
        py = py + step*(k1 + c3*x)*y

    X[0::6] = x[:] #vec[0][:]
    X[1::6] = px[:] #vec[1][:]
    X[2::6] = y[:] #vec[2][:]
    X[3::6] = py[:] #vec[3][:]
    X[4::6] = sigma[:] #vec[4][:]
    return X

def cython_test(X, step, h,k1, N, cc, ps_beta, beta):
    x =     X[0]
    px =    X[1]
    y =     X[2]
    py =    X[3]
    sigma = X[4]
    ps =    X[5]
    c1, c2, c3, c4 = cc
    for z in xrange(N - 1):
        #vec = verlet(vec, step, h, k1, k2, beta=beta, g_inv=g_inv)
        px2_py2 = px*px + py*py
        x = (x + step*px*(1. - ps_beta))/(1. - step*h*px)
        y = y + step*py*(1. + h*x - ps_beta)
        sigma = sigma + step*(-h*x/beta - px2_py2/(2.*beta) - c4)

        px = px + step*(h*ps_beta + (-h*px2_py2 + c3*y*y)/2. - (c1 + c2*x)*x)
        py = py + step*(k1 + c3*x)*y
    X[0] = x
    X[1] = px
    X[2] = y
    X[3] = py
    X[4] = sigma
    X[5] = ps
    return X