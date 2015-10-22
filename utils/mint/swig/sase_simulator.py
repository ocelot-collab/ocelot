'''
sase (GMD measurement) mock-up
'''
#import ocelot.utils.mint.mint as mint
#import ocelot.utils.mint.swig.dcs as dcs
import os, sys

from pylab import *
from scipy.optimize import *
from time import sleep
from pickle import dump, load


from time import *
from numpy.fft import *

from ocelot.common.math_op import *

#print localtime(time()).tm_min, localtime(time()).tm_sec, localtime(time())


class Settings:
    def __init__(self):
        pass

s = Settings()

def init_settings(s):
    s.sase_fluct_ampl_1 = 0.1
    s.sase_fluct_freq_1 = 0.02
    s.sase_fluct_phase_1 = 0.00

    s.sase_fluct_ampl_2 = 0.1
    s.sase_fluct_freq_2 = 0.00512
    s.sase_fluct_phase_2 = 0.00
    s.t = 0.0



def stable_sase(cor_val):
    return exp( -(cor_val - 0.55)**2 / 0.05**2)


def get_sase(cor_val, dt = 1.0):
    global s
    s.sase_fluct_ampl_1 += np.random.randn()*0.01
    s.sase_fluct_phase_1 += np.random.randn()*0.01

    s.sase_fluct_ampl_2 += np.random.randn()*0.01
    s.sase_fluct_phase_2 += np.random.randn()*0.01
    
    sase = stable_sase(cor_val)
     
    if sase > 0.1: sase *=  (1. + s.sase_fluct_ampl_1 * sin(2.*pi*s.sase_fluct_freq_1 * s.t + s.sase_fluct_phase_1) ) 
    if sase > 0.1: sase *=  (1. + s.sase_fluct_ampl_2 * sin(2.*pi*s.sase_fluct_freq_2 * s.t + s.sase_fluct_phase_2) )
    
    s.t += dt
    
    #print 'getting sase', s.t
    
    return sase


def scan_cor(minv, maxv, n, t_scan = 1):
    print s.sase_fluct_ampl_1

    gmd = np.zeros(n)
    cor_val = np.zeros(n)
    for i in range(len(gmd)):
        
        if i < n:
            cor_val[i] = minv + (maxv - minv) * i / float(n)
        else:
            cor_val[i] = maxv - (maxv - minv)*(i - n) / float(n)
        
        gmd[i] = get_sase(cor_val[i], dt=t_scan)

    return cor_val, gmd


def max_sase(correctors):

    def error_func(x):
        
        pen_max = 100.0

        for i in xrange(13):
            sase = get_sase(x[0])

        print 'sase:', sase
        print 'x[0]=', x[0] 

        pen = 0.0

        pen -= sase

        return pen
    
        
    x = [0.05]
    res  = fmin(error_func,x,xtol=1e-3, maxiter=200, maxfun=200)


def max_sase_2(correctors):
    #rough_scan
    init_settings(s)
    cor_val, gmd = scan_cor(0.0, 1.9, 20, t_scan = 20)    
    plt.plot(cor_val, gmd, '.-')
    
    idx = (-gmd).argsort()[:5] # 5 largest values
    #TODO: this can be replaced by e.g. taking all data above threshold level
    
    print idx 
    i1, i2 = sorted(idx)[0], sorted(idx)[-1]
    print cor_val[i1]
    print cor_val[i2]
    
    cor_val, gmd = scan_cor(cor_val[i1], cor_val[i2], 50, t_scan = 20) 
    plt.plot(cor_val, gmd, '.-')


    idx = (-gmd).argsort()[:5] # 5 largest values
    i1, i2 = sorted(idx)[0], sorted(idx)[-1]
    print cor_val[i1]
    print cor_val[i2]
    
    cor_val, gmd = scan_cor(cor_val[i1], cor_val[i2], 50, t_scan = 20) 
    plt.plot(cor_val, gmd, '.-')

    cor_mean = mean(cor_val)
    print 'cor_val=', cor_mean
    
    plt.show()


    
if __name__ == "__main__":
    
    '''
    init_settings(s)
    cor_val, gmd = scan_cor(0.4, 0.7, 20, t_scan = 20)    
    plt.plot(cor_val, gmd, '.-')
    
    init_settings(s)
    cor_val, gmd = scan_cor(0.4, 0.7, 20, t_scan = 20)    
    plt.plot(cor_val, gmd, '.-')
    
    init_settings(s)
    cor_val, gmd = scan_cor(0.4, 0.7, 20, t_scan = 20)    
    plt.plot(cor_val, gmd, '.-')
    plt.show()
    '''
    
    max_sase_2(['cor1'])