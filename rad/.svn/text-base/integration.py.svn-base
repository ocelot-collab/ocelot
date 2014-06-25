'''
compute oscillating integrals
'''

import sys
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter
import copy
import string


def f(x,w):
    return np.sin(w*x) * np.sin(w*x)

def I0(L,w):
    return w * (L/2.0 - np.sin(2*L*w)/(4*w))

def I1(L,w):
    dx = 0.01
    x = np.arange(0,L,dx)
    
    return w*sum(f(x, w)*dx)

def I2(L,w):
    dx = 0.01
    n  = int( L / dx)
    x = L * np.random.rand(n)
    
    return w*sum(f(x, w)*dx)



if __name__=="__main__":
     
    ws = np.arange(0.0001, 1000, 0.1)
    
    i0 = np.zeros(len(ws))
    i1 = np.zeros(len(ws))
    i2 = np.zeros(len(ws))
    
    for i in range(0,len(ws)):
        i0[i] = I0(1.0,ws[i])
        i1[i] = I1(1.0,ws[i])
        i2[i] = I2(1.0,ws[i])
    
    line, = plt.plot(np.arange(len(i0)),i0, '-')
    line, = plt.plot(np.arange(len(i1)),i1, '-')
    line, = plt.plot(np.arange(len(i2)),i2, '-')

    plt.show()