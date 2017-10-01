from __future__ import print_function
'''
basic fel calculations
'''

#from pylab import *
import numpy as np
import numpy.fft as fft
import scipy.special as sf
#from matplotlib.figure import Figure
#from mpl_toolkits.mplot3d import Axes3D

#import fel

class FelParameters:
    def __init__(self):
        pass

def calculateFelParameters(input):
    p = FelParameters()

    p.gamma0 = input.gamma0  
    p.delgam = input.delgam  
    p.xlamd = input.xlamd    # undulator period
    p.ex = input.emitx
    p.ey = input.emity
    p.rxbeam = input.rxbeam
    p.rybeam = input.rybeam
    p.aw0 = input.aw0
    p.Ip = input.curpeak
    
    p.deta =  p.delgam / p.gamma0
    
    p.lambda0 = p.xlamd / (2.0 * p.gamma0**2) *(1.0 + p.aw0**2)
    p.k0 = 2 * np.pi / p.lambda0 

    p.Ia = 17000.0
        
    a = p.aw0**2 / (2*(1+p.aw0**2))
    p.fc = sf.j0(a) - sf.j1(a)
    p.N = p.Ip * p.lambda0 / 1.4399644850445153e-10
    p.sigb = 0.5 * (p.rxbeam + p.rybeam)
    p.rho = (1.0 / p.gamma0) * np.power( (p.aw0 * p.fc * p.xlamd / (8.0 * np.pi * p.sigb) )**2 * p.Ip / p.Ia, 1.0/3.0)
    p.Pb = p.gamma0 * p.Ip * 510998.927
    p.power = 6.0 * np.sqrt(np.pi) * p.rho**2 * p.Pb / (p.N * np.log(p.N / p.rho) )
    p.lg = p.xlamd / (4*np.pi * np.sqrt(3) * p.rho)
    p.zr = 4 * np.pi * p.sigb**2 / p.lambda0
  
    xie_a = [0.45, 0.57, 0.55, 1.6, 3.0, 2.0, 0.35, 2.9, 2.4, 51.0, 0.95, 3.0, 5.4, 0.7, 1.9, 1140.0, 2.2, 2.9, 3.2]
      
    p.xie_etad = p.lg / (2 * p.k0 * p.sigb**2)
    p.xie_etae = 0
    p.xie_etagamma = p.deta / (p.rho * np.sqrt(3))
    p.xie_lscale = xie_a[0] * p.xie_etad ** xie_a[1] + xie_a[2] * p.xie_etae ** xie_a[3] + xie_a[4] * p.xie_etagamma ** xie_a[5] 
                    
    return p


def printFelParameters(input):
    
    #print (input.parameters)
    
    p = calculateFelParameters(input)
    
    print ('********    FEL Parameters    ********')
    print ('ex=', p.ex)
    print ('ey=', p.ey)
    print ('rxbeam=', p.rxbeam, ' [m]')
    print ('rybeam=', p.rybeam, ' [m]')
    print ('rel energy spread deta=', p.deta, ' [m]')
    print ('xlamd=', p.xlamd)
    print ('aw0=', p.aw0)
    print ('gamma0=', p.gamma0)
    print ('Ip=', p.Ip, ' beam peak current [A]')
    print ('lambda0=', p.lambda0)
    print ('Pb=', p.Pb, ' beam power [W]')
    print ('N=', p.N)
    print ('rho=', p.rho)
    print ('power=', p.power, ' equivalent shot noise power [W]')
    print ('coupling parameter fc=', p.fc)
    print ('gain length estimate lg=', p.lg)
    print ('Rayleigh length estimate zr=', p.zr)
    print ('')
    print ('Ming Xie gain reduction estimates:')
    print ('diffraction parameter etad=', p.xie_etad)
    print ('energy spread parameter etad=', p.xie_etagamma)
    print ('gain length degradation lscale=', p.xie_lscale)
    print ('**************************************')
    