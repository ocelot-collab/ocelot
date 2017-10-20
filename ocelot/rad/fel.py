from __future__ import print_function
'''
basic fel calculations
'''

#from pylab import *
import numpy as np
import numpy.fft as fft
import scipy.special as sf
from ocelot.common.globals import m_e_eV, epsilon_0, m_e_kg, speed_of_light, q_e
#from matplotlib.figure import Figure
#from mpl_toolkits.mplot3d import Axes3D

#import fel

class FelParameters:
    def __init__(self):
        pass
    
    def P(self, z=None, dims=3):
        '''
        unfinished
        '''
        if dims == 1:
            rho = self.rho1
            lg = self.lg1
        elif dims == 3:
            rho = self.rho3
            lg = self.lg3
        else:
            raise ValueError('dims argument should be either 1 or 3')
        
        Nc = self.Ip / (q_e * rho * self.k0 * speed_of_light)
        z_sat = 3 + 1/np.sqrt(3) * np.log(Nc)
        Psn = (3 * rho * self.Pb) / (Nc * np.sqrt(np.pi * np.log(Nc)))
        
        if z is None:
            zn = np.amin(z_sat)
        else:
            zn = z / (np.sqrt(3) * lg)
        
        Pz = Psn * (1 + 1/9 * np.exp(np.sqrt(3) * zn) / np.sqrt(np.pi * zn))
        
        return Pz

def calculateFelParameters(input):
    p = FelParameters()
    
    p.iwityp = input.iwityp # undulator type: 0 == planar, other == helical
    
    p.gamma0 = input.gamma0  
    p.delgam = input.delgam  
    p.xlamd = input.xlamd    # undulator period
    
    p.betax = input.betax
    p.betay = input.betay
    p.ex = input.emitx #normalized emittance
    p.ey = input.emity
    p.rxbeam = input.rxbeam
    p.rybeam = input.rybeam
    # p.rxbeam = np.sqrt(p.betax * p.ex / p.gamma0)
    # p.rybeam = np.sqrt(p.betay * p.ey / p.gamma0)
    
    p.aw0 = input.aw0 # rms undulator parameter K
    p.Ip = input.curpeak
    
    p.deta =  p.delgam / p.gamma0
    
    p.lambda0 = p.xlamd / (2.0 * p.gamma0**2) * (1.0 + p.aw0**2) # resonant wavelength
    p.k0 = 2 * np.pi / p.lambda0 

    p.Ia = 4 * np.pi * epsilon_0 * m_e_kg * speed_of_light**3 / q_e # Alfven (Budker) current (~17kA)
#    p.Ia = 17000
        
    if p.iwityp == 0:
        ja = p.aw0**2 / (2*(1+p.aw0**2))
        p.fc = sf.j0(ja) - sf.j1(ja)
    else:
        p.fc = 1.0
    
    # import first, ro_e * m_e_eV = 1.4399643147059695e-09
    
#    p.N = p.Ip * p.lambda0 / 1.4399644850445153e-10
    p.sigb = 0.5 * (p.rxbeam + p.rybeam) # average beam size
    
    p.rho1 = (0.5 / p.gamma0) * np.power( (p.aw0 * p.fc * p.xlamd / (2.0 * np.pi * p.sigb) )**2 * p.Ip / p.Ia, 1.0/3.0) ## check 8 in denominator
    p.Pb = p.gamma0 * p.Ip * m_e_eV# beam power [Reiche]
    
    #p.power = 6.0 * np.sqrt(np.pi) * p.rho1**2 * p.Pb / (p.N * np.log(p.N / p.rho1) ) # shot noise power [W] [Reiche]
    p.lg1 = p.xlamd / (4*np.pi * np.sqrt(3) * p.rho1) #[Xie]
    p.zr = 4 * np.pi * p.sigb**2 / p.lambda0
  
    a = [None, 0.45, 0.57, 0.55, 1.6, 3.0, 2.0, 0.35, 2.9, 2.4, 51.0, 0.95, 3.0, 5.4, 0.7, 1.9, 1140.0, 2.2, 2.9, 3.2]
    
    p.xie_etad = p.lg1 / (2 * p.k0 * p.sigb**2)
    #p.xie_etae = 4 * pi * p.lg1 / (p.betax*2*pi) * p.k0 * (p.ex / p.gamma0)
    p.xie_etae = 4 * np.pi * p.lg1 / p.lambda0 * (p.ex / p.gamma0 / p.sigb)**2 # expressed via average x-y beam size
    p.xie_etagamma = p.deta / (p.rho1 * np.sqrt(3))
    p.xie_lscale = (a[1] * p.xie_etad ** a[2] + a[3] * p.xie_etae ** a[4] + a[5] * p.xie_etagamma ** a[6] 
    + a[7] * p.xie_etae ** a[8] * p.xie_etagamma ** a[9] + a[10] * p.xie_etad ** a[11] * p.xie_etagamma ** a[12] + a[13] * p.xie_etad ** a[14] * p.xie_etae ** a[15]
    + a[16] * p.xie_etad ** a[17] * p.xie_etae ** a[18] * p.xie_etagamma ** a[19])

    p.lg3 = p.lg1 * (1 + p.xie_lscale)
    p.rho3 = p.xlamd / (4*np.pi * np.sqrt(3) * p.lg3)
    

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
    print ('coupling parameter fc=', p.fc)
    print ('gamma0=', p.gamma0)
    print ('Ip=', p.Ip, ' beam peak current [A]')
    print ('lambda0=', p.lambda0)
    print ('Pb=', p.Pb, ' beam power [W]')
    # print ('N=', p.N)
    print ('rho (1D)=', p.rho1)
    print ('gain length estimate lg (1D)=', p.lg1)
    # print ('power=', p.power, ' equivalent shot noise power [W]')
    print ('Rayleigh length estimate zr=', p.zr)
    print ('')
    print ('Ming Xie gain reduction estimates:')
    print ('diffraction parameter eta_d=', p.xie_etad)
    print ('emittance/focusing parameter eta_e=', p.xie_etae)
    print ('energy spread parameter eta_gamma=', p.xie_etagamma)
    print ('gain length degradation lscale=', p.xie_lscale)
    print ('scaled gain length lg (3D)=', p.lg3)
    print ('scaled rho (3D)=', p.rho3)
    print ('**************************************')
    