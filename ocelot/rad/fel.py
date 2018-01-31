from __future__ import print_function
'''
basic fel calculations
'''

#from pylab import *
import numpy as np
import numpy.fft as fft
import scipy.special as sf
from ocelot.common.globals import m_e_eV, epsilon_0, speed_of_light, q_e, hr_eV_s
#from matplotlib.figure import Figure
#from mpl_toolkits.mplot3d import Axes3D

#import fel

class FelParameters:
    def __init__(self):
        pass
    
    def P(self, z=None):
        '''
        unfinished
        '''
        # if dims == 1:
            # rho = self.rho1
            # lg = self.lg1
        # elif dims == 3:
            # rho = self.rho3
            # lg = self.lg3
        # else:
            # raise ValueError('dims argument should be either 1 or 3')
        
        # Nc = self.Ip / (q_e * rho * self.k0 * speed_of_light)
        # z_sat = 3 + 1/np.sqrt(3) * np.log(Nc)
        # Psn = (3 * rho * self.Pb) / (Nc * np.sqrt(np.pi * np.log(Nc)))
        
        if z is None:
            zn = self.z_sat_min / (np.sqrt(3) * self.lg3)
        elif z == 0:
            return np.array(np.size(self.P_sn)*(np.NaN,))
        else:
            if np.size(z) > 1:
                z = z[:,np.newaxis]
                if (z > self.z_sat_min).any():
                    print('Warning, estimation applicable up to z_sat_min=%.2fm' %(self.z_sat_min))
            else: 
                if (z > self.z_sat_min):
                    print('Warning, estimation applicable up to z_sat_min=%.2fm, while z=%.2fm requested' %(self.z_sat_min, z))
            
            zn = z / (np.sqrt(3) * self.lg3)
            
        
        Pz = self.P_sn * (1 + 1/9 * np.exp(np.sqrt(3) * zn) / np.sqrt(np.pi * zn))
        # Pz = self.P_sn * (1 + 1/9 * np.exp(np.sqrt(3) * zn))
        #Pz = p.P_sn * (1 + 1/9 * np.exp(np.sqrt(3) * zn))
        return Pz
        
        
    def E(self, z=None):
        P = self.P(z)
        P[np.isnan(P)] = 0
        return np.trapz(P, self.s / speed_of_light)
        
    def tcoh(self,z=None):
        #check
        if z is None:
            z = self.z_sat_min
        elif z > self.z_sat_min:
            print('Warning, estimation applicable up to z_sat_min=%.2fm, while z=%.2fm requested' %(z_sat_min, z))
        tcoh = self.lambda0 / (6 * self.rho3 * speed_of_light ) * np.sqrt(z / (2 * np.pi * self.lg3))
        return tcoh
    
    def P_sat(self):
        return self.P(self.z_sat_min)
        
    @property
    def phen0(self):
        return 2 * np.pi / self.lambda0 * hr_eV_s * speed_of_light

class FelParametersArray(FelParameters):
    def __init__(self):
        super().__init__()
    
    @property
    def idx(self):
        try:
            idx = self.I.argmax()
        except AttributeError: 
            idx = None
        return idx
    
    

def calculateFelParameters(input, array=False):
    
    if array:
        p = FelParametersArray()
    else:
        p = FelParameters()
    
    p.iwityp = input.iwityp # undulator type: 0 == planar, other == helical
    
    p.gamma0 = input.gamma0  
    p.delgam = input.delgam  
    p.xlamd = input.xlamd    # undulator period
    
    p.betax = input.betax
    p.betay = input.betay
    p.emitx = input.emitx #normalized emittance
    p.emity = input.emity
    #    p.rxbeam = input.rxbeam
    #    p.rybeam = input.rybeam
    p.rxbeam = np.sqrt(p.betax * p.emitx / p.gamma0)
    p.rybeam = np.sqrt(p.betay * p.emity / p.gamma0)
    
    p.aw0 = input.aw0 # rms undulator parameter K
    p.I = input.curpeak
    
    p.deta =  p.delgam / p.gamma0
    
    p.lambda0 = p.xlamd / (2.0 * p.gamma0**2) * (1.0 + p.aw0**2) # resonant wavelength
    p.k0 = 2 * np.pi / p.lambda0 
    
    p.Ia = 4 * np.pi * epsilon_0 * m_e_eV * speed_of_light # Alfven (Budker) current (~17kA)
    #    p.Ia = 17000
        
    if p.iwityp == 0:
        ja = p.aw0**2 / (2*(1 + p.aw0**2))
        p.fc = sf.j0(ja) - sf.j1(ja)
    else:
        p.fc = 1.0
    
    # import first, ro_e * m_e_eV = 1.4399643147059695e-09
    
    # p.N = p.I * p.lambda0 / 1.4399644850445153e-10
    # p.sigb = 0.5 * (p.rxbeam + p.rybeam) # average beam size
    
    p.rho1 = (0.5 / p.gamma0) * np.power( (p.aw0 * p.fc * p.xlamd / 2 / np.pi )**2 / (p.rxbeam * p.rybeam) * p.I / p.Ia, 1.0/3.0) ## check 8 in denominator
    p.Pb = p.gamma0 * p.I * m_e_eV# beam power [Reiche]
    
    #p.power = 6.0 * np.sqrt(np.pi) * p.rho1**2 * p.Pb / (p.N * np.log(p.N / p.rho1) ) # shot noise power [W] [Reiche]
    p.lg1 = p.xlamd / (4*np.pi * np.sqrt(3) * p.rho1) #[Xie]
    p.zr = 4 * np.pi * p.rxbeam * p.rybeam / p.lambda0
      
    a = [None, 0.45, 0.57, 0.55, 1.6, 3.0, 2.0, 0.35, 2.9, 2.4, 51.0, 0.95, 3.0, 5.4, 0.7, 1.9, 1140.0, 2.2, 2.9, 3.2]
    
    p.xie_etad = p.lg1 / (2 * p.k0 * p.rxbeam * p.rybeam)
    #p.xie_etae = 4 * pi * p.lg1 / (p.betax*2*pi) * p.k0 * (p.emitx / p.gamma0)
    p.xie_etae = 4 * np.pi * p.lg1 * (p.emitx * p.emity) / p.lambda0 / (p.rxbeam * p.rybeam) / p.gamma0**2 # expressed via average x-y beam size
    p.xie_etagamma = p.deta / (p.rho1 * np.sqrt(3))
    p.xie_lscale = (a[1] * p.xie_etad ** a[2] + a[3] * p.xie_etae ** a[4] + a[5] * p.xie_etagamma ** a[6] 
    + a[7] * p.xie_etae ** a[8] * p.xie_etagamma ** a[9] + a[10] * p.xie_etad ** a[11] * p.xie_etagamma ** a[12] + a[13] * p.xie_etad ** a[14] * p.xie_etae ** a[15]
    + a[16] * p.xie_etad ** a[17] * p.xie_etae ** a[18] * p.xie_etagamma ** a[19])
    
    p.lg3 = p.lg1 * (1 + p.xie_lscale)
    p.rho3 = p.xlamd / (4*np.pi * np.sqrt(3) * p.lg3)
    
    p.Nc = p.I / (q_e * p.rho3 * p.k0 * speed_of_light)
    p.P_sn = (3 * p.rho1 * p.Pb) / (p.Nc * np.sqrt(np.pi * np.log(p.Nc))) # shot noise power [W]
    
    p.z_sat_norm = 3 + 1/np.sqrt(3) * np.log(p.Nc) # normalized saturation length for slices
    p.z_sat_magn = p.z_sat_norm * np.sqrt(3) * p.lg3 # magnetic length to reach saturation
    
    p.z_sat_min = np.nanmin(p.z_sat_magn)
    # try:
        # p.idx = p.I.argmax()
    # except AttributeError: 
        # p.idx = 0
    
    return p


def beam2fel(beam, lu, K_peak, iwityp=0):
    '''
    tmp function to estimate fel parameters slice-wise
    '''
    if beam.len() is 0:
        raise ValueError('Beam length shoul not be zero')
    
    class tmp():
        pass
    tmp.gamma0 = beam.g
    tmp.delgam = beam.dg
    tmp.xlamd = lu # undulator period
    tmp.iwityp = iwityp
    tmp.emitx = beam.emit_xn
    tmp.emity = beam.emit_yn
    if hasattr(beam,'beta_x_eff') and hasattr(beam,'beta_y_eff'):
        tmp.betax = beam.beta_x_eff
        tmp.betay = beam.beta_y_eff
    else:
        # print('use update_effective_beta() to increase estimation accuracy')
        tmp.betax = beam.beta_x
        tmp.betay = beam.beta_y
    if K_peak == 0:
        print('Warning, undulator K=0')
    
    if iwityp == 0: #planar
        tmp.aw0 = K_peak / np.sqrt(2)
    else: #other
        tmp.aw0 = K_peak
    
    tmp.curpeak = beam.I
    
    fel=calculateFelParameters(tmp, array=True)
    fel.s = beam.s
    return(fel)


def printFelParameters(p):
    
    #print (input.parameters)
    
    print ('********    FEL Parameters    ********')
    print ('ex=', p.emitx)
    print ('ey=', p.emity)
    print ('rxbeam=', p.rxbeam, ' [m]')
    print ('rybeam=', p.rybeam, ' [m]')
    print ('rel energy spread deta=', p.deta, ' [m]')
    print ('xlamd=', p.xlamd)
    print ('aw0=', p.aw0)
    print ('coupling parameter fc=', p.fc)
    print ('gamma0=', p.gamma0)
    print ('Ip=', p.I, ' beam peak current [A]')
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
    print ('')
    print ('Saturation magn. length=', p.z_sat_min)
    print ('**************************************')
    
    
# CHECK with Xie paper parameters
#inp = GenesisInput()
#inp.curpeak = 3400
#inp.xlamd = 0.03
#inp.iwityp = 0
#inp.gamma0 = 28000
#inp.delgam = inp.gamma0 * 2e-4
#inp.betax = 18
#inp.betay = 18
#inp.emitx=1.5e-6
#inp.emity=1.5e-6
#inp.xlamd=0.03
#inp.aw0 = 3.7/sqrt(2)
#
#p = calculateFelParameters(inp)
#print(p.xie_lscale,'new')
#p.lg1
#p.rho1
#print(p.xie_etad, 0.0367)
#print(p.xie_etae, 0.739)
#print(p.xie_etagamma, 0.248)