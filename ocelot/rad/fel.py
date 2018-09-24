from __future__ import print_function
'''
basic fel calculations
'''

#from pylab import *
import numpy as np
import numpy.fft as fft
import scipy.special as sf
from ocelot.common.globals import m_e_eV, epsilon_0, speed_of_light, q_e, h_eV_s
import logging
from scipy.optimize import fmin
from copy import deepcopy

_logger = logging.getLogger('ocelot.fel')
#from matplotlib.figure import Figure
#from mpl_toolkits.mplot3d import Axes3D

#import fel

class FelParameters:
    def __init__(self):
        pass
        
    def eval(self, method='mxie'):
        _logger.debug('Calculating FEL parameters')
        
        if self.betax <=0 or self.betay <=0:
            _logger.warning('betax or betay <= 0, returning lg3=np.nan')
            self.lg3 = np.nan
            return
        
        self.rxbeam = np.sqrt(self.betax * self.emitx / self.gamma0)
        self.rybeam = np.sqrt(self.betay * self.emity / self.gamma0)
        
        if self.rxbeam <=0 or self.rybeam <=0:
            _logger.warning('rxbeam or rybeam <= 0, returning lg3=np.nan')
            self.lg3 = np.nan
            return
        
        self.deta =  self.delgam / self.gamma0
        
        self.lambda0 = self.xlamd / (2.0 * self.gamma0**2) * (1.0 + self.aw0**2) # resonant wavelength
        self.k0 = 2 * np.pi / self.lambda0
        self.Ia = 4 * np.pi * epsilon_0 * m_e_eV * speed_of_light # Alfven (Budker) current (~17kA)
            
        if self.iwityp == 0: #planar undulator
            ja = self.aw0**2 / (2*(1 + self.aw0**2))
            self.fc = sf.j0(ja) - sf.j1(ja)
        else: #helical undulator
            self.fc = 1.0 
        
        self.Pb = self.gamma0 * self.I * m_e_eV# beam power [Reiche]
        
        # import first, ro_e * m_e_eV = 1.4399643147059695e-09
        
        # self.N = self.I * self.lambda0 / 1.4399644850445153e-10
        # self.sigb = 0.5 * (self.rxbeam + self.rybeam) # average beam size
        

        # h_eV_s * speed_of_light / self.lambda0
        
        self.rho1 = (0.5 / self.gamma0) * np.power( (self.aw0 * self.fc * self.xlamd / 2 / np.pi )**2 / (self.rxbeam * self.rybeam) * self.I / self.Ia, 1.0/3.0)
        
        #self.power = 6.0 * np.sqrt(np.pi) * self.rho1**2 * self.Pb / (self.N * np.log(self.N / self.rho1) ) # shot noise power [W] [Reiche]
        self.lg1 = self.xlamd / (4*np.pi * np.sqrt(3) * self.rho1) #[Xie]
        
        self.zr = 4 * np.pi * self.rxbeam * self.rybeam / self.lambda0
        
        a = [None, 0.45, 0.57, 0.55, 1.6, 3.0, 2.0, 0.35, 2.9, 2.4, 51.0, 0.95, 3.0, 5.4, 0.7, 1.9, 1140.0, 2.2, 2.9, 3.2]
        
        self.xie_etad = self.lg1 / (2 * self.k0 * self.rxbeam * self.rybeam)
        #self.xie_etae = 4 * pi * self.lg1 / (self.betax*2*pi) * self.k0 * (self.emitx / self.gamma0)
        self.xie_etae = 4 * np.pi * self.lg1 * (self.emitx * self.emity) / self.lambda0 / (self.rxbeam * self.rybeam) / self.gamma0**2 # expressed via average x-y beam size
        self.xie_etagamma = self.deta / (self.rho1 * np.sqrt(3))
        self.xie_lscale = (a[1] * self.xie_etad ** a[2] + a[3] * self.xie_etae ** a[4] + a[5] * self.xie_etagamma ** a[6] 
        + a[7] * self.xie_etae ** a[8] * self.xie_etagamma ** a[9] + a[10] * self.xie_etad ** a[11] * self.xie_etagamma ** a[12] + a[13] * self.xie_etad ** a[14] * self.xie_etae ** a[15]
        + a[16] * self.xie_etad ** a[17] * self.xie_etae ** a[18] * self.xie_etagamma ** a[19])
        
        self.lg3 = self.lg1 * (1 + self.xie_lscale)
        self.lg3 *= self.lg_mult
        if self.lg_mult != 1:
            _logger.info('lg3 multiplied by lg_mult ({})'.format(self.lg_mult))
        self.rho3 = self.xlamd / (4*np.pi * np.sqrt(3) * self.lg3)
        
        self.Nc = self.I / (q_e * self.rho3 * self.k0 * speed_of_light)
        # self.P_sn = (3 * self.rho1 * self.Pb) / (self.Nc * np.sqrt(np.pi * np.log(self.Nc))) # shot noise power [W]
        self.P_sn = (3 * self.rho3 * self.Pb) / (self.Nc * np.sqrt(np.pi * np.log(self.Nc))) # shot noise power [W]
        
        self.z_sat_norm = 3 + 1/np.sqrt(3) * np.log(self.Nc) # normalized saturation length for slices
        self.z_sat_magn = self.z_sat_norm * np.sqrt(3) * self.lg3 # magnetic length to reach saturation
        
        self.theta_c = np.sqrt(self.lambda0 / self.lg3) #critical angle
        # _logger.debug('L_sat_norm = {}'.format(self.z_sat_norm))
        
        self.z_sat_min = np.nanmin(self.z_sat_magn)
        
    def beta_opt(self, method='mxie', apply=False, **kwargs):
        beta_orig_x, beta_orig_y = self.betax, self.betay
        beta_orig = np.mean([beta_orig_x, beta_orig_y])
        
        fel_copy = deepcopy(self)
        def f(x, method=method):
            fel_copy.betax = fel_copy.betay = x
            fel_copy.eval(method=method)
            return fel_copy.lg3
        
        err_dict = np.geterr()
        np.seterr(all='ignore')
        beta_opt = fmin(f, beta_orig, disp=0, **kwargs)
        np.seterr(**err_dict)
        
        if apply:
            self.betax = beta_opt
            self.betay = beta_opt
            self.eval()
        else:
            return beta_opt[0]
        
    
    def log(self, type='debug'):
    
        if type is 'debug':
            _log_func = _logger.debug
        elif type is 'info':
            _log_func = _logger.info
        elif type is 'log':
            _log_func = _logger.log
        
        _log_func('undulator period = {}'.format(self.xlamd))
        _log_func('undulator K (rms) = {}'.format(self.aw0))
        if self.iwityp == 0:
            _log_func('undulator type - planar')
        else:
            _log_func('undulator type - helical')
        # _log_func('beam E GeV = {}'.format(beam.E))
        _log_func('beam gamma = {}'.format(self.gamma0))
        _log_func('beam dgamma= {}'.format(self.delgam))
        _log_func('beam current = {}'.format(self.I))
        _log_func('beam power = {}'.format(self.Pb))
        # _log_func('beam alphax = {}'.format(self.alphax))
        # _log_func('beam alphay = {}'.format(self.alphay))
        _log_func('beam betax = {}'.format(self.betax))
        _log_func('beam betay = {}'.format(self.betay))
        _log_func('beam emitx = {}'.format(self.emitx))
        _log_func('beam emity = {}'.format(self.emity))
        # _log_func('beam x = {}'.format(self.xbeam))
        # _log_func('beam y = {}'.format(self.ybeam))
        # _log_func('beam px = {}'.format(self.pxbeam))
        # _log_func('beam py = {}'.format(self.pybeam))
        _log_func('beam rx = {}'.format(self.rxbeam))
        _log_func('beam ry = {}'.format(self.rybeam))
        _log_func('')
        _log_func('Estimation results')
        _log_func('Rho 1D = {}'.format(self.rho1))
        _log_func('FEL_wavelength = {:.5e} m'.format(self.lambda0))
        _log_func('FEL_E_photon   = {} eV'.format(h_eV_s * speed_of_light / self.lambda0))
        _log_func('Lg  1D = {} m'.format(self.lg1))
        _log_func('Z_Rayl = {} m'.format(self.zr))
        _log_func('xie_eta_d = {}'.format(self.xie_etad))
        _log_func('xie_eta_e = {}'.format(self.xie_etae))
        _log_func('xie_eta_gamma = {}'.format(self.xie_etagamma))
        _log_func('xie_scaling_tot = {}'.format(self.xie_lscale))
        _log_func('Lg  3D = {}'.format(self.lg3))
        _log_func('Rho 3D = {}'.format(self.rho3))
        _log_func('P_shnoise = {}'.format(self.P_sn))
        _log_func('L_sat_magn = {}'.format(self.z_sat_magn))
        _log_func('L_sat_min = {}'.format(self.z_sat_min))
        _log_func('Theta_critical = {:.5e} rad'.format(self.theta_c))
    
    def P(self, z=None):
        '''
        returns sase power at distance z
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
                    _logger.warning('Estimation applicable up to z_sat_min=%.2fm, limiting power to saturation level' %(self.z_sat_min))
                    idx = z > self.z_sat_min[:,np.newaxis]
                    z[idx] = self.z_sat_min[:,np.newaxis][idx]
            else: 
                if (z > self.z_sat_min):
                    _logger.warning('Estimation applicable up to z_sat_min=%.2fm, while z=%.2fm requested, returning saturation power' %(self.z_sat_min, z))
                    z = self.z_sat_min
            
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
            _logger.warning('estimation applicable up to z_sat_min=%.2fm, while z=%.2fm requested' %(z_sat_min, z))
        tcoh = self.lambda0 / (6 * self.rho3 * speed_of_light ) * np.sqrt(z / (2 * np.pi * self.lg3))
        return tcoh
    
    def P_sat(self):
        return self.P(self.z_sat_min)
        
    @property
    def phen0(self):
        return h_eV_s * speed_of_light / self.lambda0

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
    
    

def calculateFelParameters(input, array=False, method='mxie'):
    
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
    
    p.lg_mult = 1
    if hasattr(input,'lg_mult'):
        if p.lg_mult is None:
            p.lg_mult = input.lg_mult
    #    p.rxbeam = input.rxbeam
    #    p.rybeam = input.rybeam

    p.aw0 = input.aw0 # rms undulator parameter K
    p.I = input.curpeak
    
    p.eval(method)
    
    if array:
        p.log('log')
    else:
        pass
        # p.log('debug')
    
    
    # if not array:
        
    # try:
        # p.idx = p.I.argmax()
    # except AttributeError: 
        # p.idx = 0
    
    return p


def beam2fel(beam, lu, K_peak, iwityp=0):
    '''
    tmp function to estimate fel parameters slice-wise
    '''
    if beam.len() == 0:
        raise ValueError('Beam length should not be zero')
    
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
    return (fel)


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
    print ('Pb= %.3e beam power [W]'%(p.Pb))
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