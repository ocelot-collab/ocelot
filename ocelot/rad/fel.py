from __future__ import print_function
'''
basic fel calculations
'''

#from pylab import *
import numpy as np
import numpy.fft as fft
import scipy.special as sf
from ocelot.common.globals import m_e_eV, epsilon_0, speed_of_light, q_e, h_eV_s, lambda_C_r, I_Alfven, ro_e
import logging
from scipy.optimize import fmin
from copy import deepcopy

_logger = logging.getLogger(__name__)
#from matplotlib.figure import Figure
#from mpl_toolkits.mplot3d import Axes3D

#import fel

class FelParameters:
    def __init__(self):
        self.qf = 0
        self.inaccurate = False # True if fitting formulas do not promise good accuracy
        pass
        
    def eval(self, method='mxie'):
        _logger.debug('Calculating FEL parameters')
        
        if np.size(self.I) > 1:
            tdp=True
        else:
            tdp=False
        
        if not hasattr(self, 'hn'):
            self.hn=1 #harmonic number
        
        if np.any(self.betax <= 0) or np.any(self.betay <= 0):
            _logger.warning('betax or betay <= 0, returning lg3=np.nan')
            self.lg3 = np.nan
            return
        
        self.rxbeam = np.sqrt(self.betax * self.emitx / self.gamma0)
        self.rybeam = np.sqrt(self.betay * self.emity / self.gamma0)
        
        if np.any(self.rxbeam <= 0) or np.any(self.rybeam <= 0):
            _logger.warning('rxbeam or rybeam <= 0, returning lg3=np.nan')
            self.lg3 = np.nan
            return
        
        self.deta =  self.delgam / self.gamma0
        
        if np.isnan(self.aw0):
            _logger.warning('aw0 is nan')
            self.inaccurate = True 
            
        self.lambda0 = self.xlamd / (2.0 * self.gamma0**2) * (1.0 + self.aw0**2) # resonant wavelength
        if np.any(self.lambda0 < 0):
            _logger.error('wavelength is not reachable with und_period {} gamma {} and K {}'.format(self.xlamd,self.gamma0,self.aw0))
            self.inaccurate = True 
        self.lambdah = self.lambda0 / self.hn
        self.k0 = 2 * np.pi / self.lambda0
        # self.Ia = I_Alfven #remove
            
        if self.iwityp == 0: #planar undulator
            ja = self.aw0**2 / (2*(1 + self.aw0**2))
            self.fc = sf.j0(ja) - sf.j1(ja)
            # if self.hn != 1:
            jah = self.hn * self.aw0**2 / (2*(1 + self.aw0**2))
            self.fch = sf.jv((self.hn-1)/2, jah) - sf.jv((self.hn+1)/2, jah)
        else: #helical undulator
            self.fc = 1
            if self.hn !=1:
                _logger.warning('harmonic number != 1 and undulator is helical. Not implemented! Retunrning zero coupling at harmonic!')
                self.inaccurate = True 
                self.fch = 0
            else:
                self.fch = 1
        
        self.Pb = self.gamma0 * self.I * m_e_eV# beam power [Reiche]
        
        # import first, ro_e * m_e_eV = 1.4399643147059695e-09
        
        # self.N = self.I * self.lambda0 / 1.4399644850445153e-10
        # self.sigb = 0.5 * (self.rxbeam + self.rybeam) # average beam size
        
        emit_n = np.sqrt(self.emitx * self.emity)
        # h_eV_s * speed_of_light / self.lambda0
        self.emit_nn = 2 * np.pi * emit_n / self.lambdah / self.gamma0 ## emittance normalized as in Eq.6, 10.1103/PhysRevSTAB.15.080702
        
        if (np.any(self.emit_nn < 1) or np.any(self.emit_nn) > 5):
            self.inaccurate = True 
            if tdp:
                _logger.warning('1 <! min(emittance) {} <! 5, SSY approx. might be incorrect'.format(np.nanmin(self.emit_nn)))
            else:
                _logger.warning('1 <! emittance {} <! 5, SSY approx. might be incorrect'.format(self.emit_nn))
            #Eq.6, DOI:10.1103/PhysRevSTAB.15.080702
        
        if self.qf == 1: #account for quantum fluctuations
            if self.iwityp == 0: #planar undulator
                F_aw = 1.7 * self.aw0 + 1 / (1 + 1.88 * self.aw0 + 0.8 * self.aw0**2) 
                #eq.B2, DOI:10.1103/PhysRevSTAB.15.080702, 
                #eq.11 DOI:10.1016/j.optcom.2004.02.071
            else: #helical undulator
                F_aw = 1.42 * self.aw0 + 1 / (1 + 1.5 * self.aw0 + 0.95 * self.aw0**2)
        
        if method == 'mxie':
            '''
            M. Xie, “Exact and variational solutions of 3D eigenmodes in high gain FELs,” Nucl. Instruments Methods Phys. Res. Sect. A Accel. Spectrometers, Detect. Assoc. Equip., vol. 445, no. 1–3, pp. 59–66, 2000.
            '''
            
            # if self.hn != 1:
                # _logger.warning('MXie estimation not implemented for harmonic radaition')
            
            self.rho1 = (0.5 / self.gamma0) * np.power( (self.aw0 * self.fc * self.xlamd / 2 / np.pi )**2 / (self.rxbeam * self.rybeam) * self.I / I_Alfven, 1.0/3.0) 
            
            #self.power = 6.0 * np.sqrt(np.pi) * self.rho1**2 * self.Pb / (self.N * np.log(self.N / self.rho1) ) # shot noise power [W] [Reiche]
            self.lg1 = self.xlamd / (4*np.pi * np.sqrt(3) * self.rho1) #power gain length [Xie]
            
            self.zr = 4 * np.pi * self.rxbeam * self.rybeam / self.lambda0
            
            a = [None, 0.45, 0.57, 0.55, 1.6, 3.0, 2.0, 0.35, 2.9, 2.4, 51.0, 0.95, 3.0, 5.4, 0.7, 1.9, 1140.0, 2.2, 2.9, 3.2]
            
            self.xie_etad = self.lg1 / (2 * self.k0 * self.rxbeam * self.rybeam)
            
            #self.xie_etae = 4 * pi * self.lg1 / (self.betax*2*pi) * self.k0 * (self.emitx / self.gamma0)
            self.xie_etae = 4 * np.pi * self.lg1 * (self.emitx * self.emity) / self.lambda0 / (self.rxbeam * self.rybeam) / self.gamma0**2 * ((self.fc/self.fch)**2 / self.hn)**(1/3) / self.hn # expressed via average x-y beam size
            self.xie_etagamma = self.deta / (self.rho1 * np.sqrt(3))
            
            if self.hn !=1:
                self.xie_etad *= ((self.fc/self.fch)**2 / self.hn)**(1/3) / self.hn 
                self.xie_etae *= ((self.fc/self.fch)**2 / self.hn)**(1/3) * self.hn
                self.xie_etagamma *= ((self.fc/self.fch)**2 / self.hn)**(1/3) * self.hn #eq C2+ DOI:10.1103/PhysRevSTAB.15.080702
            
            self.delta = (a[1] * self.xie_etad ** a[2] + a[3] * self.xie_etae ** a[4] + a[5] * self.xie_etagamma ** a[6] 
            + a[7] * self.xie_etae ** a[8] * self.xie_etagamma ** a[9] + a[10] * self.xie_etad ** a[11] * self.xie_etagamma ** a[12] + a[13] * self.xie_etad ** a[14] * self.xie_etae ** a[15]
            + a[16] * self.xie_etad ** a[17] * self.xie_etae ** a[18] * self.xie_etagamma ** a[19])
            
            # self.lg3 = self.lg1 * (1 + self.xie_lscale)
            self.method = 'mxie'
            
        elif method == 'ssy_opt':
            '''
            E. L. Saldin, E. A. Schneidmiller, and M. V. Yurkov, “Design formulas for short-wavelength FELs,” Opt. Commun., vol. 235, no. 4–6, pp. 415–420, May 2004.
            '''
            
            
            self.lg1 = 0.5 * 1.67 * np.sqrt(I_Alfven / self.I) * (emit_n * self.xlamd)**(5/6) / self.lambdah**(2/3) * (1 + self.aw0**2)**(1/3) / (self.hn**(5/6) * self.aw0 * self.fch) 
            #eq.4, DOI:10.1103/PhysRevSTAB.15.080702
            # it is power gain length = 0.5 * field gain length
            
            self.delta = 131 * (I_Alfven / self.I) * emit_n**(5/4) / (self.lambdah * self.xlamd**9)**(1/8) * self.hn**(9/8) * self.delgam**2 / (self.aw0 * self.fch)**2 / (1 + self.aw0**2)**(1/8) #eq.5, DOI:10.1103/PhysRevSTAB.15.080702
            
                
            
            # if hasattr(self, 'qf'):
                # if self.qf==1:
                    # self.lg3 = self.lg1 * (1 + self.delta_eff)
            
            self.method = 'ssy_opt'
        else:
            _logger.error('method should be in ["mxie", "ssy_opt"]')
            raise ValueError('method should be in ["mxie", "ssy_opt"]')

        if self.qf == 1:
            self.delta_q = 5.5e4 * (I_Alfven / self.I)**(3/2) * lambda_C_r * ro_e * emit_n**2 / self.lambda0**(11/4) / self.xlamd**(5/4) * (1 + self.aw0**2)**(9/4) * F_aw / (self.aw0 * self.fch**3 * self.hn**(5/3))
            
            if np.any(self.delta_q >= 1):
                _logger.warning('quantum fluctuation effect exceeds 1, estimation not applicable anymore')
                self.delta_q = 0.999
                self.inaccurate = True
        else:
            self.delta_q = 0
            
        self.delta_eff = (self.delta + self.delta_q) / (1 - self.delta_q)
        
        self.delta_criterion = 2.5 * (1 - np.exp(-0.5 * self.emit_nn**2))
        if np.any(self.delta_eff > self.delta_criterion):
            if tdp:
                _logger.warning('delta_eff > delta_criterion; SSY approx. might be incorrect')
            else:
                _logger.warning('delta_eff {} > {}; SSY approx. might be incorrect'.format(self.delta_eff, self.delta_criterion))
            self.inaccurate = True
                #Eq.7, DOI:10.1103/PhysRevSTAB.15.080702
                #Eq.14+text, DOI:10.1016/j.optcom.2004.02.071
        
        self.beta_opt_calc = 11.2 * (I_Alfven / self.I)**(1/2) * (emit_n**3 * self.xlamd)**(1/2) / (self.lambdah* self.hn**(1/2) * self.aw0 * self.fch) / (1 + 8 * self.delta_eff)**(1/3)
        
        self.lg3 = self.lg1 * (1 + self.delta_eff)
        
        self.lg3 *= self.Lg_mult
        if self.Lg_mult != 1:
            _logger.info('lg3 multiplied by Lg_mult ({})'.format(self.Lg_mult))
        self.rho3 = self.xlamd / (4*np.pi * np.sqrt(3) * self.lg3)
        
        self.Nc = self.I / (q_e * self.rho3 * self.k0 * speed_of_light)
        # self.P_sn = (3 * self.rho1 * self.Pb) / (self.Nc * np.sqrt(np.pi * np.log(self.Nc))) # shot noise power [W]
        self.P_sn = (3 * self.rho3 * self.Pb) / (self.Nc * np.sqrt(np.pi * np.log(self.Nc))) # shot noise power [W]
        
        self.z_sat_norm = 3 + 1/np.sqrt(3) * np.log(self.Nc) # normalized saturation length for slices
        self.z_sat_magn = self.z_sat_norm * np.sqrt(3) * self.lg3 # magnetic length to reach saturation
        
        self.theta_c = np.sqrt(self.lambdah / self.lg3) #critical angle
        # _logger.debug('L_sat_norm = {}'.format(self.z_sat_norm))
        
        self.z_sat_min = np.nanmin(self.z_sat_magn)
        
    def beta_opt(self, method='mxie', apply=False, **kwargs):
        if method == 'mxie':
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
            
        elif method == 'ssy_opt':
            beta_opt = self.beta_opt_calc
        
        else:
            _logger.error('method should be in ["mxie", "ssy_opt"]')
            raise ValueError('method should be in ["mxie", "ssy_opt"]')
            
        if apply:
            self.betax = beta_opt
            self.betay = beta_opt
            self.eval(method)
        else:
            return beta_opt[0]
            
    
    def log(self, type='debug'):
    
        if type == 'debug':
            _log_func = _logger.debug
        elif type == 'info':
            _log_func = _logger.info
        elif type == 'log':
            _log_func = _logger.log
        elif type == 'print':
            _log_func = print
        
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
        _log_func('beam emitx_norm = {}'.format(self.emitx))
        _log_func('beam emity_norm = {}'.format(self.emity))
        # _log_func('beam x = {}'.format(self.xbeam))
        # _log_func('beam y = {}'.format(self.ybeam))
        # _log_func('beam px = {}'.format(self.pxbeam))
        # _log_func('beam py = {}'.format(self.pybeam))
        _log_func('beam rx = {}'.format(self.rxbeam))
        _log_func('beam ry = {}'.format(self.rybeam))
        _log_func('')
        _log_func('Estimation results')
        _log_func('Rho 1D = {}'.format(self.rho1))
        _log_func('FEL_wavelength = {} m'.format(self.lambda0))
        _log_func('FEL_E_photon   = {} eV'.format(h_eV_s * speed_of_light / self.lambda0))
        _log_func('Lg  1D = {} m'.format(self.lg1))
        _log_func('Z_Rayl = {} m'.format(self.zr))
        _log_func('xie_eta_d = {}'.format(self.xie_etad))
        _log_func('xie_eta_e = {}'.format(self.xie_etae))
        _log_func('xie_eta_gamma = {}'.format(self.xie_etagamma))
        # _log_func('xie_scaling_tot = {}'.format(self.xie_lscale))
        _log_func('Lg  3D = {}'.format(self.lg3))
        _log_func('Rho 3D = {}'.format(self.rho3))
        _log_func('P_shnoise = {}'.format(self.P_sn))
        _log_func('L_sat_magn = {}'.format(self.z_sat_magn))
        _log_func('L_sat_min = {}'.format(self.z_sat_min))
        _log_func('Theta_critical = {} rad'.format(self.theta_c))
    
    def P(self, z=None):
        '''
        returns sase power at distance z
        unfinished
        '''
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
        if hasattr(self,'P_mult'):
            if self.P_mult is not None:
                Pz *= self.P_mult
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
    
    @property
    def phenh(self):
        return h_eV_s * speed_of_light / self.lambdah
    
    def spectrogram(self, z=None):
        #fast spectrogram evaluation
        if z is None:
            z = self.z_sat_min
        Psat = self.P(z)
        Psat[np.isnan(Psat)]=0
        idx = self.idx

        phen0 = self.phen0
        dphen = phen0 * self.rho3
        dp = dphen[idx] / 10
        s_arr = self.s
        
        phen_arr = np.arange(np.amin(phen0 - 3 * dphen), np.amax(phen0 + 3 * dphen), dp)
        spec = np.zeros((s_arr.size, phen_arr.size))
        for i in range(s_arr.size):
            if dphen[i] != 0:
                spec[i] = np.exp(-(phen_arr - phen0[i])**2 / 2 / dphen[i]**2) / np.sqrt(2 * np.pi * dphen[i]**2)
        spec = spec * Psat[:, np.newaxis]
        
        return (s_arr, phen_arr, spec.T)
        
    def spectrum(self, z=None):
        #fast total spectrum evaluation
        s_arr, phen_arr, spectrogram = self.spectrogram(z = z)
        spectrum = np.sum(spectrogram, axis=1)
        return phen_arr, spectrum

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
    
    if hasattr(input,'hn'):
        p.hn = input.hn
    
    if hasattr(input,'qf'):
        p.qf = input.qf
    
    p.Lg_mult = 1
    if hasattr(input,'Lg_mult'):
        if input.Lg_mult is not None:
            p.Lg_mult = input.Lg_mult
    p.P_mult = 1
    if hasattr(input,'P_mult'):
        if input.P_mult is not None:
            p.P_mult = input.P_mult
    #    p.rxbeam = input.rxbeam
    #    p.rybeam = input.rybeam

    p.aw0 = input.aw0 # rms undulator parameter K
    p.I = input.curpeak
    
    p.eval(method)
    
    # if array:
        # p.log('log')
    # else:
        # pass
        # p.log('debug')
    
    
    # if not array:
        
    # try:
        # p.idx = p.I.argmax()
    # except AttributeError: 
        # p.idx = 0
    
    return p


def beam2fel(beam, lu, K_peak, iwityp=0, method='mxie', hn=1, qf=0):
    '''
    tmp function to estimate fel parameters slice-wise
    hn = harmonic number
    qf = account for quantum fluctuations
    '''
    if beam.len() == 0:
        raise ValueError('Beam length should not be zero')
    
    class tmp():
        pass
    tmp.hn=hn
    tmp.qf=qf
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
    
    fel=calculateFelParameters(tmp, array=True, method=method)
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
    # print ('gain length degradation lscale=', p.xie_lscale)
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