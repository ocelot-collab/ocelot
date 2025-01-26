"""
simple functions to estimate radiation properties.

Sergey Tomin.
"""

import scipy.special as sf
from scipy.special import jn
from ocelot.common.globals import *
from ocelot.common.ocelog import logging
_logger = logging.getLogger(__name__)

def lambda2eV(Lambda):
    Eph = h_eV_s*speed_of_light/Lambda
    return Eph

def eV2lambda(Ephoton):
    Lambda = h_eV_s*speed_of_light/Ephoton
    return Lambda
    
def THz2eV(THz):
    return np.array(THz) * h_eV_s / 1e-12

def eV2THz(Ephoton):
    return np.array(Ephoton) / h_eV_s * 1e-12

def k2lambda(k):
    Lambda = 2 * np.pi / k
    return Lambda
    
def lambda2k(Lambda):
    k = 2 * np.pi / Lambda
    return k

def k2angle(k, Lambda):
    angle = k * Lambda / 2 / np.pi
    return angle

def angle2k(angle, Lambda):
    k =  2 * np.pi * angle / Lambda
    return k

def Ephoton2K(Eph, lu=0.04, Eeb=14):
    gamma = Eeb/m_e_GeV
    K = np.sqrt(4.*gamma*gamma/lu*h_eV_s*speed_of_light/Eph - 2)
    return K

def K2Ephoton(K, lu = 0.04, E=14):
    gamma = E/m_e_GeV
    Eph = 4.*gamma*gamma*h_eV_s*speed_of_light/((K*K + 2)*lu)
    return Eph

def K2Lambda(K, lu = 0.04, E=14):
    gamma = E/m_e_GeV
    Eph = 4.*gamma*gamma*h_eV_s*speed_of_light/((K*K + 2)*lu)
    Lambda = eV2lambda(Eph)
    return Lambda

def field2K(field, lu=0.04):
    """
    function to translate magnetic field in [T] to K parameter
    :param field: in [T]
    :param lu: undulator period in [m]
    :return: K
    """
    K = field*lu*speed_of_light/(m_e_eV*2.*pi)
    return K

def K2field(K, lu = 0.04):
    field = K*m_e_eV*2.*pi/(lu*speed_of_light)
    return field
# print (K2field(27, lu = 0.25) )
# print(K2Lambda(18.8, lu = 0.25, E=2.4))
def field2Ephoton(field, lu = 0.04, E=14):
    gamma = E/m_e_GeV
    K = field2K(field, lu)
    l = lu/(2.*gamma*gamma)*(1.+K*K/2.)
    Eph = h_eV_s*speed_of_light/l
    return Eph

def Ephoton2field(energy, lu = 0.04, Eeb = 14):
    K = Ephoton2K(energy, lu, Eeb)
    field = K*2.*pi*m_e_eV/(lu*speed_of_light)
    return field

def lambda2Ebeam(Lambda = 10.4e-9, lu=0.0272, K=1.2392):
    gamma = np.sqrt(lu/(2.*Lambda)*(1. + K*K/2.))
    return gamma*m_e_GeV

def Br2B0(Br=1.1, gap=0.011, lu=0.04):
    '''
    On-axis field from remanent field calculation
    Based on eq 1 of https://accelconf.web.cern.ch/ipac2016/papers/tupmb011.pdf
    '''
    a = 0.55 * Br + 2.835
    b=-1.95 * Br + 7.225
    c=-1.3 * Br + 2.978
    B0 = a * np.exp(-b*gap/lu + c*(gap/lu)**2)
    if (gap/lu) > 0.7 or (gap/lu) < 0.07:
        _logger.warning('gap/period {:.4f} ourside of fit applicability [0.7 : 0.07]'.format(gap/lu))
    return B0

class ID_radiation:
    def __init__(self, beam, undulator):
        if beam.E == 0:
            exit("electron beam energy must be non zero!")
        if beam.I == 0:
            exit("electron beam current must be non zero!")
        try:
            if beam.sigma_x == 0 or beam.sigma_y == 0:
                beam.sizes()
        except:
            beam.sizes()
        self.beam = beam
        self.undul = undulator

        #self.distance = distance

    def f_n(self, nharm, K):
        v1 = ((nharm-1)/2.)
        v2 = ((nharm+1)/2.)
        x = nharm*K*K/(4.+2.*K*K)
        return nharm*nharm*K*K/((1+K*K/2.)**2)*(jn(v1,x) - jn(v2,x))**2

    def flux(self,current, K, nharm, energy, L, lperiod):
        alpha = 1/137.036
        gm = energy/m_e_GeV
        Nu = L/lperiod
        e = 1.602e-19
        BW = 0.001 # band width 0.1%
        k_rad_mrad = 1e-6 # coef to converse rad to mrad
        F = alpha*gm*gm*current/e*Nu*Nu*self.f_n(nharm, K)*BW*k_rad_mrad
        return F

    def eff_sizes(self, K):
        #self.beam.sizes()
        Lambda = K2Lambda(K, self.undul.lperiod, self.beam.E)
        L = self.undul.l
        self.sigma_r = np.sqrt(Lambda*L/(2*4.*pi*pi))
        self.sigma_r1 = np.sqrt(Lambda/L/2.)
        self.Sigma_x = np.sqrt(self.beam.sigma_x**2 + self.sigma_r**2)
        self.Sigma_y = np.sqrt(self.beam.sigma_y**2 + self.sigma_r**2)
        self.Sigma_x1 = np.sqrt(self.beam.sigma_xp**2 + self.sigma_r1**2)
        self.Sigma_y1 = np.sqrt(self.beam.sigma_yp**2 + self.sigma_r1**2)
        #self.size_x = sqrt(self.Sigma_x**2 + (self.Sigma_x1*self.distance)**2)
        #self.size_y = sqrt(self.Sigma_y**2 + (self.Sigma_y1*self.distance)**2)

    def Flux(self, K, nharm = 1):
        current = self.beam.I
        energy = self.beam.E
        L = self.undul.l
        lperiod = self.undul.lperiod

        return self.flux(current, K, nharm, energy, L, lperiod)

    def Flux_tot(self, K, nharm):
        current = self.beam.I
        N = self.undul.nperiods
        flux_tot = 1.431e14*current*N*self.f_n(nharm, K)*(1.+K*K/2.)/1./2.
        return flux_tot

    def Brightness(self, K,nharm):
        flux_tot = self.Flux_tot(K, nharm)
        self.eff_sizes(K)
        brightness = flux_tot/(4*pi*pi*self.Sigma_x*self.Sigma_y*self.Sigma_x1*self.Sigma_y1)*1e-12
        return brightness

def print_rad_props(beam, K, lu, L, distance):
    """
    Simple function to estimate radiation properties and print them

    :param beam: Beam
    :param K: undulator deflection parameter
    :param lu: undulator period in [m]
    :param L: undulator length in [m]
    :param distance: distance to the screen im [m]
    :return:
    """
    beam.sizes()


    def f_n(n, Ku):
        v1 = ((n-1)/2.)
        v2 = ((n+1)/2.)
        x = n*Ku*Ku/(4.+2.*Ku*Ku)
        return n*n*Ku*Ku/((1+Ku*Ku/2.)**2)*(jn(v1,x) - jn(v2,x))**2

    def flux(I, K, m):
        alpha = 1/137.036
        gm = beam.E/m_e_GeV
        Nu = L/lu
        e = 1.602e-19
        BW = 0.001 # band width 0.1%
        k_rad_mrad = 1e-6 # coef to converse rad to mrad
        F = alpha*gm*gm*I/e*Nu*Nu*f_n(m, K)*BW*k_rad_mrad
        return F
    gamma = beam.E/m_e_GeV
    Lambda = K2Lambda(K, lu, beam.E)
    sigma_r = np.sqrt(Lambda*L/(2*4.*pi*pi))
    sigma_r1 = np.sqrt(Lambda/L/2.)
    Sigma_x = np.sqrt(beam.sigma_x**2 + sigma_r**2)
    Sigma_y = np.sqrt(beam.sigma_y**2 + sigma_r**2)
    Sigma_x1 = np.sqrt(beam.sigma_xp**2 + sigma_r1**2)
    Sigma_y1 = np.sqrt(beam.sigma_yp**2 + sigma_r1**2)
    size_x = np.sqrt(Sigma_x**2 + (Sigma_x1*distance)**2)
    size_y = np.sqrt(Sigma_y**2 + (Sigma_y1*distance)**2)
    B = K2field(K, lu = lu)
    F = flux(beam.I, K, m = 1)
    N = L/lu
    flux_tot = 1.431e14*beam.I*N*f_n(1, K)*(1.+K*K/2.)/1./2.

    brightness = flux_tot/(4*pi*pi*Sigma_x*Sigma_y*Sigma_x1*Sigma_y1)*1e-12

    print ("********* ph beam ***********")
    print ("Ebeam        : ", beam.E, " GeV")
    print ("K            : ", K)
    print ("B            : ", np.round(B, 4), " T")
    print ("lambda       : ", "{0:.5E}".format(Lambda), " m ")
    print ("Eph          : ", "{0:.5E}".format(K2Ephoton(K, lu, beam.E)), " eV")
    print ("1/gamma      : ", np.round(1./gamma *1e6, 4), " um")
    print ("sigma_r      : ", np.round(sigma_r *1e6, 4), " um")
    print ("sigma_r'     : ", np.round(sigma_r1*1e6, 4), " urad")
    print ("Sigma_x      : ", np.round(Sigma_x *1e6, 4), " um")
    print ("Sigma_y      : ", np.round(Sigma_y *1e6, 4), " um")
    print ("Sigma_x'     : ", np.round(Sigma_x1*1e6, 4), "urad")
    print ("Sigma_y'     : ", np.round(Sigma_y1*1e6, 4), "urad")
    print ("H. spot size : ", np.round(size_x*1000., 4), "/", np.round(size_x/distance*1000., 4), " mm/mrad")
    print ("V. spot size : ", np.round(size_y*1000., 4), "/", np.round(size_y/distance*1000., 4), " mm/mrad")
    print ("I            : ", beam.I, " A")
    print ("Nperiods     : ", L/lu)
    print ("distance     : ", distance, " m")
    print ("flux tot     : ", "{0:.2E}".format(flux_tot), " ph/sec/0.1%BW")
    print ("flux density : ", "{0:.2E}".format(F), " ph/sec/mrad^2/0.1%BW;   ", "{0:.2E}".format(F/distance/distance), " ph/sec/mm^2/0.1%BW")
    #print "flux density : ", F/distance/distance, " ph/sec/mm^2"
    print ("brilliance   : ", "{0:.2E}".format(brightness), " ph/sec/mrad^2/mm^2/0.1%BW")



class UndulatorParameters:

    def __init__(self, und = None,el_E=1.0):

        if und == None:
            self.Nw = 100   # number of periods
            self.lw = 0.045 # period length, m
            self.E = 1.0 # electron energy, in GeV
            self.K = 1.0    # undulator parameter (peak)

        else:

            self.Nw = und.nperiods   # number of periods
            self.lw = und.lperiod # period length, m
            self.E = el_E # electron energy, in GeV
            self.K = und.Kx    # undulator parameter (peak)

        self.recalculate()


    def recalculate(self):

        self.gamma = 1.0e+9 * self.E / m_e_eV
        self.L = self.Nw*self.lw
        self.B = 2.0*pi * self.K * m_e_eV / (self.lw * speed_of_light)
        self.kw = 2.0*pi / self.lw
        self.beta = np.sqrt(1.0 - 1.0 / (self.gamma*self.gamma) )
        self.beta_z = self.beta * (1.0 - (self.K*self.K)/ (4.0*self.gamma*self.gamma))
        #self.w1, _, _ = self.computeRadiationAnalytical(200.0, 0.0, 0)
        self.w1 = 1.0 / ( (1.0 + self.K*self.K/2.0 ) / (2*self.kw*speed_of_light*self.gamma**2) )
        self.lamd_ev = h_eV_s*self.w1/(2*np.pi)


    def computeRadiationAnalytical(self,z0, theta, phi, dwRel = 0.1, npoints=100):
        w1 = 1.0 / ( (1.0 + self.K*self.K/2.0 + self.gamma*self.gamma*theta) / (2*self.kw*speed_of_light*self.gamma**2) ) # first harmonic
        dw = w1 * dwRel
        step_w = dw / npoints
        w = np.arange(w1-dw, w1+dw, step_w)
        phis = theta * theta * w * z0 / (2.0 * speed_of_light)
        delta_w = w - w1

        u = w * self.K*self.K / (8*self.kw*speed_of_light*self.gamma*self.gamma)
        ajj = sf.jv(0,u) - sf.jv(1,u)
        I = self.Nw * self.lw * self.K * w *  np.exp(1j * phis) / ( z0 * self.gamma)  * np.sinc(pi*self.Nw*(w-w1)/w1) * ajj
        I = np.real(I)

        self.fundamentalWavelength = self.lw / (2*self.gamma**2) * (1+ self.K**2 / 2.0 + self.gamma**2 * theta)

        print( 'test', h_eV_s*speed_of_light / self.fundamentalWavelength)

        return w1, w, I

    def get_k(self, E):
        # w1 in ev
        w1 = 2. * np.pi * E / h_eV_s
        return np.sqrt( (4.*self.kw*speed_of_light*self.gamma**2)/w1 - 2. )


    def printParameters(self):

        self.recalculate()

        print( "Undulator parameters:")
        print( "L=", self.L)
        print( "gamma(electron)=", self.gamma)
        print( "K=", self.K)
        print( "B[T]=", self.B)

        w1, _, _ = self.computeRadiationAnalytical(200.0, 0.0, 0)

        print ("Radiation parameters:")
        print ("w1(first harmonic, zero angle)=", w1, "Hz/2pi", h_eV_s*w1/(2*np.pi), "[eV]", self.fundamentalWavelength, "m")

        t = self.L / speed_of_light
        cb = 379.35 # 1/(gev sec Tesla^2)

        eloss = self.E - self.E / (1.0 + 0.5*self.B*self.B*cb*self.E*t)

        print( "Total energy loss [Gev]", eloss)
