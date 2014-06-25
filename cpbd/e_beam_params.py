__author__ = 'Sergey'

from numpy import linspace, array, sqrt, pi
from scipy.integrate import simps
from ocelot.common.globals import *
def radiation_integral(lattice, twiss_0, nsuperperiod = 1):

    tws_elem = twiss_0
    (I1, I2, I3,I4, I5) = (0., 0., 0., 0., 0.)

    for elem in lattice.sequence:
        if elem.type == "rbend" or elem.type == "sbend" or elem.type == "bend" or elem.type == "quadrupole":
            Dx = []
            Hinvariant = []
            k = []
            h = []
            Z = []
            for z in linspace(0, elem.l,num = 10, endpoint=True):
                tws_z = elem.transfer_map(z)*tws_elem
                Dx.append(tws_z.Dx)
                k.append(elem.k1)
                if  elem.type != "quadrupole" or elem.l == 0:
                    h.append(elem.angle/elem.l)
                else:
                    h.append(0.)
                Z.append(z)
                Hinvariant.append(tws_z.gamma_x*tws_z.Dx*tws_z.Dx + 2.*tws_z.alpha_x*tws_z.Dxp*tws_z.Dx
                                        + tws_z.beta_x*tws_z.Dxp*tws_z.Dxp)
            H = array(h)
            H2 = H*H
            H3 = abs(H*H*H)
            I1 += simps(array(Dx)*H, Z)*nsuperperiod
            I2 += simps(H2, Z)*nsuperperiod
            I3 += simps(H3, Z)*nsuperperiod
            I4 += simps(array(Dx)*H*(2*array(k)+H2), Z)*nsuperperiod
            I5 += simps(array(Hinvariant)*H3, Z)*nsuperperiod
        tws_elem = elem.transfer_map*tws_elem

    return (I1,I2,I3, I4, I5)

class EbeamParams:
    def __init__(self, lattice, twiss_0, coupling = 0.01, nsuperperiod = 1):
        self.tws0 = twiss_0
        (I1,I2,I3, I4, I5) = radiation_integral(lattice, twiss_0, nsuperperiod)
        print "I2 = ", I2
        print "I3 = ", I3
        print "I4 = ", I4
        print "I5 = ", I5
        self.Je = 2 + I4/I2
        self.Jx = 1 - I4/I2
        self.Jy = 1
        self.gamma = twiss_0.E/m_e_GeV
        self.sigma_e = self.gamma*sqrt(Cq*I3/(self.Je*I2))
        self.emittance = Cq*self.gamma*self.gamma*I5/(self.Jx*I2)
        self.U0 = Cgamma*(twiss_0.E*1000)**4*I2/(2*pi)
        #print "*********  ", twiss_0.Energy
        self.Tperiod = nsuperperiod*lattice.totalLen/speed_of_light

        self.tau0 = 2*twiss_0.E*1000*self.Tperiod/self.U0
        self.tau_e = self.tau0/self.Je
        self.tau_x = self.tau0/self.Jx
        self.tau_y = self.tau0/self.Jy
        self.alpha = I1/(speed_of_light*self.Tperiod)
        self.coupl = coupling
        self.emitt_x = self.emittance/(1 + self.coupl)
        self.emitt_y = self.emittance*self.coupl/(1 + self.coupl)
        self.sigma_x = sqrt((self.sigma_e*self.tws0.Dx)**2 + self.emitt_x*self.tws0.beta_x)
        self.sigma_y = sqrt((self.sigma_e*self.tws0.Dy)**2 + self.emitt_y*self.tws0.beta_y)
        self.sigma_xp = sqrt((self.sigma_e*self.tws0.Dxp)**2 + self.emitt_x*self.tws0.gamma_x)
        self.sigma_yp = sqrt((self.sigma_e*self.tws0.Dyp)**2 + self.emitt_y*self.tws0.gamma_y)


    def print_params(self):
        print "Je =        ", self.Je
        print "Jx =        ", self.Jx
        print "Jy =        ", self.Jy
        print "gamma =     ", self.gamma
        print "sigma_e =   ", self.sigma_e
        print "emittance = ", self.emittance*1e9, " nm*rad"
        print "U0 =        ", self.U0, "  MeV"
        print "Tperiod =   ", self.Tperiod*1e9, " nsec"
        print "alpha =     ", self.alpha
        print "tau0 =      ", self.tau0*1e3, " msec"
        print "tau_e =     ", self.tau_e*1e3, " msec"
        print "tau_x =     ", self.tau_x*1e3, " msec"
        print "tau_y =     ", self.tau_y*1e3, " msec"
        print "beta_x =    ", self.tws0.beta_x, " m"
        print "beta_y =    ", self.tws0.beta_y, " m"
        print "alpha_x =   ", self.tws0.alpha_x
        print "alpha_y =   ", self.tws0.alpha_y
        print "Dx =        ", self.tws0.Dx, " m"
        print "Dy =        ", self.tws0.Dy, " m"
        print "sigma_x =   ", self.sigma_x*1e6, " um"
        print "sigma_y =   ", self.sigma_y*1e6, " um"
        print "sigma_x' =   ", self.sigma_xp*1e6, " urad"
        print "sigma_y' =   ", self.sigma_yp*1e6, " urad"