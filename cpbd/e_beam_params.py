__author__ = 'Sergey'

from scipy.integrate import simps

from ocelot.common.globals import *
from ocelot.cpbd.optics import trace_z, twiss
from ocelot.cpbd.beam import *
from ocelot.cpbd.elements import *
from ocelot.rad.undulator_params import *


def I2_ID(L, h0):
    return L/2.*h0*h0


def I3_ID(L, h0):
    return 4.*L/(3*pi)*h0**3


def I4_ID(L, h0, lu):
    return -3.*L*lu*lu*h0**4/(32*pi)


def I5_ID(L, h0, lu, beta_xc, Dx0, Dxp0):
    # it is the same as I5_exact and I5_exact2
    # beta_c - beta_x at the center of ID
    # Dx0, Dxp0 - Dx and Dxp at the beginning of ID
    nn = int(L/lu)
    I = ((h0**3 *L)/(108000 *pi**5 *beta_xc) *(144000* pi**4* (Dx0**2 + beta_xc**2 *Dxp0**2) +
        13500* (-1)**nn* h0* pi**3* Dx0* lu**2 + 15656 *h0**2* lu**4 +
        15* (-76 + 225* (-1)**nn) *h0**2 *pi* lu**4 +
        150 *h0* pi**2* lu**2* (480* Dx0 + h0* (4* L**2 + 48 *beta_xc**2 - lu**2))))
    return I

def radiation_integrals(lattice, twiss_0, nsuperperiod = 1):
    #TODO: add I4 for rectangular magnets I4 = Integrate(2 Dx(z)*k(z)*h(z), Z)
    
    n_points_element = 20
    
    tws_elem = twiss_0
    (I1, I2, I3,I4, I5) = (0., 0., 0., 0., 0.)
    h = 0.
    for elem in lattice.sequence:
        if elem.__class__ in (SBend, RBend, Bend) and elem.l != 0:
            Dx = []
            Hinvariant = []
            Z = []
            h = elem.angle/elem.l

            for z in linspace(0, elem.l,num = n_points_element, endpoint=True):
                tws_z = elem.transfer_map(z)*tws_elem
                Dx.append(tws_z.Dx)
                Z.append(z)
                Hx = (tws_z.gamma_x*tws_z.Dx*tws_z.Dx + 2.*tws_z.alpha_x*tws_z.Dxp*tws_z.Dx
                                        + tws_z.beta_x*tws_z.Dxp*tws_z.Dxp)
                Hinvariant.append(Hx)
            #H = array(h)
            H2 = h*h
            H3 = abs(h*h*h)
            I1 += h*simps(array(Dx), Z)
            I2 += H2*elem.l  #simps(H2, Z)*nsuperperiod
            I3 += H3*elem.l  #simps(H3, Z)*nsuperperiod
            I4 += h*(2*elem.k1 + H2)*simps(array(Dx), Z)
            I5 += H3*simps(array(Hinvariant), Z)
        tws_elem = elem.transfer_map*tws_elem
    #if abs(tws_elem.beta_x - twiss_0.beta_x)>1e-7 or abs(tws_elem.beta_y - twiss_0.beta_y)>1e-7:
    #    print( "WARNING! Results may be wrong! radiation_integral() -> beta functions are not matching. ")
        #return None
    return (I1*nsuperperiod,I2*nsuperperiod,I3*nsuperperiod, I4*nsuperperiod, I5*nsuperperiod)

class EbeamParams:
    def __init__(self, lattice, beam,  coupling = 0.01, nsuperperiod = 1, tws0 = None):
        if beam.E == 0:
            exit("beam.E must be non zero!")
        self.E = beam.E
        if tws0 == None: 
            tws = twiss(lattice, Twiss(beam))
            self.tws0 = tws[0]
        else:
            #tws0.E = lattice.energy
            self.tws0 = tws0 
            tws = twiss(lattice, tws0)
            
        self.lat = lattice
        (I1,I2,I3, I4, I5) = radiation_integrals(lattice, self.tws0 , nsuperperiod)
        self.I1 = I1
        self.I2 = I2
        self.I3 = I3
        self.I4 = I4
        self.I5 = I5
        #print "I2 = ", I2
        #print "I3 = ", I3
        #print "I4 = ", I4
        #print "I5 = ", I5
        self.Je = 2 + I4/I2
        self.Jx = 1 - I4/I2
        self.Jy = 1
        self.gamma = self.E/m_e_GeV
        self.sigma_e = self.gamma*sqrt(Cq * self.I3/(self.Je*I2))
        self.emittance = Cq*self.gamma*self.gamma * self.I5/(self.Jx* self.I2)
        self.U0 = Cgamma*(beam.E*1000)**4*self.I2/(2*pi)
        #print "*********  ", twiss_0.Energy
        self.Tperiod = nsuperperiod*lattice.totalLen/speed_of_light
        self.Length = nsuperperiod*lattice.totalLen
        self.tau0 = 2*self.E*1000*self.Tperiod/self.U0
        self.tau_e = self.tau0/self.Je
        self.tau_x = self.tau0/self.Jx
        self.tau_y = self.tau0/self.Jy
        self.alpha = self.I1/(speed_of_light*self.Tperiod)
        self.coupl = coupling
        self.emitt_x = self.emittance/(1 + self.coupl)
        self.emitt_y = self.emittance*self.coupl/(1 + self.coupl)
        self.sigma_x = sqrt((self.sigma_e*self.tws0.Dx)**2 + self.emitt_x*self.tws0.beta_x)
        self.sigma_y = sqrt((self.sigma_e*self.tws0.Dy)**2 + self.emitt_y*self.tws0.beta_y)
        self.sigma_xp = sqrt((self.sigma_e*self.tws0.Dxp)**2 + self.emitt_x*self.tws0.gamma_x)
        self.sigma_yp = sqrt((self.sigma_e*self.tws0.Dyp)**2 + self.emitt_y*self.tws0.gamma_y)


    def integrals_id(self):
        L = 0.
        self.I2_IDs = 0.
        self.I3_IDs = 0.
        self.I4_IDs = 0.
        self.I5_IDs = 0.
        for elem in self.lat.sequence:
            if elem.type == "undulator":
                B = K2field(elem.Kx, lu = elem.lperiod)
                h0 = B*speed_of_light/self.E*1e-9
                #print h0, B
                tws = trace_z(self.lat, self.tws0, [L, L + elem.l/2.])
                i2 = I2_ID(elem.l,h0)
                i3 = I3_ID(elem.l,h0)
                i4 = I4_ID(elem.l,h0,elem.lperiod)
                i5 = I5_ID(elem.l,h0,elem.lperiod,tws[1].beta_x,tws[0].Dx, tws[0].Dxp)
                self.I2_IDs += i2
                self.I3_IDs += i3
                self.I4_IDs += i4
                self.I5_IDs += i5
                #print elem.type, elem.id, "B0 =  ", B, " T"
                #print elem.type, elem.id, "rho = ", 1./h0, " m"
                #print elem.type, elem.id, "L =   ", elem.l, " m"
                #print elem.type, elem.id, "beta_x cntr: ", tws[1].beta_x
                #print elem.type, elem.id, "Dx0 / Dxp0:  ", tws[0].Dx, "/", tws[0].Dxp
                #print elem.type, elem.id, "I2_ID = ", i2
                #print elem.type, elem.id, "I3_ID = ", i3
                #print elem.type, elem.id, "I4_ID = ", i4
                #print elem.type, elem.id, "I5_ID = ", i5
            L += elem.l
        self.emit_ID = self.emittance * (1.+self.I5_IDs/self.I5)/(1+(self.I2_IDs  - self.I4_IDs)/(self.I2 - self.I4))
        self.sigma_e_ID = self.sigma_e * sqrt((1.+ self.I3_IDs / self.I3)/(1 + (2*self.I2_IDs + self.I4_IDs)/(2.*self.I2 + self.I4) ) )
        self.U0_ID = Cgamma*(self.E*1000)**4.*self.I2_IDs/(2.*pi)
        print("emittance with IDs = ", self.emit_ID*1e9, " nm*rad")
        print("sigma_e with IDs =   ", self.sigma_e_ID)
        print("U0 from IDs =        ", self.U0_ID,  "  MeV")

    def __str__(self):
        val = ""
        val += ( "I1 =        " + str(self.I1) )
        val += ( "I2 =        " + str(self.I2) )
        val += ( "\nI3 =        " + str(self.I3) )
        val += ( "\nI4 =        " + str(self.I4) )
        val += ( "\nI5 =        " + str(self.I5) )
        val += ( "\nJe =        " + str(self.Je) )
        val += ( "\nJx =        " + str(self.Jx) )
        val += ( "\nJy =        " + str(self.Jy) )
        val += ( "\nenergy =    " + str(self.E) +"GeV")
        val += ( "\ngamma =     " + str(self.gamma) )
        val += ( "\nsigma_e =   " + str(self.sigma_e) )
        val += ( "\nemittance = " + str(self.emittance*1e9) +" nm*rad")
        val += ( "\nLength =    " + str(self.Length) + " m")
        val += ( "\nU0 =        " + str(self.U0) + "  MeV")
        val += ( "\nTperiod =   " + str(self.Tperiod*1e9) + " nsec")
        val += ( "\nalpha =     " + str(self.alpha) )
        val += ( "\ntau0 =      " + str(self.tau0*1e3) + " msec")
        val += ( "\ntau_e =     " + str(self.tau_e*1e3) + " msec")
        val += ( "\ntau_x =     " + str(self.tau_x*1e3) + " msec")
        val += ( "\ntau_y =     " + str(self.tau_y*1e3) + " msec")
        val += ( "\nbeta_x =    " + str(self.tws0.beta_x) + " m")
        val += ( "\nbeta_y =    " + str(self.tws0.beta_y) +" m")
        val += ( "\nalpha_x =   " + str(self.tws0.alpha_x))
        val += ( "\nalpha_y =   " + str(self.tws0.alpha_y))
        val += ( "\nDx =        " + str(self.tws0.Dx) + " m")
        val += ( "\nDy =        " + str(self.tws0.Dy) + " m")
        val += ( "\nsigma_x =   " + str(self.sigma_x*1e6) + " um")
        val += ( "\nsigma_y =   " + str(self.sigma_y*1e6) + " um")
        val += ( "\nsigma_x' =  " + str(self.sigma_xp*1e6) + " urad")
        val += ( "\nsigma_y' =  " + str(self.sigma_yp*1e6) + " urad\n")
        return val
        
    