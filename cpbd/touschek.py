from ocelot import *
from pylab import *
from scipy.special import kn


def B1():
    pass

def B2():
    pass


def F(tm, B1, B2):
    #print B1, B2
    Fi = 2. * sqrt(pi * (B1**2 - B2**2)) 
    
    km = arctan( sqrt(tm) )
    
    I2 = 0.0
    Nk = 5000
    ks = linspace(km, pi/2, Nk)
    dk = ks[1] - ks[0]
    for k in ks:
        t = tan(k)**2
        dI = (2*t+1)**2 * (t/tm/(1+t)-1.)/t + t - sqrt(t*tm*(1+t)) - (2 + 0.5/t) * log(t/tm/(1+t))
        #print t, kn(0, B2*t), B1, B2, B2*t, dI 
        #print -B1*t
        #print t, B1, ':', exp(-B1*t)
        dI *= exp(-B1*t) * kn(0,B2*t) * sqrt(1 + t)
        I2 += dI * dk
    
    return Fi * I2
    
def t_touschek(tws, beam):
    '''
    Touschek lifetime
    Piwinski DESY 98-179
    '''
    rp2c = 2.380588258494943e-21    
    gam  = beam.E / 0.000511
    beta = sqrt(1 - 1./gam**2)
    c = 299792458.0
    
    tm = beta**2 * beam.sigma_e_cut **2
    #print 'tm:', tm
    
    ex = beam.emit_x
    ey = beam.emit_y
    
    sig_e = beam.sigma_e
    
    I1 = 0.0
    
    #print sig_e, ex, ey
    #sys.exit(0)
    
    for i in range(len(tws)-1):
        
        sig_x = sqrt(ex * tws[i].beta_x)
        sig_y = sqrt(ey * tws[i].beta_y)
        
        sig_x2 = ex * tws[i].beta_x
        sig_y2 = ey * tws[i].beta_y
        
                
        Dx = tws[i].Dx
        Dy = tws[i].Dy

        betx = tws[i].beta_x
        bety = tws[i].beta_y

        
        Dxh = tws[i].alpha_x * tws[i].Dx +  tws[i].beta_x * tws[i].Dxp
        Dyh = tws[i].alpha_y * tws[i].Dy +  tws[i].beta_y * tws[i].Dyp
        
        #print 1./ sig_e**2 + (Dx**2 + Dxh**2) / sig_x**2 + (Dy**2 + Dyh**2) / sig_y**2
        
        sigh2 = 1. / (1./ sig_e**2 + (Dx**2 + Dxh**2) / sig_x2 + (Dy**2 + Dyh**2) / sig_y2) 
        
        B1 = betx**2 / (2.*beta**2 * gam**2 * sig_x2) * (1. - sigh2 * Dxh**2 / sig_x2) 
        B1 +=  bety**2 / (2.*beta**2 * gam**2 * sig_y2) * (1. - sigh2 * Dyh**2 / sig_y2)
        
        #B2 = B1**2 - tws[i].beta_x**2 * tws[i].beta_y**2 * sigh2 / (beta**4 * gam**4 * sig_x**4 * sig_y**4 * sig_e**2)* (sigp_x**2*sigp_y**2 - sig_e**4 * Dx**2 * Dy**2)
        
        B2 = ( (1. - sigh2 * Dxh**2 /sig_x2) * betx**2 / sig_x2 - (1. - sigh2 * Dyh**2 /sig_y2) * bety**2/ sig_y2)**2 / (4 * beta**4 * gam**4)
        B2 += (sigh2**2 * betx**2 * bety**2 * Dxh**2 * Dyh**2) / (beta**4 * gam**4 * sig_x2 * sig_y2)
        
        #print 'B2', B2
        if B2 < 0: print 'B2 error' 
        B2 = sqrt(B2)
        '''
        print 'B1-B2',B1**2 - B2**2
        if B1**2 - B2**2 < 0:
            print tws[i]
        '''
        
        I1 += F(tm, B1, B2) / ( sqrt(sig_x2 * sig_y2 - sig_e**4 * Dx**2 * Dy**2) ) * (tws[i+1].s - tws[i].s)
    
    
    return 1. / (  rp2c * beam.Np * I1 / (8. * pi * gam**2 * beam.sigma_s) / tws[-1].s ) 


