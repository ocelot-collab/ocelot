'''
Bending magnet
'''

from pylab import *

c = 3.e+8
hc = 1.97e-7 # ev m 
gam = 6.08/0.000511
e = 1.6e-19

sigx = 169.e-6
sigy = 10.e-6
sigxp = 5.e-4
sigyp = 1.e-4


def bm_e_a(x, y, xp, yp, x0, y0, z0, sgx, sgxp, sgy, sgyp, p_en, R):
    '''
    Bending magnet field averaged over ...
    '''
    th_x = x0/z0
    th_y = y0/z0
    
    Ex = 0j
    Ey = 0j
    
    phi_s = (p_en*z0/(2*hc))*(th_x**2 + th_y**2)
    Ex0, Ey0 = bm_e0(x0,y0,z0,p_en,R)

    for i in range(len(x)):
        for j in range(len(xp)):
            for k in range(len(y)):
                for l in range(len(yp)):
                    phi_0 = (-R*p_en*(th_x - xp[j])/(2.*hc)) * (1/ gam**2 + (th_y - yp[l])**2 + ((th_x - xp[j])**2)/3) - (p_en/hc)*(xp[j]*th_x + yp[l]*th_y)
                    Ex += (1j*p_en)/(hc*z0) * np.exp(1j*phi_s) * np.exp(1j*phi_0) * Ex0 * exp(-0.5*(x[i]/sgx)**2) * exp(-0.5*(xp[j]/sgxp)**2) * exp(-0.5*(y[k]/sgy)**2)/(sgx * sgxp * sgy * (sqrt(2*np.pi))**2)
                    Ey += (1j*p_en)/(hc*z0) * np.exp(1j*phi_s) * np.exp(1j*phi_0) * Ey0 * exp(-0.5*(x[i]/sgx)**2) * exp(-0.5*(xp[j]/sgxp)**2) * exp(-0.5*(y[k]/sgy)**2)/(sgx * sgxp * sgy * (sqrt(2*np.pi))**2)
    return Ex/(len(x) * len(xp) * len(y) * len(yp)), Ey/(len(x) * len(xp) * len(y) * len(yp))
    

def bm_e0(x0,y0,z0,p_en,R):

    th_y = y0/z0

    Ex = 0.0
    Ey = 0.0
        
    zmin = -5.378/2
    zmax = 5.378/2
    nz = 15000
    dz  = (zmax - zmin)/nz 
    
    for i in range(nz+1):
        z = zmin + i*dz
        phi = (p_en*z/(2*hc))*(1/gam**2 + z**2/(3*(R**2)) + th_y**2)
        Ex += (z/R)*np.exp(1j*phi)
        Ey += th_y*np.exp(1j*phi)    
  
    return Ex, Ey


