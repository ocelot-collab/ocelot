'''
Bending magnet radiation
for arc motion
'''

from pylab import *

c = 3.e8
hc = 1.97e-7 # ev m 
gam = 6. / 0.000511

sigx = 1.e-4
sigxp = 1.e-5


def bm_e(x,y,xp,yp,x0,y0,z,p_en, R):
    
    th_x = x0 / z
    th_y = y0 / z
        
    Ex0 = bm_e0(x0,y0,z,p_en, R)
    
    phi_s = p_en / hc * z  *  (th_x**2 + th_y**2) / (2. )
    phi_0 = -p_en / (2.*hc) * R * (th_x - xp) * (  1./gam**2 + (th_y - yp)**2 + (th_x - xp)**2 / 3.) - p_en / hc * (x*th_x + y*th_y)
    Ex = np.exp(1j*phi_s) * np.exp(1j*phi_0) * Ex0  
            
    return Ex 


def bm_e_a(x,y,xp,yp,x0,y0,z,p_en, R):
    '''
    Bending magnet field averaged over ...
    '''
    th_x = x0 / z
    th_y = y0 / z
    
    
    Ex0 = bm_e0(x0,y0,z,p_en, R)
        
    #Ex = np.zeros([len(x), len(xp)])
    Ex = 0j

    phi_s = p_en / hc * z  *  (th_x**2 + th_y**2) / (2. )
        
    for i in range(len(x)):
        for j in range(len(xp)):
            phi_0 = -p_en / (2.*hc) * R * (th_x - xp[j]) * (  1./gam**2 + (th_y - yp)**2 + (th_x - xp[j])**2 / 3.) - p_en / hc * (x[i]*xp[j] + y*yp)
            Ex += np.exp(1j*phi_s) * np.exp(1j*phi_0) * Ex0 * exp(-0.5*(x[i]/sigx)**2 - 0.5*(xp[j]/sigxp)**2) / (sigx * sigxp) 
            
    return Ex / (len(x) * len(xp)) 
    

def bm_e0(x0,y0,z,p_en, R):

    th_y = y0 / z
    
    Exr = 0.0
    Exi = 0.0
    
    zmin = -2.6
    zmax = 2.6
    nz = 15000
    dz  = (zmax - zmin)/nz 
    
    for i in range(nz):
        z = zmin + (i+0.5)*dz
        phi = (p_en/hc)*(0.5*z* th_y**2 + z**3 / (6*R**2))
        #print phi / (2*pi)
        Exr += z*cos(phi)
        Exi += z*sin(phi)
        z = zmin + i*dz
    
         
    return Exr + 1j*Exi


