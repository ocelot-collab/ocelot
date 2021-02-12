'''
crystal optics
'''

from ocelot.optics.elements import Crystal
from ocelot.optics.wave import *
from ocelot.optics.ray import Ray, trace as trace_ray
from ocelot.common.ocelog import *
import logging
_logger = logging.getLogger(__name__) 
# from ocelot.gui.optics import *

import numpy as np
# from pylab import *
import sys, os

m = 1.0
cm = 1.e-2
mm = 1.e-3
mum = 1.e-6
nm = 1.e-9
A = 1.e-10
K = 1.0

hbar = 6.58211928e-16
c = 299792458.0 

r_el        = 2.8179403267e-15       #Electron classical radius in:  m

def calc_angle(a,b):
    '''
    yields angle between vectors "a"=(ax,ay,az) and "b"=(bx,by,bz) in degrees
    '''
    a_ = np.array(a)
    b_ = np.array(b)
    return np.arccos(np.dot(a_,b_)/np.linalg.norm(a_)/np.linalg.norm(b_)) / np.pi * 180

'''
element database
'''
class Element():
    def __init__(self):
        pass

elements = {}

elements['H'] = Element()
elements['H'].Z = 1
elements['H'].A = 1

elements['C'] = Element()
elements['C'].Z = 6
elements['C'].A = 12

lattice_unit_cells = {}
cells = {}

lattice_unit_cells['NaCl'] = (np.array([5.6404*A, 0, 0]), np.array([0, 5.6404*A, 0]), np.array([0, 0, 5.6404*A]))

lattice_unit_cells['Si'] = (np.array([5.43*A, 0, 0]), np.array([0, 5.43*A, 0]), np.array([0, 0, 5.43*A]))
#error?
#cells['Si'] =(('C',[0,0,0]),('C',[0.5,0.5,0]),('C',[0.5,0,0.5]),('C',[0,0.5,0.5]),
#              ('C',[0.25,0.25,0.25]),('C',[0.75,0.25,0.75]),('C',[0.25,0.75,0.75]),('C',[0.75,0.75,0.25])) 

#diamond
lattice_unit_cells['C'] = (np.array([3.567*A, 0, 0]), np.array([0, 3.567*A, 0]), np.array([0, 0, 3.567*A]))
cells['C'] =(('C',[0,0,0]),('C',[0.5,0.5,0]),('C',[0.5,0,0.5]),('C',[0,0.5,0.5]),
                     ('C',[0.25,0.25,0.25]),('C',[0.75,0.25,0.75]),('C',[0.25,0.75,0.75]),('C',[0.75,0.75,0.25]))

lattice_unit_cells['Ge'] = (np.array([5.658*A, 0, 0]), np.array([0, 5.658*A, 0]), np.array([0, 0, 5.658*A]))
                     
class CrystalLattice():
    def __init__(self, element_name):
        self.element_name = element_name
        self.a1 = lattice_unit_cells[element_name][0]
        self.a2 = lattice_unit_cells[element_name][1]
        self.a3 = lattice_unit_cells[element_name][2]
        b_  = np.dot(self.a1, np.cross(self.a2, self.a3))
        self.b1 = np.cross(self.a2, self.a3) / b_
        self.b2 = np.cross(self.a3, self.a1) / b_
        self.b3 = np.cross(self.a1, self.a2) / b_

class CrystalCell():
    def __init__(self, element_name):
        pass

g1 = Crystal(r=[0,0,0*cm], size=[13*cm,20*cm,10*cm], no=[0,0,-1], id='cr1')
g1.lattice =  CrystalLattice('C')
g1.temp = 100 * K


s0 = [0,0,1] # beam direction

class StructureFactorFactory():
    def __init__(self):
        pass

    def atomic_structure_factor(self, element, model='spherical_elastic'):
        if model == 'spherical_elastic':
            return self.f_sperical_elastic(element)
        else: return None

    def f_spher(self, a, k):
        '''
        atomic scattering factor for hydrogen-like charge density of radius a
        k is the momentum transfer sin(theta)/lambda
        '''
        return 1. / (1. + (2.*pi*a*k)**2)**2

    def f_sperical_elastic(self, element):
        if element == 'Li': ak = 12.0
        else: ak = 13.0
        def f(k):
            return ak*k
        return f


'''
atomic scattering factor
'''
def f_scat(phi, lambd):
    '''
    atomic scattering factor
    '''
    pass

class CrystalStructureFactors():
    def __init__(self):
        pass

def load_stucture_factors(file_name):
    _logger.debug('loading structure factors {}'.format(file_name))   
    import pickle
    #sys.modules[]
    #print ('im in module', __name__)
    #print ('file ', sys.modules[__name__].__file__)
    dir_name = os.path.dirname(sys.modules[__name__].__file__)
    #print 'dir ', dir_name 
    abs_file_name = os.path.join(dir_name, 'data', file_name)
    _logger.debug(ind_str + 'full path {}'.format(abs_file_name))   
#    print ('abs path ', abs_file_name)
    cdata = pickle.load(open(abs_file_name, 'rb'))
    return cdata

def save_stucture_factors(cdata, file_name):
    import pickle
    print ('saving structure factors to', file_name)
    # pickle.dump(cdata, open(file_name, 'wb')) #original
    pickle.dump(cdata, open(file_name, 'wb'), pickle.HIGHEST_PROTOCOL) # fix for "No module named 'ocelot.optics.bragg\r'" error in Windows


def F_hkl(cryst, ref_idx, lamb, temp):
    
    _logger.debug('calculating Fhkl {} {} {}'.format(lamb, cryst.lattice.element_name, ref_idx))
    
    file_name = cryst.lattice.element_name + str(ref_idx[2]) + str(ref_idx[1]) + str(ref_idx[0]) + '.dat'
    
    target_ev = 2*pi * hbar * c / lamb
    
    #print ('reading file_name', file_name)
    #_logger.debug(ind_str + 'loading structure factors from {}'.format(file_name))    
#    cdata = load_stucture_factors(file_name)
    
    try:
        cdata = load_stucture_factors(file_name)
    except:
#        print ('form factor data not found!!!')
        _logger.error(ind_str + 'form factor data not found')
        sys.exit(0)
        return 0.0, 0.0, 0.0
    
#    print ('searching ', target_ev)
    _logger.debug(ind_str + 'searching Bragg condition for photon energy {} eV'.format(target_ev))
    
    if target_ev < cdata.ev[0] or target_ev > cdata.ev[-1]:
        _logger.error(2*ind_str + 'photon wavelength not covered in data (target_ev > cdata.ev)')
#        print ('photon wavelength not covered in data')
        return 0.0, 0.0, 0.0
        
    de = cdata.ev[1] - cdata.ev[0]
    i_e = int( (target_ev - cdata.ev[0]) / de )

    # return cdata.f000[i_e], cdata.fh[i_e], cdata.fhbar[i_e]
    f000 = np.interp(target_ev, cdata.ev, np.real(cdata.f000)) + 1j * np.interp(target_ev, cdata.ev, np.imag(cdata.f000))
    fh = np.interp(target_ev, cdata.ev, np.real(cdata.fh)) + 1j * np.interp(target_ev, cdata.ev, np.imag(cdata.fh))
    fhbar = np.interp(target_ev, cdata.ev, np.real(cdata.fhbar)) + 1j * np.interp(target_ev, cdata.ev, np.imag(cdata.fhbar))
    _logger.debug(ind_str + 'f000 {}'.format(f000))
    _logger.debug(ind_str + 'fh {}'.format(fh))
    _logger.debug(ind_str + 'hbar {}'.format(fhbar))
    return f000, fh, fhbar



def find_bragg(lambd, lattice, ord_max):
    H = {}
    d = {}
    phi = {}
    
    for h in range(0, ord_max+1):
        for k in range(0, ord_max+1):
            for l in range(0, ord_max+1):
                
                if h == k == l == 0:
                    continue
                            
                hkl = (h,k,l)
    
                H_hkl = h*lattice.b1 + k*lattice.b2 + l*lattice.b3
                d_hkl = 1. / np.linalg.norm(H_hkl)
                                        
                sin_phi = lambd / (2. * d_hkl)
                
                if np.abs(sin_phi) <= 1.0:
                    H[hkl] = H_hkl
                    d[hkl] = d_hkl
                    phi[hkl] = np.arcsin(sin_phi) * 180.0 / np.pi

    return H, d, phi
   

   
def plot_bragg_reflections(idc = [(0,0,1), (1,1,1), (2,1,1), (3,1,1), (1,2,3), (4,0,0)]):
    HH = {}
    import matplotlib.pyplot as plt
    lambds = np.linspace(0.5e-10 * m, 1.3e-9 * m, 50)
    
    for lambd in lambds:
        H, d, phi = find_bragg(lambd = lambd, lattice=g1.lattice, ord_max = 5)
        HH[lambd] = (H, d, phi) 
    
    plt.figure()
    plt.grid()
    lines = []; labels=[]
    
    for idx in idc:
        y = []
        lambds = []
        
        for l in HH.keys(): #sort(HH.keys()):
            try:    
                y.append(HH[l][2][idx])
                lambds.append(l)
            except:
                pass
                
        l1, = plt.plot(lambds, y)
        lines.append(l1), labels.append(str(idx))
    
    plt.legend(lines,labels)
    
    plt.show()
    
def plot_scattering_factors():
    f = StructureFactorFactory().atomic_structure_factor('C')
    print (f(2))
    

def eta(Dtheta, cryst):
    '''
    Description   : Calculates the eta factor as from Authier, p.136 Eq. 5.31
    '''

    etaval = (
               ( Dtheta * np.sin(2*cryst.thetaB) + cryst.chi0*(1-cryst.gamma)/2 )
               /
               ( np.sqrt(abs(cryst.gamma))*abs(cryst.c_pol)*np.sqrt(cryst.chih*cryst.chimh) )
            )
    return etaval




def xij(Dtheta, cryst):
    '''
        Description   : Calculates the xijp factor as from Authier:
                        Laue  : Eq. 5.36
                        Bragg : Eq. 5.38
    '''

    etav = eta(Dtheta, cryst)
 
    if cryst.bragg :
        xijpval = (-np.sign(cryst.c_pol) / np.sqrt(abs(cryst.gamma)) * np.sqrt(cryst.chih*cryst.chimh)/cryst.chimh * (etav + np.sqrt(etav**2 - 1)))
        xijppval = (-np.sign(cryst.c_pol) / np.sqrt(abs(cryst.gamma)) * np.sqrt(cryst.chih*cryst.chimh)/cryst.chimh * (etav - np.sqrt(etav**2 - 1)) )

    else :
        xijpval = (-np.sign(cryst.c_pol) / np.sqrt(cryst.gamma) * np.sqrt(cryst.chih*cryst.chimh)/cryst.chimh * (etav + np.sqrt(etav**2 + 1)) )
        xijppval = (-np.sign(cryst.c_pol) / np.sqrt(cryst.gamma) * np.sqrt(cryst.chih*cryst.chimh)/cryst.chimh * (etav - np.sqrt(etav**2 + 1)) )
              
    return xijpval, xijppval





def MP(Dtheta, cryst):
    '''
            Description   : Calculates the MPp factor as from Authier p185:
                        Laue  : Eq. 6.15, Bragg Eq. 7.9
    '''

    etav = eta(Dtheta, cryst)
    
    k = 1./ cryst.lamb
 
    MPpval = (
                k * cryst.chi0 / ( 2*cryst.gamma0 )
                + np.sign(cryst.gammah) * etav / ( 2*cryst.pl )
                - np.sqrt( etav**2 + np.sign(cryst.gammah) ) / ( 2*cryst.pl )
             )
    MPppval = (
                k * cryst.chi0 / ( 2*cryst.gamma0 )
                + np.sign(cryst.gammah) * etav / ( 2*cryst.pl )
                + np.sqrt( etav**2 + np.sign(cryst.gammah) ) / ( 2*cryst.pl )
             )

    return MPpval, MPppval



'''
formerly
D0dfracD0a
DhdfracD0a
TODO: add a reference to Authier
'''
def D0_Dh(Dtheta, cryst):
    
    thick = cryst.size[2]
 
    if cryst.gamma<0:
        BRAGG = True
    else:
        BRAGG = False


    cryst.bragg = BRAGG
    
    xijpv,  xijppv = xij(Dtheta, cryst)
    #xijppv = xijpp(Dtheta)
    MPpv, MPppv   = MP(Dtheta, cryst)
     
 
    if BRAGG :
         D0dfracD0aval = (    
              (xijpv - xijppv) * np.exp( -2*np.pi*1j*(MPpv+MPppv)*thick ) /
              (  xijpv * np.exp( -2*np.pi*1j*MPpv *thick )
               - xijppv* np.exp( -2*np.pi*1j*MPppv*thick )
              )
         )
         DhdfracD0aval = (    
              xijpv * xijppv
              * ( np.exp( -2*np.pi*1j*MPpv*thick ) - np.exp( -2*np.pi*1j*MPppv*thick ) )
            /
              (  xijpv * np.exp( -2*np.pi*1j*MPpv *thick )
               - xijppv* np.exp( -2*np.pi*1j*MPppv*thick )
              )
         )

    else :
         D0dfracD0aval = (
              ( xijpv * np.exp( -2*np.pi*1j*MPpv *thick )
               -xijppv* np.exp( -2*np.pi*1j*MPppv*thick ) ) / (xijpv - xijppv) 
        )
         DhdfracD0aval = (
              ( xijpv * xijppv / (xijppv - xijpv) 
                * (   np.exp( -2*np.pi*1j*MPpv *thick )
                    - np.exp( -2*np.pi*1j*MPppv*thick ) )
               )
        )
              
    return D0dfracD0aval, DhdfracD0aval





def transmissivity_reflectivity(klist, cryst):
    _logger.debug('calculating transmissivity and reflectivity')
    t  = np.zeros(len(klist),'complex')
    r  = np.zeros(len(klist),'complex')
    for i in range(len(klist)):
        t[i], r[i] = D0_Dh( -(cryst.kb - klist[i]) / ( cryst.kb * (1/np.tan(cryst.thetaB)) ) , cryst)
    abs_r = abs(r)
    abs_t = abs(t)
    if np.nan in abs_t or np.nan in abs_r:
        _logger.error(ind_str + 'abs(tr) = {}-{}'.format(abs(t[0]), abs(t[-1])))
        _logger.error(ind_str + 'abs(ref) = {}-{}'.format(abs(r[0]), abs(r[-1])))
    else:
        _logger.debug(ind_str + 'abs(tr) = {}-{}'.format(abs(t[0]), abs(t[-1])))
        _logger.debug(ind_str + 'abs(ref) = {}-{}'.format(abs(r[0]), abs(r[-1])))
    return t, r


def get_crystal_filter(cryst, ev_seed, nk=10000, k = None, n_width = 100):
    _logger.debug('getting crystal filter')
    #import crystal as cry
    
    #n_width - number of Darwin widths
    lamb=h_eV_s*speed_of_light/ev_seed
    ref_idx=cryst.ref_idx
    n = cryst.cut
    _logger.debug(ind_str + 'reflection index = {}'.format(ref_idx))
    kb = 2*np.pi/lamb
    
    H, d, phi = find_bragg(lambd = lamb, lattice=cryst.lattice, ord_max = 15)
    
    polarization = 'sigma'

    abs_ref_idx = tuple(abs(np.array(ref_idx)))
    dhkl = d[abs_ref_idx]
    
    '''
    print 'lamb=', lamb
    print 'thick=', cryst.size[2]
    print 'd_hkl',d[ref_idx]
    print 'theta bragg', phi[ref_idx]
    '''
    
    thetaB = phi[abs_ref_idx]* np.pi / 180.0
    _logger.debug(ind_str + 'Bragg angle {} [rad]'.format(thetaB))
    _logger.debug(ind_str + 'Bragg angle {} [deg]'.format(thetaB/np.pi*180))
    
    if polarization == 'sigma':
        c_pol = 1
    else :
        c_pol = - np.cos(2.0 * thetaB)
    
    ref_idx_length = np.sqrt(np.sum(np.array(ref_idx)**2))
    h_norm = tuple([i/ref_idx_length for i in ref_idx])
    h = np.array([i*2*np.pi/dhkl for i in h_norm]) # normalized q-vector of crystalline plane

    #quadratic equation coefficients to find incidence angle
    a = - 0.5 * np.sin(thetaB)**2 + 0.5*h_norm[2]**2
    b = 2/np.sqrt(2) * (h_norm[0] + h_norm[1]) * h_norm[2]
    c = - 2 * np.sin(thetaB)**2 + (h_norm[0] + h_norm[1])**2

    ctg2_tilt = np.polynomial.polynomial.polyroots((a,b,c))
    tilt = np.arctan(ctg2_tilt*2)/np.pi*180
    tilt_min = tilt[np.argmin(abs(tilt))]
    _logger.debug(ind_str + 'normal angle = {} deg'.format(calc_angle((0,0,1), n)))
    _logger.debug(ind_str + 'tilt solutions = {} deg'.format(tilt))

    #s0_0 = np.array([-1*sqrt(2)*np.tan((tilt[0] - 90)/180*np.pi)/2, -1*sqrt(2)*np.tan((tilt[0] - 90)/180*np.pi)/2, -1])
    s0_0 = np.array([-1*np.sqrt(2)*ctg2_tilt[0], -1*np.sqrt(2)*ctg2_tilt[0], -1])
    s0_0 = s0_0 / np.linalg.norm(s0_0)

    #s0_1 = np.array([-1*sqrt(2)*np.tan((tilt[1] - 90)/180*np.pi)/2, -1*sqrt(2)*np.tan((tilt[1] - 90)/180*np.pi)/2, -1])
    s0_1 = np.array([-1*np.sqrt(2)*ctg2_tilt[1], -1*np.sqrt(2)*ctg2_tilt[1], -1])
    s0_1 = s0_1 / np.linalg.norm(s0_1)

    gamma0_0 = np.dot(n, s0_0)/np.linalg.norm(s0_0)/np.linalg.norm(n)
    gamma0_1 = np.dot(n, s0_1)/np.linalg.norm(s0_1)/np.linalg.norm(n)

    if gamma0_0 < gamma0_1:
        s0 = s0_0 * 2 * np.pi / lamb
        gamma0 = gamma0_0
    else:
        s0 = s0_1 * 2 * np.pi / lamb
        gamma0 = gamma0_1

    #h = np.array(cryst.ref_idx) / np.linalg.norm(cryst.ref_idx) * 2 * np.pi / dhkl
    sh = s0 + h
    _logger.debug(ind_str + 'n   = {}'.format(n))
    _logger.debug(ind_str + 's0 = {}'.format(s0))
    _logger.debug(ind_str + 'h  = {}'.format(h))
    _logger.debug(ind_str + 'sh  = {}'.format(sh))
    
    
    gammah = np.dot(n, sh)/np.linalg.norm(sh)/np.linalg.norm(n)
    
    gamma0, gammah = -gamma0, -gammah ####bugfix?
    
    gamma = gamma0 / gammah

    _logger.debug(ind_str + '2 theta = {} deg'.format(calc_angle(s0,-sh)))
    _logger.debug(ind_str + 'angle s0 to -n (0==best): {}'.format(np.arccos(gamma0)/np.pi*180))
    _logger.debug(ind_str + 'angle sh to -n (90==worst?): {}'.format(180 - np.arccos(gammah)/np.pi*180))

    #################
    # cryst.psi_n=-np.pi/2
    # psih   = cryst.psi_n - thetaB
    # psi0   = cryst.psi_n + thetaB
    # gamma0 = np.cos(psi0)
    # gammah = np.cos(psih)
    # gamma  = gammah/gamma0
    ###################
    
    _logger.debug(ind_str + 'gamma0: {}'.format(gamma0))
    _logger.debug(ind_str + 'gammah: {}'.format(gammah))
    _logger.debug(ind_str + 'gamma: {}'.format(gamma))
    if gamma<0:
        _logger.debug(2*ind_str + 'BRAGG geometry')
    else:
        _logger.warning(2*ind_str + 'LAUE geometry')


    
    f0, fh, fmh =  F_hkl(cryst = cryst, ref_idx = abs_ref_idx, lamb=lamb, temp = 300*K)
    
    
#    _logger.debug(ind_str + 'structure factors {} {} {}'.format(f0, fh, fmh))
    #plt.figure()

     
        
    vcell = ( dhkl * np.sqrt( np.dot(ref_idx,ref_idx) ) ) **3
        
    delta    = r_el * np.abs(c_pol) * lamb**2 * np.sqrt( np.abs(gamma)*fh*fmh ) / ( np.pi * vcell * np.sin(2*thetaB) )
    mid_TH = np.real( r_el * lamb**2 * f0 * (1-gamma) / (2*np.pi * vcell * np.sin(2*thetaB)) )
    mid_k  = kb / (- mid_TH * 1/np.tan(thetaB) + 1 )
    dk  =  ( -2*kb*np.sin(thetaB) + np.sqrt(16*mid_k**2 * np.real(delta)**2 * np.cos(thetaB)**2 + 4*kb**2 * np.sin(thetaB)**2) ) / ( 2 * np.real(delta) * np.cos(thetaB) )
    
    print('vcell', vcell)
    print('mid_k', mid_k)
    
    ki, kf = mid_k - n_width*dk, mid_k + n_width*dk    #final   k    
    
    if k == None:    
        k = np.linspace(ki, kf, nk)    
        
    
    cryst.lamb = lamb
    cryst.kb = kb
    cryst.kb = 2 * kb - mid_k # bugfix
    cryst.thetaB = thetaB
    cryst.gamma = gamma
    cryst.gamma0 = gamma0
    cryst.gammah = gammah
    cryst.c_pol = c_pol
    cryst.chi0  = - ( r_el * lamb**2 * f0  ) / ( np.pi * vcell )
    cryst.chih  = - ( r_el * lamb**2 * fh  ) / ( np.pi * vcell )
    cryst.chimh = - ( r_el * lamb**2 * fmh ) / ( np.pi * vcell )
    cryst.pl = np.pi * vcell * np.sqrt( gamma0 * np.abs(gammah) ) / ( r_el * lamb * np.abs(c_pol) * np.sqrt( fh*fmh ) )
    
    
    from ocelot.optics.wave import TransferFunction #hotfix here
    f = TransferFunction()
    
    f.tr, f.ref = transmissivity_reflectivity(k, cryst)   
    
    f.k = k
    f.thetaB = thetaB
    f.mid_k = mid_k
    f.dk = dk
    f.xlamds = lamb
    f.cryst = cryst
    return f



def unfold_angles(Phlist):
    '''
    
    '''

    NK = len(Phlist)
    Phout  = np.zeros(NK)
    
    Phout[0] = Phlist[0] 
    count  = 0.0
    for i in range(1,NK):
        count = count + np.floor( (Phlist[i]-Phlist[i-1])/(2*np.pi) + 0.5 )
        Phout[i] = Phlist[i] - 2*np.pi*count

    return Phout


if __name__ == "__main__":
    plot_bragg_reflections()
    #plot_scattering_factors()
    