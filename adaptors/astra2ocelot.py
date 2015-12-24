import numpy as np
#from numpy.core.umath import sqrt
from ocelot.common.globals import m_e_eV
from ocelot.cpbd.beam import *

def exact_xp_2_xxstg(xp, gamref):
    N = xp.shape[0]
    xxstg = np.zeros((N, 6))
    pref = m_e_eV*sqrt(gamref**2-1)
    u = np.c_[xp[:, 3], xp[:, 4], xp[:, 5]+pref]
    gamma = np.sqrt(1 + np.sum(u*u, 1)/m_e_eV**2)
    beta = np.sqrt(1-gamma**-2)
    if np.__version__ > "1.8":
        p0 = np.linalg.norm(u, 2, 1).reshape((N, 1))
    else:
        p0 = np.sqrt(u[:, 0]**2 + u[:, 1]**2 + u[:, 2]**2).reshape((N, 1))
    u = u/p0
    cdt = -xp[:, 2]/(beta*u[:, 2])
    xxstg[:, 0] = xp[:, 0]+beta*u[:, 0]*cdt
    xxstg[:, 2] = xp[:, 1]+beta*u[:, 1]*cdt
    xxstg[:, 4] = cdt
    xxstg[:, 1] = u[:, 0]/u[:, 2]
    xxstg[:, 3] = u[:, 1]/u[:, 2]
    xxstg[:, 5] = gamma/gamref-1
    return xxstg


def exact_xxstg_2_xp(xxstg, gamref):
    N = xxstg.shape[0]
    xp = np.zeros((N, 6))
    pref = m_e_eV*sqrt(gamref**2-1)
    gamma = gamref*(1+xxstg[:, 5])
    beta = np.sqrt(1-gamma**-2)
    u = np.c_[xxstg[:, 1], xxstg[:, 3], np.ones(N)]
    if np.__version__ > "1.8":
        norm = np.linalg.norm(u, 2, 1).reshape((N, 1))
    else:
        norm = np.sqrt(u[:, 0]**2 + u[:, 1]**2 + u[:, 2]**2).reshape((N, 1))
    u = u/norm
    xp[:, 0] = xxstg[:, 0]-u[:, 0]*beta*xxstg[:, 4]
    xp[:, 1] = xxstg[:, 2]-u[:, 1]*beta*xxstg[:, 4]
    xp[:, 2] = -u[:, 2]*beta*xxstg[:, 4]
    xp[:, 3] = u[:, 0]*gamma*beta*m_e_eV
    xp[:, 4] = u[:, 1]*gamma*beta*m_e_eV
    xp[:, 5] = u[:, 2]*gamma*beta*m_e_eV-pref
    return xp


def astraBeam2particleArray(filename):
    P0 = np.loadtxt(filename)
    charge_array = -P0[:, 7]*1e-9  #charge in nC -> in C
    print "charge = ", sum(charge_array)
    xp = P0[:, :6]
    Pref = xp[0, 5]
    xp[0,5] = 0
    gamref = sqrt((Pref/m_e_eV)**2+1)
    xxstg = exact_xp_2_xxstg(xp, gamref)

    p_array = ParticleArray(len(charge_array))
    p_array.E =  sqrt((Pref/m_e_eV)**2+1)*m_e_GeV
    p_array.particles[0::6] = xxstg[:,0]
    p_array.particles[1::6] = xxstg[:,1]
    p_array.particles[2::6] = xxstg[:,2]
    p_array.particles[3::6] = xxstg[:,3]
    p_array.particles[4::6] = xxstg[:,4]
    p_array.particles[5::6] = xxstg[:,5]
    return p_array, charge_array

def particleArray2astraBeam(p_array):
    gamref = p_array.E/m_e_GeV
    print gamref
    P = p_array.particles.view()
    P.shape = len(P)/6,6
    xp = exact_xxstg_2_xp(P, gamref)
    xp[0, 5] = xp[0, 5] + p_array.E*1e9
    np.savetxt('pytest.ast',xp)
