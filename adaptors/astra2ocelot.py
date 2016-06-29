
#from numpy.core.umath import sqrt
from ocelot.common.globals import m_e_eV
from ocelot.cpbd.beam import *
import numpy as np

def exact_xp_2_xxstg(xp, gamref):
    N = xp.shape[0]
    xxstg = np.zeros((N, 6))
    pref = m_e_eV*np.sqrt(gamref**2-1)
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
    N = len(xxstg)/6
    xp = np.zeros((N, 6))
    pref = m_e_eV*np.sqrt(gamref**2-1)
    gamma = gamref*(1+xxstg[5::6])
    beta = np.sqrt(1-gamma**-2)
    u = np.c_[xxstg[1::6], xxstg[3::6], np.ones(N)]
    if np.__version__ > "1.8":
        norm = np.linalg.norm(u, 2, 1).reshape((N, 1))
    else:
        norm = np.sqrt(u[:, 0]**2 + u[:, 1]**2 + u[:, 2]**2).reshape((N, 1))
    u = u/norm
    xp[:, 0] = xxstg[0::6]-u[:, 0]*beta*xxstg[4::6]
    xp[:, 1] = xxstg[2::6]-u[:, 1]*beta*xxstg[4::6]
    xp[:, 2] = -u[:, 2]*beta*xxstg[4::6]
    xp[:, 3] = u[:, 0]*gamma*beta*m_e_eV
    xp[:, 4] = u[:, 1]*gamma*beta*m_e_eV
    xp[:, 5] = u[:, 2]*gamma*beta*m_e_eV-pref
    return xp


def astraBeam2particleArray(filename):
    P0 = np.loadtxt(filename)
    charge_array = -P0[:, 7]*1e-9  #charge in nC -> in C
    print( "charge = ", sum(charge_array))
    xp = P0[:, :6]
    Pref = xp[0, 5]
    s_ref = xp[0, 2]
    xp[0, 5] = 0
    xp[0, 2] = 0.
    gamref = np.sqrt((Pref/m_e_eV)**2+1)
    xxstg = exact_xp_2_xxstg(xp, gamref)

    p_array = ParticleArray(len(charge_array))
    p_array.s = s_ref
    p_array.E =  np.sqrt((Pref/m_e_eV)**2+1)*m_e_GeV
    p_array.particles[0::6] = xxstg[:,0]
    p_array.particles[1::6] = xxstg[:,1]
    p_array.particles[2::6] = xxstg[:,2]
    p_array.particles[3::6] = xxstg[:,3]
    p_array.particles[4::6] = xxstg[:,4]
    p_array.particles[5::6] = xxstg[:,5]
    return p_array, charge_array

def particleArray2astraBeam(p_array, charge_array, filename="tytest.ast"):
    gamref = p_array.E/m_e_GeV
    #print gamref
    P = p_array.particles.view()
    #energy = P[5]
    xp = exact_xxstg_2_xp(P, gamref)
    Pref = np.sqrt(p_array.E**2/m_e_GeV**2 - 1)*m_e_eV
    xp[0, 5] = xp[0, 5] + Pref # xp[0, 5] + p_array.E*1e9
    #np.savetxt('D:/pytest.ast',xp)
    charge_array = -charge_array.reshape(len(charge_array), 1)*1e+9 #charge in C -> in nC
    flag = np.zeros((len(charge_array),1))
    #print(np.shape(xp), np.shape(charge_array))
    astra = np.append(xp, flag, axis=1)
    astra = np.append(astra, charge_array, axis=1)
    astra = np.append(astra, flag, axis=1)
    astra = np.append(astra, flag, axis=1)
    np.savetxt(filename, astra)
