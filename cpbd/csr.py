"""
@ author Martin Dohlus DESY, 2015
@ Python realization Sergey Tomin XFEL, 2016
"""
import numpy as np
from scipy import interpolate
from ocelot.common.globals import *
def interp1(x, y, xnew, k=1):
    tck = interpolate.splrep(x, y, k=k)
    ynew = interpolate.splev(xnew, tck, der=0)
    return ynew


def K0_inf_anf( i, traj, wmin):
    #%function [ w,KS ] = K0_inf_anf( i,traj,wmin )

    i1 = i-1 #% ignore points i1+1:i on linear path to observer
    ra = np.arange(0, i1+1)
    s = traj[0, ra] - traj[0, i]
    n = np.array([traj[1, i] - traj[1, ra],
                  traj[2, i] - traj[2, ra],
                  traj[3, i] - traj[3, ra]])

    R = np.sqrt(np.sum(n**2, axis=0))
    n = np.array([n[0, :]/R, n[1, :]/R, n[2, :]/R])

    w = s + R
    j = np.where(w<=wmin)[0]
    if len(j) > 0:
        j = j[-1]
        ra = np.arange(j, i1+1)
        w = w[ra]
        s = s[ra]

    #% kernel
    K = (n[0, ra]*(traj[4, ra] - traj[4, i]) +
         n[1, ra]*(traj[5, ra] - traj[5, i]) +
         n[2, ra]*(traj[6, ra] - traj[6, i]) -
       (1. - traj[4, ra]*traj[4, i] - traj[5, ra]*traj[5, i] - traj[6, ra]*traj[6, i]))/R[ra]

    #% integrated kernel: KS=int_s^0{K(u)*du}=int_0^{-s}{K(-u)*du}
    if np.shape(K)[0] > 1:
        a = np.append(0.5*(K[0:-1] + K[1:])*np.diff(s), 0.5*K[-1]*s[-1])
        KS = np.cumsum(a[::-1])[::-1]
        #KS = np.fliplr(np.cumsum(np.fliplr([0.5*(K[1:-1] + K[2:])*np.diff(s), 0.5*K[-1]*s[-1]])))
    else:
        KS = 0.5*K[-1]*s[-1]

    return w, KS

def K0_fin_anf(i, traj, wmin, gamma):
    #%function [ w,KS ] = K0_inf_anf( i,traj,wmin,gamma )

    g2i = 1./gamma**2
    b2 = 1. - g2i
    beta = np.sqrt(b2)
    i1 = i-1 #% ignore points i1+1:i on linear path to observer
    ra = np.arange(0, i1+1)
    s = traj[0, ra]-traj[0, i]
    n=np.array([traj[1, i] - traj[1, ra],
                traj[2, i] - traj[2, ra],
                traj[3, i] - traj[3, ra]])
    R = np.sqrt(np.sum(n**2, axis=0))
    n = np.array([n[0, :]/R, n[1, :]/R, n[2, :]/R])
    w = s + beta*R
    j = np.where(w<=wmin)[0]

    if len(j) > 0:
        j = j[-1]
        ra = np.arange(j, i1+1)
        w = w[ra]
        s = s[ra]

    #% kernel
    K = ((beta*(n[0, ra]*(traj[4, ra] - traj[4, i]) +
                n[1, ra]*(traj[5, ra] - traj[5, i]) +
                n[2, ra]*(traj[6, ra] - traj[6, i])) -
        b2*(1. - traj[4, ra]*traj[4, i] - traj[5, ra]*traj[5, i] - traj[6, ra]*traj[6, i]) - g2i)/R[ra] -
        (1. - beta*(n[0, ra]*traj[4, ra] + n[1, ra]*traj[5, ra] + n[2, ra]*traj[6, ra]))/w*g2i)

    #% integrated kernel: KS=int_s^0{K(u)*du}=int_0^{-s}{K(-u)*du}
    if len(K) > 1:
        a = np.append(0.5*(K[0:-1] + K[1:])*np.diff(s), 0.5*K[-1]*s[-1])
        KS = np.cumsum(a[::-1])[::-1]
    else:
        KS = 0.5*K[-1]*s[-1]

    return w, KS


def K0_fin_inf(i, traj, w_range, gamma):
    #%function [ KS ] = K0_inf_inf( i,traj,w_range,gamma )

    g2 = gamma**2
    g2i = 1./g2
    b2 = 1. - g2i
    beta = np.sqrt(b2)
    #% winf
    Rv1 = traj[1:4, i] - traj[1:4, 0]
    s1 =  traj[0, 0] - traj[0, i]
    ev1 = traj[4:, 0]
    evo = traj[4:, i]
    winfms1 = np.dot(Rv1, ev1)

    aup = -Rv1 + winfms1*ev1
    a2 = np.dot(aup, aup)
    a = np.sqrt(a2)
    uup = aup/a
    winf = s1 + winfms1
    s = winf + gamma*(gamma*(w_range - winf)-beta*np.sqrt(g2*(w_range-winf)**2+a2))
    R = (w_range-s)/beta
    if a2/R[1]**2 > 1e-7:
        KS = (beta*(1. - np.dot(ev1, evo))*np.log(R[0]/R) - beta*np.dot(uup, evo)*(np.arctan((s[0] - winf)/a) - np.arctan((s-winf)/a))
           - (b2*np.dot(ev1, evo) - 1)*np.log((winf - s + R)/(winf-s[0] + R[0]))
           + g2i*np.log(w_range[0]/w_range))
    else:
        KS = (beta*(1. - np.dot(ev1, evo))*np.log(R[0]/R)
           - (b2*np.dot(ev1, evo) - 1.)*np.log((winf - s + R)/(winf - s[0] + R[0]))
           + g2i*np.log(w_range[0]/w_range))

    return KS

def K0_inf_inf(i, traj, w_range):
    #%function [ KS ] = K0_inf_inf( i,traj,w_range )

    #% winf
    Rv1 = traj[1:4, i] - traj[1:4, 0]
    s1 =  traj[0, 0] - traj[0, i]
    ev1 = traj[4:, 0]
    evo = traj[4:, i]
    winfms1 = np.dot(Rv1, ev1)
    aup = -Rv1 + winfms1*ev1
    a2 = np.dot(aup, aup)
    a = np.sqrt(a2)
    uup = aup/a
    winf = s1 + winfms1

    Nvalid = np.where(winf<w_range)[0]
    if len(Nvalid) >0:
        Nvalid = Nvalid[0]
        w = w_range[Nvalid:]
        s = (winf+w)/2. + a2/2./(winf-w)
        #%K =(1-ev1'*evo)*((winf-s)./(w-s)-1)./(w-s)+aup'*evo./(w-s).^2;
        KS = (1. - np.dot(ev1, evo))*np.log(0.5*a2/(w-winf)/(w-s)) + np.dot(uup, evo)*(np.arctan((s-winf)/a) + np.pi/2.)
        KS = np.append(np.zeros(Nvalid), KS)
    else:
        KS = np.zeros(len(w_range))
    return KS


def CSR_K1( i, traj, NdW, gamma=None):
    """
    :param i: index of the trajectories points for the convolution kernel is calculated;
    :param traj: trajectory. traj[0,:] - longitudinal coordinate,
                             traj[1,:],traj[2,:],traj[3,:] - rectangular coordinates, \
                             traj[4,:],traj[5,:],traj[6,:] - tangential unit vectors
    :param NdW: list N[0] 0 number of mesh points, N[1] = dW> 0 - increment, Mesh = Mesh = (N: 0) * dW
    :param gamma:
    :return:
    """
    #%function [ K1 ] = CSR_K1( i,traj,NdW,gamma )

    #L_fin=nargin==4 && ~isempty(gamma) && gamma>1
    if gamma != None:
        L_fin = True
    else:
        L_fin = False
    w_range = np.arange(-NdW[0]-1, 0)*NdW[1]
    if L_fin:
        w, KS = K0_fin_anf(i, traj, w_range[0], gamma)
    else:
        w, KS = K0_inf_anf(i, traj, w_range[0])
    KS1 = KS[0]
    idx = np.argsort(w)
    w = w[idx]
    KS = KS[idx]
    w, idx = np.unique(w, return_index=True)
    KS = KS[idx]
    #% sort and unique takes time, but is required to avoid numerical trouble
    if w_range[0] < w[0]:
        m = np.where(w_range < w[0])[0][-1]
        if L_fin:
            KS2 = K0_fin_inf(i, traj, np.append(w_range[0:m+1], w[0]), gamma)
        else:
            KS2 = K0_inf_inf(i, traj, np.append(w_range[0:m+1], w[0]))
        KS2 = (KS2[-1] - KS2) + KS1
        KS = np.append(KS2[0:-1], interp1(w, KS, w_range[m+1:]))
    else:
        KS = interp1(w, KS, w_range)
    four_pi_eps0 = 1./(1e-7*speed_of_light**2)
    K1 = np.diff(np.append(np.diff(np.append(KS, 0)), 0))/NdW[1]/four_pi_eps0
    return K1











