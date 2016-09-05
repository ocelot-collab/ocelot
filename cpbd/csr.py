"""
@ authors Martin Dohlus DESY, 2015, Sergey Tomin XFEL, 2016
"""

import numpy as np
from scipy import interpolate
from ocelot.common.globals import *
from ocelot.common import math_op
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
from ocelot.cpbd.beam import *
from ocelot.cpbd.high_order import *
from ocelot.cpbd.magnetic_lattice import *


from matplotlib import pyplot as plt
from ocelot.cpbd.wake3D import *
import copy


def interp1(x, y, xnew, k=1):
    if len(xnew)>0:
        tck = interpolate.splrep(x, y, k=k)
        ynew = interpolate.splev(xnew, tck, der=0)
    else:
        ynew = []
    return ynew


class CSR:

    def __init__(self):
        self.mesh_params = None # [n_mesh, mesh_step]. list N[0] 0 number of mesh points, N[1] = dW> 0 - increment, Mesh = Mesh = (N: 0) * dW
        #self.csr_traj = np.transpose([[0, 0, 0, 0, 0, 0, 0]])
        self.start_elem = None
        self.end_elem = None
        self.z_csr_start = 0.  # z [m] position of the start_elem
        self.z0 = 0.           # self.z0 = navigator.z0 in track.track()

    def K0_inf_anf(self, i, traj, wmin):
        #%function [ w,KS ] = K0_inf_anf( i,traj,wmin )

        i1 = i-1 #% ignore points i1+1:i on linear path to observer
        ra = np.arange(0, i1+1)
        s = traj[0, ra] - traj[0, i]
        n = np.array([traj[1, i] - traj[1, ra],
                      traj[2, i] - traj[2, ra],
                      traj[3, i] - traj[3, ra]])

        R = np.sqrt(np.sum(n**2, axis=0))
        n = np.array([n[0, :]/R, n[1, :]/R, n[2, :]/R])
        #print(ra)
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

    def K0_fin_anf(self, i, traj, wmin, gamma):
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


    def K0_fin_inf(self, i, traj, w_range, gamma):
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
        #print(a, aup, a2)
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

    def K0_inf_inf(self, i, traj, w_range):
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
        #print(a, aup, a2)
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


    def CSR_K1(self, i, traj, NdW, gamma=None):
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
            w, KS = self.K0_fin_anf(i, traj, w_range[0], gamma)

        else:
            w, KS = self.K0_inf_anf(i, traj, w_range[0])
            #print("w=", len(w))
        #print("w=", len(w))
        KS1 = KS[0]
        idx = np.argsort(w)

        w = w[idx]
        #print("w=", len(w))
        KS = KS[idx]
        w, idx = np.unique(w, return_index=True)
        #print("w=", len(w))
        KS = KS[idx]
        #% sort and unique takes time, but is required to avoid numerical trouble
        if w_range[0] < w[0]:
            m = np.where(w_range < w[0])[0][-1]
            if L_fin:
                KS2 = self.K0_fin_inf(i, traj, np.append(w_range[0:m+1], w[0]), gamma)
            else:
                KS2 = self.K0_inf_inf(i, traj, np.append(w_range[0:m+1], w[0]))
            KS2 = (KS2[-1] - KS2) + KS1
            #print("lens", len(w), m, len(KS), len(w_range), len(w_range[m+1:]))
            KS = np.append(KS2[0:-1], interp1(w, KS, w_range[m+1:]))
        else:
            KS = interp1(w, KS, w_range)
        four_pi_eps0 = 1./(1e-7*speed_of_light**2)
        K1 = np.diff(np.append(np.diff(np.append(KS, 0)), 0))/NdW[1]/four_pi_eps0
        return K1

    def prepare(self, lat):
        """
        calculation of trajectory in rectangular coordinates
        calculation of the z_csr_start
        :param lat: Magnetic Lattice
        :return: self.csr_traj: trajectory. traj[0,:] - longitudinal coordinate,
                                 traj[1,:], traj[2,:], traj[3,:] - rectangular coordinates, \
                                 traj[4,:], traj[5,:], traj[6,:] - tangential unit vectors
        """
        csr_lat = MagneticLattice(copy.deepcopy(lat.sequence), start=self.start_elem, stop=self.end_elem)

        self.z_csr_start = sum([p.l for p in lat.sequence[:lat.sequence.index(self.start_elem)]])

        #Particle(x=0.0, y=0.0, px=0.0, py=0.0, s=0.0, p=0.0,  tau=0.0, E=0.0)
        p = Particle()
        #energy = p.E

        #beta = 1.#sqrt(1. - 1./gamma**2)
        self.csr_traj = np.transpose([[0, p.x, p.y, p.s, p.px, p.py, 1.-p.px*p.px/2. - p.py*p.py/2.]])

        #angle = 0.
        for elem in csr_lat.sequence:
            #print"elem.l = ", elem.l
            if elem.l == 0 :
                continue
            delta_s = elem.l
            #L = elem.l
            step = 0.0002
            if elem.__class__ in [Bend, RBend, SBend]:
                R = elem.l/elem.angle
                #B = energy*1e9*beta/(R*speed_of_light)
                R_vect = [0, -R, 0.]

                #L = R*np.sin(elem.angle)
                #angle += elem.angle
                #print("B = ", B)
            else:
                #B = 0.
                R_vect = [0, 0, 0.]
                #L = elem.l*np.cos(angle)
            self.csr_traj = arcline(self.csr_traj, delta_s, step, R_vect )

        return self.csr_traj

    #def apply(self, p_array, delta_s):
    #    pass

    def apply(self, p_array, delta_s):
        #print("APPLY CSR")
        s_cur = self.z0 - self.z_csr_start
        z = p_array.particles[4::6]
        #bins_start, hist_start = get_current(p_array, charge=p_array.q_array[0], num_bins=200)
        #plt.plot(bins_start, hist_start)
        #plt.show()
        bunch_size = max(z) - min(z)
        if self.mesh_params == None:
            n_mesh = 2000
            mesh_step = bunch_size/n_mesh*10.
            Ndw = [n_mesh, mesh_step]
        else:
            Ndw = self.mesh_params

        s_array = self.csr_traj[0,:]
        indx = (np.abs(s_array-s_cur)).argmin()
        #print(s_cur, indx, Ndw)
        gamma = p_array.E/m_e_GeV
        K1 = self.CSR_K1( indx, self.csr_traj, Ndw, gamma=gamma)
        #plt.plot(K1)
        #plt.show()

        # filtering
        Ns = np.ceil(bunch_size/2./Ndw[1])
        #s = np.arange(-Ns, Ns+1)*Ndw[1]
        #hist, bin_edges = np.histogram(z, bins=len(s))
        #b_distr = hist*p_array.q_array[0]/(bin_edges[1] - bin_edges[0])
        #b_distr_filt = gaussian_filter(b_distr, sigma=2)


        I = s2current(z, p_array.q_array, n_points=2*Ns+1, filter_order=5, mean_vel=speed_of_light)
        b_distr_filt = I[:, 1]/speed_of_light
        #print("charge1 = ", sum(b_distr_filt*charge_array[0]))
        #plt.plot(s, b_distr_filt, "r")
        #plt.plot(I[:, 0], I[:, 1]/speed_of_light, "b")
        #plt.show()
        #bins_start, hist_start = get_current(particle_list, charge=charge_array[0], num_bins=1000)
        #print("rms = ", np.std(z))
        lam_K1 = np.convolve(b_distr_filt, K1)

        Nend = len(lam_K1)
        x = np.arange(Ns+1-Nend, Ns+1)*Ndw[1]
        tck = interpolate.splrep(x, lam_K1, k=1)
        dE = interpolate.splev(z, tck, der=0)

        p_array.E = p_array.E * (1. + p_array.particles[5])
        p_array.particles[5::6] -= p_array.particles[5]
        #print("dE = ", particle_list.E)
        #plt.plot(x, lam_K1, "b")
        #print("E = ", particle_list.E)
        p_array.particles[5::6] += dE * 1e-9 / p_array.E * delta_s
        #plt.plot(z, particle_list.particles[5::6], "r.")
        #plt.plot(z, dE*1e-9/particle_list.E, "b.")
        #plt.show()







