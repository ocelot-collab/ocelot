'''
@author: Igor Zagorodnov @ Martin Dohlus
Created on 27.03.2015
Revision on 01.06.2017: coordinate transform to the velocity direction
'''

import scipy.ndimage as ndimage
import numpy as np
from time import time
from ocelot.common.globals import *
from ocelot.cpbd.coord_transform import *
try:
    from pyfftw.interfaces.numpy_fft import fftn
    from pyfftw.interfaces.numpy_fft import ifftn
except:
    print("cs.py: module PYFFTW is not install. Install it if you want speed up your calculation")
    from numpy.fft import ifftn
    from numpy.fft import fftn

try:
    import numexpr as ne
    ne_flag = True
except:
    print("sc.py: module NUMEXPR is not installed. Install it if you want higher speed calculation.")
    ne_flag = False


def smooth_z(Zin, mslice):
    def myfunc(x, A):
        if x < 2*A:
            y = x-x*x/(4*A)
        else:
            y = A
        return y
    inds = np.argsort(Zin, axis=0)
    Zout = np.sort(Zin, axis=0)
    N = Zin.shape[0]
    S = np.zeros(N+1)
    S[N] = 0
    S[0] = 0
    for i in range(N):
        S[i+1] = S[i] + Zout[i]
    Zout2 = np.zeros(N)
    Zout2[N-1] = Zout[N-1]
    Zout2[0] = Zout[0]
    for i in range(1, N-1):
        m = min(i, N-i+1)
        m = np.floor(myfunc(0.5*m, 0.5*mslice)+0.500001)
        Zout2[i] = (S[i+m+1]-S[i-m])/(2*m+1)
    Zout[inds] = Zout2
    return Zout


class SpaceCharge():
    """
    The space charge forces are calculated by solving the Poisson equation in the bunch frame.
    Then the Lorentz transformed electromagnetic field is applied as a kick in the laboratory frame.
    For the solution of the Poisson equation we use an integral representation of the electrostatic potential
    by convolution of the free-space Green's function with the charge distribution.
    The convolution equation is solved with the help of the Fast Fourier Transform (FFT). The same algorithm for
    solution of the 3D Poisson equation is used, for example, in ASTRA
    """
    def __init__(self):
        #PhysicsProcess.__init__(self)
        self.step = 1 # in unit step
        self.nmesh_xyz = [31, 31, 31]
        self.low_order_kick = True

        self.start_elem = None
        self.end_elem = None
        self.debug = False

        self.random_seed = 1 # random seeding number. if None seeding is random

    def prepare(self, lat):
        pass

    def sym_kernel(self, ijk2, hxyz):
        i2 = ijk2[0]
        j2 = ijk2[1]
        k2 = ijk2[2]
        hx = hxyz[0]
        hy = hxyz[1]
        hz = hxyz[2]
        x = hx*np.r_[0:i2+1] - hx/2
        y = hy*np.r_[0:j2+1] - hy/2
        z = hz*np.r_[0:k2+1] - hz/2
        x, y, z = np.ix_(x, y, z)
        r = np.sqrt(x*x + y*y + z*z)
        if ne_flag:
            IG = ne.evaluate("(-x*x*0.5*arctan(y*z/(x*r)) + y*z*log(x+r) - y*y*0.5*arctan(z*x/(y*r)) + z*x*log(y+r) - z*z*0.5*arctan(x*y/(z*r)) + x*y*log(z+r))")
        else:
            IG = (-x*x*0.5*np.arctan(y*z/(x*r)) + y*z*np.log(x+r)
                  -y*y*0.5*np.arctan(z*x/(y*r)) + z*x*np.log(y+r)
                  -z*z*0.5*np.arctan(x*y/(z*r)) + x*y*np.log(z+r))

        kern = (IG[1:i2+1, 1:j2+1, 1:k2+1] - IG[0:i2, 1:j2+1, 1:k2+1]
               -IG[1:i2+1, 0:j2, 1:k2+1] + IG[0:i2, 0:j2, 1:k2+1]
               -IG[1:i2+1, 1:j2+1, 0:k2] + IG[0:i2, 1:j2+1, 0:k2]
               +IG[1:i2+1, 0:j2, 0:k2] - IG[0:i2, 0:j2, 0:k2])
        return kern

    def potential(self, q, steps):
        hx = steps[0]
        hy = steps[1]
        hz = steps[2]
        Nx = q.shape[0]
        Ny = q.shape[1]
        Nz = q.shape[2]
        out = np.zeros((2*Nx-1, 2*Ny-1, 2*Nz-1))
        out[:Nx, :Ny, :Nz] = q
        K1 = self.sym_kernel(q.shape, steps)
        K2 = np.zeros((2*Nx-1, 2*Ny-1, 2*Nz-1))
        K2[0:Nx, 0:Ny, 0:Nz] = K1
        K2[0:Nx, 0:Ny, Nz:2*Nz-1] = K2[0:Nx, 0:Ny, Nz-1:0:-1] #z-mirror
        K2[0:Nx, Ny:2*Ny-1,:] = K2[0:Nx, Ny-1:0:-1, :]        #y-mirror
        K2[Nx:2*Nx-1, :, :] = K2[Nx-1:0:-1, :, :]             #x-mirror
        t0 = time()
        out = np.real(ifftn(fftn(out)*fftn(K2)))
        t1 = time()
        if self.debug: print( 'fft time:', t1-t0, ' sec')
        out[:Nx, :Ny, :Nz] = out[:Nx,:Ny,:Nz]/(4*pi*epsilon_0*hx*hy*hz)
        return out[:Nx, :Ny, :Nz]

    def el_field(self, X, Q, gamma, nxyz):
        if self.random_seed != None:
            np.random.seed(self.random_seed)
        N = X.shape[0]
        X[:, 2] = X[:, 2]*gamma
        XX = np.max(X, axis=0)-np.min(X, axis=0)
        #XX = XX*np.random.uniform(low=1.0, high=1.1)
        if self.debug: print( 'mesh steps:', XX)
        # here we use a fast 3D "near-point" interpolation
        # we need a stand-alone module with 1D,2D,3D parricles-to-grid functions
        steps = XX/(nxyz-3)
        X = X/steps
        X_min = np.min(X, axis=0)
        X_mid = np.dot(Q, X)/np.sum(Q)
        X_off = np.floor(X_min-X_mid) + X_mid
        X = X - X_off
        nx = nxyz[0]
        ny = nxyz[1]
        nz = nxyz[2]
        nzny = nz*ny
        Xi = np.int_(np.floor(X)+1)
        inds = np.int_(Xi[:, 0]*nzny+Xi[:, 1]*nz+Xi[:, 2])  # 3d -> 1d
        q = np.bincount(inds, Q, nzny*nx).reshape(nxyz)
        p = self.potential(q, steps)
        Ex = np.zeros(p.shape)
        Ey = np.zeros(p.shape)
        Ez = np.zeros(p.shape)
        Ex[:nx-1, :, :] = (p[:nx-1, :, :] - p[1:nx, :, :])/steps[0]
        Ey[:, :ny-1, :] = (p[:, :ny-1, :] - p[:, 1:ny, :])/steps[1]
        Ez[:, :, :nz-1] = (p[:, :, :nz-1] - p[:, :, 1:nz])/steps[2]
        Exyz = np.zeros((N, 3))
        Exyz[:, 0] = ndimage.map_coordinates(Ex, np.c_[X[:, 0], X[:, 1]+0.5, X[:, 2]+0.5].T, order=1)*gamma
        Exyz[:, 1] = ndimage.map_coordinates(Ey, np.c_[X[:, 0]+0.5, X[:, 1], X[:, 2]+0.5].T, order=1)*gamma
        Exyz[:, 2] = ndimage.map_coordinates(Ez, np.c_[X[:, 0]+0.5, X[:, 1]+0.5, X[:, 2]].T, order=1)
        return Exyz


    def apply(self, p_array, zstep):

        nmesh_xyz = np.array(self.nmesh_xyz)
        gamref = p_array.E / m_e_GeV
        betref2 = 1 - gamref ** -2
        betref = np.sqrt(betref2)

        # MAD coordinates!!!
        # Lorentz transformation with V-axis and gamma_av
        xp = np.zeros(p_array.rparticles.shape)
        xp = xxstg_2_xp_mad(p_array.rparticles, xp, gamref)

        # coordinate transformation to the velocity direction
        t3 = np.mean(xp[3:6], axis=1)
        Pav = np.linalg.norm(t3)
        t3 = t3 / Pav
        ey = np.array([0, 1, 0])
        t1 = np.cross(ey, t3)
        t1 = t1 / np.linalg.norm(t1)
        t2 = np.cross(t3, t1)
        T = np.c_[t1, t2, t3]
        xyz = np.dot(xp[0:3].T, T)
        xp[3:6] = np.dot(xp[3:6].T, T).T

        # electric field in the rest frame of bunch
        gamma0 = sqrt((Pav / m_e_eV) ** 2 + 1)
        beta02 = 1 - gamma0 ** -2
        beta0 = np.sqrt(beta02)

        Exyz = self.el_field(xyz, p_array.q_array, gamma0, nmesh_xyz)

        # equations of motion in the lab system
        gamma = gamref * (1 + p_array.rparticles[5] * betref)
        Energy = m_e_eV * gamma
        betax = xp[3] / Energy
        betay = xp[4] / Energy
        betaz = xp[5] / Energy
        cdT = zstep / betref
        xp[3] = xp[3] + cdT * (1 - beta0 * betaz) * Exyz[:, 0]
        xp[4] = xp[4] + cdT * (1 - beta0 * betaz) * Exyz[:, 1]
        xp[5] = xp[5] + cdT * (Exyz[:, 2] + beta0 * (betax * Exyz[:, 0] + betay * Exyz[:, 1]))
        T = np.transpose(T)
        xp[3:6] = np.dot(xp[3:6].T, T).T
        xp_2_xxstg_mad(xp, p_array.rparticles, gamref)

    def SC_xp_update(self, xp, Q, gamref, dS, nxyz):
        #Lorentz transformation with z-axis and gamref
        betref2 = 1 - gamref**-2
        betref = np.sqrt(betref2)
        Eref = gamref*m_e_eV
        pref = Eref*betref
        Exyz = self.el_field(np.c_[xp[:, 0], xp[:, 1], xp[:, 2]], Q, gamref, nxyz)
        u = np.c_[xp[:, 3], xp[:, 4], xp[:, 5] + pref]
        gamma = np.sqrt(1+np.sum(u*u, 1)/m_e_eV**2).reshape((xp.shape[0], 1))
        cdT = dS/betref
        u = u/(gamma*m_e_eV)
        xp[:,3] = xp[:, 3] + cdT*(1-betref*u[:, 2])*Exyz[:, 0]
        xp[:,4] = xp[:, 4] + cdT*(1-betref*u[:, 2])*Exyz[:, 1]
        xp[:,5] = xp[:, 5] + cdT*(Exyz[:, 2] + betref*(u[:, 0]*Exyz[:, 0] + u[:, 1]*Exyz[:, 1]))


"""
def sc_track(lattice):
    navi = Navigator(lattice=lattice)
    #start = time.time()

    for i, zi in enumerate(Z[1:]):
        print zi
        dz = zi - Z[i]
        step(lat=lat, particle_list=p_array, dz=dz, navi=navi, order=order)
        if SC:
            #SC_xxstg_update(P, Q, p_array.E / 0.000511, dz, True, nxnynz)
            #sc.sc_apply(p_array, Q, dz, True, nxnynz)
            sc.sc_apply(p_array, q_array=Q, zstep=dz, nmesh_xyz=[63, 63, 63], low_order_kick=True)
        tw = get_envelope(p_array)
        tw.s = navi.z0
        tws_track.append(tw)
        #f.add_subplot(211)
        #plt.plot(p_array.particles[::6], p_array.particles[2::6], '.')
        #f.add_subplot(212)
        #plt.plot(p_array.particles[4::6],p_array.particles[5::6],'.')
        #plt.draw()
        #plt.pause(0.1)
    print "time exc = ", time.time() - start
"""
        

    
    
