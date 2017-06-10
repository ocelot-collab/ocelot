'''
Created on 17.05.2016
@author: Igor Zagorodnov
'''

from ocelot.adaptors import *
from ocelot.adaptors.astra2ocelot import *

def triang_filter(x, filter_order):
    Ns = x.shape[0]
    for i in range(filter_order):
        x[1:Ns] = (x[1:Ns] + x[0:Ns-1])*0.5
        x[0:Ns-1] = (x[1:Ns] + x[0:Ns-1])*0.5
    return x

def Der(x, y):
    #numerical derivative
    n=x.shape[0]
    dy = np.zeros(n)
    dy[1:n-1] = (y[2:n]-y[0:n-2])/(x[2:n]-x[0:n-2])
    dy[0] = (y[1]-y[0])/(x[1]-x[0])
    dy[n-1] = (y[n-1]-y[n-2])/(x[n-1]-x[n-2])
    return dy

def Int1(x, y):
    n = x.shape[0]
    Y = np.zeros(n)
    for i in range (1,n):
        Y[i] = Y[i-1]+0.5*(y(i)+y(i-1))*(x(i)-x(i-1))
    return Y

def Int1h(h, y):
    n = y.shape[0]
    Y = np.zeros(n)
    # slow, switch to vector operations to be done
    for i in range(1,n):
        Y[i] = Y[i-1] + 0.5*(y[i]+y[i-1])
    Y = Y*h
    return Y


def s2current(s_array, q_array, n_points, filter_order, mean_vel):
    """
    I = s2current(P0,q,Ns,NF)
    :param s_array: s-vector, coordinates in longitudinal direction
    :param q_array: charge-vector
    :param n_points: number of sampling points
    :param filter_order: filter order
    :param mean_vel: mean velocity
    :return:
    """
    s0 = np.min(s_array)
    s1 = np.max(s_array)
    NF2 = int(np.floor(filter_order / 2.))
    n_points = n_points + 2 * NF2
    Np = s_array.shape[0]
    ds = (s1 - s0) / (n_points - 2 - 2 * NF2)
    s = s0 + np.arange(-NF2, n_points - NF2) * ds
    # here we need a fast 1D linear interpolation of charges on the grid
    # in sc.py we use a fast 3D "near-point" interpolation
    # we need a stand-alone module with 1D,2D,3D particles-to-grid functions
    Ip = (s_array - s0) / ds
    I0 = np.floor(Ip)
    dI0 = Ip - I0
    I0 = I0 + NF2
    Ro = np.zeros(n_points)
    # slow, switch to vector operations to be done
    for i in range(Np):
        i0 = int(I0[i])
        di0 = dI0[i]
        Ro[i0] = Ro[i0] + (1 - di0) * q_array[i]
        Ro[i0 + 1] = Ro[i0 + 1] + di0 * q_array[i]
    if filter_order > 0:
        triang_filter(Ro, filter_order)
    I = np.zeros([n_points, 2])
    I[:, 0] = s
    I[:, 1] = Ro * mean_vel / ds
    return I


class WakeTable:
    """
    WakeTable(wake_file) - load and prepare wake table
    wake_file - path to the wake table
    """
    def __init__(self, wake_file):
        self.TH = self.load_wake_table(wake_file)

    def load_wake_table(self, wake_file):
        """
        :param wake_file: file name
        :return: (T, H): T- table of wakes coefs, H- matrix of the coefs place in T
        """
        W = np.loadtxt(wake_file)
        # head format %Nt 0 %N0 N1 %R L %C nm
        H = np.zeros([5, 5])
        Nt = int(W[0, 0])
        T = []
        ind = 0
        for i in range(Nt):
            ind = ind+1
            N0 = int(W[ind, 0])
            N1 = int(W[ind, 1])
            R = W[ind + 1, 0]
            L = W[ind + 1, 1]
            Cinv = W[ind + 2, 0]
            nm = int(W[ind+2, 1])
            n = int(np.floor(nm/10))
            m = int(nm - n*10)
            H[n, m] = i
            ind = ind + 2
            if N0 > 0:
                W0 = np.zeros([N0, 2])
                W0[0:N0, :] = W[ind+1:ind+N0+1, :]
                ind = ind + N0
            else:
                W0 = 0
            if N1 > 0:
                W1 = np.zeros([N1, 2])
                W1[0:N1, :]=W[ind+1:ind+N1+1, :]
                ind = ind + N1
            else:
                W1 = 0
            T = T + [(R, L, Cinv, nm, W0, N0, W1, N1)]
        return (T, H)


class Wake():
    """
    The wake field impact on the beam is included as series of kicks.
    In order to take into account the impact of the wake field on the beam the longitudinal wake function
    of point charge through the second order Taylor expansion is used.
    In general case it uses 13 one-dimensional functions to represent the  longitudinal component of the wake
    function for arbitrary sets of the source and the wittness particles near to the reference axis.

    parameters:
    -----------
    w_sampling = 500 -  defines the number of the equidistant sampling points for the one-dimensional
                        wake coefficients in the Taylor expansion of the 3D wake function.
    filter_order = 20 - smoothing filter order
    wake_table = None - wake table [WakeTable()]
    factor = 1. - scaling coefficient
    """
    def __init__(self):
        self.w_sampling = 500  # wake sampling
        self.filter_order = 20   # smoothing filter order
        #self.wake_file = ""
        self.wake_table = None
        self.factor = 1.
        self.step = 1

    def convolution(self, xu, u, xw, w):
        #convolution of equally spaced functions
        hx = xu[1] - xu[0]
        wc = np.convolve(u, w)*hx
        nw = w.shape[0]
        nu = u.shape[0]
        x0 = xu[0] + xw[0]
        xc = x0 + np.arange(nw + nu)*hx
        return xc, wc

    def wake_convolution(self, xb, bunch, xw, wake):
        #convolution of unequally spaced functions
        #bunch defines the parameters
        nb = xb.shape[0]
        xwi = np.zeros(nb)
        xwi = xb - xb[0]
        wake1 = np.interp(xwi, xw, wake, 0, 0)
        wake1[0] = wake1[0]*0.5
        xW, Wake = self.convolution(xb, bunch, xwi, wake1)
        return xW[0:nb], Wake[0:nb]

    def add_wake(self, I, T):
        """
        [x, W] = AddWake(I, T)
        :param I: wake table in V/C, W in V (R,L,Cinv,nm,W0,N0,W1,N1)
        :param T:
        :return:
        """
        """[x, W] =AddWake (I,T)
            T - wake table in V/C, W in V
            (R,L,Cinv,nm,W0,N0,W1,N1)"""
        R, L, Cinv, nm, W0, N0, W1, N1 = T
        c = speed_of_light
        x = I[:, 0]
        bunch = I[:, 1]
        if L != 0 or N1 > 0:
            d1_bunch=Der(x,bunch)
        nb=x.shape[0]
        W=np.zeros(nb)
        if N0 > 0:
            x, ww = self.wake_convolution(x, bunch, W0[:, 0], W0[:, 1])
            W = W-ww[0:nb]/c
        if N1>0:
            x, ww = self.wake_convolution(x, d1_bunch, W1[:, 0], W1[:, 1])
            #W = W - ww[0:nb]
            W = W + ww[0:nb]
        if R != 0:
            W = W-bunch*R
        if L != 0:
            W = W-d1_bunch*L*c
        if Cinv != 0:
          int_bunch = Int1(x, bunch)
          W = W - int_bunch*Cinv/c
        return x, W

    def add_total_wake(self, X, Y, Z, q, TH, Ns, NF):
        #function [Px Py Pz I00]=AddTotalWake (P,q,wakeFile,Ns,NF)
        T, H = TH
        c = speed_of_light
        #Z=-Z
        Np=X.shape[0]
        X2 = X**2
        Y2 = Y**2
        XY = X*Y
        #generalized currents;
        I00 = s2current(Z, q, Ns, NF, c)
        Nw=I00.shape[0]
        if (H[0,2]>0)or(H[2,3]>0)or(H[2,4]>0):
            qn=q*Y
            I01 = s2current(Z,qn,Ns,NF,c)
        if (H[0, 1] > 0)or(H[1, 3] > 0) or (H[1, 4] > 0):
            qn=q*X
            I10 = s2current(Z,qn,Ns,NF,c)
        if H[1,2]>0:
            qn=q*XY
            I11 = s2current(Z,qn,Ns,NF,c)
        if H[1,1]>0:
            qn=q*(X2-Y2)
            I20_02 = s2current(Z, qn, Ns, NF, c)
        #longitudinal wake
        #mn=0
        x, Wz = self.add_wake (I00, T[int(H[0, 0])])
        if H[0, 1] > 0:
            x, w = self.add_wake(I10, T[int(H[0, 1])])
            Wz = Wz+w
        if H[0,2]>0:
            x, w = self.add_wake(I01, T[int(H[0, 2])])
            Wz = Wz+w
        if H[1,1]>0:
            x, w = self.add_wake(I20_02, T[int(H[1, 1])])
            Wz = Wz+w
        if H[1,2]>0:
            x, w = self.add_wake(I11, T[int(H[1, 2])])
            Wz = Wz+2*w
        Pz = np.interp(Z, x, Wz, 0, 0)
        Py = np.zeros(Np)
        Px = np.zeros(Np)
        #mn=01
        Wz[0:Nw] = 0
        Wy = np.zeros(Nw)
        if H[0, 4] > 0:
            x, w = self.add_wake(I00, T[int(H[0, 4])])
            Wz=Wz+w
            Wy=Wy+w
        if H[1,4]>0:
            x, w = self.add_wake(I10, T[int(H[1, 4])])
            Wz = Wz + 2*w
            Wy = Wy + 2*w
        if H[2,4]>0:
            x, w = self.add_wake(I01, T[int(H[2, 4])])
            Wz = Wz + 2*w
            Wy = Wy + 2*w
        Pz = Pz + np.interp(Z, x, Wz, 0, 0)*Y
        h = x[1] - x[0]
        Wy = -Int1h(h, Wy)
        Py = Py + np.interp(Z, x, Wy, 0, 0)
        #mn=10
        Wz[0:Nw] = 0
        Wx = np.zeros(Nw)
        if H[0, 3] > 0:
            x, w = self.add_wake(I00, T[int(H[0, 3])])
            Wz = Wz + w
            Wx = Wx + w
        if H[1,3]>0:
            x, w = self.add_wake(I10, T[int(H[1, 3])])
            Wz = Wz + 2*w
            Wx = Wx + 2*w
        if H[2,3]>0:
            x, w = self.add_wake(I01, T[int(H[2, 3])])
            Wz = Wz + 2*w
            Wx = Wx + 2*w
        Wx=-Int1h(h,Wx)
        Pz = Pz + np.interp(Z, x, Wz, 0, 0)*X
        Px = Px + np.interp(Z, x, Wx, 0, 0)
        #mn=11
        if H[3,4]>0:
            x, w = self.add_wake(I00, T[int(H[3, 4])])
            Wx=-2*Int1h(h,w)
            p=np.interp(Z,x,Wx,0,0)
            Px = Px + p*Y
            Py = Py + p*X
            Pz = Pz + 2*np.interp(Z, x, w, 0, 0)*XY
        #mn=02,20
        if H[3,3]>0:
            x, w = self.add_wake(I00, T[int(H[3, 3])])
            Pz = Pz+np.interp(Z,x,w,0,0)*(X2-Y2)
            Wx = -2*Int1h(h,w)
            p = np.interp(Z,x,Wx,0,0)
            Px = Px + p*X
            Py = Py - p*Y
        I00[:,0]=-I00[:,0]
        #Z=-Z
        return Px, Py, Pz, I00

    def prepare(self, lat):
        #pass
        #self.TH = self.load_wake_table(self.wake_file)
        if self.wake_table == None:
            print("Wake.wake_table is None! Please specify the WakeTable()")
        else:
            self.TH = self.wake_table.TH

    def apply(self, p_array, dz):
        #print("apply: WAKE")
        #Px = 0
        #Py = 0
        #Pz = 0
        #ziw = zi - dz * 0.5
        #if (1.0 < ziw <= 3.0) or (5.0 < ziw <= 7.0):  # or(10.0<ziw<=12.0):
        ps = p_array.rparticles
        Px, Py, Pz, I00 = self.add_total_wake(ps[0], ps[2], ps[4], p_array.q_array, self.TH, self.w_sampling, self.filter_order)
        #if (3.0 < ziw <= 5.0):  # or(8.0<ziw<=10.0)or(12.0<ziw<=14.0):
        #    Px, Py, Pz, I00 = self.add_total_wake(Ps[:, 0], Ps[:, 2], Ps[:, 4], p_array.q_array, THh, Ns, NF)
        #print(zi, dz, ziw)

        p_array.rparticles[5] = p_array.rparticles[5] + Pz * dz*self.factor / (p_array.E * 1e9)
        p_array.rparticles[3] = p_array.rparticles[3] + Py * dz*self.factor / (p_array.E * 1e9)
        p_array.rparticles[1] = p_array.rparticles[1] + Px * dz*self.factor / (p_array.E * 1e9)


class WakeKick(Wake):
    def __init__(self):
        Wake.__init__(self)

    def apply(self, p_array, dz):
        #print("Apply WakeKick")
        ps = p_array.rparticles
        Px, Py, Pz, I00 = self.add_total_wake(ps[0], ps[2], ps[4], p_array.q_array, self.TH, self.w_sampling,
                                              self.filter_order)

        p_array.rparticles[5] = p_array.rparticles[5] + self.factor * Pz / (p_array.E * 1e9)
        p_array.rparticles[3] = p_array.rparticles[3] + self.factor * Py / (p_array.E * 1e9)
        p_array.rparticles[1] = p_array.rparticles[1] + self.factor * Px / (p_array.E * 1e9)