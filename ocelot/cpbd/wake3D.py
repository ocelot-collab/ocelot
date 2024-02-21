"""
Created on 17.05.2016
@author: Igor Zagorodnov
Added wake table WakeTableDechirperOffAxis on 11.2019
@authors: S. Tomin and I. Zagorodnov
"""

from ocelot.adaptors import *
from ocelot.adaptors.astra2ocelot import *
from ocelot.cpbd.physics_proc import PhysProc

import logging

_logger = logging.getLogger(__name__)

try:
    import numba as nb

    nb_flag = True
except:
    _logger.info("wake3D.py: module NUMBA is not installed. Install it to speed up calculation")
    nb_flag = False


def triang_filter(x, filter_order):
    Ns = x.shape[0]
    for i in range(filter_order):
        x[1:Ns] = (x[1:Ns] + x[0:Ns - 1]) * 0.5
        x[0:Ns - 1] = (x[1:Ns] + x[0:Ns - 1]) * 0.5
    return x


def Der(x, y):
    # numerical derivative
    n = x.shape[0]
    dy = np.zeros(n)
    dy[1:n - 1] = (y[2:n] - y[0:n - 2]) / (x[2:n] - x[0:n - 2])
    dy[0] = (y[1] - y[0]) / (x[1] - x[0])
    dy[n - 1] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2])
    return dy


def Int1(x, y):
    n = x.shape[0]
    Y = np.zeros(n)
    for i in range(1, n):
        Y[i] = Y[i - 1] + 0.5 * (y(i) + y(i - 1)) * (x(i) - x(i - 1))
    return Y


def Int1h(h, y):
    n = y.shape[0]
    Y = np.zeros(n)
    # slow, switch to vector operations to be done
    for i in range(1, n):
        Y[i] = Y[i - 1] + 0.5 * (y[i] + y[i - 1])
    Y = Y * h
    return Y


def project_on_grid_py(Ro, I0, dI0, q_array):
    """
    Simple function to project particles charge on grid

    :param Ro: grid
    :param I0: grid index for each particle
    :param dI0: coefficient how particle close to Ro[i+1]. Example: particle i, dI0=0.7 -> Ro[i]---------dI0[i]--Ro[i+1]
    :param q_array: charge array in [C]
    :return: Ro
    """
    Np = len(q_array)
    for i in range(Np):
        i0 = int(I0[i])
        di0 = dI0[i]
        Ro[i0] = Ro[i0] + (1 - di0) * q_array[i]
        Ro[i0 + 1] = Ro[i0 + 1] + di0 * q_array[i]
    return Ro


project_on_grid = project_on_grid_py if not nb_flag else nb.jit(project_on_grid_py)


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
    # with numba project charge on grid
    Ro = project_on_grid(Ro, I0, dI0, q_array)

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

    def __init__(self, wake_file=None):
        if wake_file is not None:
            self.TH = self.load_table(wake_file)

    def load_table(self, wake_file):
        wake_table = self.read_file(wake_file)
        TH = self.process_wake_table(wake_table)
        return TH

    def read_file(self, wake_file):
        W = np.loadtxt(wake_file)
        return W

    def process_wake_table(self, wake_table):
        """
        :param wake_file: file name
        :return: (T, H): T- table of wakes coefs, H - matrix of the coefs place in T
        """
        W = wake_table
        # head format %Nt 0 %N0 N1 %R L %C nm
        H = np.zeros([5, 5])
        Nt = int(W[0, 0])
        T = []
        ind = 0
        for i in range(Nt):
            ind = ind + 1
            N0 = int(W[ind, 0])
            N1 = int(W[ind, 1])
            R = W[ind + 1, 0]
            L = W[ind + 1, 1]
            Cinv = W[ind + 2, 0]
            nm = int(W[ind + 2, 1])
            n = int(np.floor(nm / 10))
            m = int(nm - n * 10)
            H[n, m] = i
            ind = ind + 2
            if N0 > 0:
                W0 = np.zeros([N0, 2])
                W0[0:N0, :] = W[ind + 1:ind + N0 + 1, :]
                ind = ind + N0
            else:
                W0 = 0
            if N1 > 0:
                W1 = np.zeros([N1, 2])
                W1[0:N1, :] = W[ind + 1:ind + N1 + 1, :]
                ind = ind + N1
            else:
                W1 = 0
            T = T + [(R, L, Cinv, nm, W0, N0, W1, N1)]
        return (T, H)



class WakeTable3(WakeTable):
    """
    WakeTable3(wake_file) - load and prepare wake table with 3rd order coefficents
    wake_file - path to the wake table
    """

    def __init__(self, wake_file=None):
        if wake_file is not None:
            self.TH = self.load_table(wake_file)

    def process_wake_table(self, wake_table):
        """
        :param wake_file: file name
        :return: (T, H): T- table of wakes coefs, H - matrix of the coefs place in T
        """
        W = wake_table
        # head format %Nt 0 %N0 N1 %R L %C nmk
        H = np.zeros([5, 5, 5])
        Nt = int(W[0, 0])
        T = []
        ind = 0
        for i in range(Nt):
            ind = ind + 1
            N0 = int(W[ind, 0])
            N1 = int(W[ind, 1])
            R = W[ind + 1, 0]
            L = W[ind + 1, 1]
            Cinv = W[ind + 2, 0]
            nmk = int(W[ind + 2, 1])
            n = int(np.floor(nmk / 100))
            nmk = nmk - 100 * n
            m = int(np.floor(nmk / 10))
            k = int(nmk - 10 * m)
            H[n, m, k] = i
            ind = ind + 2
            if N0 > 0:
                W0 = np.zeros([N0, 2])
                W0[0:N0, :] = W[ind + 1:ind + N0 + 1, :]
                ind = ind + N0
            else:
                W0 = 0
            if N1 > 0:
                W1 = np.zeros([N1, 2])
                W1[0:N1, :] = W[ind + 1:ind + N1 + 1, :]
                ind = ind + N1
            else:
                W1 = 0
            T = T + [(R, L, Cinv, nmk, W0, N0, W1, N1)]
        return (T, H)

class WakeTableDechirperOffAxis(WakeTable):
    """
    WakeTableDechirperOffAxis() - creates two wake tables for horizontal and vertical corrugated plates.
    Based on https://doi.org/10.1016/j.nima.2016.09.001 and SLAC-PUB-16881

    :param b: distance from the plate in [m]
    :param a: half gap between plates in [m]
    :param width: width of the corrugated structure in [m]
    :param t: longitudinal gap in [m]
    :param p: period of corrugation in [m]
    :param length: length of the corrugated structure in [m]
    :param sigma: characteristic (rms) longitudinal beam size in [m]
    :param orient: "horz" or "vert" plate orientation
    :return: hor_wake_table, vert_wake_table
    """

    def __init__(self, b=500 * 1e-6, a=0.01, width=0.02, t=0.25 * 1e-3, p=0.5 * 1e-3, length=1, sigma=30e-6,
                 orient="horz"):
        WakeTable.__init__(self)
        weke_horz, wake_vert = self.calculate_wake_tables(b=b, a=a, width=width, t=t, p=p, length=length, sigma=sigma)
        if orient == "horz":
            self.TH = self.process_wake_table(weke_horz)
        else:
            self.TH = self.process_wake_table(wake_vert)

    def calculate_wake_tables(self, b, a, width, t, p, length, sigma):
        """
        Function creates two wake tables for horizontal and vertical corrugated plates

        :param b: distance from the plate in [m]
        :param a: half gap between plates in [m]
        :param width: width of the corrugated structure in [m]
        :param t: longitudinal gap in [m]
        :param p: period of corrugation in [m]
        :param length: length of the corrugated structure in [m]
        :param sigma: characteristic longitudinal beam size in [m]]
        :param filename: save to files if filename is not None
        :return: hor_wake_table, vert_wake_table
        """
        p = p * 1e3  # m -> mm
        t = t * 1e3  # m -> mm
        L = length  # in m
        D = width * 1e3  # m -> mm
        a = a * 1e3  # m -> mm, in mm half gap
        b = b * 1e3  # m -> mm, distance from the plate
        sigma = sigma * 1e3  # m -> mm, characteristic beam size in mm
        y0 = a - b  # in mm POSITION of the charge CHANGE ONLY THIS PARAMETER
        # now distance from the plate 500um = a-y0
        c = speed_of_light
        Z0 = 3.767303134695850e+02
        y = y0
        x0 = D / 2.
        x = x0
        Nm = 300
        s = np.arange(0, 50 + 0.01, 0.01) * sigma
        ns = len(s)

        t2p = t / p
        alpha = 1 - 0.465 * np.sqrt(t2p) - 0.07 * t2p
        # s0r_fit=0.41*a**1.8*0.25^1.6/0.5^2.4;
        s0r_bane = a * a * t / (2 * np.pi * alpha ** 2 * p ** 2)
        s0r_igor = s0r_bane * np.pi / 4
        s0 = 4 * s0r_igor

        W = np.zeros(ns)

        dWdx0 = np.zeros(ns)
        dWdy0 = np.zeros(ns)
        dWdx = np.zeros(ns)
        dWdy = np.zeros(ns)
        ddWdx0dx0 = np.zeros(ns)  # h11
        ddWdy0dx0 = np.zeros(ns)  # h12=0
        ddWdxdx0 = np.zeros(ns)
        ddWdxdy0 = np.zeros(ns)  # h13 h23=0
        ddWdydx0 = np.zeros(ns)
        ddWdydy0 = np.zeros(ns)
        ddWdydx = np.zeros(ns)  # %h14=0 h24 h34=0

        Fz = np.zeros((Nm, ns))
        Fyd = np.zeros((Nm, ns))
        A = Z0 * c / (2 * a) * L
        for i in np.arange(Nm):
            m = i + 1
            M = np.pi / D * m
            X = M * a
            dx = np.sin(M * x0) * np.sin(M * x)
            # to avoid overflow
            coeff = X / (np.cosh(X) * np.sinh(X)) if X < 350. else 0.
            Wcc = A * coeff * np.exp(-(s / s0) ** 0.5 * (X / np.tanh(X)))
            Wss = A * coeff * np.exp(-(s / s0) ** 0.5 * (X * np.tanh(X)))

            Fz[i, :] = Wcc * np.cosh(M * y) * np.cosh(M * y0) + Wss * np.sinh(M * y) * np.sinh(M * y0)

            dW = Fz[i, :] * dx
            W = W + dW
            ddx0 = np.cos(M * x0) * np.sin(M * x)
            dWdx0 = dWdx0 + M * Fz[i, :] * ddx0
            ddy0 = Wcc * np.cosh(M * y) * np.sinh(M * y0) + Wss * np.sinh(M * y) * np.cosh(M * y0)
            dWdy0 = dWdy0 + M * ddy0 * dx
            ddx = np.sin(M * x0) * np.cos(M * x)
            dWdx = dWdx + M * Fz[i, :] * ddx
            ddy = Wcc * np.sinh(M * y) * np.cosh(M * y0) + Wss * np.cosh(M * y) * np.sinh(M * y0)
            dWdy = dWdy + M * ddy * dx
            Fyd[i, :] = M * ddy * dx
            ddWdx0dx0 = ddWdx0dx0 - M ** 2 * dW
            ddWdy0dx0 = ddWdy0dx0 + M ** 2 * (ddy0 * ddx0)
            ddWdxdx0 = ddWdxdx0 + M ** 2 * Fz[i, :] * np.cos(M * x0) * np.cos(M * x)
            ddWdxdy0 = ddWdxdy0 + M ** 2 * ddy0 * ddx
            ddWdydx0 = ddWdydx0 + M ** 2 * ddy * ddx0
            ddWdydy0 = ddWdydy0 + M ** 2 * (
                    Wcc * np.sinh(M * y) * np.sinh(M * y0) + Wss * np.cosh(M * y) * np.cosh(M * y0)) * dx
            ddWdydx = ddWdydx + M ** 2 * ddy * ddx

        W = W * 2 / D * 1e6
        h00 = W
        dWdx0 = dWdx0 * 2 / D * 1e9
        h01 = dWdx0
        dWdy0 = dWdy0 * 2 / D * 1e9
        h02 = dWdy0
        dWdx = dWdx * 2 / D * 1e9
        h03 = dWdx
        dWdy = dWdy * 2 / D * 1e9
        h04 = dWdy
        ddWdx0dx0 = ddWdx0dx0 * 2 / D * 1e12
        h11 = ddWdx0dx0 * 0.5
        ddWdy0dx0 = ddWdy0dx0 * 2 / D * 1e12
        h12 = ddWdy0dx0 * 0.5
        ddWdxdx0 = ddWdxdx0 * 2 / D * 1e12
        h13 = ddWdxdx0 * 0.5
        ddWdydx0 = ddWdydx0 * 2 / D * 1e12
        h14 = ddWdydx0 * 0.5
        ddWdxdy0 = ddWdxdy0 * 2 / D * 1e12
        h23 = ddWdxdy0 * 0.5
        ddWdydy0 = ddWdydy0 * 2 / D * 1e12
        h24 = ddWdydy0 * 0.5
        ddWdydx = ddWdydx * 2 / D * 1e12
        h34 = ddWdydx * 0.5
        h33 = h11

        s = s * 1e-3
        N = len(s)
        out00 = np.append([[N, 0], [0, 0], [0, 0]], np.vstack((s, h00)).T, axis=0)
        out02 = np.append([[N, 0], [0, 0], [0, 2]], np.vstack((s, h02)).T, axis=0)
        out04 = np.append([[N, 0], [0, 0], [0, 4]], np.vstack((s, h04)).T, axis=0)
        out11 = np.append([[N, 0], [0, 0], [0, 11]], np.vstack((s, h11)).T, axis=0)
        out13 = np.append([[N, 0], [0, 0], [0, 13]], np.vstack((s, h13)).T, axis=0)
        out24 = np.append([[N, 0], [0, 0], [0, 24]], np.vstack((s, h24)).T, axis=0)
        out33 = np.append([[N, 0], [0, 0], [0, 33]], np.vstack((s, h33)).T, axis=0)
        wake_horz = np.vstack(([[7, 0]], out00, out02, out04, out11, out13, out24, out33))

        out11 = np.append([[N, 0], [0, 0], [0, 11]], np.vstack((s, -h11)).T, axis=0)
        out33 = np.append([[N, 0], [0, 0], [0, 33]], np.vstack((s, -h33)).T, axis=0)
        out01 = np.append([[N, 0], [0, 0], [0, 1]], np.vstack((s, h02)).T, axis=0)
        out03 = np.append([[N, 0], [0, 0], [0, 3]], np.vstack((s, h04)).T, axis=0)
        out13 = np.append([[N, 0], [0, 0], [0, 13]], np.vstack((s, h24)).T, axis=0)
        out24 = np.append([[N, 0], [0, 0], [0, 24]], np.vstack((s, h13)).T, axis=0)
        wake_vert = np.vstack(([[7, 0]], out00, out01, out03, out11, out13, out24, out33))
        return wake_horz, wake_vert

class WakeTableParallelPlate_origin(WakeTable):
    """
    WakeTableParallelPlate_origin() - create wake tables for parallel plate corrugated 
    strcutres based on 0th order analytical results.
    """
    def __init__(self, b=500e-6, a=500e-6, t=250e-6, p=500e-6, length=1, sigma=30e-6, orient="horz"):
        WakeTable.__init__(self)
        wake_horz, wake_vert = self.calculate_wake_table(b, a, t, p, length, sigma)
        if orient == "horz":
            self.TH = self.process_wake_table(wake_horz)
        else:
            self.TH = self.process_wake_table(wake_vert)
       
    def calculate_wake_table(self, b, a, t, p, length, sigma):
        """
        :param b: distance from the +y/+x plates [m], offset of beam from center = a - b
        :param a: half gap between plates [m],   0 < b < 2a
        :param t: longitudinal gap of corrugations [m]
        :param p: period of corrugations [m]
        :param length: length of the strucutre [m]
        :param sigma: longitudinal beam size in [m], to calculate s window length
        :param orient: structure orientation: "horz" kick in y, or "vert" kick in x
        :return: hor_wake_table, vert_wake_table
        """
        offset = a - b
        if np.abs(offset) >= a:
            raise ValueError("distance to +y/+x plate must be smaller than 2a")
        c = speed_of_light
        Z0 = 3.767303134695850e+02
        s = np.arange(0, 50 + 0.01, 0.01) * sigma
        ns = len(s)
        
        t2p = t / p
        alpha = 1 - 0.465 * np.sqrt(t2p) - 0.07 * t2p
        s0r = a * a * t / (2 * np.pi * alpha * alpha * p * p )   # Bane s0r      
        
        # convert coefficient from cgs unit to mks unit, scale by length 
        mks = Z0 * c / (4 * np.pi) *length

        # initialize
        h00 = np.zeros(ns)
        h01 = np.zeros(ns)
        h02 = np.zeros(ns)
        h03 = np.zeros(ns)
        h04 = np.zeros(ns)
        h11 = np.zeros(ns)
        h12 = np.zeros(ns)
        h13 = np.zeros(ns)
        h14 = np.zeros(ns)
        h23 = np.zeros(ns)
        h24 = np.zeros(ns)
        h33 = np.zeros(ns)
        h34 = np.zeros(ns)
        
        # trig function
        csc = lambda x : 1 / np.sin(x)
        sec = lambda x : 1 / np.cos(x)
        cot = lambda x : 1 / np.tan(x)
        
        
        # scale factors, hab
        Y = pi * offset / (2 * a)
        
        if Y == 0:
            sl = 9/4 * s0r
            sd = (15/14)**2 * s0r
            sq = (15/16)**2 * s0r
            
            h00 = mks * np.pi**2 / (4 * a**2)  * np.ones(s.shape)
            h11 = mks * (-1) * np.pi**4 / (64 * a**4)  * np.ones(s.shape)
            h13 = -h11
            h24 = mks * np.pi**4 / (64 * a**4)  * np.ones(s.shape)
            h33 = h11
        else:            
            sl = 4 * s0r * (1 + np.cos(Y)**2 / 3 + Y*np.tan(Y)) ** (-2) 
            sm = 4 * s0r * (1.5 - Y*cot(2*Y) + Y*csc(Y)*sec(Y) ) ** (-2)
            sd = 4 * s0r * ( (64 + np.cos(2*Y))/30 + 2*Y*np.tan(Y) + (0.3 - Y*np.sin(2*Y))/(np.cos(2*Y) - 2) ) ** (-2)
            sq = 4 * s0r * ( (56 - np.cos(2*Y))/30 + 2*Y*np.tan(Y) - (0.3 + Y*np.sin(2*Y))/(np.cos(2*Y) - 2) ) ** (-2)
            
            h00 = mks * np.pi**2 / (4 * a**2) * sec(Y)**2 *np.ones(s.shape)
            h02 = mks * np.pi**3 / (16 * a**3) * np.sin(2*Y) * sec(Y)**4  * np.ones(s.shape)
            h04 = h02
            h11 = mks * np.pi**4 / (64 * a**4) * ( np.cos(2*Y) - 2 ) * sec(Y)**4  * np.ones(s.shape)
            h13 = -h11
            h24 = mks * np.pi**4 / (64 * a**4) * ( 2 - np.cos(2*Y) ) * sec(Y)**4  * np.ones(s.shape)
            h33 = h11
        
        N = len(s)
        # coefficients, horizontal, kick is in y 
        out00 = np.append([[N, 0], [0, 0], [0, 0]], np.vstack((s, h00)).T, axis=0)
        out02 = np.append([[N, 0], [0, 0], [0, 2]], np.vstack((s, h02)).T, axis=0)
        out04 = np.append([[N, 0], [0, 0], [0, 4]], np.vstack((s, h04)).T, axis=0)
        out11 = np.append([[N, 0], [0, 0], [0, 11]], np.vstack((s, h11)).T, axis=0)
        out13 = np.append([[N, 0], [0, 0], [0, 13]], np.vstack((s, h13)).T, axis=0)
        out24 = np.append([[N, 0], [0, 0], [0, 24]], np.vstack((s, h24)).T, axis=0)
        out33 = np.append([[N, 0], [0, 0], [0, 33]], np.vstack((s, h33)).T, axis=0)
        wake_horz = np.vstack(([[7, 0]], out00, out02,  out04, out11, out13, out24, out33))
        
        # coefficients, vertical, swap x -> y, x0 -> y0, i.e., 1 -> 2, 3 -> 4, keep 0
        # note h22 = -h11, h44 = -h33
        out00 = np.append([[N, 0], [0, 0], [0, 0]], np.vstack((s, h00)).T, axis=0)
        out01 = np.append([[N, 0], [0, 0], [0, 1]], np.vstack((s, h02)).T, axis=0)
        out03 = np.append([[N, 0], [0, 0], [0, 3]], np.vstack((s, h04)).T, axis=0)
        out11 = np.append([[N, 0], [0, 0], [0, 11]], np.vstack((s, -h11)).T, axis=0)
        out13 = np.append([[N, 0], [0, 0], [0, 13]], np.vstack((s, h24)).T, axis=0)
        out24 = np.append([[N, 0], [0, 0], [0, 24]], np.vstack((s, h13)).T, axis=0)
        out33 = np.append([[N, 0], [0, 0], [0, 33]], np.vstack((s, -h33)).T, axis=0)
        wake_vert = np.vstack(([[7, 0]], out00, out01, out03, out11, out13, out24, out33))
        return wake_horz, wake_vert

class WakeTableParallelPlate(WakeTable):
    """
    WakeTableParallelPlate() - create wake tables for parallel plate corrugated 
    strcutres based on 1st order analytical results.
    """
    def __init__(self, b=500e-6, a=500e-6, t=250e-6, p=500e-6, length=1, sigma=30e-6, orient="horz"):
        WakeTable.__init__(self)
        wake_horz, wake_vert = self.calculate_wake_table(b, a, t, p, length, sigma)
        if orient == "horz":
            self.TH = self.process_wake_table(wake_horz)
        else:
            self.TH = self.process_wake_table(wake_vert)
       
    def calculate_wake_table(self, b, a, t, p, length, sigma):
        """
        :param b: distance from the +y/+x plates [m], offset of beam from center = a - b
        :param a: half gap between plates [m],   0<b < 2a
        :param t: longitudinal gap of corrugations [m]
        :param p: period of corrugations [m]
        :param length: length of the strucutre [m]
        :param sigma: longitudinal beam size in [m], to calculate s window length
        :param orient: structure orientation: "horz" kick in y, or "vert" kick in x
        :return: hor_wake_table, vert_wake_table
        """
        offset = a - b
        if np.abs(offset) >= a:
            raise ValueError("distance to +y/+x plate must be smaller than 2a")
        c = speed_of_light
        Z0 = 3.767303134695850e+02
        s = np.arange(0, 50 + 0.01, 0.01) * sigma
        ns = len(s)
        
        t2p = t / p
        alpha = 1 - 0.465 * np.sqrt(t2p) - 0.07 * t2p
        s0r = a * a * t / (2 * np.pi * alpha * alpha * p * p )   # Bane s0r      
        
        # convert coefficient from cgs unit to mks unit, scale by length 
        mks = Z0 * c / (4 * np.pi) *length

        # initialize
        h00 = np.zeros(ns)
        h01 = np.zeros(ns)
        h02 = np.zeros(ns)
        h03 = np.zeros(ns)
        h04 = np.zeros(ns)
        h11 = np.zeros(ns)
        h12 = np.zeros(ns)
        h13 = np.zeros(ns)
        h14 = np.zeros(ns)
        h23 = np.zeros(ns)
        h24 = np.zeros(ns)
        h33 = np.zeros(ns)
        h34 = np.zeros(ns)
        
        # trig function
        csc = lambda x : 1 / np.sin(x)
        sec = lambda x : 1 / np.cos(x)
        cot = lambda x : 1 / np.tan(x)
        
        
        # scale factors, hab
        Y = pi * offset / (2 * a)
        
        if Y == 0:
            sl = 9/4 * s0r
            sd = (15/14)**2 * s0r
            sq = (15/16)**2 * s0r
            
            h00 = mks * np.pi**2 / (4 * a**2) *  np.exp( - np.sqrt(s/sl) )
            h11 = mks * (-1) * np.pi**4 / (64 * a**4) * np.exp( - np.sqrt(s/sq) )
            h13 = -h11
            h24 = mks * np.pi**4 / (64 * a**4) * np.exp( - np.sqrt(s/sd) )
            h33 = h11
        else:            
            sl = 4 * s0r * (1 + np.cos(Y)**2 / 3 + Y*np.tan(Y)) ** (-2)
            sm = 4 * s0r * (1.5 - Y*cot(2*Y) + Y*csc(Y)*sec(Y) ) ** (-2)
            sd = 4 * s0r * ( (64 + np.cos(2*Y))/30 + 2*Y*np.tan(Y) + (0.3 - Y*np.sin(2*Y))/(np.cos(2*Y) - 2) ) ** (-2)
            sq = 4 * s0r * ( (56 - np.cos(2*Y))/30 + 2*Y*np.tan(Y) - (0.3 + Y*np.sin(2*Y))/(np.cos(2*Y) - 2) ) ** (-2)
            
            h00 = mks * np.pi**2 / (4 * a**2) * sec(Y)**2 * np.exp( - np.sqrt(s/sl) )
            h02 = mks * np.pi**3 / (16 * a**3) * np.sin(2*Y) * sec(Y)**4 * np.exp( - np.sqrt(s/sm) )
            h04 = h02
            h11 = mks * np.pi**4 / (64 * a**4) * ( np.cos(2*Y) - 2 ) * sec(Y)**4 * np.exp( - np.sqrt(s/sq) )
            h13 = -h11
            h24 = mks * np.pi**4 / (64 * a**4) * ( 2 - np.cos(2*Y) ) * sec(Y)**4 * np.exp( - np.sqrt(s/sd) )
            h33 = h11
        
        
        N = len(s)
        # coefficients, horizontal, kick is in y 
        out00 = np.append([[N, 0], [0, 0], [0, 0]], np.vstack((s, h00)).T, axis=0)
        out02 = np.append([[N, 0], [0, 0], [0, 2]], np.vstack((s, h02)).T, axis=0)
        out04 = np.append([[N, 0], [0, 0], [0, 4]], np.vstack((s, h04)).T, axis=0)
        out11 = np.append([[N, 0], [0, 0], [0, 11]], np.vstack((s, h11)).T, axis=0)
        out13 = np.append([[N, 0], [0, 0], [0, 13]], np.vstack((s, h13)).T, axis=0)
        out24 = np.append([[N, 0], [0, 0], [0, 24]], np.vstack((s, h24)).T, axis=0)
        out33 = np.append([[N, 0], [0, 0], [0, 33]], np.vstack((s, h33)).T, axis=0)
        wake_horz = np.vstack(([[7, 0]], out00, out02,  out04, out11, out13, out24, out33))
        
        # coefficients, vertical, swap x -> y, x0 -> y0, i.e., 1 -> 2, 3 -> 4, keep 0
        # note h22 = -h11, h44 = -h33
        out00 = np.append([[N, 0], [0, 0], [0, 0]], np.vstack((s, h00)).T, axis=0)
        out01 = np.append([[N, 0], [0, 0], [0, 1]], np.vstack((s, h02)).T, axis=0)
        out03 = np.append([[N, 0], [0, 0], [0, 3]], np.vstack((s, h04)).T, axis=0)
        out11 = np.append([[N, 0], [0, 0], [0, 11]], np.vstack((s, -h11)).T, axis=0)
        out13 = np.append([[N, 0], [0, 0], [0, 13]], np.vstack((s, h24)).T, axis=0)
        out24 = np.append([[N, 0], [0, 0], [0, 24]], np.vstack((s, h13)).T, axis=0)
        out33 = np.append([[N, 0], [0, 0], [0, 33]], np.vstack((s, -h33)).T, axis=0)
        wake_vert = np.vstack(([[7, 0]], out00, out01, out03, out11, out13, out24, out33))
        return wake_horz, wake_vert
    
class WakeTableParallelPlate3_origin(WakeTable3):
    """
    WakeTableParallelPlate() - create 3rd order wake tables for parallel plate corrugated 
    strcutres based on 0th order analytical results.
    """
    def __init__(self, b=500e-6, a=500e-6, t=250e-6, p=500e-6, length=1, sigma=30e-6, orient="horz"):
        WakeTable3.__init__(self)
        wake_horz, wake_vert = self.calculate_wake_table(b, a, t, p, length, sigma)
        if orient == "horz":
            self.TH = self.process_wake_table(wake_horz)
        else:
            self.TH = self.process_wake_table(wake_vert)
       
    def calculate_wake_table(self, b, a, t, p, length, sigma):
        """
        :param b: distance from the +y/+x plates [m], offset of beam from center = a - b
        :param a: half gap between plates [m],   0<b < 2a
        :param t: longitudinal gap of corrugations [m]
        :param p: period of corrugations [m]
        :param length: length of the strucutre [m]
        :param sigma: longitudinal beam size in [m], to calculate s window length
        :param orient: structure orientation: "horz" kick in y, or "vert" kick in x
        :return: hor_wake_table, vert_wake_table
        """
        offset = a - b
        if np.abs(offset) >= a:
            raise ValueError("distance to +y/+x plate must be smaller than 2a")
        c = speed_of_light
        Z0 = 3.767303134695850e+02
        s = np.arange(0, 50 + 0.01, 0.01) * sigma
        ns = len(s)
        
        t2p = t / p
        alpha = 1 - 0.465 * np.sqrt(t2p) - 0.07 * t2p
        s0r = a * a * t / (2 * np.pi * alpha * alpha * p * p )   # Bane s0r      
        
        # convert coefficient from cgs unit to mks unit, scale by length 
        mks = Z0 * c / (4 * np.pi) *length

        # initialize
        h000 = np.zeros(ns)
        h001 = np.zeros(ns)
        h002 = np.zeros(ns)
        h003 = np.zeros(ns)
        h004 = np.zeros(ns)
        h011 = np.zeros(ns)
        h012 = np.zeros(ns)
        h013 = np.zeros(ns)
        h014 = np.zeros(ns)
        h023 = np.zeros(ns)
        h024 = np.zeros(ns)
        h033 = np.zeros(ns)
        h034 = np.zeros(ns)
        
        # higher order
        h111 = np.zeros(ns)
        h112 = np.zeros(ns)
        h113 = np.zeros(ns)
        h114 = np.zeros(ns)
        h122 = np.zeros(ns)
        h123 = np.zeros(ns)
        h124 = np.zeros(ns)
        h133 = np.zeros(ns)
        h134 = np.zeros(ns)
        h144 = np.zeros(ns)
        h222 = np.zeros(ns)
        h223 = np.zeros(ns)
        h224 = np.zeros(ns)
        h233 = np.zeros(ns)
        h234 = np.zeros(ns)
        h244 = np.zeros(ns)
        h333 = np.zeros(ns)
        h334 = np.zeros(ns)
        h344 = np.zeros(ns)
        h444 = np.zeros(ns)
        
        # trig function
        csc = lambda x : 1 / np.sin(x)
        sec = lambda x : 1 / np.cos(x)
        cot = lambda x : 1 / np.tan(x)
        
        
        # scale factors, hab
        Y = pi * offset / (2 * a)
        
        if Y == 0:
            sl = 9/4 * s0r
            sd = (15/14)**2 * s0r
            sq = (15/16)**2 * s0r
            
            h000 = mks * np.pi**2 / (4 * a**2) * np.ones(s.shape)
            h011 = mks * (-1) * np.pi**4 / (64 * a**4) *  np.ones(s.shape)
            h013 = -h011
            h024 = mks * np.pi**4 / (64 * a**4) * np.ones(s.shape)
            h033 = h011
        else:            
            sl = 4 * s0r * (1 + np.cos(Y)**2 / 3 + Y*np.tan(Y)) ** (-2)
            sm = 4 * s0r * (1.5 - Y*cot(2*Y) + Y*csc(Y)*sec(Y) ) ** (-2)
            sd = 4 * s0r * ( (64 + np.cos(2*Y))/30 + 2*Y*np.tan(Y) + (0.3 - Y*np.sin(2*Y))/(np.cos(2*Y) - 2) ) ** (-2)
            sq = 4 * s0r * ( (56 - np.cos(2*Y))/30 + 2*Y*np.tan(Y) - (0.3 + Y*np.sin(2*Y))/(np.cos(2*Y) - 2) ) ** (-2)
            
            h000 = mks * np.pi**2 / (4 * a**2) * sec(Y)**2 *  np.ones(s.shape)
            h002 = mks * np.pi**3 / (16 * a**3) * np.sin(2*Y) * sec(Y)**4 * np.ones(s.shape)
            h004 = h002
            h011 = mks * np.pi**4 / (64 * a**4) * ( np.cos(2*Y) - 2 ) * sec(Y)**4 * np.ones(s.shape)
            h013 = -h011
            h024 = mks * np.pi**4 / (64 * a**4) * ( 2 - np.cos(2*Y) ) * sec(Y)**4 * np.ones(s.shape)
            h033 = h011
        
            # higher order
            st = 16 * s0r * ( 5*Y*np.tan(Y) + Y*cot(Y) + ( 2*Y*np.sin(2*Y) )/(5-np.cos(2*Y)) + 5 ) ** (-2)
            h112 = mks * np.pi**5 / (768 * a**5) * sec(Y)**5 * ( np.sin(3*Y) - 11*np.sin(Y) ) * np.ones(s.shape)
            h114 = h112
            h233 = h112
            h334 = h112
            h123 = -h112
            h134 = -h112
            h222 = -h112
            h224 = -h112
            h244 = -h112
            h444 = -h112

        
        N = len(s)
        # full coefficients, horizontal, kick is in y 
        out000 = np.append([[N, 0], [0, 0], [0, 0]], np.vstack((s, h000)).T, axis=0)
        out002 = np.append([[N, 0], [0, 0], [0, 2]], np.vstack((s, h002)).T, axis=0)
        out004 = np.append([[N, 0], [0, 0], [0, 4]], np.vstack((s, h004)).T, axis=0)
        out011 = np.append([[N, 0], [0, 0], [0, 11]], np.vstack((s, h011)).T, axis=0)
        out013 = np.append([[N, 0], [0, 0], [0, 13]], np.vstack((s, h013)).T, axis=0)
        out024 = np.append([[N, 0], [0, 0], [0, 24]], np.vstack((s, h024)).T, axis=0)
        out033 = np.append([[N, 0], [0, 0], [0, 33]], np.vstack((s, h033)).T, axis=0)
        
        out112 = np.append([[N, 0], [0, 0], [0, 112]], np.vstack((s, h112)).T, axis=0)
        out114 = np.append([[N, 0], [0, 0], [0, 114]], np.vstack((s, h114)).T, axis=0)
        out233 = np.append([[N, 0], [0, 0], [0, 233]], np.vstack((s, h233)).T, axis=0)
        out334 = np.append([[N, 0], [0, 0], [0, 334]], np.vstack((s, h334)).T, axis=0)
        
        out123 = np.append([[N, 0], [0, 0], [0, 123]], np.vstack((s, h123)).T, axis=0)
        out134 = np.append([[N, 0], [0, 0], [0, 134]], np.vstack((s, h134)).T, axis=0)
        out222 = np.append([[N, 0], [0, 0], [0, 222]], np.vstack((s, h222)).T, axis=0)
        out224 = np.append([[N, 0], [0, 0], [0, 224]], np.vstack((s, h224)).T, axis=0)
        out244 = np.append([[N, 0], [0, 0], [0, 244]], np.vstack((s, h244)).T, axis=0)
        out444 = np.append([[N, 0], [0, 0], [0, 444]], np.vstack((s, h444)).T, axis=0)
        

        wake_horz = np.vstack(([[17, 0]], out000,  out002,  out004, out011,  out013, out024, out033,
                               out112, out114, out233, out334, out123, out134, out222, out224, out244, out444
                               ))
        # full coefficients, vertical, swap x -> y, x0 -> y0, i.e., 1 -> 2, 3 -> 4, keep 0
        # note h22 = -h11, h44 = -h33
        out000 = np.append([[N, 0], [0, 0], [0, 0]], np.vstack((s, h000)).T, axis=0)
        out001 = np.append([[N, 0], [0, 0], [0, 1]], np.vstack((s, h002)).T, axis=0)
        out003 = np.append([[N, 0], [0, 0], [0, 3]], np.vstack((s, h004)).T, axis=0)
        out011 = np.append([[N, 0], [0, 0], [0, 11]], np.vstack((s, -h011)).T, axis=0)
        out024 = np.append([[N, 0], [0, 0], [0, 24]], np.vstack((s, h013)).T, axis=0)
        out013 = np.append([[N, 0], [0, 0], [0, 13]], np.vstack((s, h024)).T, axis=0)
        out033 = np.append([[N, 0], [0, 0], [0, 33]], np.vstack((s, -h033)).T, axis=0)
        
        out122 = np.append([[N, 0], [0, 0], [0, 122]], np.vstack((s, h112)).T, axis=0)
        out223 = np.append([[N, 0], [0, 0], [0, 223]], np.vstack((s, h114)).T, axis=0)
        out144 = np.append([[N, 0], [0, 0], [0, 144]], np.vstack((s, h233)).T, axis=0)
        out344 = np.append([[N, 0], [0, 0], [0, 344]], np.vstack((s, h334)).T, axis=0)
        
        out124 = np.append([[N, 0], [0, 0], [0, 124]], np.vstack((s, h123)).T, axis=0)
        out234 = np.append([[N, 0], [0, 0], [0, 234]], np.vstack((s, h134)).T, axis=0)
        out111 = np.append([[N, 0], [0, 0], [0, 111]], np.vstack((s, h222)).T, axis=0)
        out113 = np.append([[N, 0], [0, 0], [0, 113]], np.vstack((s, h224)).T, axis=0)
        out133 = np.append([[N, 0], [0, 0], [0, 133]], np.vstack((s, h244)).T, axis=0)
        out333 = np.append([[N, 0], [0, 0], [0, 333]], np.vstack((s, h444)).T, axis=0)
        

        wake_vert = np.vstack(([[17, 0]], out000,  out001,  out003, out011,  out013, out024, out033,
                               out122, out223, out144, out344, out124, out234, out111, out113, out133, out333
                               ))
        
        return wake_horz, wake_vert
    
    
class WakeTableParallelPlate3(WakeTable3):
    """
    WakeTableParallelPlate() - create 3rd order wake tables for parallel plate corrugated 
    strcutres based on 1st order analytical results.
    """
    def __init__(self, b=500e-6, a=500e-6, t=250e-6, p=500e-6, length=1, sigma=30e-6, orient="horz"):
        WakeTable3.__init__(self)
        wake_horz, wake_vert = self.calculate_wake_table(b, a, t, p, length, sigma)
        if orient == "horz":
            self.TH = self.process_wake_table(wake_horz)
        else:
            self.TH = self.process_wake_table(wake_vert)
       
    def calculate_wake_table(self, b, a, t, p, length, sigma):
        """
        :param b: distance from the +y/+x plates [m], offset of beam from center = a - b
        :param a: half gap between plates [m],   0<b < 2a
        :param t: longitudinal gap of corrugations [m]
        :param p: period of corrugations [m]
        :param length: length of the strucutre [m]
        :param sigma: longitudinal beam size in [m], to calculate s window length
        :param orient: structure orientation: "horz" kick in y, or "vert" kick in x
        :return: hor_wake_table, vert_wake_table
        """
        offset = a - b
        if np.abs(offset) >= a:
            raise ValueError("distance to +y/+x plate must be smaller than 2a")
        c = speed_of_light
        Z0 = 3.767303134695850e+02
        s = np.arange(0, 50 + 0.01, 0.01) * sigma
        ns = len(s)
        
        t2p = t / p
        alpha = 1 - 0.465 * np.sqrt(t2p) - 0.07 * t2p
        s0r = a * a * t / (2 * np.pi * alpha * alpha * p * p )   # Bane s0r      
        
        # convert coefficient from cgs unit to mks unit, scale by length 
        mks = Z0 * c / (4 * np.pi) *length

        # initialize
        h000 = np.zeros(ns)
        h001 = np.zeros(ns)
        h002 = np.zeros(ns)
        h003 = np.zeros(ns)
        h004 = np.zeros(ns)
        h011 = np.zeros(ns)
        h012 = np.zeros(ns)
        h013 = np.zeros(ns)
        h014 = np.zeros(ns)
        h023 = np.zeros(ns)
        h024 = np.zeros(ns)
        h033 = np.zeros(ns)
        h034 = np.zeros(ns)
        
        # higher order
        h111 = np.zeros(ns)
        h112 = np.zeros(ns)
        h113 = np.zeros(ns)
        h114 = np.zeros(ns)
        h122 = np.zeros(ns)
        h123 = np.zeros(ns)
        h124 = np.zeros(ns)
        h133 = np.zeros(ns)
        h134 = np.zeros(ns)
        h144 = np.zeros(ns)
        h222 = np.zeros(ns)
        h223 = np.zeros(ns)
        h224 = np.zeros(ns)
        h233 = np.zeros(ns)
        h234 = np.zeros(ns)
        h244 = np.zeros(ns)
        h333 = np.zeros(ns)
        h334 = np.zeros(ns)
        h344 = np.zeros(ns)
        h444 = np.zeros(ns)
        
        # trig function
        csc = lambda x : 1 / np.sin(x)
        sec = lambda x : 1 / np.cos(x)
        cot = lambda x : 1 / np.tan(x)
        
        
        # scale factors, hab
        Y = pi * offset / (2 * a)
        
        if Y == 0:
            sl = 9/4 * s0r
            sd = (15/14)**2 * s0r
            sq = (15/16)**2 * s0r
            
            h000 = mks * np.pi**2 / (4 * a**2) *  np.exp( - np.sqrt(s/sl) )
            h011 = mks * (-1) * np.pi**4 / (64 * a**4) * np.exp( - np.sqrt(s/sq) )
            h013 = -h011
            h024 = mks * np.pi**4 / (64 * a**4) * np.exp( - np.sqrt(s/sd) )
            h033 = h011
        else:            
            sl = 4 * s0r * (1 + np.cos(Y)**2 / 3 + Y*np.tan(Y)) ** (-2)
            sm = 4 * s0r * (1.5 - Y*cot(2*Y) + Y*csc(Y)*sec(Y) ) ** (-2)
            sd = 4 * s0r * ( (64 + np.cos(2*Y))/30 + 2*Y*np.tan(Y) + (0.3 - Y*np.sin(2*Y))/(np.cos(2*Y) - 2) ) ** (-2)
            sq = 4 * s0r * ( (56 - np.cos(2*Y))/30 + 2*Y*np.tan(Y) - (0.3 + Y*np.sin(2*Y))/(np.cos(2*Y) - 2) ) ** (-2)
            
            h000 = mks * np.pi**2 / (4 * a**2) * sec(Y)**2 * np.exp( - np.sqrt(s/sl) )
            h002 = mks * np.pi**3 / (16 * a**3) * np.sin(2*Y) * sec(Y)**4 * np.exp( - np.sqrt(s/sm) )
            h004 = h002
            h011 = mks * np.pi**4 / (64 * a**4) * ( np.cos(2*Y) - 2 ) * sec(Y)**4 * np.exp( - np.sqrt(s/sq) )
            h013 = -h011
            h024 = mks * np.pi**4 / (64 * a**4) * ( 2 - np.cos(2*Y) ) * sec(Y)**4 * np.exp( - np.sqrt(s/sd) )
            h033 = h011
        
            # higher order
            st = 16 * s0r * ( 5*Y*np.tan(Y) + Y*cot(Y) + ( 2*Y*np.sin(2*Y) )/(5-np.cos(2*Y)) + 5 ) ** (-2)
            h112 = mks * np.pi**5 / (768 * a**5) * sec(Y)**5 * ( np.sin(3*Y) - 11*np.sin(Y) ) *  np.exp( - np.sqrt(s/st) )
            h114 = h112
            h233 = h112
            h334 = h112
            h123 = -h112
            h134 = -h112
            h222 = -h112
            h224 = -h112
            h244 = -h112
            h444 = -h112
        
        N = len(s)
        # full coefficients, horizontal, kick is in y 
        out000 = np.append([[N, 0], [0, 0], [0, 0]], np.vstack((s, h000)).T, axis=0)
        out002 = np.append([[N, 0], [0, 0], [0, 2]], np.vstack((s, h002)).T, axis=0)
        out004 = np.append([[N, 0], [0, 0], [0, 4]], np.vstack((s, h004)).T, axis=0)
        out011 = np.append([[N, 0], [0, 0], [0, 11]], np.vstack((s, h011)).T, axis=0)
        out013 = np.append([[N, 0], [0, 0], [0, 13]], np.vstack((s, h013)).T, axis=0)
        out024 = np.append([[N, 0], [0, 0], [0, 24]], np.vstack((s, h024)).T, axis=0)
        out033 = np.append([[N, 0], [0, 0], [0, 33]], np.vstack((s, h033)).T, axis=0)
        
        out112 = np.append([[N, 0], [0, 0], [0, 112]], np.vstack((s, h112)).T, axis=0)
        out114 = np.append([[N, 0], [0, 0], [0, 114]], np.vstack((s, h114)).T, axis=0)
        out233 = np.append([[N, 0], [0, 0], [0, 233]], np.vstack((s, h233)).T, axis=0)
        out334 = np.append([[N, 0], [0, 0], [0, 334]], np.vstack((s, h334)).T, axis=0)
        
        out123 = np.append([[N, 0], [0, 0], [0, 123]], np.vstack((s, h123)).T, axis=0)
        out134 = np.append([[N, 0], [0, 0], [0, 134]], np.vstack((s, h134)).T, axis=0)
        out222 = np.append([[N, 0], [0, 0], [0, 222]], np.vstack((s, h222)).T, axis=0)
        out224 = np.append([[N, 0], [0, 0], [0, 224]], np.vstack((s, h224)).T, axis=0)
        out244 = np.append([[N, 0], [0, 0], [0, 244]], np.vstack((s, h244)).T, axis=0)
        out444 = np.append([[N, 0], [0, 0], [0, 444]], np.vstack((s, h444)).T, axis=0)
        

        wake_horz = np.vstack(([[17, 0]], out000,  out002,  out004, out011,  out013, out024, out033,
                               out112, out114, out233, out334, out123, out134, out222, out224, out244, out444
                               ))
        # full coefficients, vertical, swap x -> y, x0 -> y0, i.e., 1 -> 2, 3 -> 4, keep 0
        # note h22 = -h11, h44 = -h33
        out000 = np.append([[N, 0], [0, 0], [0, 0]], np.vstack((s, h000)).T, axis=0)
        out001 = np.append([[N, 0], [0, 0], [0, 1]], np.vstack((s, h002)).T, axis=0)
        out003 = np.append([[N, 0], [0, 0], [0, 3]], np.vstack((s, h004)).T, axis=0)
        out011 = np.append([[N, 0], [0, 0], [0, 11]], np.vstack((s, -h011)).T, axis=0)
        out024 = np.append([[N, 0], [0, 0], [0, 24]], np.vstack((s, h013)).T, axis=0)
        out013 = np.append([[N, 0], [0, 0], [0, 13]], np.vstack((s, h024)).T, axis=0)
        out033 = np.append([[N, 0], [0, 0], [0, 33]], np.vstack((s, -h033)).T, axis=0)
        
        out122 = np.append([[N, 0], [0, 0], [0, 122]], np.vstack((s, h112)).T, axis=0)
        out223 = np.append([[N, 0], [0, 0], [0, 223]], np.vstack((s, h114)).T, axis=0)
        out144 = np.append([[N, 0], [0, 0], [0, 144]], np.vstack((s, h233)).T, axis=0)
        out344 = np.append([[N, 0], [0, 0], [0, 344]], np.vstack((s, h334)).T, axis=0)
        
        out124 = np.append([[N, 0], [0, 0], [0, 124]], np.vstack((s, h123)).T, axis=0)
        out234 = np.append([[N, 0], [0, 0], [0, 234]], np.vstack((s, h134)).T, axis=0)
        out111 = np.append([[N, 0], [0, 0], [0, 111]], np.vstack((s, h222)).T, axis=0)
        out113 = np.append([[N, 0], [0, 0], [0, 113]], np.vstack((s, h224)).T, axis=0)
        out133 = np.append([[N, 0], [0, 0], [0, 133]], np.vstack((s, h244)).T, axis=0)
        out333 = np.append([[N, 0], [0, 0], [0, 333]], np.vstack((s, h444)).T, axis=0)
        

        wake_vert = np.vstack(([[17, 0]], out000,  out001,  out003, out011,  out013, out024, out033,
                               out122, out223, out144, out344, out124, out234, out111, out113, out133, out333
                               ))
        
        return wake_horz, wake_vert
    


class Wake(PhysProc):
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
    TH - list from WakeTable, (T, H): T- table of wakes coefs, H - matrix of the coefs place in T
    """

    def __init__(self, step=1):
        PhysProc.__init__(self)
        self.w_sampling = 500  # wake sampling
        self.filter_order = 20  # smoothing filter order
        # self.wake_file = ""
        self.wake_table = None
        self.factor = 1.
        self.step = step
        self.TH = None

    def convolution(self, xu, u, xw, w):
        # convolution of equally spaced functions
        hx = xu[1] - xu[0]
        wc = np.convolve(u, w) * hx
        nw = w.shape[0]
        nu = u.shape[0]
        x0 = xu[0] + xw[0]
        xc = x0 + np.arange(nw + nu) * hx
        return xc, wc

    def wake_convolution(self, xb, bunch, xw, wake):
        # convolution of unequally spaced functions
        # bunch defines the parameters
        nb = xb.shape[0]
        xwi = xb - xb[0]
        wake1 = np.interp(xwi, xw, wake, 0, 0)
        wake1[0] = wake1[0] * 0.5
        xW, Wake = self.convolution(xb, bunch, xwi, wake1)
        return xW[0:nb], Wake[0:nb]

    def add_wake(self, I, T):
        """
        [x, W] = AddWake(I, T)
        :param I: wake table in V/C, W in V (R, L, Cinv, nm, W0, N0, W1, N1)
        :param T: wake table in V/C, W in V
        :return:
        """
        R, L, Cinv, nm, W0, N0, W1, N1 = T
        c = speed_of_light
        x = I[:, 0]
        bunch = I[:, 1]
        if L != 0 or N1 > 0:
            d1_bunch = Der(x, bunch)
        nb = x.shape[0]
        W = np.zeros(nb)
        if N0 > 0:
            x, ww = self.wake_convolution(x, bunch, W0[:, 0], W0[:, 1])
            W = W - ww[0:nb] / c
        if N1 > 0:
            x, ww = self.wake_convolution(x, d1_bunch, W1[:, 0], W1[:, 1])
            # W = W - ww[0:nb]
            W = W + ww[0:nb]
        if R != 0:
            W = W - bunch * R
        if L != 0:
            # W = W - d1_bunch*L*c
            W = W + d1_bunch * L * c
        if Cinv != 0:
            int_bunch = Int1(x, bunch)
            W = W - int_bunch * Cinv / c
        return x, W

    def add_total_wake(self, X, Y, Z, q, TH, Ns, NF):
        T, H = TH
        c = speed_of_light
        Np = X.shape[0]
        X2 = X ** 2
        Y2 = Y ** 2
        XY = X * Y
        # generalized currents;
        I00 = s2current(Z, q, Ns, NF, c)
        Nw = I00.shape[0]
        if (H[0, 2] > 0) or (H[2, 3] > 0) or (H[2, 4] > 0):
            qn = q * Y
            I01 = s2current(Z, qn, Ns, NF, c)
        if (H[0, 1] > 0) or (H[1, 3] > 0) or (H[1, 4] > 0):
            qn = q * X
            I10 = s2current(Z, qn, Ns, NF, c)
        if H[1, 2] > 0:
            qn = q * XY
            I11 = s2current(Z, qn, Ns, NF, c)
        if H[1, 1] > 0:
            qn = q * (X2 - Y2)
            I20_02 = s2current(Z, qn, Ns, NF, c)
        # longitudinal wake
        # mn=0
        x, Wz = self.add_wake(I00, T[int(H[0, 0])])
        if H[0, 1] > 0:
            x, w = self.add_wake(I10, T[int(H[0, 1])])
            Wz = Wz + w
        if H[0, 2] > 0:
            x, w = self.add_wake(I01, T[int(H[0, 2])])
            Wz = Wz + w
        if H[1, 1] > 0:
            x, w = self.add_wake(I20_02, T[int(H[1, 1])])
            Wz = Wz + w
        if H[1, 2] > 0:
            x, w = self.add_wake(I11, T[int(H[1, 2])])
            Wz = Wz + 2 * w
        Pz = np.interp(Z, x, Wz, 0, 0)
        Py = np.zeros(Np)
        Px = np.zeros(Np)
        # mn=01
        Wz[0:Nw] = 0
        Wy = np.zeros(Nw)
        if H[0, 4] > 0:
            x, w = self.add_wake(I00, T[int(H[0, 4])])
            Wz = Wz + w
            Wy = Wy + w
        if H[1, 4] > 0:
            x, w = self.add_wake(I10, T[int(H[1, 4])])
            Wz = Wz + 2 * w
            Wy = Wy + 2 * w
        if H[2, 4] > 0:
            x, w = self.add_wake(I01, T[int(H[2, 4])])
            Wz = Wz + 2 * w
            Wy = Wy + 2 * w
        Pz = Pz + np.interp(Z, x, Wz, 0, 0) * Y
        h = x[1] - x[0]
        Wy = -Int1h(h, Wy)
        Py = Py + np.interp(Z, x, Wy, 0, 0)
        # mn=10
        Wz[0:Nw] = 0
        Wx = np.zeros(Nw)
        if H[0, 3] > 0:
            x, w = self.add_wake(I00, T[int(H[0, 3])])
            Wz = Wz + w
            Wx = Wx + w
        if H[1, 3] > 0:
            x, w = self.add_wake(I10, T[int(H[1, 3])])
            Wz = Wz + 2 * w
            Wx = Wx + 2 * w
        if H[2, 3] > 0:
            x, w = self.add_wake(I01, T[int(H[2, 3])])
            Wz = Wz + 2 * w
            Wx = Wx + 2 * w
        Wx = -Int1h(h, Wx)
        Pz = Pz + np.interp(Z, x, Wz, 0, 0) * X
        Px = Px + np.interp(Z, x, Wx, 0, 0)
        # mn=11
        if H[3, 4] > 0:
            x, w = self.add_wake(I00, T[int(H[3, 4])])
            Wx = -2 * Int1h(h, w)
            p = np.interp(Z, x, Wx, 0, 0)
            Px = Px + p * Y
            Py = Py + p * X
            Pz = Pz + 2 * np.interp(Z, x, w, 0, 0) * XY
        # mn=02,20
        if H[3, 3] > 0:
            x, w = self.add_wake(I00, T[int(H[3, 3])])
            Pz = Pz + np.interp(Z, x, w, 0, 0) * (X2 - Y2)
            Wx = -2 * Int1h(h, w)
            p = np.interp(Z, x, Wx, 0, 0)
            Px = Px + p * X
            Py = Py - p * Y
        I00[:, 0] = - I00[:, 0]
        # Z=-Z
        return Px, Py, Pz, I00

    def prepare(self, lat):
        if self.wake_table is None:
            _logger.info("Wake.wake_table is None! Please specify the WakeTable()")
        else:
            self.TH = self.wake_table.TH

    def get_long_wake(self, current_profile):
        """
        method to extract a longitudinal wake from the Table for specific current profile

        :param current_profile: 2D array with shape (n, 2) where first column is position and second is a beam current
        :return: wake
        """
        T, H = self.TH
        x, Wz = self.add_wake(current_profile, T[int(H[0, 0])])
        return x, Wz * self.factor

    def apply(self, p_array, dz):
        _logger.debug(" Wake: apply fraction: dz/L = " + str(dz))

        ps = p_array.rparticles
        Px, Py, Pz, I00 = self.add_total_wake(ps[0], ps[2], ps[4], p_array.q_array, self.TH,
                                              self.w_sampling, self.filter_order)

        L = self.s_stop - self.s_start
        if L == 0:
            dz = 1.0
        else:
            dz = dz / L

        p_array.rparticles[5] = p_array.rparticles[5] + Pz * dz * self.factor / (p_array.E * 1e9)
        p_array.rparticles[3] = p_array.rparticles[3] + Py * dz * self.factor / (p_array.E * 1e9)
        p_array.rparticles[1] = p_array.rparticles[1] + Px * dz * self.factor / (p_array.E * 1e9)



class Wake3(Wake):
    """
    The wake field impact on the beam is included as series of kicks.
    In order to take into account the impact of the wake field on the beam the longitudinal wake function
    of point charge through the third order Taylor expansion is used.
    In general case it uses 13 + 20 one-dimensional functions to represent the  longitudinal component of the wake
    function for arbitrary sets of the source and the wittness particles near to the reference axis.

    parameters:
    -----------
    w_sampling = 500 -  defines the number of the equidistant sampling points for the one-dimensional
                        wake coefficients in the Taylor expansion of the 3D wake function.
    filter_order = 20 - smoothing filter order
    wake_table = None - wake table [WakeTable3()]
    factor = 1. - scaling coefficient
    TH - list from WakeTable, (T, H): T- table of wakes coefs, H - matrix of the coefs place in T
    """

    def __init__(self, step=1):
        PhysProc.__init__(self)
        self.w_sampling = 500  # wake sampling
        self.filter_order = 20  # smoothing filter order
        # self.wake_file = ""
        self.wake_table = None
        self.factor = 1.
        self.step = step
        self.TH = None

    def add_total_wake(self, X, Y, Z, q, TH, Ns, NF):
        T, H = TH
        c = speed_of_light
        Np = X.shape[0]
        X2 = X ** 2
        Y2 = Y ** 2
        XY = X * Y
        X3 = X ** 3
        Y3 = Y ** 3
        X2Y = X**2 * Y
        XY2 = X * Y**2
        # generalized currents;
        I00 = s2current(Z, q, Ns, NF, c)
        Nw = I00.shape[0]
        if (H[0, 0, 2] > 0) or (H[0, 2, 3] > 0) or (H[0, 2, 4] > 0) or (H[2, 3, 3] > 0) or (H[2, 3, 4] > 0) or (H[2, 4, 4] > 0):
            qn = q * Y
            I01 = s2current(Z, qn, Ns, NF, c)
        if (H[0, 0, 1] > 0) or (H[0, 1, 3] > 0) or (H[0, 1, 4] > 0) or (H[1, 3, 3] > 0) or (H[1, 3, 4] > 0) or (H[1, 4, 4] > 0):
            qn = q * X
            I10 = s2current(Z, qn, Ns, NF, c)
        if (H[0, 1, 2] > 0) or (H[1, 2, 3] > 0) or (H[1, 2, 4] > 0):
            qn = q * XY
            I11 = s2current(Z, qn, Ns, NF, c)
        if H[0, 1, 1] > 0:
            qn = q * (X2 - Y2)
            I20_02 = s2current(Z, qn, Ns, NF, c)
        if (H[1, 1, 3] > 0) or (H[1, 1, 4] > 0):
            qn = q * X2
            I20 = s2current(Z, qn, Ns, NF, c)
        if (H[2, 2, 3] > 0) or (H[2, 2, 4] > 0):
            qn = q * Y2
            I02 = s2current(Z, qn, Ns, NF, c)
        if H[1, 1, 1] > 0:
            qn = q * X3
            I30 = s2current(Z, qn, Ns, NF, c)
        if H[1, 1, 2] > 0:
            qn = q * X2Y
            I21 = s2current(Z, qn, Ns, NF, c)           
        if H[1, 2, 2] > 0:
            qn = q * XY2
            I12 = s2current(Z, qn, Ns, NF, c)
        if H[2, 2, 2] > 0:
            qn = q * Y3
            I03 = s2current(Z, qn, Ns, NF, c)
            
        # longitudinal wake
        # mn=0
        x, Wz = self.add_wake(I00, T[int(H[0, 0, 0])])
        if H[0, 0, 1] > 0:
            x, w = self.add_wake(I10, T[int(H[0, 0, 1])])
            Wz = Wz + w
        if H[0, 0, 2] > 0:
            x, w = self.add_wake(I01, T[int(H[0, 0, 2])])
            Wz = Wz + w
        if H[0, 1, 1] > 0:
            x, w = self.add_wake(I20_02, T[int(H[0, 1, 1])])
            Wz = Wz + w
        if H[0, 1, 2] > 0:
            x, w = self.add_wake(I11, T[int(H[0, 1, 2])])
            Wz = Wz + 2 * w
        if H[1, 1, 1] > 0:
            x, w = self.add_wake(I30, T[int(H[1, 1, 1])])
            Wz = Wz + w
        if H[1, 1, 2] > 0:
            x, w = self.add_wake(I21, T[int(H[1, 1, 2])])
            Wz = Wz + 3 * w
        if H[1, 2, 2] > 0:
            x, w = self.add_wake(I12, T[int(H[1, 2, 2])])
            Wz = Wz + 3 * w
        if H[2, 2, 2] > 0:
            x, w = self.add_wake(I03, T[int(H[2, 2, 2])])
            Wz = Wz + w     
        Pz = np.interp(Z, x, Wz, 0, 0)
        Py = np.zeros(Np)
        Px = np.zeros(Np)
        # mn=01
        Wz[0:Nw] = 0
        Wy = np.zeros(Nw)
        if H[0, 0, 4] > 0:
            x, w = self.add_wake(I00, T[int(H[0, 0, 4])])
            Wz = Wz + w
            Wy = Wy + w
        if H[0, 1, 4] > 0:
            x, w = self.add_wake(I10, T[int(H[0, 1, 4])])
            Wz = Wz + 2 * w
            Wy = Wy + 2 * w
        if H[0, 2, 4] > 0:
            x, w = self.add_wake(I01, T[int(H[0, 2, 4])])
            Wz = Wz + 2 * w
            Wy = Wy + 2 * w
        if H[1, 1, 4] > 0:
            x, w = self.add_wake(I20, T[int(H[1, 1, 4])])
            Wz = Wz + 3 * w
            Wy = Wy + 3 * w
        if H[1, 2, 4] > 0:
            x, w = self.add_wake(I11, T[int(H[1, 2, 4])])
            Wz = Wz + 6 * w
            Wy = Wy + 6 * w            
        if H[2, 2, 4] > 0:
            x, w = self.add_wake(I02, T[int(H[2, 2, 4])])
            Wz = Wz + 3 * w
            Wy = Wy + 3 * w            
        Pz = Pz + np.interp(Z, x, Wz, 0, 0) * Y
        h = x[1] - x[0]
        Wy = -Int1h(h, Wy)
        Py = Py + np.interp(Z, x, Wy, 0, 0)
        # mn=10
        Wz[0:Nw] = 0
        Wx = np.zeros(Nw)
        if H[0, 0, 3] > 0:
            x, w = self.add_wake(I00, T[int(H[0, 0, 3])])
            Wz = Wz + w
            Wx = Wx + w
        if H[0, 1, 3] > 0:
            x, w = self.add_wake(I10, T[int(H[0, 1, 3])])
            Wz = Wz + 2 * w
            Wx = Wx + 2 * w
        if H[0, 2, 3] > 0:
            x, w = self.add_wake(I01, T[int(H[0, 2, 3])])
            Wz = Wz + 2 * w
            Wx = Wx + 2 * w
        if H[1, 1, 3] > 0:
            x, w = self.add_wake(I20, T[int(H[1, 1, 3])])
            Wz = Wz + 3 * w
            Wx = Wx + 3 * w            
        if H[1, 2, 3] > 0:
            x, w = self.add_wake(I11, T[int(H[1, 2, 3])])
            Wz = Wz + 6 * w
            Wx = Wx + 6 * w                     
        if H[2, 2, 3] > 0:
            x, w = self.add_wake(I02, T[int(H[2, 2, 3])])
            Wz = Wz + 3 * w
            Wx = Wx + 3 * w     
        Wx = -Int1h(h, Wx)
        Pz = Pz + np.interp(Z, x, Wz, 0, 0) * X
        Px = Px + np.interp(Z, x, Wx, 0, 0)
        # mn=11
        Wz[0:Nw] = 0
        Wx[0:Nw] = 0
        Wy[0:Nw] = 0
        if H[0, 3, 4] > 0:
            x, w = self.add_wake(I00, T[int(H[0, 3, 4])])
            Wz = Wz + 2 * w
            Wx = Wx + 2 * w
            Wy = Wy + 2 * w
        if H[1, 3, 4] > 0:
            x, w = self.add_wake(I10, T[int(H[1, 3, 4])])
            Wz = Wz + 6 * w
            Wx = Wx + 6 * w
            Wy = Wy + 6 * w
        if H[2, 3, 4] > 0:
            x, w = self.add_wake(I01, T[int(H[2, 3, 4])])
            Wz = Wz + 6 * w
            Wx = Wx + 6 * w
            Wy = Wy + 6 * w            
        Wx = -Int1h(h, Wx)
        Wy = -Int1h(h, Wy)
        Px = Px + np.interp(Z, x, Wx, 0, 0) * Y
        Py = Py + np.interp(Z, x, Wy, 0, 0) * X
        Pz = Pz + np.interp(Z, x, Wz, 0, 0) * XY
        # mn=02,20
        if H[0, 3, 3] > 0:
            x, w = self.add_wake(I00, T[int(H[0, 3, 3])])
            Pz = Pz + np.interp(Z, x, w, 0, 0) * (X2 - Y2)
            Wx = -2 * Int1h(h, w)
            p = np.interp(Z, x, Wx, 0, 0)
            Px = Px + p * X
            Py = Py - p * Y
        # other terms for 02
        Wz[0:Nw] = 0
        Wy[0:Nw] = 0
        if H[1, 4, 4] > 0:
            x, w = self.add_wake(I10, T[int(H[1, 4, 4])])
            Wz = Wz + 3 * w
            Wy = Wy + 6 * w
        if H[2, 4, 4] > 0:
            x, w = self.add_wake(I01, T[int(H[2, 4, 4])])
            Wz = Wz + 3 * w
            Wy = Wy + 6 * w
        Wy = -Int1h(h, Wy)
        Py = Py + np.interp(Z, x, Wy, 0, 0) * Y
        Pz = Pz + np.interp(Z, x, Wz, 0, 0) * Y2
        # other terms for 20
        Wz[0:Nw] = 0
        Wx[0:Nw] = 0
        if H[1, 3, 3] > 0:
            x, w = self.add_wake(I10, T[int(H[1, 3, 3])])
            Wz = Wz + 3 * w
            Wx = Wx + 6 * w
        if H[2, 3, 3] > 0:
            x, w = self.add_wake(I01, T[int(H[2, 3, 3])])
            Wz = Wz + 3 * w
            Wx = Wx + 6 * w
        Wx = -Int1h(h, Wx)
        Px = Px + np.interp(Z, x, Wx, 0, 0) * X
        Pz = Pz + np.interp(Z, x, Wz, 0, 0) * X2        
        # mn=12
        if H[3, 4, 4] > 0:
            x, w = self.add_wake(I00, T[int(H[3, 4, 4])])
            Wz = 3 * w
            Wx = 3 * w
            Wy = 6 * w            
            Wx = -Int1h(h, Wx)
            Wy = -Int1h(h, Wy)
            Px = Px + np.interp(Z, x, Wx, 0, 0) * Y2
            Py = Py + np.interp(Z, x, Wy, 0, 0) * XY
            Pz = Pz + np.interp(Z, x, Wz, 0, 0) * XY2
        # mn=21
        if H[3, 3, 4] > 0:
            x, w = self.add_wake(I00, T[int(H[3, 3, 4])])
            Wz = 3 * w
            Wx = 6 * w
            Wy = 3 * w            
            Wx = -Int1h(h, Wx)
            Wy = -Int1h(h, Wy)
            Px = Px + np.interp(Z, x, Wx, 0, 0) * XY
            Py = Py + np.interp(Z, x, Wy, 0, 0) * X2
            Pz = Pz + np.interp(Z, x, Wz, 0, 0) * X2Y        
        # mn=30
        if H[3, 3, 3] > 0:
            x, w = self.add_wake(I00, T[int(H[3, 3, 3])])
            Wz = w
            Wx = 3 * w         
            Wx = -Int1h(h, Wx)
            Px = Px + np.interp(Z, x, Wx, 0, 0) * X2
            Pz = Pz + np.interp(Z, x, Wz, 0, 0) * X3
        # mn=03
        if H[4, 4, 4] > 0:
            x, w = self.add_wake(I00, T[int(H[4, 4, 4])])
            Wz =  w
            Wy = 3 * w            
            Wy = -Int1h(h, Wy)
            Py = Py + np.interp(Z, x, Wy, 0, 0) * Y2
            Pz = Pz + np.interp(Z, x, Wz, 0, 0) * Y3        

        I00[:, 0] = - I00[:, 0]
        # Z=-Z
        return Px, Py, Pz, I00

    def get_long_wake(self, current_profile):
        """
        method to extract a longitudinal wake from the Table for specific current profile

        :param current_profile: 2D array with shape (n, 2) where first column is position and second is a beam current
        :return: wake
        """
        T, H = self.TH
        x, Wz = self.add_wake(current_profile, T[int(H[0, 0, 0])])
        return x, Wz * self.factor



class WakeKick(Wake):
    def __init__(self, factor=1):
        print("WakeKick physics process is obsolete. Use Wake.")
        Wake.__init__(self)
        self.factor = factor

    def apply(self, p_array, dz):
        _logger.debug(" WakeKick: apply")
        ps = p_array.rparticles
        Px, Py, Pz, I00 = self.add_total_wake(ps[0], ps[2], ps[4], p_array.q_array, self.TH, self.w_sampling,
                                              self.filter_order)

        p_array.rparticles[5] = p_array.rparticles[5] + self.factor * Pz / (p_array.E * 1e9)
        p_array.rparticles[3] = p_array.rparticles[3] + self.factor * Py / (p_array.E * 1e9)
        p_array.rparticles[1] = p_array.rparticles[1] + self.factor * Px / (p_array.E * 1e9)
