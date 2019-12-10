"""
Planar wakefield structure
S.Tomin, 11.2018

added function which generates wake tables
S.Tomin and I.Zagorodnov, 11.2019
"""
import numpy as np
from ocelot.cpbd.beam import *
from ocelot.cpbd.physics_proc import *
from ocelot.gui import *
import os

def wake_tables(b=500*1e-6, a=0.01, width=0.02, t=0.25*1e-3, p=0.5*1e-3, length=1, sigma=30e-6, filename=None):
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
    p = p * 1e3 # m -> mm
    t = t * 1e3 # m -> mm
    L = length  # in m
    D = width * 1e3 # m -> mm
    a = a*1e3  # m -> mm, in mm half gap
    b = b * 1e3 # m -> mm, distance from the plate
    sigma = sigma*1e3  # m -> mm, characteristic beam size in mm
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
    f = []
    for i in np.arange(Nm):
        m = i + 1
        # f[i] = np.pi/D*m
        M = np.pi / D * m
        X = M * a
        dx = np.sin(M * x0) * np.sin(M * x)
        Wcc = A * X / (np.cosh(X) * np.sinh(X)) * np.exp(-(s / s0) ** 0.5 * (X / np.tanh(X)))
        Wss = A * X / (np.cosh(X) * np.sinh(X)) * np.exp(-(s / s0) ** 0.5 * (X * np.tanh(X)))

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
        # Fyq[i, :] = M ** 2 * Fz[i,:]
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

    if filename is not None:
        part, ext = os.path.splitext(filename)
        np.savetxt(part + "_horz_" + str(length)+"m" + ext, wake_horz, fmt='%1.7e')
        np.savetxt(part + "_vert_" + str(length)+"m" + ext, wake_vert, fmt='%1.7e')
    return wake_horz, wake_vert

def wake_kick(p_array, b=500*1e-6, t=0.25*1e-3, p=0.5*1e-3):
    """
    Function to calculate transverse kick by corrugated structure [SLAC-PUB-16881]

    :param p_array: ParticleArray
    :param b: distance to the dechirper wall
    :param t: longitudinal gap
    :param p: period
    :return: (wake, current) - wake[:, 0] - s in [m],  wake[:, 1] - kick in [V]
                            - current[:, 0] - s in [m], current[:, 1] - current in [A]
    """
    def dipole_wake(s, b, t, p):
        """
        dipole wake
        :param s: position along a bunch
        :param b: distance to the dechirper wall
        :param t: gap
        :param p: period
        :return:
        """
        alpha = 1 - 0.465*np.sqrt(t/p) - 0.07*(t/p)
        #s0l = 2*b*b*t/(np.pi*alpha**2*p**2)

        s0yd = 8*b**2*t/(9*np.pi * alpha**2 * p**2)
        #print(s0yd, alpha, s0l)
        w = 2/b**3 * s0yd * (1 - (1 + np.sqrt(s/s0yd))*np.exp(- np.sqrt(s/s0yd)))
        return w


    I = s_to_cur(p_array.tau(), sigma=0.03*np.std(p_array.tau()), q0=np.sum(p_array.q_array), v=speed_of_light)
    I[:, 0] -= I[0, 0]
    s = I[:, 0]


    step = (s[-1] - s[0])/(len(s) - 1)
    #s = np.append(s, s+s[-1]+step)
    q = I[:, 1]/speed_of_light
    #I = np.append(B[:, 1], np.zeros(len(s) - len(B[:, 1])))/speed_of_light
    #print("charge = ", np.sum(q)*step)

    w = np.array([dipole_wake(si, b, t, p) for si in s])*377*speed_of_light/(4*np.pi)
    wake = np.convolve(q, w)*step
    s_new = np.cumsum(np.ones(len(wake)))*step
    kick_array = np.vstack((s_new, wake))
    return kick_array.T, I


if __name__ == "__main__":
    from matplotlib import pyplot as plt
    from ocelot import *

    data_dir = "/Users/tomins/ownCloud/DESY/repository/ocelot/projects/dechirper/dataxxl/ocelot"
    p_array = load_particle_array(data_dir + "/particles/section_SASE1.npz")
    print(p_array)
    kick, I = wake_kick(p_array, b=100 * 1e-6, t=25 * 1e-3, p=0.5 * 1e-3)



    fig, ax1 = plt.subplots()
    color = "C0"
    ax1.plot(I[:, 0], I[:, 1], color=color)
    ax1.set_ylabel("I, A", color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()
    color = "C1"
    ax2.plot(kick[:, 0], kick[:, 1]*1e-6, color=color)
    ax2.set_ylabel("W, MV", color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    plt.show()
    p_array = load_particle_array(data_dir + "/particles/section_SASE1.npz")
    ws = SPDKick(b=100 * 1e-6, gap=25 * 1e-3, period=0.5 * 1e-3)
    ws.alpha = 0
    ws.apply(p_array, 1)
    #show_e_beam(p_array)
    show_density(p_array.tau(), p_array.py())
    plt.show()

