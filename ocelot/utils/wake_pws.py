"""
Planar wakefield structure

S.Tomin, 11.2018
"""
import numpy as np
from ocelot.cpbd.beam import *


def wake_kick(p_array, b=500*1e-6, t=0.25*1e-3, p=0.5*1e-3):
    """
    Function to calculate transverse kick by corrugated structure [SLAC-PUB-16881]

    :param p_array: ParticleArray
    :param b: distance to the dechirper wall
    :param t: gap
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

    kick, I = wake_kick(p_array, b=500 * 1e-6, t=0.25 * 1e-3, p=0.5 * 1e-3)

    fig, ax1 = plt.subplots()
    color="C0"
    ax1.plot(I[:, 0], I[:, 1], color=color)
    ax1.set_ylabel("I, A", color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()
    color = "C1"
    ax2.plot(kick[:, 0], kick[:, 1]*1e-6, color=color)
    ax2.set_ylabel("W, MV", color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    plt.show()