"""
writen by I. Zagorodnov, DESY and S.Tomin XFEL, 2015.
"""
from scipy.integrate import simps
from numpy import arange, sqrt, append, zeros, conj, dot, linspace
from numpy.fft import fft, irfft, ifft
from ocelot.common.globals import *

def wake2impedance(s, w):
    """
    %Fourier transform with exp(iwt)
    %             s    - Meter
    %             w -    V/C
    %             f -    Hz
    %             y -    Om
    """
    ds = s[1] - s[0]
    dt = ds/speed_of_light
    n = len(s)
    shift = 1.
    y = dt*fft(w, n)*shift
    return y


def impedance2wake(f, y):
    """
    % Fourier transform with exp(-iwt)
    %             f -    Hz
    %             y -    Om
    %             s    - Meter
    %             w -    V/C
    """
    df = f[1] - f[0]
    n = len(f)
    s = 1./df*arange(n)/n*speed_of_light
    w = n*df*irfft(y, n)
    return s, w


def imp_resistiveAC_SI(f, cond, a, tau, L):
    """
    # resistive impedance of round pipe (in SI Units)
    # dimensions: f - Hertz
    #             cond - in 1/Second
    #             a - pipe radius in m
    #             tau - the relaxation time in s
    #             L-inductive for dielectric layer
    """
    n = len(f)
    f2w = 2.*pi
    koef = a*0.5*1j/(speed_of_light*Z0)
    Z = zeros(n, dtype=complex)
    for i in range(0, n):
        w = f[i]*f2w
        kw = cond/(1. + 1j*w*tau)
        Zs = sqrt(1j*w*mu_0/kw) + 1j*w*L
        Z[i] = Zs/(f2w*a*(1. + w*Zs*koef))
    return Z


def ResistiveZaZb(xb, bunch, a, conductivity, tau, Ind):
    nb = len(xb)
    ds = xb[1] - xb[0]
    n = 2*nb

    dt = ds/speed_of_light
    f = 1./dt*arange(n)/n
    Za = 1e-12*imp_resistiveAC_SI(f[:nb], conductivity, a, tau, Ind) # -> v/pC/m

    xb1 = linspace(xb[0], xb[0]+ds*(n-1), num=n)

    bunch1 = append(bunch, zeros(nb))
    Zb = wake2impedance(xb1, bunch1*speed_of_light)

    Z = zeros(n, dtype=complex)

    Z[:nb] = Za[:]*Zb[:nb]
    #print len(Z[nb-1::-1]), len(Z[nb:])
    Z[nb:] = Z[nb-1::-1]
    #plt.plot(Z.real)
    #plt.show()
    xa, wa = impedance2wake(f, Z)
    res = zeros(nb, dtype=complex)
    res[:] = -wa[:nb]

    return res


def LossShape(bunch, wake):
    """
    % loss, spread, peak
    % dimensions:
    %             wake - m , Volt/pC
    %             out - V/pC;
    """
    w = wake[1]
    bi2 = bunch[1]
    s = wake[0]
    loss = simps(-bi2*w, s)
    spread = sqrt(simps(bi2*(w + loss)**2, s))
    peak = max(abs(w))
    return loss, spread, peak


def pipe_wake(z, current, tube_radius, tube_len, conductivity, tau, roughness, d_oxid):

    Q = simps(current, z)/speed_of_light
    print ("Charge = ", Q*1e12, "pC")
    xb = -z[::-1]
    yb = current[::-1]/(speed_of_light*Q)
    Q = Q*1e12 #C->pC

    ds=xb[3]-xb[0]
    xb = append(xb, arange(1, 100001)*ds + xb[-1])
    yb = append(yb, arange(1, 100001)*0)

    # roughness and axid layer are here
    eps_r = 2.
    Ind = mu_0*((eps_r-1.)/eps_r*d_oxid + 0.01035*roughness)

    # the result is in V
    W = ResistiveZaZb(xb, yb, tube_radius, conductivity, tau, Ind)#*Q*L

    W = W.real*Q*tube_len
    n = len(current)
    bunch = [xb[:n], yb[:n]]
    wake = [xb[:n], W[:n]]
    # postprocessing
    L, S, P = LossShape(bunch, wake)
    print ('Loss [V]:  ', L)
    print ('Spread [V]:', S)
    print ('Peak [V]:  ', P)
    return bunch, wake

def xfel_pipe_wake(s, current):
    """
    :param s: smaller number (usually negative) is the head
    :param current:
    :return: s, current, wake
    """
    conductivity=3.66e+7   # Ohm aluminium Dohlus
    tau=7.1e-15            # s relaxation time - aluminium Dohlus
    tube_radius = 5e-3     # m radius
    tube_len = 1.          # m length
    roughness = 600e-9
    d_oxid = 5e-9          # m thickness of oxide layer

    bunch, wake = pipe_wake(s, current, tube_radius, tube_len, conductivity, tau, roughness, d_oxid)
    return bunch[0], bunch[1], wake[1]


if __name__ == "__main__":
    from numpy import loadtxt
    import matplotlib.pyplot as plt
    tube_radius = 5e-3 #m radius
    tube_len = 1.   #m length

    #conductivity=1.9e+6;    t=2.4e-15; %kicker
    conductivity=3.66e+7  #aluminium Dohlus
    tau=7.1e-15            #relaxation time - aluminium Dohlus
    #conductivity=5.8e+7;    t=2.46e-14; %copper Dohlus
    #conductivity=2.78e+7;    t=2.46e-14; %BeCu
    #conductivity=1.4e+6;    t=2.4e-15; %stainless steel 304
    #conductivity=0.6e+6; t=0; %titanium
    roughness=600e-9
    d_oxid=5e-9

    # bunch shape
    # sigma=25e-6; %m
    # q=1000; %in pC!!!
    # dx=0.01;xb(:,1)=[-5:dx:100]*sigma; yb(:,1)=gauss(xb,sigma);

    # (can be arbitrary)
    current=loadtxt('current.txt')
    s = current[:,0]
    current = current[:,1]
    bunch, wake = pipe_wake(s, current, tube_radius, tube_len, conductivity, tau, roughness, d_oxid)
    nk = max(abs(wake[1]))/max(abs(bunch[1]))
    plt.plot(bunch[0], bunch[1]*nk, wake[0], wake[1])
    plt.grid(True)
    plt.show()

