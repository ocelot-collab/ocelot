__author__ = 'tomins'
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation
from matplotlib import animation
import numpy as np
from ocelot.rad import *
from ocelot import *
from ocelot.gui import *

font = {'size'   : 10}
matplotlib.rc('font', **font)

# Set up formatting for the movie files
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

beam = Beam()
beam.E = 17.5
beam.I = 0.1

#beam.beta_x = 12.84
#beam.beta_y = 6.11
#beam.Dx = 0.526
def phase_shifter(i):
    und = Undulator(Kx = 4., nperiods=125, lperiod=0.04, eid= "und")
    D = Drift(l=0.5, eid="D")
    b1 = Hcor(l=0.1, angle = 0.1*i*-0.00001, eid="b1")
    b2 = Hcor(l=0.2, angle = 0.1*i*0.00002, eid="b2")
    b3 = Hcor(l=0.1, angle = 0.1*i*-0.00001, eid="b3")
    phase_shift =  (b1, b2, b3)
    cell = (und, D, phase_shift, D, und)
    lat = MagneticLattice(cell)

    screen = Screen()
    screen.z = 100.0
    screen.size_x = 0.0
    screen.size_y = 0.0
    screen.nx = 1
    screen.ny = 1

    screen.start_energy = 7900 #eV
    screen.end_energy = 8200 #eV
    screen.num_energy = 1000

    #print_rad_props(beam, K=und.Kx, lu=und.lperiod, L=und.l, distance=screen.z)
    screen = calculate_radiation(lat, screen, beam)

    # trajectory
    Z = np.array([])
    X = np.array([])
    for u in screen.motion:
        Z = np.append(Z, u[4::9])
        X = np.append(X, u[0::9])
        #print(X)
        #plt.plot(u[4::9], u[0::9], "r")
    #X = np.array(X)
    #X.flatten()
    #Z = np.array(Z)
    #Z.flatten()
    #plt.show()
    #print(X)
    return Z, X, screen.Eph, screen.Total
    #plt.plot(screen.Eph, screen.Total)
    #plt.show()

#show_flux(screen, unit="mrad")



def init_animation():

    global line_traj, line_spectrum

    line_traj, = ax0.plot(z, x)
    ax0.set_ylim(-6, 0.2)

    line_spectrum, = ax1.plot(eph, total)
    ax1.set_ylim(0,5e15)


def animate(i):
    print(i)
    z, x, eph, total = phase_shifter(i)
    line_traj.set_ydata( x*1e6)
    line_spectrum.set_ydata(total)
    return line_traj,line_spectrum

#fig = plt.figure()
#ax = fig.add_subplot(111)
fig, (ax0, ax1) = plt.subplots(nrows=2)
ax0.grid(True)
ax0.set_ylabel(r"X, $\mu m$")
ax0.set_xlabel(r"Z, $m$")
ax1.set_ylabel(r"$I$, $\frac{ph}{sec \cdot mm^2 10^{-3}BW}$")
ax1.set_xlabel(r"$E_{ph}$, eV")
#x = np.linspace(0, 2*np.pi, 200)
z, x, eph, total = phase_shifter(0)
#print(x, z)
#plt.plot(z, x)
#plt.show()

ani = matplotlib.animation.FuncAnimation(fig, animate, init_func=init_animation, frames=60)
# worked for windows
#ani.save('animation.gif',writer='imagemagick', fps=30)
ani.save('animation.html', fps=30)
#ani.save('animation.mp4', writer=writer)