'''
tuning fel bandwidth for flash
'''
from ocelot.adaptors.genesis import *
import matplotlib.animation as anim
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import *

g = readGenesisOutput('/home/iagapov/tmp/workshop/run_1/run.1.gout')

npoints = g('ncar')
zstop = g('zstop')
delz = g('delz')
xlamd = g('xlamd')
xlamds = g('xlamds')
nslice = g('nslice')
zsep = g('zsep')


smax = nslice * zsep * xlamds 
                 
print 'zstop=', zstop
print 'delz=', delz
print 'xlamd=', xlamd
print 'xlamds=', xlamds
print 'npoints', npoints
print 'nslice', nslice
print 'zsep', zsep

nZ = len(g.sliceValues[1]['power']) 


def get_u(iz = 0):
    global g
    power = np.zeros( len(g.sliceValues.keys()))
    phi = np.zeros( len(g.sliceValues.keys()))
    power_mid = np.zeros( len(g.sliceValues.keys()))
    E = np.zeros( len(g.sliceValues.keys()), dtype=complex)
        
    for i in g.sliceValues.keys():
        power[i-1] = g.sliceValues[i]['power'][iz]
        power_mid[i-1] = g.sliceValues[i]['p_mid'][iz]
        phi[i-1] = g.sliceValues[i]['phi_mid'][iz]
        #E[i-1] = power_mid[i-1] * np.exp(1j*phi[i-1])
        E[i-1] = power[i-1] * np.exp(1j*phi[i-1])


    pad = np.zeros(len(E))
    E = np.concatenate((pad,E))
    E = np.concatenate((E,pad))
    spec = fft(E)
    spec = np.sqrt( spec * np.conj(spec) )
    spec = np.real( np.roll(spec, len(spec)/2) )

    return power, power_mid, spec


fig = plt.figure()

p, p_mid, s = get_u(1)

ax = fig.add_subplot(211)
ax.grid(True)
#m, = plt.plot(p, lw = 3)
I = np.array(g.I)
I = I * np.max(p) / np.max(I)
m, = plt.plot(p, 'b', lw = 3)
m1, = plt.plot(I, 'g--', lw = 3)


ax.set_ylim(0, 100.0)

ax2 = fig.add_subplot(212)
ax2.grid(True)
m2, = plt.plot(s, lw = 3)
ax2.set_ylim(0, 100.0)


def plot_brightness():
    
    ss = []
    pp = []
    
    for i in arange(1, nZ):
        p,p_mid, s = get_u(i)
        ss.append(np.max(s))
        pp.append(np.mean(p))
    
    fig = plt.figure()
    ax3 = fig.add_subplot(111)
    print len(ss), len(arange(1, nZ))
    m3, = plt.plot(delz*xlamd*2*arange(1, nZ), ss, lw = 3)
    ax4 = ax3.twinx()
    m4, = plt.plot(delz*xlamd*2*arange(1, nZ), pp, 'r--', lw = 3)
    ax3.grid(True)
    ax4.grid(True)
    ax3.legend([m3,m4],['max spectral density','power'])

plt.grid(True)
i = 0

def updatefig(*args):
    global i, I
    
    i = i+1
    if i >= nZ-1:
        i = i - 1

    p,p_mid, s = get_u(i)
    
    m.set_data(np.arange(len(g.sliceValues.keys())), p)
    
    I = I * np.max(p) / np.max(I)
    m1.set_data(np.arange(len(g.sliceValues.keys())), I)
    
    m2.set_data(np.arange(len(g.sliceValues.keys())), s)
    
    
    print 2*i*delz*xlamd, np.max(s)
    
    ymin, ymax = ax.get_ylim()
    
    if np.max(p) >= ymax:
        ax.set_ylim(ymin, 2*ymax)
        ax.figure.canvas.draw()

    ymin, ymax = ax2.get_ylim()
    
    if np.max(s) >= ymax:
        ax2.set_ylim(ymin, 2*ymax)
        ax2.figure.canvas.draw()


    
    return m,
# Set up formatting for the movie files
#Writer = anim.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

#ani = anim.FuncAnimation(fig, updatefig, interval = 1)
#ani.save('im.mp4', fps=7)

plot_brightness()

fig = plt.figure()
ax = fig.add_subplot(221)
p,p_mid, s = get_u(int(nZ/3))
plt.plot(s, 'r--', lw=3)
ax = fig.add_subplot(222)
plt.plot(p, 'r',lw=3)


ax = fig.add_subplot(223)
p,p_mid, s = get_u(int(nZ-1))
plt.plot(s, 'g--', lw=3)
ax = fig.add_subplot(224)

fig = plt.figure()
plt.grid(True)
plt.plot(p, 'g', lw=3)


plt.show()