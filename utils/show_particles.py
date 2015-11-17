'''

'''
from ocelot.adaptors.genesis import *
try:
    import matplotlib.animation as anim
except:
    print 'animation not installed'
import numpy as np
import matplotlib.pyplot as plt

def show_plots(displays, fig):
    n1 = (len(displays) -1 )/2 + 1
    n2 = (len(displays) -1) / n1 +1
    #print n1, n2
    fmt = str(n1)+ str(n2)
    print fmt
    
    for i in xrange(len(displays)):
        ax = fig.add_subplot(fmt + str(i+1))
        ax.grid(True)
        for f in displays[i].data:
            #x,y = f(x = np.linspace(-10, 10, 100))
            ax.plot(f[0], f[1], '.')

    plt.show()

class Display:
    def __init__(self, data=lambda x: (x, 0*x) , xlabel='',ylabel=''):
        self.data = (data,)
        self.xlabel = xlabel
        self.ylabel = ylabel


sys.path.append('../utils/')


if len(sys.argv)>1:
    outf = sys.argv[1]
else:
    outf = '/home/iagapov/tmp/run_25/run.25.gout'

g = readGenesisOutput(outf)


npoints = g('ncar')
zstop = g('zstop')
delz = g('delz')
xlamd = g('xlamd')
xlamds = g('xlamds')
#nslice = int(g('nslice'))
nslice = len(g.sliceValues.keys())
zsep = g('zsep')
npart = int(g('npart'))

I = np.array(g.I)

smax = nslice * zsep * xlamds 
                 
print 'zstop=', zstop
print 'delz=', delz
print 'xlamd=', xlamd
print 'xlamds=', xlamds
print 'npoints', npoints
print 'nslice', nslice
print 'zsep', zsep
print 'npart', npart


#particles = readParticleFile(fileName=outf + '.dpa', npart=npart, nslice = nslice)

particles = read_particle_file(file_name=outf + '.dpa', npart=npart, nslice = nslice)

'''
#plt.plot(particles[90][2],particles[90][3], '.')
f = plt.figure()
#f = plt.figure('x/px')
ax=f.add_subplot(221)
ax.set_xlabel('x')
ax.set_ylabel('px')
#n, bins = np.histogram(particles[1][3], 50, normed=False)
#plt.bar( bins[1:] - (bins[1]-bins[0])/2.0 , n, width = (bins[1]-bins[0]), alpha=0.5)
plt.plot(particles[islice][2],particles[islice][3],'.')
plt.grid(True)

ax=f.add_subplot(222)
ax.set_xlabel('y')
ax.set_ylabel('py')
#f = plt.figure('y/py')
plt.plot(particles[islice][4],particles[islice][5],'.')
plt.grid(True)

ax=f.add_subplot(223)
ax.set_xlabel('E')
ax.set_ylabel('phi')
#f = plt.figure('E/pz')
plt.plot(particles[islice][0],particles[islice][1],'.')
plt.grid(True)

ax=f.add_subplot(224)
ax.set_xlabel('E')
ax.set_ylabel('phi')
#f = plt.figure('E/pz')
plt.plot(particles[islice2][0],particles[islice][1],'.')
plt.grid(True)
'''


E = np.zeros(nslice)
phi = np.zeros(nslice)
phi_std = np.zeros(nslice)

for i in xrange(nslice):
    E[i] = np.mean( np.array(particles[i][0]) )
    phi[i] = np.mean( np.array(particles[i][1]) )
    phi_std[i] = np.std( np.array(particles[i][1]) )

f = plt.figure()

ax=f.add_subplot(311)
ax.set_xlabel('E')
ax.plot(E)
ax=f.add_subplot(312)
ax.set_xlabel('phi')
ax.plot(phi)
ax=f.add_subplot(313)
ax.set_xlabel('phi_std')
ax.plot(phi_std)


fig = plt.figure()

displays = []

i_slices = [100,200,300]

if len(sys.argv)>2: i_slices = map(int, sys.argv[2:])

for i in i_slices:
    d = Display()
    #d.data = (lambda x: (particles[i][0],particles[i][1]),)
    d.data = ((particles[i][1],particles[i][0]),)  # phase horizontal axis
    displays.append(d)


show_plots(displays, fig)

plt.show()

