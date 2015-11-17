'''
extract genesis 'beam' file from genesis 'particle' file
usage python particle2beam.py in.gout out.beam  
'''
from xframework.adaptors.genesis import *
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
    beamf = sys.argv[2]
else:
    outf = '/home/iagapov/tmp/run_2/run.2.gout'
    beamf = '/home/iagapov/tmp/run_2/tmp.beam'

g = readGenesisOutput(outf)
b = read_beam_file(beamf)
        
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


print 'number of slices', len(b.z)


#particles = readParticleFile(fileName=outf + '.dpa', npart=npart, nslice = nslice)

particles = read_particle_file(file_name=outf + '.dpa', npart=npart, nslice = nslice)


beam = Beam()

beam.g0 = np.zeros(nslice)   # energy
beam.dg = np.zeros(nslice)
beam.eloss = b.eloss
beam.phi = np.zeros(nslice)
beam.phi_std = np.zeros(nslice)

beam.z = b.z
beam.I = I

beam.betax = np.zeros(nslice)
beam.alphax = np.zeros(nslice)
beam.betay = np.zeros(nslice)
beam.alphay = np.zeros(nslice)
beam.ex = np.zeros(nslice)
beam.ey = np.zeros(nslice)

beam.x = np.zeros(nslice)
beam.y = np.zeros(nslice)

beam.px = np.zeros(nslice)
beam.py = np.zeros(nslice)


beam.columns = ['ZPOS','CURPEAK','GAMMA0', 'DELGAM', 'BETAX', 'BETAY', 'ALPHAX', 
                'ALPHAY', 'EMITX','EMITY','XBEAM','YBEAM','PXBEAM','PYBEAM','ELOSS']
beam.column_values  = {}

beam.idx_max = 0


for i in xrange(nslice):
    
    print 'processing slice ', i, '/', nslice
    
    beam.g0[i] = np.mean( np.array(particles[i][0]) )
    beam.dg[i] = np.std( np.array(particles[i][0]) )
    beam.phi[i] = np.mean( np.array(particles[i][1]) )
    beam.phi_std[i] = np.std( np.array(particles[i][1]) )
    
    beam.x[i] = np.mean( np.array(particles[i][2]) )
    beam.px[i] = np.mean( np.array(particles[i][3]) )

    beam.y[i] = np.mean( np.array(particles[i][4]) )
    beam.py[i] = np.mean( np.array(particles[i][5]) )
    
    #beam.z[i] = i
    
    gx = np.matrix([[0,0],[0,0]])
    gy = np.matrix([[0,0],[0,0]])
    
    q = np.array(particles[i][2])
    p = np.array(particles[i][3])
    
    gx = np.matrix([[np.mean(q*q), np.mean(q*p)],
                    [np.mean(q*p), np.mean(p*p)]])

    q = np.array(particles[i][4])
    p = np.array(particles[i][5])
    
    gy = np.matrix([[np.mean(q*q), np.mean(q*p)],
                    [np.mean(q*p), np.mean(p*p)]])
    
            
    beam.ex[i] =  np.sqrt( np.linalg.det(gx) )
    beam.ey[i] =  np.sqrt( np.linalg.det(gy) )
    
    beam.betax[i] = gx[0,0] / ( beam.ex[i] / beam.g0[i] )
    beam.alphax[i] = gx[0,1] / ( beam.ex[i] / beam.g0[i] )

    beam.betay[i] = gy[0,0] / ( beam.ey[i] / beam.g0[i] )
    beam.alphay[i] = gy[0,1] / ( beam.ey[i] / beam.g0[i] ) 

    
        
print beam_file_str(beam)

plot_beam(plt.figure(), beam)    
plt.show()


f=open('new.beam','w')
f.write(beam_file_str(beam))





