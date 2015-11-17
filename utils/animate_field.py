from ocelot.adaptors.genesis import *
import matplotlib.animation as anim
import numpy as np
import matplotlib.pyplot as plt

#file='/home/iagapov/data/fel/genesis_runs/flash_40fsec/2900A/run_1/run.1.gout'
file='/home/iagapov/tmp/workshop/run_1/run.1.gout'
g = readGenesisOutput(file)

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


radSlices = readRadiationFile(fileName=file+'.dfl', npoints=npoints)
print 'read', len(radSlices), 'slices'
E = np.copy(radSlices[0])

for i in np.arange(1,len(radSlices)):
    #E1 += radSlices[i][0]
    E += np.abs(radSlices[i]) **2

Z = np.abs(E)**2

i1 = 1
i2 = 251


def get_z(i = 0):
    global radSlices
    #E1 = radSlices[i][0]
    #E2 = radSlices[i][1]
    #Z = E1*E1 + E2*E2
    return np.abs(radSlices[i])**2


for isl in [320,321,322,323,324,325,326,327,328,329,330]:#[300,320,340,360,380,400]:
    fig = plt.figure((str(isl)))
    z = get_z(isl)[i1:i2,i1:i2]
    m = plt.imshow(z, cmap='gist_gray')

i = 0

def updatefig(*args):
    global i
    m.set_data( get_z(i))
    i = (i+1) % len(radSlices)
    return m,

#ani = anim.FuncAnimation(fig, updatefig, interval = 1)
#ani.save('sase_pulse.mp4', fps=25, clear_temp=True)
plt.show()