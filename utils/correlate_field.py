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
ncar = g('ncar')

smax = nslice * zsep * xlamds 
                 
print ('zstop=', zstop)
print ('delz=', delz)
print ('xlamd=', xlamd)
print ('xlamds=', xlamds)
print ('npoints', npoints)
print ('nslice', nslice)
print ('zsep', zsep)
print ('ncar', ncar)


radSlices = readRadiationFile(fileName=file+'.dfl', npoints=npoints)
print ('read', len(radSlices), 'slices')
E = np.copy(radSlices[0])

for i in np.arange(1,len(radSlices)):
    #E1 += radSlices[i][0]
    E += np.abs(radSlices[i]) **2

Z = np.abs(E)**2

i1 = 1
i2 = 251


i1 = int(ncar/2)
i2 = int(ncar/2 + ncar/6)

fig = plt.figure()
x = radSlices[:,i1,i1]
y = radSlices[:,i2,i2]

ax = fig.add_subplot(221)
plot(np.abs(x))

ax = fig.add_subplot(222)
plot(np.abs(y))

ax = fig.add_subplot(223)
plot(np.angle(x))

ax = fig.add_subplot(224)
plot(np.angle(y))


plt.show()