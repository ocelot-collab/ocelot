from ocelot.optics.bragg import *

f_name = 'C444'


lines=open('data/'+ f_name + '.txt').read().split('\n')

f000 = {}
fh = {}
fhbar = {}

for l in lines:
    tks = l.split()
    if len(tks) < 2: continue
    if tks[1] == 'energy':
        energy = float(tks[3])
        #print energy
    if tks[2] == 'F(0,0,0)':
        val = map(float, tks[4][1:-1].split(',') )
        f000[energy] = val[0] + 1j*val[1]    
    if tks[2] == 'FH':
        val = map(float, tks[4][1:-1].split(',') )
        fh[energy] = val[0] + 1j*val[1]    
    if tks[2] == 'FH_BAR':
        val = map(float, tks[4][1:-1].split(',') )
        fhbar[energy] = val[0] + 1j*val[1]    
    
    
ev = sorted(f000.keys())


cdata = CrystalStructureFactors()

cdata.ev = ev
cdata.f000 = map(lambda x: f000[x], ev)
cdata.fh = map(lambda x: fh[x], ev)
cdata.fhbar = map(lambda x: fhbar[x], ev)

print(cdata.ev)
print(cdata.fh)

save_stucture_factors(cdata, 'data/'+f_name + '.dat')
bdata = load_stucture_factors(f_name + '.dat')

plt.plot(bdata.ev, np.imag(bdata.f000), 'b.--')
plt.show()
