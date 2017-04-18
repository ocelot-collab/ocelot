'''
example07 -- crystal -- dynamical diffraction (ocelot.optics.bragg)
'''

from ocelot.optics.utils import *

def save_filter(filt, f_name):
    f= open(f_name,'w')
    for i in xrange( len(filt.ev)):
        f.write(str(filt.ev[i]) + '\t' + str(np.abs(filt.tr[i])**2) + '\t' +  str(np.abs(filt.ref[i])**2) + '\n')

E_ev = 10000
ref_idx = (2,2,0)
thickness = 5 * mum


cr1 = Crystal(r=[0,0,0*cm], size=[5*cm,5*cm,thickness], no=[0,0,-1], id="cr1")
cr1.lattice =  CrystalLattice('Si')
#cr1.lattice =  CrystalLattice('Si')
#cr1.psi_n = -(pi/2. - 54.7356*(pi/180.0)) #input angle psi_n according to Authier 
cr1.psi_n = -pi/2. #input angle psi_n according to Authier (symmetric reflection, Si)


r = Ray(r0=[0,0.0,-0.5], k=[0,0.0,1]) 
r.lamb = 2 * pi * hbar * c / E_ev
print('wavelength', r.lamb)

w1 = read_signal(file_name='data/pulse_9kev_20fs.txt', npad =10, E_ref = E_ev)
plt.figure()
plot_signal(w1)

#plt.figure()
f_test = get_crystal_filter(cryst=cr1, ray=r, nk=3000, ref_idx = ref_idx)
filt = get_crystal_filter(cr1, r, ref_idx = ref_idx, k = w1.freq_k)

#save_filter(f_test, 'C400_8000ev_filter.txt')

plot_filters(filt, f_test)
plot_filters(filt, f_test, param='ref')

fig=plt.figure()

plot_spec_filt(w1, filt, ax=fig.add_subplot(111))

cr1.filter = filt

i1=np.sum(w1.sp*np.conj(w1.sp))*(w1.freq_ev[1] - w1.freq_ev[0])

def transform_field(cr, wave):
    print('transforming field')
    wave.sp = wave.sp * cr.filter.tr
    wave.sp_ref = wave.sp * cr.filter.ref
    wave.f = np.fft.ifft(wave.sp)

fig = plt.figure()
plt.grid(True)
ax = fig.add_subplot(111)
plt.plot(w1.t, np.abs(w1.f))
transform_field(cr1, w1)

i2 = np.sum(w1.sp*np.conj(w1.sp))*(w1.freq_ev[1] - w1.freq_ev[0])
i3 = np.sum(w1.sp_ref*np.conj(w1.sp_ref))*(w1.freq_ev[1] - w1.freq_ev[0])

print('transmission (%)', 100*np.real(i2/i1), 'reflection (%)', 100*np.real(i3/i1))

plt.plot(w1.t, np.abs(w1.f))
ax.set_yscale('log')    
    
plt.figure(), plt.grid(True)
plt.plot(w1.freq_ev, np.abs(w1.sp))


plt.show()
