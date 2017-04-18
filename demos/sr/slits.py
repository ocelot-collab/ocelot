'''
Bending magnet radiation
'''

from ocelot.rad.bmrad import *


def slit_func(xc,yc):
    for i in range(len(xc)):
        for j in range(len(yc)):
            is_slit = False
            if np.abs(xc[i] - 0.05) < 0.01 and np.abs(yc[j] - 0.00) < 0.01 : is_slit = True 
            if np.abs(xc[i] + 0.05) < 0.01 and np.abs(yc[j] - 0.00) < 0.01 : is_slit = True
            if np.abs(xc[i] + 0.00) < 0.01 and np.abs(yc[j] - 0.05) < 0.01 : is_slit = True
            if np.abs(xc[i] + 0.00) < 0.01 and np.abs(yc[j] + 0.05) < 0.01 : is_slit = True
            if is_slit:
                E[i,j] = bm_e0(xc[i], yc[j], 7.0, p_en = 3.0, R=170.0)
            else: 
                E[i,j] = 0.0
            
    
    return E

def interference(E,R):
    # interference pattern
    pass
    '''
    ys = np.linspace(-0.5,0.5,500)
    Ei = np.zeros([len(ys)], dtype='complex')
    d = 1.0
    p_en = 2.0
    for i in range(len(ys)):
        r1 = sqrt((ys[i] - x0)**2 + d**2)
        r2 = sqrt((ys[i] + x0)**2 + d**2)
        Ei[i] = E[0]*exp(2j*pi*r1*p_en / hc)/r1 + E2[0] * exp(2j*pi*r2*p_en/hc)/r2

    print np.abs(Ei * conjugate(Ei))

    plot(ys, np.abs(Ei * conjugate(Ei)))
    '''


# vertical cut

sigx = 1.e-4
sigxp = 1.e-5

xc = np.linspace(-5.*sigx,5.*sigx,20)
xpc = np.linspace(-5.*sigxp,5.*sigxp,20)


#yo = np.linspace(-0.2,0.2,50)
yo = np.linspace(-1.e-1,1.e-1,50)
E = np.zeros([len(yo)], dtype='complex')


#E2 = np.zeros([len(yo)], dtype='complex')

x0 = 0.05

for i in range(len(yo)):
    E[i] = bm_e_a(xc,0,xpc, 0,  x0, yo[i], 10.0, p_en = 2.0, R=170.0)
    #E2[i] = bm_e_a(xc,0,xpc, 0,  -x0, yo[i], 10.0, p_en = 2.0)


ax1 = figure().add_subplot(111)
p1, = ax1.plot(yo, np.abs(E))

sigx = 1.e-5
sigxp = 1.e-5

xc = np.linspace(-5.*sigx,5.*sigx,250)
xpc = np.linspace(-5.*sigxp,5.*sigxp,250)

for i in range(len(yo)):
    E[i] = bm_e_a(xc,0,xpc, 0,  0.05, yo[i], 10.0, p_en = 2.0, R=170.0)


p2, = ax1.plot(yo, np.abs(E))

ax1.set_xlabel('Y [m]')
ax1.set_ylabel('arb. units')

ax1.legend([p1,p2],['1.e-4 x 1.e-5','1.e-5 x 1.e-5'])

plt.show()
sys.exit(0)


# slits


sigx = 1.e-4
sigxp = 1.e-5


yc = np.linspace(-0.1,0.1,200)
xc = np.linspace(-0.1,0.1,200)

X, Y = np.meshgrid(xc,yc)
E = np.zeros([len(xc),len(yc)], dtype='complex')
Z = np.zeros([len(xc),len(yc)])

E = slit_func(xc,yc)

Z = np.abs(E)

plt.imshow(Z.T)
plt.figure()
plt.plot(yc, Z[len(xc)/2, :])

plt.show()


