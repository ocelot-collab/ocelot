'''
basic undulator trajectory and radiation calculations
'''

import scipy.special as sf
#from matplotlib.figure import Figure
#from mpl_toolkits.mplot3d import Axes3D

from ocelot.cpbd.beam import *
from ocelot.cpbd.elements import *

import numpy as np
from numpy import arange

me = 1.e+6 * 0.510998928  # electron mass, eV / c^2
c = 299792458.0
ce = -1.602176565e-19     # electron charge, C
h = 4.135667516e-15       # Planck's constant, eV*s


def total(intens):
    nf = np.shape(intens)[2]
    itot = np.zeros(nf)
    for i in xrange(0,nf):
        itot[i] = np.sum(intens[:,:,i])
    return itot

class UndulatorParameters:

    def __init__(self, und = None): 

        if und == None:
            self.Nw = 100   # number of periods
            self.lw = 0.045 # period length, m
            self.E = 1.0 # electron energy, in GeV
            self.K = 1.0    # undulator parameter
     
        else:
            
            self.Nw = und.nperiods   # number of periods
            self.lw = und.lperiod # period length, m
            self.E = 1.0 # electron energy, in GeV
            self.K = und.Kx    # undulator parameter

        self.recalculate()


    def recalculate(self):
                
        self.gamma = 1.0e+9 * self.E / me
        self.L = self.Nw*self.lw
        self.B = 2.0*pi * self.K * me / (self.lw * c)
        self.kw = 2.0*pi / self.lw
        self.beta = sqrt(1.0 - 1.0 / (self.gamma*self.gamma) )
        self.beta_z = self.beta * (1.0 - (self.K*self.K)/ (4.0*self.gamma*self.gamma))
        #self.w1, _, _ = self.computeRadiationAnalytical(200.0, 0.0, 0)
        self.w1 = 1.0 / ( (1.0 + self.K*self.K/2.0 ) / (2*self.kw*c*self.gamma**2) )
        self.lamd_ev = h*self.w1/(2*np.pi)
        
        
    def computeRadiationAnalytical(self,z0, theta, phi, dwRel = 0.1, npoints=100):
        w1 = 1.0 / ( (1.0 + self.K*self.K/2.0 + self.gamma*self.gamma*theta) / (2*self.kw*c*self.gamma**2) ) # first harmonic
        dw = w1 * dwRel
        step_w = dw / npoints
        w = arange(w1-dw, w1+dw, step_w) 
        phis = theta * theta * w * z0 / (2.0 * c)
        delta_w = w - w1
    
        u = w * self.K*self.K / (8*self.kw*c*self.gamma*self.gamma)
        ajj = sf.jv(0,u) - sf.jv(1,u)
        I = self.Nw * self.lw * self.K * w *  np.exp(1j * phis) / ( z0 * self.gamma)  * np.sinc(pi*self.Nw*(w-w1)/w1) * ajj
        I = np.real(I) 
        
        self.fundamentalWavelength = self.lw / (2*self.gamma**2) * (1+ self.K**2 / 2.0 + self.gamma**2 * theta)
        
        print 'test', h*c / self.fundamentalWavelength 
        
        return w1, w, I
    
    def get_k(self, E):
        # w1 in ev
        w1 = 2. * np.pi * E / h
        return np.sqrt( (4.*self.kw*c*self.gamma**2)/w1 - 2. )  

        
    def printParameters(self):
    
        self.recalculate()

        print "Undulator parameters:"
        print "L=", self.L
        print "gamma(electron)=", self.gamma
        print "K=", self.K
        print "B[T]=", self.B

        w1, _, _ = self.computeRadiationAnalytical(200.0, 0.0, 0)

        print "Radiation parameters:"
        print "w1(first harmonic, zero angle)=", w1, "Hz/2pi", h*w1/(2*np.pi), "[eV]", self.fundamentalWavelength, "m"

        t = self.L / c
        cb = 379.35 # 1/(gev sec Tesla^2)

        eloss = self.E - self.E / (1.0 + 0.5*self.B*self.B*cb*self.E*t)

        print "Total energy loss [Gev]", eloss

        
parameters = UndulatorParameters()

def computeTrajectory(par, z0,z1,npoints):
        
    beta = par.beta * (1 - par.K**2 / (4 * par.gamma**2))
        
    tr = Trajectory()    
        
    tr.ct = np.arange(0, z1-z0, (z1-z0)/npoints)
    tr.z = np.arange(0, z1-z0, (z1-z0)/npoints)
    tr.x = par.K  / (par.gamma* par.kw ) * cos(par.kw*beta*tr.ct)
    tr.xp = -c * beta* par.K / par.gamma * sin (par.kw*beta*tr.ct) # xp == velocity
    
    tr.y = 0*tr.ct
    tr.yp =0*tr.ct
    tr.s = tr.z * par.beta / par.beta_z - par.K**2/(8.0*par.gamma**2*par.kw) * sin(2.0*par.kw*tr.z)  

    return tr

def computeTrajectoryNum(par, z0,z1,npoints):
        
    tr = Trajectory() 
        
    beta = par.beta * (1.0 - par.K**2 / (4 * par.gamma**2))
        
    tr.ct = np.arange(0, z1-z0, (z1-z0)/npoints)
    tr.z = np.zeros(len(tr.ct))
    tr.x = np.zeros(len(tr.ct))
    tr.xp = np.zeros(len(tr.ct))
    tr.E = np.zeros(len(tr.ct))
    
    tr.xp[0] = 0.0
    tr.x[0] = par.K  / (par.gamma* par.kw )
    
    dct = (tr.ct[1] - tr.ct[0]) 
    dt = (tr.ct[1] - tr.ct[0]) / c
    
    vz = c * par.beta
        
    tr.E[0] = me / np.sqrt(1 - par.beta**2)
    
    for i in range(1,len(tr.ct)):
        tr.z[i] = tr.ct[i] 
        tr.xp[i] = tr.xp[i-1] - dct * beta* par.K * c * beta * par.kw / par.gamma * np.cos(par.kw * 0.5 * (tr.ct[i-1] + tr.ct[i])*beta)
        tr.x[i] = tr.x[i-1] + tr.xp[i]*dt
        
        cb = 379.35
        dE = 1.e-9 * (tr.E[i-1] - tr.E[i-1] / (1.0 + 0.5*par.B * par.B*cb*tr.E[i-1]*dt))
        
        d2E = np.sqrt(dE) * np.random.randn() * np.sqrt(dt)
        #print dE, d2E
        
        tr.E[i] = tr.E[i-1] - dE + d2E    
        
        par.beta = np.sqrt( 1.0 - me**2 / tr.E[i]**2 )
        par.gamma = tr.E[i] / me
        beta = par.beta * (1.0 - par.K**2 / (4 * par.gamma**2))
        
        #print dt,tr.E[i-1], 0.5*par.B * par.B*cb*tr.E[i-1]*dt, tr.E[i-1] - tr.E[i-1] / (1.0 + 0.5*par.B * par.B*cb*tr.E[i-1]*dt), dE, d2E, par.gamma
        
    tr.y = 0*tr.ct
    tr.yp =0*tr.ct
    tr.s = tr.z * par.beta / par.beta_z - par.K**2/(8.0*par.gamma**2*par.kw) * sin(2.0*par.kw*tr.z)  

    return tr

'''
def computeRadiation1(undParam, h5File = 'my_sampler.h5'):

    from common.xio import *

    io = XIO(h5File)
    
    undParam.computeRadiationAnalytical(100.0, 0.0, 0.0)
    
    intens = srwutil.Intensity()
    
    theta = arange(0.0,1.e-5,1.e-7)

    Itheta = {}
    Wtheta = {}
    Itotal = {}

    for i in range(0,len(theta)):
        _, w, I = computeRadiationAnalytical(200.0, theta[i], 0,  0.05)
        Itheta[theta[i]] = I
        Wtheta[theta[i]] = w
        Itotal[theta[i]] = sum(I)

        if 'Intensities' in self.file.keys():
            gi = self.file['Intensities']
            intens.xmin = gi.attrs['xmin']
            intens.xmax = gi.attrs['xmax']
            intens.ymin = gi.attrs['ymin']
            intens.ymax = gi.attrs['ymax']
            intens.nx = gi.attrs['nx']
            intens.ny = gi.attrs['ny']

            for e in gi.keys():
                
                #print gi[e].attrs['energy']
                #print gi[e].value
                
                intens.intensity[gi[e].attrs['energy']] = gi[e].value
                
                    
    io.writeIntensity(intens)
    
    io.file.close()
'''

class IntegrationLog:
    def __init__(self):
        self.phis = []
        self.z = []
        self.s = []
        self.ct = []
        

def computeRadiationNumerical3(z0, theta, phi, dwRel = 1, npoints=400):
    
    xo = z0*theta*cos(phi)
    yo = z0*theta*sin(phi)
    
    ct,x,vx,y,vy,z,s = computeTrajectory(-L/2.0,L/2.0,Nw * 100)
    
    w1 = 1.0 / ( (1 + K*K/2.0 + gamma*gamma*theta) / (2*kw*c*gamma) ) # first harmonic
    dw = w1 * dwRel 
    step_w = dw / npoints
    w = arange(w1-dw, w1+dw, step_w) 
    I = np.zeros(len(w))

    ss = 0
    
    log = IntegrationLog()
    log.z = z
    log.ct = ct

    for i in range(0,len(ct)-1):
        dz = z[i+1] - z[i]
        
        ds = (x[i+1] - x[i]) / (ct[i+1] - ct[i])
                
        zobs = z0 - 0.5*(z[i+1] + z[i])
        #phi = w * ( 0.5* (s[i+1] + s[i])/(beta*c) - 0.5* (z[i+1] + z[i])/ c + ( pow( xo - 0.5*(x[i] + x[i+1]), 2 ) + pow(yo - 0.5*(y[i] + y[i+1]), 2 ) )/(2*c*zobs))
        #I += -1j * w * dz * np.exp(1j*phi)/ zobs * (vx[i] / c - (xo - x[i])/zobs)

        phi = (w/c) * ( s[i]/beta - z[i] + ( x[i]*x[i] ) / (2*zobs) )
        I += -1j * w * dz * np.exp(1j*phi)/ zobs * (vx[i] / c + x[i]/zobs)

        ss += sqrt(ds*ds+1) * (ct[i+1] - ct[i])
        
        log.s.append(ss)
        log.phis.append(phi[0])
    
        
    I = np.real(I) 

    return w1, w, I, log


def computeRadiationNumerical2(par, z0, theta, phi, dwRel = 1, npoints=100):
    
    xo = z0*theta*cos(phi)
    yo = z0*theta*sin(phi)
    
    ct,x,vx,y,vy,z,s = computeTrajectory(par, -par.L/2.0,par.L/2.0,par.Nw * 10)
    
    w1 = 1.0 / ( (1 + par.K**2/2.0 + par.gamma**2*theta) / (2*par.kw*c*par.gamma**2) ) # first harmonic
    
    dw = w1 * dwRel 
    step_w = dw / npoints
    w = arange(w1-dw, w1+dw, step_w) 
    I = np.zeros(len(w))

    ss = 0
    
    log = IntegrationLog()
    log.z = z
    log.ct = ct

    for i in range(0,len(ct)-1):
        dz = z[i+1] - z[i]
        ds = (x[i+1] - x[i]) / (ct[i+1] - ct[i])
        zobs = z0 - 0.5*(z[i+1] + z[i])

        phi = (w/w1) * ( par.kw * z[i] )
        I += -1j * w / z0 *  dz * par.K / par.gamma * sin(par.kw * z[i]) * np.exp(1j*phi)

        log.phis.append(phi[0])
    
        
    I = np.real(I) 

    return w1, w, I, log


def computeRadiationNumerical(par, z0, theta, phi, dwRel = 1, npoints=100):
    
    xo = z0*theta*cos(phi)
    yo = z0*theta*sin(phi)
    
    ct,x,vx,y,vy,z,s = computeTrajectory(par, -par.L/2.0,par.L/2.0,par.Nw * 10)
    
    w1 = 1.0 / ( (1 + par.K**2/2.0 + par.gamma**2*theta) / (2*par.kw*c*par.gamma**2) ) # first harmonic
    
    dw = w1 * dwRel 
    step_w = dw / npoints
    w = arange(w1-dw, w1+dw, step_w) 
    I = np.zeros(len(w))

    ss = 0
    
    log = IntegrationLog()
    log.z = z
    log.ct = ct

    for i in range(0,len(ct)-1):
        dz = z[i+1] - z[i]
        ds = (x[i+1] - x[i]) / (ct[i+1] - ct[i])
        zobs = z0 - 0.5*(z[i+1] + z[i])

        phi = (w/w1) * ( par.kw * z[i] )
        I += -1j * w / z0 *  dz * par.K / par.gamma * sin(par.kw * z[i]) * np.exp(1j*phi)

        log.phis.append(phi[0])
    
        
    I = np.real(I) 

    return w1, w, I, log


def debugIntegration():
    theta = 0.0
    phi = 0.0
    z0 = 100
    xo = z0*theta*cos(phi)
    yo = z0*theta*sin(phi)
    
    ct,x,vx,y,vy,z,s = computeTrajectory(-L/2.0,L/2.0,Nw * 10)
    
    w = 1.0 / ( (1 + K*K/2.0 + gamma*gamma*theta) / (2*kw*c*gamma) ) # first harmonic

    ss = 0
    
    log = IntegrationLog()
    log.z = z
    log.ct = ct

    I = 0

    for i in range(0,len(ct)-1):
                      
        dz = z[i+1] - z[i]        
        zobs = z0 - 0.5*(z[i+1] + z[i])
        phi = ( kw * z[i] )
        I += -1j * w / z0 *  dz * K / gamma * sin(kw * z[i]) * np.exp(1j*phi)

        log.phis.append(phi)
    
    
    log.phis.append(kw * z[ len(z)-1])
        
    I = np.real(I) 


    plt.figure()
    line, = plot(z,log.phis,'-')
    line, = plot(z,x*1.e8,'-')
    
    plt.draw()
        
    

def plotSpectrum(par, method='analytical'):
    
    plt.figure()
    
    theta = arange(0.0,1.e-5,1.e-7)

    Itheta = {}
    Wtheta = {}
    Itotal = {}

    for i in range(0,len(theta)):
        if method == 'analytical':
            _, w, I = par.computeRadiationAnalytical(200.0, theta[i], 0,  0.05)
        else:
            pass
            _, w, I, log = computeRadiationNumerical(par, 200.0, theta[i], 0,  0.05)
        
        Itheta[theta[i]] = I
        Wtheta[theta[i]] = w
        Itotal[theta[i]] = sum(I)

    thetas = Itheta.keys()
    thetas.sort()

    points = [thetas[0]]
    
    lines = []

    for i in range(0,len(points)):
        lines.append(plot(h/(2*np.pi) * Wtheta[points[i]],Itheta[points[i]],lw='2'))
    
    legend(lines, map(lambda x: r'$\theta =$' +str(x), points) )
    plt.xlabel('[eV]')
    plt.ylabel('[A.U]')
   
    draw()
            
    
def plotIntensity():

    plt.figure()

    theta = arange(0.0,1.e-3,1.e-5)

    Itheta = {}
    Wtheta = {}
    Itotal = {}

    Itheta2 = {}
    Wtheta2 = {}
    Itotal2 = {}


    for i in range(0,len(theta)):
        _, w, I = computeRadiationAnalytical(100.0, theta[i], 0.0)
        _, w2, I2 = computeRadiationAnalytical(110.0, theta[i], 0.0)
        Itheta[theta[i]] = I
        Wtheta[theta[i]] = w
        Itotal[theta[i]] = sum(I)
        Itheta2[theta[i]] = I2
        Wtheta2[theta[i]] = w2
        Itotal2[theta[i]] = sum(I2)


    thetas = np.array(sorted(Itotal.keys()))
    line, = plot(thetas, map(Itotal.get, thetas), lw='2')

    thetas2 = np.array(sorted(Itotal2.keys()))
    line2, = plot(thetas2, map(Itotal2.get, thetas2), lw='2')
    
    legend([line, line2], ['Total intensity (analytical)', 'Total intensity (numerical)'] )

    draw()
    
def plotIntensity3D():
    #theta = arange(0.0,1.e-3,1.e-5)
    #phi = arange(0,2*pi,1.e-03)
    
    #X = np.outer(theta,sin(phi))
    #Y = np.outer(theta,cos(phi))
    
    xs = np.arange(-1.0, 1.0, 0.1)
    ys = np.arange(-1.0, 1.0, 0.1)
    Z = np.zeros([len(xs),len(ys)]) 
    
    for i in range(0,len(xs)):
        for j in range(0,len(ys)):
            Z[i,j] = sin(xs[i])
    
    X, Y = np.meshgrid(xs, ys)
       
    fig = plt.figure()
    axes = Axes3D(fig)
    surf = axes.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False)

    fig.colorbar(surf, shrink=0.5, aspect=5)

    draw()
