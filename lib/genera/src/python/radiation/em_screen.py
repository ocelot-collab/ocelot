__author__ = 'Sergey Tomin'


#import convolution as conv
from ctypes import c_double
from numpy import zeros, empty_like, linspace, array, sin, cos
from xframework.common.screen import Screen

def Py2C(array):
    arr_type =  c_double*len(array)
    c_array = arr_type(*array)
    return c_array
        
class EMScreen():
    def __init__(self, screen=None):
        if screen !=None:
            self.screen_to_emscreen(screen)
        else:
            self.create_empty_emclass()


    def screen_to_emscreen(self, screen):
        # X
        self.x_start = (screen.x - screen.size_x)*1000.

        if(screen.nx != 1):
            self.x_step = float(screen.size_x*2.)/(screen.nx - 1)*1000.
        else:
            self.x_step = 0
            self.x_start = screen.x*1000.
        self.nx = screen.nx
        # Y
        self.y_start =  (screen.y - screen.size_y)*1000.
        if(screen.ny != 1):
            self.y_step = float(screen.size_y*2.)/(screen.ny - 1)*1000.
        else:
            self.y_step = 0
            self.y_start = screen.y*1000.
        self.ny = screen.ny
        # E
        self.e_start = screen.start_energy
        if(screen.num_energy != 1):
            self.e_step = (float(screen.end_energy) - screen.start_energy)/(screen.num_energy - 1)
        else:
            self.e_step = 0
        self.ne = screen.num_energy

        self.Zstart = 0. # this parameter is needed to calculation radiation
        # it is initial position in longitudinal direction for chain of undulators

        self.Distance = screen.z*1000. # Please note this parameter is just distance from ZERO to screen!!!

        Nscr = self.ne*self.nx*self.ny
        #self.current = 1
        self.memory_screen = zeros(Nscr*5)
        self.arReEx = self.memory_screen[0:Nscr]
        self.arImEx = self.memory_screen[Nscr:2*Nscr]
        self.arReEy = self.memory_screen[2*Nscr:3*Nscr]
        self.arImEy = self.memory_screen[3*Nscr:4*Nscr]
        self.arPhase = self.memory_screen[4*Nscr:5*Nscr]
        self.Xph = linspace(self.x_start, self.x_start + self.x_step*(self.nx -1), self.nx)
        self.Yph = linspace(self.y_start, self.y_start + self.y_step*(self.ny -1), self.ny)
        self.Eph = linspace(self.e_start, self.e_start + self.e_step*(self.ne -1), self.ne)

        # additional
        #self.fund_harm_eV = screen.fund_harm_eV
        self.start_energy = screen.start_energy

    def create_empty_emclass(self):
        self.x_step = 0
        self.x_start = 0
        self.nx = 0

        self.y_step = 0
        self.y_start = 0
        self.ny = 0

        self.e_step = 0
        self.e_start = 0
        self.ne = 0

        self.Zstart = 0

        self.Distance =0

        Nscr = self.ne*self.nx*self.ny
        self.memory_screen = zeros(Nscr)

        self.arReEx = self.memory_screen[0:Nscr]
        self.arImEx = self.memory_screen[Nscr:2*Nscr]
        self.arReEy = self.memory_screen[2*Nscr:3*Nscr]
        self.arImEy = self.memory_screen[3*Nscr:4*Nscr]
        self.arPhase = self.memory_screen[4*Nscr:5*Nscr]

    def create_like(self, em_screen):
        self.x_step = em_screen.x_step
        self.x_start = em_screen.x_start
        self.nx = em_screen.nx

        self.y_step = em_screen.y_step
        self.y_start = em_screen.y_start
        self.ny = em_screen.ny

        self.e_step = em_screen.e_step
        self.e_start = em_screen.e_start
        self.ne = em_screen.ne

        self.Zstart = em_screen.Zstart
        #self.First = 0
        self.Distance = em_screen.Distance
        #self.current = em_screen.current
        Nscr = self.ne*self.nx*self.ny
        self.memory_screen = empty_like(em_screen.memory_screen)

        self.arReEx = self.memory_screen[0:Nscr]
        self.arImEx = self.memory_screen[Nscr:2*Nscr]
        self.arReEy = self.memory_screen[2*Nscr:3*Nscr]
        self.arImEy = self.memory_screen[3*Nscr:4*Nscr]
        self.arPhase = self.memory_screen[4*Nscr:5*Nscr]

    def screenPy2C(self, lperiod, nperiods, status):

        scrPrm = zeros(14)
        
        scrPrm[0] = self.ne
        scrPrm[1] = self.e_start
        scrPrm[2] = self.e_step

        scrPrm[3] = self.nx
        scrPrm[4] = self.x_start
        scrPrm[5] = self.x_step

        scrPrm[6] = self.ny
        scrPrm[7] = self.y_start
        scrPrm[8] = self.y_step
        scrPrm[9] = self.Distance 
        scrPrm[10] = self.Zstart
        scrPrm[11] = lperiod*1000.
        scrPrm[12] = nperiods
        scrPrm[13] = status  #status is just for fast mode calculation (ideal undulator and zero initial conditions)
        c_scrPrm = Py2C(scrPrm)
        return c_scrPrm

    def screenC2Py(self, c_screen):
        Nscr = self.nx*self.ny*self.ne

        self.arReEx = c_screen[0:Nscr]
        self.arImEx = c_screen[Nscr:2*Nscr]
        self.arReEy = c_screen[2*Nscr:3*Nscr]
        self.arImEy = c_screen[3*Nscr:4*Nscr]
        self.arPhase = c_screen[4*Nscr:5*Nscr]


    def distPhoton(self, gamma, current):
        LenPntrConst = self.Distance - self.Zstart
        constQuant = 3.461090202456155e+9*current*gamma*gamma/LenPntrConst/LenPntrConst
        Ex2r = array(self.arReEx)*array(self.arReEx)
        Ex2i = array(self.arImEx)*array(self.arImEx)
        self.Sigma = (Ex2r + Ex2i)*constQuant
        Ey2r = array(self.arReEy)*array(self.arReEy)
        Ey2i = array(self.arImEy)*array(self.arImEy)
        self.Pi = (Ey2r + Ey2i)*constQuant
        self.Total = self.Pi + self.Sigma
        self.Xph = linspace(self.x_start, self.x_start + self.x_step*(self.nx -1), self.nx)
        self.Yph = linspace(self.y_start, self.y_start + self.y_step*(self.ny -1), self.ny)
        self.Eph = linspace(self.e_start, self.e_start + self.e_step*(self.ne -1), self.ne)
    """
    def convolution(self, sigma_x, sigma_y):
        data = numpy.reshape(self.Total,(len(self.Xph), len(self.Xph)))

        conv.convolution(self.Xph , self.Yph, data, sigma_x, sigma_y)
    """
        
    def zerosArray(self):
        Nscr = self.nx*self.ny*self.ne
        self.memory_screen = zeros(Nscr*5)
        self.arReEx = self.memory_screen[0:Nscr]
        self.arImEx = self.memory_screen[Nscr:2*Nscr]
        self.arReEy = self.memory_screen[2*Nscr:3*Nscr]
        self.arImEy = self.memory_screen[3*Nscr:4*Nscr]
        self.arPhase = self.memory_screen[4*Nscr:5*Nscr]


def sum_screens(screen_down, screen_up):
    """
    the function accumulates the radiations from different emitters which placed in consecutive order (u_down, d, u_up, ...).
    This means that before summation of screens we must rotate the electric field vectors of a screen_up on angle screen_down.arPhase anticlockwise,
    because the radiation calculation starts from zero phase from each emitter.
    """
    #screen = Screen()
    em_screen = EMScreen()
    em_screen.create_like(screen_down)

    sinfa = sin(screen_down.arPhase)
    cosfa = cos(screen_down.arPhase)
    em_screen.arReEx = screen_down.arReEx + screen_up.arReEx*cosfa - screen_up.arImEx*sinfa
    em_screen.arImEx = screen_down.arImEx + screen_up.arImEx*cosfa + screen_up.arReEx*sinfa
    em_screen.arReEy = screen_down.arReEy + screen_up.arReEy*cosfa - screen_up.arImEy*sinfa
    em_screen.arImEy = screen_down.arImEy + screen_up.arImEy*cosfa + screen_up.arReEy*sinfa

    em_screen.arPhase = screen_down.arPhase + screen_up.arPhase
    return em_screen

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import cm

def show_flux(screen, show = 'Total', xlim = (0,0), ylim = (0,0),  file_name = None):
    if show == 'Total':
        data = screen.Total
    elif show == 'Sigma':
        data = screen.Sigma
    else:
        data = screen.Pi

    if screen.nx == 1 or screen.ny == 1:
        if screen.nx == 1 and screen.ny == 1:
            X = screen.Eph
            xlabel = 'Eph'
        elif screen.nx == 1:
            X = screen.Yph
            xlabel = 'Yph'
        else:
            X = screen.Xph
            xlabel = 'Xph'

        D1(data, X, xlabel = xlabel, xlim = xlim, ylim = ylim,  file_name = file_name)
    else:
        if screen.ne!=1:
            print " ******** ERROR into show.screen ! *********** "
            return
        D3(screen, data, file_name = file_name)


def D1(data, X, xlabel = "Eph", xlim = (0,0), ylim = (0,0),  file_name = None):

    maxS = max(data)
    index = np.where(data== max(data))[0][0]
    energy = X[index]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(X, data)

    if xlim != (0,0):
        ax.set_xlim(xlim)
    if ylim != (0,0):
        ax.set_ylim(ylim)
    #ax.set_title()
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Intensity")
    ax.grid(True)

    ax.annotate('$\epsilon_1 = ' + str(int(energy*10)/10.) +'$', xy=(-30, -20),
               xycoords='axes points',
               horizontalalignment='right', verticalalignment='top',
               fontsize=20)

    ax.annotate('I = ' + str(int(maxS*1e-10)*1.e+10) , xy=(-25, -50),
               xycoords='axes points',
               horizontalalignment='right', verticalalignment='top',
               fontsize=15)

    if file_name != None:
        figg = plt.gcf()
        k_size = 1.4
        figg.set_size_inches( (4*k_size, 3.01*k_size) )
        figg.savefig(file_name)
    else:
        plt.show()
    plt.show()


def D3(screen,Data, file_name = None):
    #print " showme.any = ", np.shape(Data)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    X,Y = np.meshgrid(screen.Xph, screen.Yph)
    #print " showme.any = ", np.shape(X)
    #print " showme.any = ", np.shape(Y)
    data = np.zeros((screen.yNstep, screen.xNstep))
    for j in range(screen.yNstep):
        for i in range(screen.xNstep):
            data[j,i] = Data[screen.xNstep*j + i]
    ax.plot_surface(X, Y, data, rstride=1, cstride=1, cmap=cm.jet)
    #ax.set_zlim3d(0, 1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Intensity')
    ax.set_xticks([])
    if file_name != None:
        figg = plt.gcf()
        k_size = 1.7
        figg.set_size_inches( (4*k_size, 3.01*k_size) )
        figg.savefig(file_name)
    else:
        plt.show()
    #plt.show()

def plot3D_data(data, x = None, y = None):
    if x != None and y != None:
        X,Y = np.meshgrid(x,y)
    else:
        X,Y = np.meshgrid(np.arange(np.shape(data)[1]), np.arange(np.shape(data)[0]))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, data, rstride=1, cstride=1, cmap=cm.jet)
    plt.show()