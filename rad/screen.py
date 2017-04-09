__author__ = 'Sergey Tomin'


#import convolution as conv
from ctypes import c_double
from numpy import zeros, empty_like, linspace, array, sin, cos
#from ocelot.common.screen import Screen


def Py2C(array):
    arr_type =  c_double*len(array)
    c_array = arr_type(*array)
    return c_array


class Screen:
    def __init__(self):
        # position of screen center
        self.x = 0.0   # in [m]
        self.y = 0.0   # in [m]
        self.z = 100.0   # in [m]

        # Screen sizes and resolution
        self.size_x = 1.0  # half-size horizontal in [m]
        self.size_y = 1.0  # half-size vertical in [m]
        self.nx = 1  # points per horizontal direction
        self.ny = 1  # points per vertical direction

        # parameters of energy options
        self.start_energy = 100.0
        self.end_energy = 10000.0
        self.num_energy = 1000

        # half of angle aperture.  Angle relative to undulator axis for emittance influence on spectrum
        self.theta_x = 0  # rad
        self.theta_y = 0  # rad
        self.update()

    def update(self):
        self.x_start = (self.x - self.size_x)*1000.  # in mm

        if(self.nx != 1):
            self.x_step = float(self.size_x*2.)/(self.nx - 1)*1000.  # in mm
        else:
            self.x_step = 0.                # in mm
            self.x_start = self.x*1000.  # in mm
        # Y
        self.y_start = (self.y - self.size_y)*1000.  # in mm
        if self.ny != 1:
            self.y_step = float(self.size_y*2.)/(self.ny - 1)*1000.  # in mm
        else:
            self.y_step = 0                 # in mm
            self.y_start = self.y*1000.  # in mm
        # E
        self.e_start = self.start_energy
        if(self.num_energy != 1):
            self.e_step = (float(self.end_energy) - self.start_energy)/(self.num_energy - 1)
        else:
            self.e_step = 0
        self.ne = self.num_energy

        self.Zstart = 0. # this parameter is needed to calculation radiation

        self.Distance = self.z*1000. # Note this parameter is just distance from ZERO to screen # in mm!!!

        Nscr = self.ne*self.nx*self.ny

        self.memory_screen = zeros(Nscr*5)
        self.arReEx = self.memory_screen[0:Nscr]
        self.arImEx = self.memory_screen[Nscr:2*Nscr]
        self.arReEy = self.memory_screen[2*Nscr:3*Nscr]
        self.arImEy = self.memory_screen[3*Nscr:4*Nscr]
        self.arPhase = self.memory_screen[4*Nscr:5*Nscr]
        self.Xph = linspace(self.x_start, self.x_start + self.x_step*(self.nx -1), self.nx)
        self.Yph = linspace(self.y_start, self.y_start + self.y_step*(self.ny -1), self.ny)
        self.Eph = linspace(self.e_start, self.e_start + self.e_step*(self.ne -1), self.ne)


        #class EMScreen():
        #    def __init__(self, screen=None):
        #        if screen !=None:
        #            self.screen_to_emscreen(screen)
        #        else:
        #            self.create_empty_emclass()


    def screen_to_emscreen(self, screen):

        self.x = screen.x   # in [m]
        self.y = screen.y   # in [m]
        self.size_x = screen.size_x  # in [m] half-size horizontal
        self.size_y = screen.size_y  # in [m] half-size vertical

        self.start_energy = screen.start_energy # in [eV]
        self.end_energy = screen.end_energy     # in [eV]
        self.num_energy = screen.num_energy     # in [eV]
        # X
        self.x_start = (screen.x - screen.size_x)*1000.  # in mm

        if(screen.nx != 1):
            self.x_step = float(screen.size_x*2.)/(screen.nx - 1)*1000.  # in mm
        else:
            self.x_step = 0.
            self.x_start = screen.x*1000.  # in mm
        self.nx = screen.nx # points per horizontal direction
        # Y
        self.y_start =  (screen.y - screen.size_y)*1000.  # in mm
        if(screen.ny != 1):
            self.y_step = float(screen.size_y*2.)/(screen.ny - 1)*1000.  # in mm
        else:
            self.y_step = 0
            self.y_start = screen.y*1000.  # in mm
        self.ny = screen.ny # points per vertical direction
        # E
        self.e_start = screen.start_energy
        if(screen.num_energy != 1):
            self.e_step = (float(screen.end_energy) - screen.start_energy)/(screen.num_energy - 1)
        else:
            self.e_step = 0
        self.ne = screen.num_energy

        self.Zstart = 0. # this parameter is needed to calculation radiation
        # it is initial position in longitudinal direction for chain of undulators
        self.z = screen.z
        self.Distance = screen.z*1000. # Note this parameter is just distance from ZERO to screen # in mm!!!

        # half of angle aperture.  Angle relative to undulator axis for emittance influence on spectrum
        self.theta_x = screen.theta_x
        self.theta_y = screen.theta_y

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
        #self.start_energy = screen.start_energy

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

        self.Distance = 0
        # half of angle aperture.  Angle relative to undulator axis for emittance influence on spectrum
        self.theta_x = 0.
        self.theta_y = 0.

        Nscr = self.ne*self.nx*self.ny
        self.memory_screen = zeros(Nscr)

        self.arReEx = self.memory_screen[0:Nscr]
        self.arImEx = self.memory_screen[Nscr:2*Nscr]
        self.arReEy = self.memory_screen[2*Nscr:3*Nscr]
        self.arImEy = self.memory_screen[3*Nscr:4*Nscr]
        self.arPhase = self.memory_screen[4*Nscr:5*Nscr]

    def nullify(self):
        Nscr = self.ne*self.nx*self.ny
        #self.current = 1
        self.memory_screen = zeros(Nscr*5)
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
    screen = Screen()
    screen.create_like(screen_down)

    sinfa = sin(screen_down.arPhase)
    cosfa = cos(screen_down.arPhase)
    screen.arReEx = screen_down.arReEx + screen_up.arReEx*cosfa - screen_up.arImEx*sinfa
    screen.arImEx = screen_down.arImEx + screen_up.arImEx*cosfa + screen_up.arReEx*sinfa
    screen.arReEy = screen_down.arReEy + screen_up.arReEy*cosfa - screen_up.arImEy*sinfa
    screen.arImEy = screen_down.arImEy + screen_up.arImEy*cosfa + screen_up.arReEy*sinfa

    screen.arPhase = screen_down.arPhase + screen_up.arPhase
    return screen

