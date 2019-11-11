__author__ = 'Sergey Tomin'
"""
Screen class for SR module. The first version was written in 2011 - 2012. 
S.Tomin
"""

from ctypes import c_double
import numpy as np
from ocelot.common.globals import *


def Py2C(array):
    arr_type =  c_double*len(array)
    c_array = arr_type(*array)
    return c_array


class Screen:
    """
    Class to store radiation field and to provide information about screen parameters where radiation will be observed.
    Format for electric fields in arrays (arReEx, arImEx, ...) is following: ReEx[ny*nx*je + nx*jy + jx]

    self.z: 100.0 [m], distance from the beginning of the lattice to the screen
    self.size_x: 1 [m], half of screen size in horizontal plane
    self.size_y: 1 [m], half of screen size in vertical
    self.nx: 1, number of points in horizontal plane
    self.ny: 1, number of points in vertical plane
    self.start_energy: 100.0 [eV], starting photon energy
    self.end_energy: 10000.0 [eV], ending photon energy
    self.num_energy: 1000,   number of energy points
    self.arReEx = [],  Real part of horizontal component of the electric field
    self.arImEx = [],  Imaginary part of horizontal component of the electric field
    self.arReEy = [],  Real part of the vertical component of the electric field
    self.arImEy = [],  Imaginary part of the vertical component of the electric field
    self.arPhase = [], phase between Re and Im components
    """
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

        self.arReEx = []     # array, Real part of horizontal component of the electric field
        self.arImEx = []     # array, Imaginary part of horizontal component of the electric field
        self.arReEy = []     # array, Real part of the vertical component of the electric field
        self.arImEy = []     # array, Imaginary part of the vertical component of the electric field
        self.arPhase = []    # array, phase between Re and Im components
        self.Xph = []
        self.Yph = []
        self.Eph = []

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

        self.memory_screen = np.zeros(Nscr*5)
        self.arReEx = self.memory_screen[0:Nscr]
        self.arImEx = self.memory_screen[Nscr:2*Nscr]
        self.arReEy = self.memory_screen[2*Nscr:3*Nscr]
        self.arImEy = self.memory_screen[3*Nscr:4*Nscr]
        self.arPhase = self.memory_screen[4*Nscr:5*Nscr]
        self.Xph = np.linspace(self.x_start, self.x_start + self.x_step*(self.nx - 1), self.nx)
        self.Yph = np.linspace(self.y_start, self.y_start + self.y_step*(self.ny - 1), self.ny)
        self.Eph = np.linspace(self.e_start, self.e_start + self.e_step*(self.ne - 1), self.ne)

    def rebuild_efields(self, x0=0, y0=0, z0=0):
        """
        the method recalculates the field phase and electrical fields to obtain the correct values that can be used
        to propagate the wave front.

        :param x0: initial the electron coordinate
        :param y0: initial the electron coordinate
        :param z0: initial the electron coordinate
        :return:
        """
        hc = 1.239841874330e-3  # h_eV_s*speed_of_light*1000  // mm
        shape_array = [self.ne, self.ny, self.nx]

        Xscr = np.linspace(self.x_start, self.x_start + self.x_step * (self.nx - 1), num=self.nx)
        Yscr = np.linspace(self.y_start, self.y_start + self.y_step * (self.ny - 1), num=self.ny)
        Yscr = Yscr.reshape((self.ny, 1))
        Erad = np.linspace(self.e_start, self.e_start + self.e_step * (self.ne - 1), num=self.ne)
        Erad = Erad.reshape((self.ne, 1, 1))
        prXconst = Xscr - x0
        prYconst = Yscr - y0
        phaseConstIn = np.pi * Erad / hc * (prXconst * prXconst + prYconst * prYconst) / (self.Distance - z0)
        phaseConstIn = phaseConstIn.flatten()
        self.arPhase += phaseConstIn
        cosf = np.cos(phaseConstIn)
        sinf = np.sin(phaseConstIn)
        arReEx = self.arReEx * cosf - self.arImEx * sinf    # sum of cos
        arImEx = self.arImEx * cosf + self.arReEx * sinf    # sum of sin
        arReEy = self.arReEy * cosf - self.arImEy * sinf    # sum of cos
        arImEy = self.arImEy * cosf + self.arReEy * sinf    # sum of sin
        self.arReEx = arReEx
        self.arImEx = arImEx
        self.arReEy = arReEy
        self.arImEy = arImEy

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
        self.memory_screen = np.zeros(Nscr*5)
        self.arReEx = self.memory_screen[0:Nscr]
        self.arImEx = self.memory_screen[Nscr:2*Nscr]
        self.arReEy = self.memory_screen[2*Nscr:3*Nscr]
        self.arImEy = self.memory_screen[3*Nscr:4*Nscr]
        self.arPhase = self.memory_screen[4*Nscr:5*Nscr]
        self.Xph = np.linspace(self.x_start, self.x_start + self.x_step*(self.nx -1), self.nx)
        self.Yph = np.linspace(self.y_start, self.y_start + self.y_step*(self.ny -1), self.ny)
        self.Eph = np.linspace(self.e_start, self.e_start + self.e_step*(self.ne -1), self.ne)

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
        self.memory_screen = np.zeros(Nscr)

        self.arReEx = self.memory_screen[0:Nscr]
        self.arImEx = self.memory_screen[Nscr:2*Nscr]
        self.arReEy = self.memory_screen[2*Nscr:3*Nscr]
        self.arImEy = self.memory_screen[3*Nscr:4*Nscr]
        self.arPhase = self.memory_screen[4*Nscr:5*Nscr]

    def nullify(self):
        Nscr = self.ne*self.nx*self.ny
        #self.current = 1
        self.memory_screen = np.zeros(Nscr*5)
        self.arReEx = self.memory_screen[0:Nscr]
        self.arImEx = self.memory_screen[Nscr:2*Nscr]
        self.arReEy = self.memory_screen[2*Nscr:3*Nscr]
        self.arImEy = self.memory_screen[3*Nscr:4*Nscr]
        self.arPhase = self.memory_screen[4*Nscr:5*Nscr]

    def create_like(self, em_screen):
        # added new check is needed
        # self.x = em_screen.x
        # self.y = em_screen.y
        # self.z = em_screen.z
        # self.size_x = em_screen.size_x
        # self.size_y = em_screen.size_y
        # self.nx = em_screen.nx
        # self.ny = em_screen.ny
        # self.start_energy = em_screen.start_energy
        # self.end_energy = em_screen.end_energy
        # self.num_energy = em_screen.num_energy
        # added new check is needed

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
        self.memory_screen = np.zeros(len(em_screen.memory_screen))

        self.arReEx = self.memory_screen[0:Nscr]
        self.arImEx = self.memory_screen[Nscr:2*Nscr]
        self.arReEy = self.memory_screen[2*Nscr:3*Nscr]
        self.arImEy = self.memory_screen[3*Nscr:4*Nscr]
        self.arPhase = self.memory_screen[4*Nscr:5*Nscr]

    def screenPy2C(self, lperiod, nperiods, status):

        scrPrm = np.zeros(14)
        
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
        """
        On the area ds during 1 sec falls dN photons in spectral width (dlambda/lambda)
        dN = ds/Distance**2 * (dlambda/lambda) * (I/qe) * 3*alpha*gamma**2/(4*pi**2) * |Eul(lambda, Xscreen)|**2
        Eul(lambda, Xscreen) is unitless electric field: Eul(lambda, Xscreen) = - (c/qe) * D/(sqrt(3)*gamma**2) * (E(lambda, Xscreen))

        :param gamma:
        :param current: in A
        :return:
        """

        LenPntrConst = self.Distance - self.Zstart

        # old constant with current in [mA]
        # constQuant = 3.461090202456155e+9*current*gamma*gamma/LenPntrConst/LenPntrConst

        constQuant = 3*alpha/q_e/(4*pi**2)*1e-3 * current * gamma * gamma / LenPntrConst / LenPntrConst
        Ex2r = np.array(self.arReEx)*np.array(self.arReEx)
        Ex2i = np.array(self.arImEx)*np.array(self.arImEx)
        self.Sigma = (Ex2r + Ex2i)*constQuant
        Ey2r = np.array(self.arReEy)*np.array(self.arReEy)
        Ey2i = np.array(self.arImEy)*np.array(self.arImEy)
        self.Pi = (Ey2r + Ey2i)*constQuant
        self.Total = self.Pi + self.Sigma
        self.Xph = np.linspace(self.x_start, self.x_start + self.x_step*(self.nx -1), self.nx)
        self.Yph = np.linspace(self.y_start, self.y_start + self.y_step*(self.ny -1), self.ny)
        self.Eph = np.linspace(self.e_start, self.e_start + self.e_step*(self.ne -1), self.ne)

    def coherent_photon_dist(self):
        """
        On the area ds during 1 sec falls dN photons in spectral width (dlambda/lambda)
        dN = ds/Distance**2 * (dlambda/lambda) * (I/qe) * 3*alpha*gamma**2/(4*pi**2) * |Eul(lambda, Xscreen)|**2
        Eul(lambda, Xscreen) is unitless electric field: Eul(lambda, Xscreen) = - (c/qe) * D/(sqrt(3)*gamma**2) * (E(lambda, Xscreen))

        For coherent radiation calculation:
        I = qe
        dN = ds/Distance**2 * (dlambda/lambda) * 3*alpha/(4*pi**2) * |Eul(lambda, Xscreen) * gamma * n_e|**2
        |Eul(lambda, Xscreen) * gamma * n_e| is calculated in function coherent_radiation()

        :return:
        """

        LenPntrConst = self.Distance - self.Zstart

        constQuant = 3*alpha/(4*pi**2)*1e-3/LenPntrConst / LenPntrConst
        Ex2r = np.array(self.arReEx) * np.array(self.arReEx)
        Ex2i = np.array(self.arImEx) * np.array(self.arImEx)
        self.Sigma = (Ex2r + Ex2i) * constQuant
        Ey2r = np.array(self.arReEy) * np.array(self.arReEy)
        Ey2i = np.array(self.arImEy) * np.array(self.arImEy)
        self.Pi = (Ey2r + Ey2i) * constQuant
        self.Total = self.Pi + self.Sigma
        self.Xph = np.linspace(self.x_start, self.x_start + self.x_step * (self.nx - 1), self.nx)
        self.Yph = np.linspace(self.y_start, self.y_start + self.y_step * (self.ny - 1), self.ny)
        self.Eph = np.linspace(self.e_start, self.e_start + self.e_step * (self.ne - 1), self.ne)

    """
    def convolution(self, sigma_x, sigma_y):
        data = numpy.reshape(self.Total,(len(self.Xph), len(self.Xph)))

        conv.convolution(self.Xph , self.Yph, data, sigma_x, sigma_y)
    """
        
    def zerosArray(self):
        Nscr = self.nx*self.ny*self.ne
        self.memory_screen = np.zeros(Nscr*5)
        self.arReEx = self.memory_screen[0:Nscr]
        self.arImEx = self.memory_screen[Nscr:2*Nscr]
        self.arReEy = self.memory_screen[2*Nscr:3*Nscr]
        self.arImEy = self.memory_screen[3*Nscr:4*Nscr]
        self.arPhase = self.memory_screen[4*Nscr:5*Nscr]

    def screen2dict(self):
        new_dict = {}
        new_dict["memory_screen"] = list(self.memory_screen)
        new_dict["x_step"] = self.x_step
        new_dict["x_start"] = self.x_start
        new_dict["nx"] = self.nx
        new_dict["y_step"] = self.y_step
        new_dict["y_start"] = self.y_start
        new_dict["ny"] = self.ny
        new_dict["e_step"] = self.e_step
        new_dict["e_start"] = self.e_start
        new_dict["ne"] = self.ne
        new_dict["Zstart"] = self.Zstart
        new_dict["Distance"] = self.Distance

        new_dict["x"] = self.x
        new_dict["y"] = self.y
        new_dict["z"] = self.z
        new_dict["size_x"] = self.size_x
        new_dict["size_y"] = self.size_y
        new_dict["nx"] = self.nx
        new_dict["ny"] = self.ny
        new_dict["start_energy"] = self.start_energy
        new_dict["end_energy"] = self.end_energy
        new_dict["num_energy"] = self.num_energy
        return new_dict

    def dict2screen(self, dictionaty):
        for key in dictionaty.keys():
            self.__dict__[key] = dictionaty[key]
        Nscr = self.nx * self.ny * self.ne
        self.memory_screen = np.array(self.memory_screen)
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

    sinfa = np.sin(screen_down.arPhase)
    cosfa = np.cos(screen_down.arPhase)
    screen.arReEx = screen_down.arReEx + screen_up.arReEx*cosfa - screen_up.arImEx*sinfa
    screen.arImEx = screen_down.arImEx + screen_up.arImEx*cosfa + screen_up.arReEx*sinfa
    screen.arReEy = screen_down.arReEy + screen_up.arReEy*cosfa - screen_up.arImEy*sinfa
    screen.arImEy = screen_down.arImEy + screen_up.arImEy*cosfa + screen_up.arReEy*sinfa

    screen.arPhase = screen_down.arPhase + screen_up.arPhase
    return screen

