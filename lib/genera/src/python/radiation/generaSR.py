__author__ = 'Sergey Tomin'

from em_screen import EMScreen, show_flux
from ocelot.common.globals import *
#from codes.genera.src.python.trajectory.tr_solver import trajectory_body
from ocelot.lib.genera.src.python.trajectory.lat_trajectory import trace4radiation
from emitt_spread import  change_sizes_screen, convolution_all
from ocelot.cpbd.beam import *
from ctypes import CDLL, c_double, c_int, POINTER
from numpy import array, zeros, abs
from sys import path
import os
from scipy.special import *
from ocelot.common.globals import *
from time import time
from os import name as os_name
from sys import platform
import socket

"""
if os_name == "nt":
    tail = r"codes/genera/build/windows/Radiation.dll"
else:
    if socket.gethostname() == "exflwgs03":
        tail = "codes/genera/build/exflwgs03/radiation.so"
    else:
        tail = "codes/genera/build/ubuntu/radiation.so"

    if platform == "darwin":
        tail = "codes/genera/build/mac/radiation.so"
"""

tail = "ocelot/lib/genera/build/genera_libs/radiation.so"
home_dir = path[0]
index =  path[0].find("siberia2")
pathToDll = path[0][:index]+ tail
try:
    my_rad = CDLL(pathToDll)
except:
    exec(open(path[0][:index]+ "ocelot/lib/genera/src/cpp/compile.py"))
    os.chdir(home_dir)
    my_rad = CDLL(pathToDll)


class Trajectory:
    """Charged Particle Trajectory"""
    arX = array('d') #arrays of horizontal, vertical and longitudinal positions [m] and relative velocities
    arXp = array('d')
    arY = array('d')
    arYp = array('d')
    arZ = array('d')
    arZp = array('d')
    np = 0 #number of trajectory points
    ctStart = 0 #start and end values of independent variable (c*t) for which the trajectory should be (/is) calculated (is constant step enough?)
    ctEnd = 0


def list_of_lists_to_LIST(list_lists):
    LIST = []
    for list in list_lists:
        LIST.extend(list)
    return LIST


def pointer_to_list(list_motion):
    number_elements = len(list_motion)

    dPtr = POINTER(c_double)

    nPoints_elements = []#zeros((number_elements) , dtype = int32) # dtype = int32 it can cause PROBLEMS!!!
    B = []
    M = [0.]*number_elements
    for j in range(number_elements):
        motion = list_motion[j]
        nPoints_elements.append(len(motion.Z))
        #print "npoints on traject =",  nPoints_elements
        M[j] = motion.memory_motion
        B.append(M[j].ctypes.data_as(dPtr))
    pointer = dPtr*number_elements
    ppMotion = pointer(*B)
    arr_type =  c_int*number_elements
    ptr_nPoints_elements = arr_type(*nPoints_elements)
    return ppMotion, ptr_nPoints_elements


def radiation(em_screen, list_motion, gamma, beam_current,undulator, mode_proc):
    em_screen.zerosArray()
    c_scrPrm = em_screen.screenPy2C(undulator.lperiod, undulator.nperiods, undulator.status)
    """  undulator.status is just for fast mode calculation (ideal undulator and zero initial conditions)  """
    ppMotion, ptr_nPoints_elements = pointer_to_list(list_motion)

    dPtr = POINTER(c_double)
    ret = my_rad.emsolver(c_int(len(list_motion)),
            ptr_nPoints_elements,
            ppMotion,
            c_double(gamma),
            c_scrPrm,
            em_screen.memory_screen.ctypes.data_as(dPtr))

    if ret != 1:
        print "radiation return = ", ret

    em_screen.distPhoton(gamma, current = beam_current)
    return em_screen


def motion_to_trj(list_motions):
    trj = Trajectory()
    arX = []
    arY = []
    arZ = []
    arXp = []
    arYp = []
    arZp = []
    for motion in list_motions:
        arX.extend(array(motion.X))
        arY.extend(array(motion.Y))
        arZ.extend(array(motion.Z))
        arXp.extend(array(motion.Xbeta))
        arYp.extend(array(motion.Ybeta))
        arZp.extend(array(motion.Zbeta))
        trj.np += len(motion.Z)
    trj.arX = array(arX)
    trj.arY = array(arY)
    trj.arZ = array(arZ)
    trj.arXp = array(arXp)
    trj.arYp = array(arYp)
    trj.arZp = array(arZp)
    return trj


def define_status(cell, beam, mode_traj):
    """
    new feature - the function analyses the parameters of task and define one of two possible modes:
    a. fast calculation ( without ending poles)
    b. classic calculation ( with ending poles)
    """

    undul_ok = 0
    beam_ok = 0

    for elem in cell.sequence:

        if len(cell.sequence) == 1 and elem.type == "undulator" and elem.ax < 0 and mode_traj != "trajectory":
            if elem.field_map == None :
                undul_ok = 1
        #print beam.x,beam.y,beam.xp,beam.yp
        if beam.x == 0 and beam.y == 0 and beam.xp == 0 and beam.yp == 0:
            beam_ok = 1
            #print "####"
        if undul_ok == 1 and beam_ok:
            elem.status = 13
            #beam.xp = - elem.Kx/beam.E*m_e_GeV
            #beam.yp = - elem.Ky/beam.E*m_e_GeV
            #print beam.xp
            #print "fast calculation"
        else:
            elem.status = 0

def data_format(emscreen):
    #ypoint*xpoint*je + xpoint*jy + jx -->> array(arI1).reshape(screen.nx, screen.ny,screen.num_energy)

    intens = zeros((emscreen.nx, emscreen.ny, emscreen.ne))
    total = emscreen.Total.reshape(emscreen.ne, emscreen.ny,emscreen.nx)
    for i in xrange(emscreen.ne):
        intens[:,:,i] = total[i,:,:]
    return intens


def checking_step(lat, screen, beam, list_motions):

    Q = 0.5866740802042227 #mm^-1*T^-1
    gamma = beam.E/m_e_GeV
    g2 = gamma*gamma
    z = screen.z - lat.totalLen
    xc_max = screen.size_x + screen.x
    yc_max = screen.size_y + screen.y

    lph = h_eV_s*speed_of_light/screen.end_energy
    Bmax = []
    Rx = []
    Ry = []
    L = 0
    #import matplotlib.pyplot as plt
    for motion in list_motions:
        if len(motion.Z)>3:
            z = screen.z - motion.Z
            rx = max(abs(xc_max/z*0 - motion.X/z - motion.Xbeta))
            ry = max(abs(yc_max/z*0 - motion.Y/z - motion.Ybeta))
            B = max(motion.Bx**2 + motion.By**2)
            Bmax.append(B)
            L = motion.Z[-1] - motion.Z[0]
            Rx.append(rx)
            Ry.append(ry)
    N1 = L*Q*max(Bmax)
    rx = max(Rx)
    ry = max(Ry)

    df = pi/(lph*g2)*(1 + g2*(rx*rx + ry*ry))
    dz = pi/df/2.
    N2 = L/dz/1000
    #print "N2/N1 = ", N2, "/",N1
    #for i in xrange(int(len(list_motions)/11)):
    #    df_dz = pi/(lph*g2)*(1 + g2*())
    #print "N1+N2 = ", 1.3*(N1+N2 + 100)+100
"""
def calculateSR_py(lat, beam, screen, runParameters = None):

    #1. find trajectory and mag field on the trajectory. system of unit is [mm,rad]
    #2. display undulator parameter, for the controlling of calculation
    #3. creation new class EMScreen. transform units of Screen [m] to units of EMScreen [mm]
    #    New system of definition em_screen.x_start, x_step and etc
    #4. calculation of radiation.
    #    a. Current must be in mA.
    #    b. mode_proc can be "CPU" or "GPU". code can choose automatically mode_proc if OS has not pyopencl or dimentions of task are very large
    #    c. Code can choose right library (dll / so) automatically depend on OS
    #    d. list_lists_motion will be transform to list_motiona

    accuracy = 1
    for elem in lat.sequence:
        if elem.type == "undulator":
            undulator = elem


    list_lists_motion = trajectory_body(lat, beam, accuracy = accuracy, mode_traj = "radiation")
    list_motions = list_of_lists_to_LIST(list_lists_motion)

    checking_step(lat, screen, beam, list_motions)

    trj = motion_to_trj(list_motions)
    #showMeTrajectory(list_motions)
    #display_undulator_param(screen, cell, beam, list_motion, show_param)
    em_screen = EMScreen(screen)
    beam_current = beam.I*1000

    #print "before x:", em_screen.nx, em_screen.x_start, em_screen.x_step
    #print "before y:", em_screen.ny, em_screen.y_start, em_screen.y_step
    #print "before e:", em_screen.ne, em_screen.e_start, em_screen.e_step
    change_sizes_screen(em_screen, beam)
    #print "after x:", em_screen.nx, em_screen.x_start, em_screen.x_step, em_screen.nx_add
    #print "after y:", em_screen.ny, em_screen.y_start, em_screen.y_step, em_screen.ny_add
    #print "after e:", em_screen.ne, em_screen.e_start, em_screen.e_step, em_screen.ne_add

    em_screen = radiation(em_screen, list_motions, beam.gamma, beam_current, undulator, mode_proc = "CPU")

    convolution_all(em_screen)
    intens = data_format(em_screen)
    return trj, em_screen
"""

def print_rad_props(beam, K, lu, L, E, distance):
    print "********* e beam ***********"
    beam.print_sizes()


    def F_n(n, Ku):
        v1 = ((n-1)/2.)
        v2 = ((n+1)/2.)
        x = n*Ku*Ku/(4.+2.*Ku*Ku)
        return n*n*Ku*Ku/((1+Ku*Ku/2.)**2)*(jn(v1,x) - jn(v2,x))**2

    def flux(I, K, m):
        alpha = 1/137.036
        gm = E/m_e_GeV
        Nu = L/lu
        e = 1.602e-19
        BW = 0.001 # band width 0.1%
        k_rad_mrad = 1e-6 # coef to converse rad to mrad
        F = alpha*gm*gm*I/e*Nu*Nu*F_n(m, K)*BW*k_rad_mrad
        return F
    gamma = E/m_e_GeV
    Lambda = K2Lambda(K, lu, E)
    sigma_r = sqrt(Lambda*L/(2*4.*pi*pi))
    sigma_r1 = sqrt(Lambda/L/2.)
    Sigma_x = sqrt(beam.sigma_x**2 + sigma_r**2)
    Sigma_y = sqrt(beam.sigma_y**2 + sigma_r**2)
    Sigma_x1 = sqrt(beam.sigma_xp**2 + sigma_r1**2)
    Sigma_y1 = sqrt(beam.sigma_yp**2 + sigma_r1**2)
    size_x = sqrt(Sigma_x**2 + (Sigma_x1*distance)**2)
    size_y = sqrt(Sigma_y**2 + (Sigma_y1*distance)**2)
    B = K2field(K, lu = lu)
    F = flux(beam.I, K, m = 1)
    N = L/lu
    flux_tot = 1.431e14*beam.I*N*F_n(1, K)*(1.+K*K/2.)/1./2.

    brightness = flux_tot/(4*pi*pi*Sigma_x*Sigma_y*Sigma_x1*Sigma_y1)*1e-12

    print "********* ph beam ***********"
    print "Ebeam        : ", E
    print "K            : ", K
    print "B            : ", B, " T"
    print "lambda       : ", Lambda, " m "
    print "Eph          : ", K2Ephoton(K, lu, E), " eV"
    print "1/gamma      : ", 1./gamma *1e6, " um"
    print "sigma_r      : ", sigma_r*1e6, " um"
    print "sigma_r'     : ", sigma_r1*1e6, " urad"
    print "Sigma_x      : ", Sigma_x *1e6, " um"
    print "Sigma_y      : ", Sigma_y *1e6, " um"
    print "Sigma_x'     : ", Sigma_x1*1e6, "urad"
    print "Sigma_y'     : ", Sigma_y1*1e6, "urad"
    print "H. spot size : ", size_x*1000., "/", size_x/distance*1000., " mm/mrad"
    print "V. spot size : ", size_y*1000., "/", size_y/distance*1000., " mm/mrad"
    print "I            : ", beam.I, " A"
    print "Nperiods     : ", L/lu
    print "distance     : ", distance, " m"
    print "flux tot     : ", flux_tot, " ph/sec"
    print "flux density : ", F, " ph/sec/mrad^2;   ", F/distance/distance, " ph/sec/mm^2"
    #print "flux density : ", F/distance/distance, " ph/sec/mm^2"
    print "brilliance   : ", brightness, " ph/sec/mrad^2/mm^2"



def calculateSR_py(lat, beam, screen, runParameters = None):
    """
    1. find trajectory and mag field on the trajectory. system of unit is [mm,rad]
    2. display undulator parameter, for the controlling of calculation
    3. creation new class EMScreen. transform units of Screen [m] to units of EMScreen [mm]
        New system of definition em_screen.x_start, x_step and etc
    4. calculation of radiation.
        a. Current must be in mA.
        b. mode_proc can be "CPU" or "GPU". code can choose automatically mode_proc if OS has not pyopencl or dimentions of task are very large
        c. Code can choose right library (dll / so) automatically depend on OS
        d. list_lists_motion will be transform to list_motiona
    """
    accuracy = 2

    #print "in calculator ", screen.size_x, screen.x
    for elem in lat.sequence:
        if elem.type == "undulator":
            print_rad_props(beam, elem.Kx, elem.lperiod, elem.l, lat.energy, screen.z)
            undulator = elem
            undulator.status = 0
    beam.gamma = beam.E/m_e_GeV
    particle0 = Particle(x=beam.x, y=beam.y, px=beam.xp, py=beam.xp, s=0.0, p=0,  tau=0)
    list_motions = trace4radiation(lat,particle0, accuracy = accuracy)
    #TODO: include in process checking_step
    checking_step(lat, screen, beam, list_motions)

    trj = motion_to_trj(list_motions)
    #showMeTrajectory(list_motions)
    #display_undulator_param(screen, cell, beam, list_motion, show_param)
    em_screen = EMScreen(screen)

    beam_current = beam.I*1000

    print "before x:", em_screen.nx, em_screen.x_start, em_screen.x_step
    print "before y:", em_screen.ny, em_screen.y_start, em_screen.y_step
    print "before e:", em_screen.ne, em_screen.e_start, em_screen.e_step
    change_sizes_screen(em_screen, beam)
    print "after x:", em_screen.nx, em_screen.x_start, em_screen.x_step, em_screen.nx_add
    print "after y:", em_screen.ny, em_screen.y_start, em_screen.y_step, em_screen.ny_add
    print "after e:", em_screen.ne, em_screen.e_start, em_screen.e_step, em_screen.ne_add
    start = time()
    em_screen = radiation(em_screen, list_motions, beam.gamma, beam_current, undulator, mode_proc = "CPU")
    print "radiation solver: ", time() - start, " sec"
    start = time()

    convolution_all(em_screen)
    print "convolution solver: ", time() - start, " sec"
    #intens = data_format(em_screen)
    #show_flux(em_screen)
    return trj, em_screen