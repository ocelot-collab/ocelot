__author__ = 'Sergey Tomin'
from em_screen import EMScreen
from xframework.common.globals import *
from trajectory.tr_solver import trajectory_body
from trajectory.show_traject import showMeTrajectory
from convolution.convolution_gauss import convolution_2D_cpp
from show_screen import any_show as shw
from display import display_undulator_param
from save_intensity import save_intens as svr
from ctypes import CDLL, c_double, c_int, POINTER
from time import clock
from numpy import sqrt, shape, zeros, pi, exp
#from numpy import zeros, int32
from sys import path
from os import name as os_name
from os import environ




#from pylab import *
flag_pyOCL = True
try:
    import pyopencl as cl
except:
    flag_pyOCL = False

if flag_pyOCL:
    import radiation_cl as rad_cl

if os_name == "nt":
    pathToDll = path[0] + r"/radiation"
else:
    pathToDll = path[0] [:-7]+ "codes/renera/gcc_codes/radiation"
#print pathToDll
environ['PATH'] = environ['PATH'] + ";"+pathToDll
#print os.environ['PATH']

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
    #print len(list_lists_motion), len(list_lists_motion[0])
    #list_motion = list_of_lists_to_LIST(list_lists_motion)
    #print len(list_motion)
    if flag_pyOCL and mode_proc == "GPU":
        platform = cl.get_platforms()
        device = platform[0].get_devices()
        max_work_item = device[0].get_info(cl.device_info.MAX_WORK_ITEM_SIZES)

        max_size = max(max_work_item)*2
        if (em_screen.nx > max_size or em_screen.ne > max_size  or em_screen.ny > max_size) and mode_proc != "CPU":
            mode_proc = "CPU"
            print "MAX_WORK_ITEM_SIZES = ", max_work_item
            print "[xNstep, yNstep, eNstep] = ", em_screen.nx, em_screen.ny, em_screen.ne
            print "GPU-mode was changed on CPU-mode", max_work_item
    if flag_pyOCL == False and mode_proc == "GPU":
        mode_proc = "CPU"
        print "GPU-mode was changed on CPU-mode. pyOpenCL is absent "
    em_screen.zerosArray()
    if mode_proc == "CPU":
        if os_name == "nt":
            pathToDll = path[0] [:-7]+r"codes/genera/build/Radiation.dll"
        else:
            pathToDll = path[0] [:-7]+"codes/genera/build/radiation.so"

        my_rad = CDLL(pathToDll)

        c_scrPrm = em_screen.screenPy2C(undulator.lperiod, undulator.nperiods, undulator.status)

        ppMotion, ptr_nPoints_elements = pointer_to_list(list_motion)
        start = clock()
        dPtr = POINTER(c_double)
        ret = my_rad.emsolver(c_int(len(list_motion)),
                              ptr_nPoints_elements,
                              ppMotion,
                              c_double(gamma),
                              c_scrPrm,
                              em_screen.memory_screen.ctypes.data_as(dPtr))

        print "CPU time execution: ", clock() - start
        if ret != 1:
            print "radiation return = ", ret


    elif mode_proc == "GPU" and flag_pyOCL == True:
        em_screen = rad_cl.rad_cl(em_screen, list_motion, gamma, undulator)

    else:
        print " mode - ERROR !! (emsolver.py) "
    #print "beam current: ", beam_current
    em_screen.distPhoton(gamma, current = beam_current)
    return em_screen

from emitt_spread import   change_sizes_screen, convolution_all

def solver(screen, cell, beam):
    """
    1. find trajectory and mag field on the trajectory. system of unit is [mm,rad]
    2. display undulator parameter, for the controlling of calculation
    3. creation new class EMScreen. transform units of Screen [m] to units of EMScreen [mm]
        New system of definition em_screen.x_start, x_step and etc
    4. calculation of radiation.
        a. Current must be in mA.
        b. mode_proc can be "CPU" or "GPU". code can choose automatically mode_proc if OS has not pyopencl or dimentions of task are very large
        c. Code can choose right library (dll / so) automatically depend on OS
        d. list_lists_motion will be transform to list_motion
    """
    accuracy = 2
    show_param = True
    for elem in cell:
        if elem.type == "undulator":
            undulator = elem

            #elem.status = 13


    list_lists_motion = trajectory_body(cell, beam, accuracy = accuracy, mode_traj = "radiation")
    list_motion = list_of_lists_to_LIST(list_lists_motion)

    showMeTrajectory(list_motion)
    #beam.I = beam.Q

    display_undulator_param(screen, cell, beam, list_motion, show_param)
    em_screen = EMScreen(screen)
    beam_current = beam.I*1000


    #bool_es = energy_spread_analyse(em_screen, beam)
    print "before x:", em_screen.nx, em_screen.x_start, em_screen.x_step
    print "before y:", em_screen.ny, em_screen.y_start, em_screen.y_step
    print "before e:", em_screen.ne, em_screen.e_start, em_screen.e_step
    change_sizes_screen(em_screen, beam)
    print "after x:", em_screen.nx, em_screen.x_start, em_screen.x_step, em_screen.nx_add
    print "after y:", em_screen.ny, em_screen.y_start, em_screen.y_step, em_screen.ny_add
    print "after e:", em_screen.ne, em_screen.e_start, em_screen.e_step, em_screen.ne_add

    em_screen = radiation(em_screen, list_motion, beam.gamma, beam_current, undulator, mode_proc = "CPU")
    #print "Ok", shape(em_screen.Total)
    convolution_all(em_screen)
    #print "after conv x:", em_screen.nx, em_screen.x_start, em_screen.x_step, em_screen.nx_add
    #print "after conv y:", em_screen.ny, em_screen.y_start, em_screen.y_step, em_screen.ny_add
    #print "after conv e:", em_screen.ne, em_screen.e_start, em_screen.e_step, em_screen.ne_add

    #print len(em_screen.Eph), em_screen.Eph[:10]

    #print
    #print "totla 2 = ", shape(em_screen.Total)
    #emittance_effect(beam, screen)

    #energy_spread_effect(em_screen, beam, bool_es)
    #print em_screen.Xph

    #print shape(em_screen.Total)
    shw(em_screen)
    svr(em_screen, "intens.gnri")
    return em_screen



"""
def emittance_effect_analys(beam, screen):
    pass

def emittance_effect(beam, screen):
    if beam.beta_x> 0 or beam.beta_y > 0:
        beam.sigma_x = sqrt(beam.emit_x*beam.beta_x)
        beam.sigma_y = sqrt(beam.emit_y*beam.beta_y)
        beam.sigma_x1 = sqrt(beam.emit_x/beam.beta_x)
        beam.sigma_y1 = sqrt(beam.emit_y/beam.beta_y)
        screen.beam_size_x = sqrt(beam.sigma_x*beam.sigma_x + (beam.sigma_x1*screen.z)**2)
        screen.beam_size_y = sqrt(beam.sigma_y*beam.sigma_y + (beam.sigma_y1*screen.z)**2)
    else:
        screen.beam_size_x = 0
        screen.beam_size_y = 0
        return 0
    if screen.beam_size_x and screen.beam_size_y:
        print "calculkate emittance effects .... ",
        screen.Sigma = convolution_cpp(screen.Sigma, screen.Xph, screen.Yph, screen.beam_size_x*1000., screen.beam_size_y*1000.)
        screen.Pi = convolution_cpp(screen.Pi, screen.Xph, screen.Yph, screen.beam_size_x*1000., screen.beam_size_y*1000.)
        screen.Total = screen.Sigma + screen.Pi
        print "OK"
    elif (screen.beam_size_x and not screen.beam_size_y) or (not screen.beam_size_x and screen.beam_size_y):
        print "one of the emittances is equal to zero. Calculation of the emittance effect will skip"
        return 0

def energy_spread_analyse(screen, beam):
    if beam.sigma_E:
        print "energy spread effect is switched on!"
        if screen.ne == 1:
            sigma_e = 2*beam.sigma_E/beam.E*screen.fund_harm_eV
            # part of class Screen!!!
            screen.e_start = screen.e_start - 3.*sigma_e
            screen.ne  = 7
            screen.e_step = float(6*sigma_e)/(screen.ne - 1)

        elif screen.nx != 1 or screen.ny !=1:
            exit("screen.ne != 1 and screen.nx or screen.ny != 1. Conditions is wrong")
        else:
            return 10
        return 1
    else:
        print "energy spread effect is switched off!"
        return 0

def energy_spread_effect(screen, beam, bool_es):
    #ypoint*xpoint*je + xpoint*jy + jx
    if bool_es == 0:
        return 0
    se = 2*beam.sigma_E/beam.E*screen.fund_harm_eV
    if bool_es ==1:

        print "se = ", se
        k = lambda energy: exp(-(energy - screen.fund_harm_eV)**2/(2.*se*se))/(sqrt(2.*pi)*se)
        print "coeff = ",k(screen.Eph)
        nx_ny = len(screen.Xph)*len(screen.Yph)
        new_Pi = zeros(nx_ny)
        new_Sigma = zeros(nx_ny)
        for i, eph in enumerate(screen.Eph[:-1]):
            data_Pi = screen.Pi[nx_ny*i:nx_ny*(i+1)]*k(eph)*(screen.Eph[i+1]-eph)
            new_Pi += data_Pi
            data_Sigma = screen.Sigma[nx_ny*i:nx_ny*(i+1)]*k(eph)*(screen.Eph[i+1]-eph)
            new_Sigma += data_Sigma
        screen.Sigma = new_Sigma
        screen.Pi = new_Pi
        screen.Total = screen.Sigma + screen.Pi
        screen.Eph = screen.Eph[int(len(screen.Eph)/2)]
        screen.ne = 1
        print "shape = ", shape(screen.Total )
        print "Eph = ", screen.Eph
    elif bool_es ==10:
        from convolution.convolution_1D import convolution_1D_cpp
        screen.Sigma = convolution_1D_cpp(screen.Sigma, screen.Eph,se)
        screen.Pi = convolution_1D_cpp(screen.Pi, screen.Eph,se)
        screen.Total = screen.Sigma + screen.Pi

"""
