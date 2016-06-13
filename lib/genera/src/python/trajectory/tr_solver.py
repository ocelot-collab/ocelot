__author__ = 'Sergey Tomin'

from numpy import where, array, delete, split
from sys import exit, path
import motion as mtn
from xframework.common.globals import *
from ctypes import CDLL, c_double, POINTER, c_int
from os import name as os_name
from xframework.cpbd.beam import Particle
from undulator import field_map_for_dll, define_points_undul
from sys import platform

und_type = "undulator"
import socket
if os_name == "nt":
    tail = r"codes/genera/build/windows/lattice_trajectory.dll"
    tail_u = r"codes/genera/build/windows/undulator.dll"
else:
    if socket.gethostname() == "exflwgs03":
        tail_u = "codes/genera/build/exflwgs03/undulator.so"
        tail = "codes/genera/build/exflwgs03/trajectsolver.so"
    else:
        tail_u = "codes/genera/build/ubuntu/undulator.so"
        tail = "codes/genera/build/ubuntu/lattice_trajectory.so"

    if platform == "darwin":
        tail_u = "codes/genera/build/mac/undulator.so"
        tail = "codes/genera/build/mac/lattice_trajectory.so"

index =  path[0].find("scripts")
pathToDll_u = path[0][:index]+ tail_u
pathToDll = path[0][:index]+ tail

undulator_dll = CDLL(pathToDll_u)
myfunc = CDLL(pathToDll)


def errorsBlock(error):
    def printOfStars(phrases, error):
        print '**********************************************'
        print '*', phrases, error
        print '**********************************************'
    if error>0:
        phrases =  "Errors Block - OK. Error cod = "
    elif error == -100:
        phrases =  "Errors Block - Unknow type of undulator. Error code = "
    elif error == -300:
        phrases =  "Errors Block - Z must be increase. Error code = "
    elif error == -301:
        phrases =  "Errors Block - X must be increase. Error code = "
    elif error == -302:
        phrases =  "Errors Block - Y must be increase. Error code = "
    elif error == -200:
        phrases = "x or y out of bounds. spline.clspline2d. Error code = "
    else:
        phrases = "Unknown error. Error code = "
    if error<0:
        printOfStars(phrases, error)


def Py2C(array):
    arr_type =  c_double*len(array)
    c_array = arr_type(*array)
    return c_array


def divide_cell(cell):
    parts_cell = []
    types = []
    for element in cell:
        types.append(element.type)
    indexes_of_und = where(array(types)==und_type)[0]
    if len(indexes_of_und) == 0:
        parts_cell.append(cell)
        return parts_cell
    else:
        indexes_of_und = delete(indexes_of_und, 0)
        parts_cell = split(array(cell), indexes_of_und)
        return parts_cell


def define_points(element, accuracy):
    accuracy = int(accuracy)
    if element.type != und_type:
        element.npoints_traj = int(element.l*10*1.5+10)*accuracy
    else:
        element.npoints_traj = define_points_undul(element, accuracy)
        """
        if element.field_map.type == "analytic":
            element.npoints_traj = int(element.field_map.nperiods*50+10)*accuracy
        else:
            element.npoints_traj = int(((element.field_map.Z_arr[-1] - element.field_map.Z_arr[0])/0.5)+10)*accuracy
        """

def create_motions(parts_cell, mode):
    list_list_motion = []
    bRough = 0
    for part in parts_cell:
        list_motion = []
        for elem in part:
            npoints_traj = elem.npoints_traj
            if mode == "radiation":
                size_motion = npoints_traj*3-1
            else:
                size_motion = npoints_traj

            if elem.type == und_type: # structure.type == "undulator"
                #size_motion = size_motion*elem.field_file_rep - (elem.field_file_rep - 1)*(2-bRough)
                size_motion = size_motion 
                if elem.field_map:
                    size_motion = size_motion*elem.field_map.field_file_rep - (elem.field_map.field_file_rep - 1)*(2-bRough)
            motion = mtn.Motion(N = size_motion)
            #print "length memory motion:",  len(motion.memory_motion)
            motion.bRough = 0
            motion.mode = mode
            motion.type_element = elem.type
            list_motion.append(motion)
        list_list_motion.append(list_motion)
    return list_list_motion


def pointer_to_pointer(list_list_motion):
    array_pointer_pointer = []
    for list in list_list_motion:
        array_pointer = []
        size = len(list)
        dPtr = POINTER(c_double)
        M = [0.]*size
        for i, motion in enumerate(list):
            M[i] = motion.memory_motion
            array_pointer.append(M[i].ctypes.data_as(dPtr))
        pointer = dPtr*size
        ppMotion = pointer(*array_pointer)
        array_pointer_pointer.append(ppMotion)
    return array_pointer_pointer


def c_format_cell(parts_cell):
    structures = []
    for part in parts_cell:
        structure = []
        for elem in part:
            length = elem.l*1000 # [mm]
            npoints = elem.npoints_traj
            if elem.type == "drift":
                line = [10,length, 0.,0.,0.,npoints]
            elif elem.type == "sbend":
                line = [20, length, elem.k1*1e-6, elem.angle/length, 0., npoints]
            elif elem.type == "rbend":
                h = elem.angle/length # 1/ro
                phi = elem.angle/180.*pi
                line = [25, 0., phi/2, h, 0., 0] # first edge. See trajectory.accelerator_physics
                structure.extend(line)
                line = [20, length, elem.k1*1e-6, h, 0., npoints]
                structure.extend(line)
                line = [25, 0., phi/2, elem.angle/length, 0., 0] # second edge
            elif elem.type == "quadrupole":
                line = [30, length, elem.k1*1e-6, 0.,0.,npoints]
            elif elem.type == "undulator":
                line = [40, 0., 0.,0.,0.,npoints]

            else:
                exit("c_format_cell() is wrong. tr_solver.py")

            #line[0] = c_code_elem
            #line[1] = elem.L*1000 # transform m -> mm
            #line[2] = elem.k1/Bro*1e-6  #field in vertical direction creates 1/Rx [1/m^2 -> 1/mm^2]
            #line[3] = elem.angle/line[1] #field in vertical direction creates 1/Rx [1/m -> 1/mm]
            #line[4] = 0. # elem.angle/line[1] #field in horizontal direction creates Ry [1/m -> 1/mm]
            #line[5] = elem.npoints_traj   #number of points on trajectory
            structure.extend(line)
        structures.append(structure)
    return structures


def prepare_field(part_cell):
    field_file_rep = 1
    c_magField = Py2C([0])
    len_mf = c_int(0)
    cols_mf = c_int(0)
    c_undul_param = Py2C([0,0,0,0,0,0])
    for elem in part_cell:
        if elem.type == und_type:
            mag_field, len_mf, cols_mf, undul_param = field_map_for_dll(elem) # params in [T] and [mm]
            field_file_rep = elem.field_map.field_file_rep
            """
            if elem.status == 13:
                lenMF = c_int(13)
            """
            c_magField = Py2C(mag_field)
            c_undul_param = Py2C(undul_param)
    return c_magField, len_mf, cols_mf, c_undul_param, field_file_rep


def motion2particle(list_list_motion):
    list_particle = []
    for list_motion in list_list_motion:
        for motion in list_motion:

            for i in range(len(motion.Z)):
                particle = Particle()

                particle.s = motion.Z[i]
                particle.x = motion.X[i]
                particle.y = motion.XbetaI2[i]
                particle.xp = motion.Xbeta[i]
                particle.yp = motion.Ybeta[i]
                list_particle.append(particle)
    return list_particle


def define_status(cell, beam, mode_traj):
    """
    new feature - the function analyses the parameters of task and define one of two possible modes:
    a. fast calculation ( without ending poles)
    b. classic calculation ( with ending poles)
    """

    undul_ok = 0
    beam_ok = 0

    for elem in cell.sequence:

        if len(cell.sequence) == 1 and elem.type == und_type and elem.ax < 0 and mode_traj != "trajectory":
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


def trajectory_body(cell, beam, accuracy = 1, mode_traj = "radiation"):
    ACCURACY = accuracy
    define_status(cell, beam, mode_traj)   # fast or classic calculation: only for ideal undulator with infinite width
                                # zero conditions (x = y = xp = yp = 0)
                                # cell includes only undulator
    #print cell[0].status
    """
    The C++ code allows to accept only one map of magnetic field, that's why we must divide CELL on several part.
    Each part of CELL will include only one undulator. But also there is possibility that part of CELL does not include
    any undulator. in this case the C++ code will run as well.
    0. The C++ code can use two modes "trajectory" and "radiation" The first mode is for equidistant points | | | (bRough = 1).
        The second mode is designed for Gauss integration |...(\)...| bRough = 0(three points integration)
    1. defining number of points for each element of magnetic lattice. It is needed to develop.
    2. dividing CELL on parts which will include only one undulator or does not include any undulator
    3. Creation and preparation arrays of MOTION. I create the structure analog of CELL and
        I will store in the structure motion after that I will create class particle
    4(unit). prepare initial conditions: [X, Y, 1, Xp, Yp, 1] - and here we make the transform of SI UNIT [m,rad]
        to my system of UNIT [mm, rad](it is historical cause). Other big topic is dp/p = 1 in initial conditions
        Note: add to class beam gamma. please use gamma only from class beam.
    5(unit). prepare C++ format of parts of CELL (we will do this procedure the last because we need to know npoints_traj).
        Moreover, only here we make the transform of SI UNIT [m, rad] to my system of UNIT [mm, rad] (it is historical cause),
        and last procedure is adding element EDGE to C++ structure. EDGE is infinitely thin element and we dont allocate memory
        for it. So all procedures will be ignore element EDGE except matrix multiply in the C++ code
    6.  We must remember that reading magnetic map and preparation field structure from file will be done at the moment of
        initialization of undulator class. Simple structure for analytic undulator we will do in function prepare_field()
    """
    #0
    #mode_traj = "trajectory" # or "radiation"
    #mode = "radiation"
    bRough = 1
    if mode_traj == "trajectory":
        bRough = 1
    elif mode_traj == "radiation":
        bRough = 0
    else:
        exit("mode of TRACK is wrong")
        #1
    for elem in cell.sequence:
        define_points(elem, accuracy=ACCURACY)
        #2
    parts_cell = divide_cell(cell.sequence) # it will divide CELL on parts, ones will include undulator or not.
    #3
    list_list_motion = create_motions(parts_cell, mode_traj)
    array_p_p = pointer_to_pointer(list_list_motion)
    #4
    beam.gamma = beam.E/m_e_GeV
    xp = beam.xp
    yp = beam.yp
    
    if cell.sequence[0].type == und_type:
        if cell.sequence[0].status == 13:
            xp = - elem.Kx/beam.E*m_e_GeV 
            yp = - elem.Ky/beam.E*m_e_GeV 
    
    init_cond = [beam.x*1000., beam.y*1000., 0., xp, yp, beam.gamma] # [m] -> [mm]
    #print init_cond
    #5
    c_parts_cell = c_format_cell(parts_cell)

    for i, part in enumerate(parts_cell):
        c_init_cond = Py2C(init_cond)
        #6
        #search in parts_cell undulator and preparation C++ description of undulator
        c_magField,lenMF,cols,c_undul_param, field_file_rep = prepare_field(part)
        
        #print c_magField,lenMF,cols,c_undul_param, field_file_rep
        
        structure = c_parts_cell[i]
        #print structure
        size = int(len(structure)/6)
        c_structure = Py2C(structure)
        ppMotion = array_p_p[i]
        #print c_magField
        ret = myfunc.trajectory2(c_magField,    # 1D magnetic field array [X,Y,Z,Bx,By,Bz] or [Z,By] or [Z,By,Bx]
            cols,               # number of cols in magnetic data
            lenMF,              # number of rows in magnetic data
            bRough,             # if bRough = 1 then |  |  | else  | . . . (|) . . . | (gauss integration)
            c_init_cond,        # [beam.x [mm], beam.y [mm], 0 [long. pos. mm], beam.xp [rad!], beam.yp[rad!], gamma]
            c_undul_param,      # [Bx (amplitude [T]), By (amplitude [T]), phase (between Bx and By),
            # number_of_periods, lenPeriod [m], ax (width of undulator [m])]
            c_int(field_file_rep), # Number of repetition of field file. for example, undul = (end1, period*40, end2)
            ppMotion,           # pointer to array [motion.memory_motion (for und),
            # motion.memory_motion (for quad),  motion.memory_motion (for drift) etc...]
            c_structure,        # 1D array with structure of cell [drift,undulator, ...] =
            # [10, 2.0, 0.0, 0.0, 0.0, 3,  40, 0.0, 0.0, 0.0, 0.0, 9995, ...]
            # drift = [10, # - code of element (10 - drift, 20 - bend, 30 - quad, 40 - undulator)
            #       2.0, # (length [mm])
            #       0.0, # K/Bro*1e-6 [1/m^2 -> 1/mm^2] gradient, where Bro = gamma*m0*beta/c, K [T/m])
            #       0.0, # field in vertical direction creates [1/m -> 1/mm]
            #       0.0, # field in horizontal direction [1/m -> 1/mm]
            #       3]   # number of points on trajectory on element
            c_int(size))        # number of elements = len(c_structure)/6
        errorsBlock(ret)
        init_cond = list_list_motion[i][-1].defineInitCond(beam.gamma)
        """
        very important is to have possibility to send "list_list_motion" to radiation calculation
        """
    
    return list_list_motion


def track(cell, beam):

    list_list_motion = trajectory_body(cell, beam, accuracy= 1, mode_traj = "trajectory")
    list_particles = motion2particle(list_list_motion)
    return list_particles
