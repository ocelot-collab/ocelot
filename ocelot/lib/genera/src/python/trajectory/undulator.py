__author__ = 'Sergey Tomin'

from ctypes import CDLL, c_double, POINTER, c_int, cdll
from numpy import pi, array, hstack, zeros, linspace
from sys import exit, path
from os import name as os_name
import os
from ocelot.common.globals import *
from ocelot.cpbd.beam import Particle
from ocelot.lib.genera.src.python.trajectory.motion import Motion
und_type = "undulator"
import socket
from sys import platform
"""
if os_name == "nt":
    tail_u = r"codes/genera/build/windows/undulator.dll"
else:
    if socket.gethostname() == "exflwgs03":
        tail_u = "codes/genera/build/exflwgs03/undulator.so"
    else:
        tail_u = "codes/genera/build/ubuntu/undulator.so"

    if platform == "darwin":
        tail_u = "codes/genera/build/mac/undulator.so"
"""
import ocelot
#print ocelot.__file__

import os
path_to_ocelot = os.path.dirname(ocelot.__file__)
#print path
#tail = "ocelot/lib/genera/build/genera_libs/undulator.so"
home_dir = path[0]
#index =  path[0].find("siberia2")
#pathToDll = path[0][:index]+ tail
pathToDll = path_to_ocelot + "/lib/genera/build/genera_libs/undulator.so"
#print pathToDll
try:
    cundul= CDLL(pathToDll)
except:
    #os.system("python "+ path_to_ocelot + "/lib/genera/src/cpp/compile.py")
    #execfile(path_to_ocelot + "/lib/genera/src/cpp/compile.py")
    import ocelot.lib.genera.src.cpp.compile
    os.chdir(home_dir)
    cundul= CDLL(pathToDll)


def errorsBlock(error):
    def printOfStars(phrases, error):
        print ('**********************************************')
        print ('*', phrases, error)
        print ('**********************************************')
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


def define_points_undul(undulator, accuracy):
    accuracy = int(accuracy)
    #print undulator.type, undulator.id,undulator.transfer_map.order
    if undulator.field_file == None:
        npoints_traj = int(undulator.nperiods*50+10)*accuracy
    else:
        if undulator.field_map.units == "mm":
            npoints_traj = int(((undulator.field_map.z_arr[-1] - undulator.field_map.z_arr[0])*4)+10)*accuracy
        else:
            npoints_traj = int(((undulator.field_map.z_arr[-1] - undulator.field_map.z_arr[0])*2000)+10)*accuracy
    return npoints_traj


def field_map_for_dll(undulator):
    fm = undulator.field_map
    if undulator.field_file !=None:
        fm.z_arr = fm.z_arr - fm.z_arr[0]
        if fm.units == "m":
            x_arr = fm.x_arr*1000.
            y_arr = fm.y_arr*1000.
            z_arr = fm.z_arr*1000.
        else:
            x_arr = fm.x_arr
            y_arr = fm.y_arr
            z_arr = fm.z_arr
        if len(fm.x_arr) and len(fm.y_arr) and len(fm.z_arr):
            x_arr_ex = []
            y_arr_ex = []
            z_arr_ex = []
            for z in z_arr:
                for y in y_arr:
                    for x in x_arr:
                        x_arr_ex.append(x)
                        y_arr_ex.append(y)
                        z_arr_ex.append(z)
            x_arr = x_arr_ex
            y_arr = y_arr_ex
            z_arr = z_arr_ex
        mag_field = hstack((x_arr, y_arr, z_arr, fm.Bx_arr, fm.By_arr, fm.Bz_arr))
        len_mf = len(z_arr)
        cols_mf = len(mag_field)/len_mf
        if cols_mf == 2:
            fm.type = "planar"
        elif cols_mf == 3:
            fm.type = "spiral"
        elif cols_mf == 6:
            fm.type = "3D"
        else:
            exit("unknown type of undulator, file includes wrong number of cols")
        # undul_param = [Bx = 0, By = 0, phase = 0, nperiods = 0,
        #                 lperiod -(used for "Planar undulator unly), ax -(used for "Planar undulator unly)]
        undul_param = [0.,0.,0, 0, undulator.lperiod*1000., undulator.ax*1000.]
        fm.l = fm.z_arr[-1] - fm.z_arr[0]
    else:
        fm.type = "analytic"
        len_mf = 0
        cols_mf = 0
        mag_field = array([0])

        undulator.By = undulator.Kx*m_e_eV*2.*pi/(undulator.lperiod*speed_of_light) # [T]
        undulator.Bx = undulator.Ky*m_e_eV*2.*pi/(undulator.lperiod*speed_of_light) # [T}
        undul_param = [undulator.Bx, undulator.By, undulator.phase, undulator.nperiods, undulator.lperiod*1000., undulator.ax*1000.]
        fm.l = undulator.lperiod*undulator.nperiods*1000.
    return mag_field, len_mf, cols_mf, undul_param # params in [T] and [mm]


def motion2particles(aMotion, particle0):
    list_particle = []
    npoints_traj = int(len(aMotion)/11)
    for i in range(npoints_traj):
        particle = Particle()
        particle.x = aMotion[i]/1000.
        particle.px = aMotion[3*npoints_traj + i]
        particle.y = aMotion[npoints_traj + i]/1000.
        particle.py = aMotion[4*npoints_traj + i]
        particle.s =  particle0.s + aMotion[2*npoints_traj + i]/1000.
        particle.p = particle0.p
        list_particle.append(particle)
    return list_particle


def und_trace(undulator, particle0, energy, bRough, n_trajectory_points, accuracy = 1):
    # energy in [GeV]
    N = define_points_undul(undulator, accuracy = accuracy)
    n_trajectory_points = int(N/3)*3-1
    fm = undulator.field_map
    dPtr = POINTER(c_double)
    magField, lenMF, colsMF, undul_param = field_map_for_dll(undulator)

    npoints_traj = n_trajectory_points*undulator.field_map.field_file_rep - (undulator.field_map.field_file_rep - 1)*(2-bRough)
    #print "npoints_traj = ", npoints_traj,undulator.field_map.field_file_rep
    motion = Motion(N = npoints_traj)
    aMotion = motion.memory_motion
    #print "len(aMotion) = ", len(aMotion), npoints_traj*11
    #aMotion = zeros(npoints_traj*11)
    init_cond = [particle0.x*1000,particle0.y*1000, particle0.s*0., particle0.px, particle0.py, energy*(1+particle0.p)/m_e_GeV]

    #print n_trajectory_points*11
    #print len(motion.memory_motion)
    error = cundul.trajectory(Py2C(magField),       # 1D magnetic field array [X,Y,Z,Bx,By,Bz] or [Z,By] or [Z,By,Bx]
        colsMF,                         # number of cols in magnetic data
        lenMF,                          # number of rows in magnetic data
        Py2C([undulator.dx,             # misaligmnent of undulator  [dx, dy, v_angle, h_angle, tilt]
              undulator.dy, 0., 0.,
              undulator.dtilt]),
        c_int(bRough),                  # if bRough = 1 then |  |  | else  | . . . (|) . . . | (gauss integration)
        Py2C(init_cond),                # [beam.x [mm], beam.y [mm], 0 [long. pos. mm], beam.xp [rad!], beam.yp[rad!], gamma]
        Py2C(undul_param),              # [Bx (amplitude [T]), By (amplitude [T]), phase (between Bx and By),
                                        # number_of_periods, lenPeriod [m], ax (width of undulator [m])]

        n_trajectory_points,            # number of points on trajectory on element
        c_int(fm.field_file_rep),# Number of repetition of field file. for example, undul = (end1, period*40, end2)

        aMotion.ctypes.data_as(dPtr)) # pointer to array of motion

    errorsBlock(error)
    #print "motion = ", motion.X[-2]
    return motion


def trajectory_through_undul(undulator, particle0, energy, bRough, n_trajectory_points):

    aMotion = und_trace(undulator, particle0, energy, bRough, n_trajectory_points)

    list_particle = motion2particles(aMotion.memory_motion, particle0)
    return list_particle

"""
def particle2c_array(particle):
    npart = len(list_particle)
    c_particle = zeros(7*npart)
"""


def da_undul_list(undulator, list_particle, bRough = 1):
    #print undulator.type
    # energy in [GeV]
    n_trajectory_points = define_points_undul(undulator, accuracy = 1)
    #print "npoint = ", n_trajectory_points
    fm = undulator.field_map
    dPtr = POINTER(c_double)
    magField, lenMF, colsMF, undul_param = field_map_for_dll(undulator)

    npart = len(list_particle)
    c_particle = zeros(7*npart)
    for i, part in enumerate(list_particle):
        #part = pxy.p
        c_particle[i*7 + 0] = part.x
        c_particle[i*7 + 1] = part.y
        c_particle[i*7 + 2] = part.px
        c_particle[i*7 + 3] = part.py
        c_particle[i*7 + 4] = part.p
        c_particle[i*7 + 5] = part.s
        c_particle[i*7 + 6] = part.tau

    gamma = undulator.transfer_map.energy/m_e_GeV
    #print "start undulator ", [undulator.dx,             # misaligmnent of undulator  [dx, dy, v_angle, h_angle, tilt]
    #                           undulator.dy, undulator.v_angle, undulator.h_angle,
    #                           undulator.dtilt]
    error = cundul.da_undulator(Py2C(magField),       # 1D magnetic field array [X,Y,Z,Bx,By,Bz] or [Z,By] or [Z,By,Bx]
                colsMF,                         # number of cols in magnetic data
                lenMF,                          # number of rows in magnetic data
                Py2C([undulator.dx,             # misaligmnent of undulator  [dx, dy, v_angle, h_angle, tilt]
                      undulator.dy, undulator.v_angle, undulator.h_angle,
                      undulator.dtilt]),
                c_int(bRough),                  # if bRough = 1 then |  |  | else  | . . . (|) . . . | (gauss integration)
                Py2C(undul_param),              # [Bx (amplitude [T]), By (amplitude [T]), phase (between Bx and By),
                # number_of_periods, lenPeriod [m], ax (width of undulator [m])]

                n_trajectory_points,            # number of points on trajectory on element
                c_int(fm.field_file_rep),# Number of repetition of field file. for example, undul = (end1, period*40, end2)

                c_double(gamma),
                c_particle.ctypes.data_as(dPtr),
                c_int(npart))
    #print "stop undulator "
    errorsBlock(error)
    #list_particle = []
    for i, particle in enumerate(list_particle):
        #particle = pxy.p
        particle.x = c_particle[i*7 + 0]
        particle.y = c_particle[i*7 + 1]
        particle.px = c_particle[i*7 + 2]
        particle.py = c_particle[i*7 + 3]
        particle.p = c_particle[i*7 + 4]
        particle.s += c_particle[i*7 + 5]
        particle.tau = c_particle[i*7 + 5]
        #list_particle.append(particle)
    return list_particle



def track_with_IDs(lattice,particle0, ndiv_lin = 1,ndiv_und = 1000):

    particle_list = [particle0]
    particle = particle0
    L = 0.
    for element in lattice.sequence:

        if element.type == "undulator":

            energy = element.transfer_map.energy
            if energy == 0:
                print ("energy = 0 GeV: undulator as matrix")
                particle = element.transfer_map*particle
                particle_list.append(particle)
                continue
            else:

                list_particle = trajectory_through_undul(element, particle, energy, bRough = 0, n_trajectory_points = ndiv_und)
                particle_list.extend(list_particle)
                particle= list_particle[-1]
                continue
        else:
            p0 = particle
            for i in range(ndiv_lin):
                particle = element.transfer_map((i+1)*element.l/ndiv_lin)*p0
                particle_list.append(particle)

    return particle_list



"""
        part_loc_list = []
        if element.field_map == None and element.type != "sextupole":
            for z in linspace(0, element.l, num = nsteps):
                particle = element.transfer_map(z)*part_element
                particle.s = z + L
                part_loc_list.append(particle)
        elif element.field_map == None and element.type == "sextupole":
            particle = element.transfer_map*part_element
            #particle.s = z + L
            part_loc_list.append(particle)
        else:
            part_loc_list = element.transfer_map*part_element
        L += element.l
        #print L, map(lambda p:p.s, part_loc_list)
        particle_list.extend(part_loc_list)
        part_element = particle_list[-1]
"""
