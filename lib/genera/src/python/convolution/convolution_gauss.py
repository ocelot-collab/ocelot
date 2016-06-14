from ctypes import CDLL, c_double, c_int, POINTER
from numpy import ones, shape
import os
#import convolution_fft
#from convolution_py_fast import convolution_fast
"""
The C++ code completely copies approach of the py- programm.
1. analytic integration of gauss G with function F ~ Intergation(exp(-x**2/sx^2-y**2/sy**2)*(ax+by+cxy+d), dx,dy
(see function integ_cell_2D )
where F we can describe approximately as  F ~ (ax+by+cxy+d) if we know 4 points and 4 value in these points.
To close this theme I put mathematics file for understanding of this approach.
"""
from sys import path
from os import name as os_name
from sys import platform
"""
from os import environ
if os_name == "nt":
    pathToDll = path[0] + r"/convolution"
else:
    pathToDll = path[0] [:-7]+ "/convolution/"

environ['PATH'] = environ['PATH'] + ";"+pathToDll
"""


tail = "/lib/genera/build/genera_libs/convolution.so"
home_dir = path[0]
#index =  path[0].find("siberia2")
import ocelot
#print ocelot.__file__
import os
path_to_ocelot = os.path.dirname(ocelot.__file__)
pathToDll = path_to_ocelot+ tail
try:
    my_conv= CDLL(pathToDll)
except:
    exec(open(path_to_ocelot+  "/lib/genera/src/cpp/compile.py"))
    #exec(open("ocelot/lib/codes/genera/src/cpp/compile.py"))
    os.chdir(home_dir)
    my_conv= CDLL(pathToDll)

def Py2C(array):
    arr_type =  c_double*len(array)
    c_array = arr_type(*array)
    return c_array

def convolution_2D_cpp(screen,X,Y, nx_add, ny_add, sx,sy ):


    #my_conv = CDLL("convolution.dll")

    dPtr = POINTER(c_double)
    nx = len(X)
    ny = len(Y)
    new_screen = ones(shape(screen))*screen
    pscreen = new_screen.ctypes.data_as(dPtr)
    ret = my_conv.conv_2D(pscreen, Py2C(X), Py2C(Y), c_int(nx), c_int(ny),  c_int(nx_add), c_int(ny_add), c_double(sx), c_double(sy))
    return new_screen

def convolution_1D_cpp(screen,X,nx_add, sx):


    #my_conv = CDLL("convolution.dll")

    dPtr = POINTER(c_double)
    nx = len(X)
    new_screen = ones(shape(screen))*screen
    pscreen = new_screen.ctypes.data_as(dPtr)
    ret = my_conv.conv_1D(pscreen, Py2C(X), c_int(nx),c_int(nx_add),  c_double(sx))
    return new_screen