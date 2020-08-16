__author__ = 'Sergey Tomin'

import numpy as np

def length_array(xyz_array):
    inds = xyz_array[xyz_array[:-1]>xyz_array[1:]]
    if len(inds)>0:
        nxyz = np.argwhere(xyz_array == inds[0])[0][0]+1
    else:
        nxyz = len(xyz_array)
    return nxyz


def read_flat_file(field_file, sCom = "#"):
    # File must have 3 cols, else will be error
    f = open(field_file, 'r')
    f.readline() #1st line: just pass
    x_start = float(f.readline().split(sCom, 2)[1]) #2nd line: initial X position [m]; it will not actually be used
    x_step = float(f.readline().split(sCom, 2)[1]) #3rd line: step vs X [m]
    x_np = int(f.readline().split(sCom, 2)[1]) #4th line: number of points vs X
    y_start = float(f.readline().split(sCom, 2)[1]) #5th line: initial Y position [m]; it will not actually be used
    y_step = float(f.readline().split(sCom, 2)[1]) #6th line: step vs Y [m]
    y_np = int(f.readline().split(sCom, 2)[1]) #7th line: number of points vs Y
    z_start = float(f.readline().split(sCom, 2)[1]) #8th line: initial Z position [m]; it will not actually be used
    z_step = float(f.readline().split(sCom, 2)[1]) #9th line: step vs Z [m]
    z_np = int(f.readline().split(sCom, 2)[1]) #10th line: number of points vs Z
    tot_np = x_np*y_np*z_np

    Bx_array = np.zeros([tot_np])
    By_array = np.zeros([tot_np])
    Bz_array = np.zeros([tot_np])

    for i in range(tot_np):
        curLineParts = f.readline().split('\t')
        Bx_array[i] = float(curLineParts[0])
        By_array[i] = float(curLineParts[1])
        Bz_array[i] = float(curLineParts[2])

    f.close()

    x_array = np.linspace(x_start, x_start + x_step*x_np, num=x_np, endpoint=False)
    y_array = np.linspace(y_start, y_start + y_step*y_np, num=y_np, endpoint=False)
    z_array = np.linspace(z_start, z_start + z_step*z_np, num=z_np, endpoint=False)
    return x_array, y_array, z_array, Bx_array, By_array, Bz_array


def read_tabular_file(field_file):
    try:
        field_data = np.loadtxt(field_file, delimiter=' ', unpack = True)
    except:
        try:
            print( "read_map: try to use delimiter = ','")
            field_data = np.loadtxt(field_file, delimiter=',', unpack = True)
        except:
            field_data = np.loadtxt(field_file, unpack=True)

    ncols = np.shape(field_data)[0]

    Bx_array = np.array([])
    By_array = np.array([])
    Bz_array = np.array([])
    x_array = np.array([])
    y_array = np.array([])
    z_array = np.array([])

    if ncols == 2:
        # planar undulator or dipol magnet/phase shifter
        # So there is only one possibility - (z_array, By_array)
        z_array = field_data[0, :]
        By_array = field_data[1, :]
    elif ncols == 3:
    # spiral undulator or 2D map
    # So, (z_array, Bx_array, By_array) or  (x_array, z_array, By_array)
    # BUT we start with spiral undulator.
        z_array = field_data[0, :]
        Bx_array = field_data[1, :]
        By_array = field_data[2, :]

    elif ncols == 6:
        x_array = field_data[0, :]
        y_array = field_data[1, :]
        z_array = field_data[2, :]
        Bx_array = field_data[3, :]
        By_array = field_data[4, :]
        Bz_array = field_data[5, :]
        nx_t = length_array(x_array)
        ny_t = length_array(y_array)
        nz_t = length_array(z_array)
        if nz_t > ny_t and ny_t > nx_t:
            x_array = x_array[:nx_t]
            y_array = y_array[:ny_t:nx_t]
            ny = len(y_array)
            z_array = z_array[:nz_t:ny*nx_t]
        else:
            print("wrong coordinates order in the field file (field_map.py)")
    return x_array, y_array, z_array, Bx_array, By_array, Bz_array


class FieldMap:
    def __init__(self, field_file, xyz_units="mm", format="flat"):
        self.field_file = field_file
        self.format = format
        self.units = xyz_units
        self.field_file_rep = 1
        self.Bx_arr = np.array([])
        self.By_arr = np.array([])
        self.Bz_arr = np.array([])
        self.x_arr = np.array([])
        self.y_arr = np.array([])
        self.z_arr = np.array([])
        if self.field_file != None:
            try:
                self.format = "flat"
                self.read(field_file)
            except :
                self.format = "tabular"
                self.read(field_file)
        else:
            pass

    def read(self, field_file):
        if self.format == "flat":
            self.x_arr, self.y_arr, self.z_arr, self.Bx_arr, self.By_arr, self.Bz_arr = read_flat_file(field_file)

        elif self.format == "tabular":
            self.x_arr, self.y_arr, self.z_arr, self.Bx_arr, self.By_arr, self.Bz_arr = read_tabular_file(field_file)
        self.l = self.z_arr[-1] - self.z_arr[0]
            #return SRWLMagFld3D(locArBx, locArBy, locArBz, xNp, yNp, zNp, xRange, yRange, zRange, 1)

"""

class FieldMap:
    def __init__(self, filename, lperiod, nperiods, Kx , Ky , phase ,ax):
        self.ax = ax                            # width of undulator, when ax is negative undulator width is infinite
        self.field_file_rep = 1    # Number of repetition of field file.
        self.lperiod = lperiod    #
        self.nperiods = nperiods
        self.type = 0
        self.field_map_array = array([0])
        self.len_field_map = 0
        self.ncols_field_map = 0
        self.l = lperiod*nperiods
        if filename != None:
            self.read(filename)
            self.define_type()
        else:
            self.type = "analytic"
            By = Kx*m_e_MeV*1.e+6*2.*pi/(self.lperiod*speed_of_light) # [T]
            Bx = Ky*m_e_MeV*1.e+6*2.*pi/(self.lperiod*speed_of_light) # [T}
            self.undul_param = array([Bx, By, phase, self.nperiods, self.lperiod*1000, self.ax*1000])
            # [Bx = 0, By = 0, phase = 0, nperiods = 0, lperiod -(used for "Planar undulator unly), ax -(used for "Planar undulator unly)]


    def read(self, filename):
        # what do you think about this? I think format of field map will be standard for ours codes,
        # therefore we can standardize this function
        try:
            self.field_data = loadtxt(filename, delimiter=' ', unpack = True)
        except:
            print "read_map: try to use delimiter = ','"
            #print self.field_file_path
            self.field_data = loadtxt(filename, delimiter=',', unpack = True)

    def write(self, filename):
        pass

    def B(self, x, y, z):
        pass

    def define_type(self):
        # get magnetic field
        self.ncols_field_map = shape(self.field_data)[0]
        if self.ncols_field_map == 2:
            self.type = "planar"
            self.Z_arr = self.field_data[0,:]
            self.By_arr = self.field_data[1,:]
            fld_XYZ = self.Z_arr
            fld_B = self.By_arr
        elif self.ncols_field_map == 3:
            self.type = "spiral"
            self.Z_arr = self.field_data[0,:]
            self.By_arr = self.field_data[1,:]
            self.Bx_arr = self.field_data[2,:]
            fld_XYZ = self.Z_arr
            fld_B = hstack((self.By_arr, self.Bx_arr))
        elif self.ncols_field_map == 6:
            self.type = "3D"
            self.X_arr = self.field_data[0,:]
            self.Y_arr = self.field_data[1,:]
            self.Z_arr = self.field_data[2,:]
            self.Bx_arr = self.field_data[3,:]
            self.By_arr = self.field_data[4,:]
            self.Bz_arr = self.field_data[5,:]
            fld_XYZ = hstack((self.X_arr, self.Y_arr, self.Z_arr))
            fld_B = hstack((self.Bx_arr, self.By_arr, self.Bz_arr))
        else:
            exit("unknown type of undulator, file includes wrong number of cols")
        # undul_param = [Bx = 0, By = 0, phase = 0, nperiods = 0, lperiod -(used for "Planar undulator unly), ax -(used for "Planar undulator unly)]
        self.undul_param = [0.,0.,0, 0, self.lperiod, self.ax]
        self.field_map_array = hstack((fld_XYZ, fld_B))
        self.len_field_map = len(self.Z_arr)
        self.l = self.Z_arr[-1] - self.Z_arr[0]

"""
