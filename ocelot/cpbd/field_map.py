__author__ = 'Sergey Tomin'

import numpy as np
from scipy.interpolate import splrep, splev

def length_array(xyz_array):
    inds = xyz_array[xyz_array[:-1]>xyz_array[1:]]
    if len(inds)>0:
        nxyz = np.argwhere(xyz_array == inds[0])[0][0]+1
    else:
        nxyz = len(xyz_array)
    return nxyz


def read_flat_file(field_file, sCom = "#"):
    # File must have 3 cols, else will be error
    with open(field_file, 'r') as f:
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
        # Normalize tabular rows to the flat-file convention: X fastest.
        order = np.lexsort((x_array, y_array, z_array))
        points = np.column_stack((x_array, y_array, z_array))
        x_array, y_array, z_array = (np.unique(axis) for axis in (x_array, y_array, z_array))
        expected = len(x_array) * len(y_array) * len(z_array)
        if len(points) != expected or len(np.unique(points, axis=0)) != expected:
            raise ValueError("Field map coordinates must form a complete rectangular grid")
        Bx_array, By_array, Bz_array = (component[order] for component in (Bx_array, By_array, Bz_array))
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


def field_map2field_func(z, By):
    tck = splrep(z, By, k=3)

    def func(x, y, z):
        return 0., splev(z, tck, der=0), 0.

    return func


def field_map_to_field_func(field_map):
    """Return all field components in local Cartesian coordinates (metres).

    Tabular line maps are interpolated along Z; rectangular 3D maps use their
    X/Y/Z grid. Z=0 denotes the start of the map. Missing components are zero.
    """
    unit_coef = 0.001 if field_map.units == "mm" else 1.
    z_array = (field_map.z_arr - field_map.z_arr[0]) * unit_coef
    components = (field_map.Bx_arr, field_map.By_arr, field_map.Bz_arr)
    nx, ny, nz = max(1, len(field_map.x_arr)), max(1, len(field_map.y_arr)), len(z_array)
    if nz < 2 or np.any(np.diff(z_array) <= 0):
        raise ValueError("Field map requires at least two strictly increasing Z coordinates")
    if nx == ny == 1:
        splines = [splrep(z_array, values, k=min(3, nz - 1)) if len(values) else None
                   for values in components]

        def interpolate(x, y, z):
            return tuple(splev(z, spline) if spline is not None else np.zeros_like(z)
                         for spline in splines)
    else:
        from scipy.interpolate import RegularGridInterpolator

        x_array = field_map.x_arr * unit_coef
        y_array = field_map.y_arr * unit_coef
        values = np.column_stack([v if len(v) else np.zeros(nx * ny * nz) for v in components])
        # Flat files enumerate X fastest, then Y, then Z.
        values = values.reshape(nz, ny, nx, 3).transpose(2, 1, 0, 3)
        grid = RegularGridInterpolator((x_array, y_array, z_array), values,
                                       bounds_error=False, fill_value=None)

        def interpolate(x, y, z):
            points = np.stack(np.broadcast_arrays(x, y, z), axis=-1)
            values = grid(points).reshape(points.shape[:-1] + (3,))
            return tuple(values[..., i] for i in range(3))

    def field(x, y, z):
        if field_map.field_file_rep > 1:
            z = np.where(np.asarray(z) == z_array[-1] * field_map.field_file_rep,
                         z_array[-1], np.mod(z, z_array[-1]))
        return interpolate(x, y, z)
    return field

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
