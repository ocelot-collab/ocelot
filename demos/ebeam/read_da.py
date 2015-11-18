import numpy as np
from ocelot.gui import *

da = np.loadtxt("da.txt")
xy = np.loadtxt("da_axis.txt")
x_array = xy[:,0]
y_array = xy[:,1]
show_da(da, x_array, y_array)
