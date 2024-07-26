"""
Adaptors between WPG and Ocelot
Utilizes WPG branch: https://github.com/twguest/WPG

author: twguest (trey.guest@xfel.eu)
"""

from scipy.constants import c, h, e
from wpg.srw import srwlpy
from wpg.wavefront import Wavefront
import numpy as np

def complex_to_wpg(arr):
    """
    Converts a complex wavefield array into a WPG-style electric field array.
    
    The WPG-style array separates the real and imaginary parts of the input 
    complex wavefield.

    Parameters
    ----------
    arr : np.ndarray
        Complex wavefield array of shape (x, y, t) and dtype complex128.

    Returns
    -------
    new_arr : np.ndarray
        WPG-style electric field array of shape (nx, ny, nz, 2) and dtype float64.
        The last dimension separates the real and imaginary components.
    """
    new_arr = np.zeros([arr.shape[0], arr.shape[1], arr.shape[2], 2])
    new_arr[:,:,:,0] = arr.real
    new_arr[:,:,:,1] = arr.imag
    return new_arr

def dfl2wpg(dfl):
    """
    Converts an Ocelot dfl to a WPG wavefront.
    
    The function converts the Ocelot dfl into a real-space, frequency-domain 
    representation suitable for use with WPG.

    Parameters
    ----------
    dfl : ocelot.rad.RadiationField
        Ocelot RadiationField object containing the field distribution.

    Returns
    -------
    wfr : wpg.wavefront.Wavefront
        WPG Wavefront object with the converted field distribution.
    """ 

    ny, nx, nt = dfl.fld.T.shape

    wfr = Wavefront()
        
    # Setup E-field.
    wfr.data.arrEhor = np.zeros(shape=(nx, ny, nt, 2))
    wfr.data.arrEver = np.zeros(shape=(nx, ny, nt, 2))

    wfr.params.wEFieldUnit = 'sqrt(W/mm^2)'
    wfr.params.photonEnergy = (h * c) / (dfl.xlamds * e)
    
    wfr.params.Mesh.nSlices = nt
    wfr.params.Mesh.nx = nx
    wfr.params.Mesh.ny = ny      
    
    wfr.params.Mesh.sliceMin = -nt / 2 * dfl.dz
    wfr.params.Mesh.sliceMax = +nt / 2 * dfl.dz
    
    if dfl.domain_z == 't':
        wfr.params.wDomain = 'time'
    elif dfl.domain_z == 'f': 
        wfr.params.wDomain = 'frequency'
    
    if dfl.domain_xy == 's':
        wfr.params.wSpace = 'R-space'
    elif dfl.domain_xy == 'k':
        wfr.params.wSpace = 'Q-space'

    wfr.params.Mesh.xMin = -nx / 2 * dfl.dx
    wfr.params.Mesh.xMax = +nx / 2 * dfl.dx
    wfr.params.Mesh.yMin = -ny / 2 * dfl.dy
    wfr.params.Mesh.yMax = +ny / 2 * dfl.dy

    wfr.params.Rx = 1
    wfr.params.Ry = 1

    wfr.data.arrEhor = complex_to_wpg(dfl.fld.T)
    
    return wfr

def wpg2dfl(wfr):
    """
    Converts a WPG wavefront to an Ocelot dfl.
    
    The function ensures the output is in a real-space, frequency-domain 
    representation.

    Parameters
    ----------
    wfr : wpg.wavefront.Wavefront
        WPG Wavefront object containing the electric field distribution.

    Returns
    -------
    dfl : ocelot.rad.RadiationField
        Ocelot RadiationField object with the converted field distribution.
    """
    
    wfr.set_electric_field_representation('f')
    wfr.set_electric_field_representation('c')
    
    dfl = RadiationField(shape=wfr.get_efield().T.shape)
    
    dfl.fld = wfr.get_efield().T
    
    dfl.dx = wfr.get_pixel_size('x')
    dfl.dy = wfr.get_pixel_size('y')
    dfl.dz = wfr.get_pixel_size('f')
    
    dfl.xlamds = (h * c) / (e * wfr.params.photonEnergy)
    dfl.domain_z = 'f'
    dfl.domain_xy = 's'
    
    return dfl
