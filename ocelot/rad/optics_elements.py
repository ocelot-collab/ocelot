import numpy as np

from ocelot.common.ocelog import *
_logger = logging.getLogger(__name__)


class OpticsElement:
    """
    Parent class Optics element
    :param eid: element id, (for example 'KB')  
    """
    def __init__(self, eid=None):
        self.eid = eid
        
#    def apply(self, dfl): #is this method need?
#        """
#        TODO
#        write documentation
#        """
#        get_transfer_function(self)
#        self.mask.apply(self.mask, dfl)
           
class FreeSpace(OpticsElement):
    """
    Class for Drift element 
    :param OpticsElement(): optics element parent class with following parameters
        :param eid: element id, (for example 'KB')  
    :param l: propagation distance, [m]
    :param mx: is the output x mesh size in terms of input mesh size (mx = Lx_out/Lx_inp)
    :param my: is the output y mesh size in terms of input mesh size (my = Ly_out/Ly_inp)
    """
    def __init__(self, l=0., mx=1, my=1, method='PropMask_kf', eid=None):
        OpticsElement.__init__(self, eid=eid)
        self.l = l
        self.mx = mx
        self.my = my
        self.method = method  #method of propagation, also may be 
                                # 'Fraunhofer_Propagator'
                                # 'Fresnel_Propagator' . . .

class ThinLens(OpticsElement):
    """
    Class for Lens element
    :param OpticsElement(): optics element parent class with following parameters
        :param eid: element id, (for example 'KB') 
    :param fx: focus length in x direction, [m]
    :param fy: focus length in y direction, [m]
    """
    def __init__(self, fx=np.inf, fy=np.inf, eid=None):
        OpticsElement.__init__(self, eid=eid)
        self.fx = fx
        self.fy = fy
        
class Aperture(OpticsElement):
    """
    Aperture
    :param eid: id of the optical element
    """
    def __init__(self, eid=None):
        OpticsElement.__init__(self, eid=eid)
        
class ApertureRect(Aperture):
    """
    Rectangular aperture
    :param Aperture(): optics element parent class with following parameters
        :param eid: element id, (for example 'KB', 'Aperture') 
    :param lx: horizonatal size of the aperture, [m]
    :param ly: vertical size of the aperture, [m]
    :param cx: x coprdinate coordinate of the aperture center, [m] 
    :param cy: y coprdinate coordinate of the aperture center, [m] 
    :param eid: id of the optical element
    """

    def __init__(self, lx=np.inf, ly=np.inf, cx=0., cy=0., eid=None):
        Aperture.__init__(self, eid=eid)
        self.lx = lx
        self.ly = ly
        self.cx = cx
        self.cy = cy

class ApertureEllips(Aperture):
    """
    Elliptical Aperture
    :param ax: ellipse x main axis, [m]
    :param ay: ellipse y main axis, [m]
    :param cx: ellipse x coordinate of center, [m] 
    :param cy: ellipse y coordinate of center, [m] 
    :param eid: id of the optical element
    """

    def __init__(self, ax=np.inf, ay=np.inf, cx=0., cy=0., eid=None):
        Aperture.__init__(self, eid=eid)
        self.ax = ax
        self.ay = ay
        self.cx = cx
        self.cy = cy
             
      
class DispersiveSection(OpticsElement):
    """
    Dispersive Section
        :param coeff:  
            coeff[0] =: measured in [rad]      --- phase
            coeff[1] =: measured in [fm s ^ 1] --- group delay
            coeff[2] =: measured in [fm s ^ 2] --- group delay dispersion (GDD)
            coeff[3] =: measured in [fm s ^ 3] --- third-order dispersion (TOD)
            ...
        :param E_ph0: energy with respect to which the phase shift is calculated
    """  
    def __init__(self, coeff=[0], E_ph0=None, eid=None):
        OpticsElement.__init__(self, eid=None)
        self.coeff = coeff
        self.E_ph0 = E_ph0


class ImperfectMirrorSurface(OpticsElement):
    """
    TODO
    write documentation
    """
    def __init__(self, height_profile=None, hrms=0, lx=np.inf, ly=np.inf, angle=np.pi * 2 / 180, plane='x', eid=None):
        OpticsElement.__init__(self, eid=eid)
        self.height_profile = height_profile
        self.hrms = hrms
        self.lx = lx
        self.ly = ly
        self.angle=angle
        self.plane=plane
