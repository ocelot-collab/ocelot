from copy import deepcopy
import numpy as np

from ocelot.common.globals import *
from ocelot.rad.transfer_function import *
from ocelot.rad.optics_elements import *

from ocelot.common.ocelog import *
_logger = logging.getLogger(__name__)

flatten = lambda *n: (e for a in n
                      for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))

class OpticsLine():
    def __init__(self, sequence, start=None, stop=None):
        self.sequence = list(flatten(sequence))
        self.stop = stop
        self.start = start
        
        _logger.info(ind_str + 'beamline sequence {}'.format(self.sequence))

        try:
            if start is not None:
                id1 = self.sequence.index(self.start)

            else:
                id1 = 0 
            if stop is not None:
                id2 = self.sequence.index(self.stop) + 1
                self.sequence = self.sequence[id1:id2]
            else: 
                self.sequence = self.sequence[id1:]
            _logger.info(ind_str + 'sequence the first element {} and its index {}'.format(self.start, id1))
            _logger.info(ind_str + 'sequence the last element {} and its index {}'.format(self.start, id2))
            _logger.info(ind_str + 'beamline sequence {}'.format(self.sequence))
        except:
            raise ValueError(ind_str + 'cannot construct sequence, element not found')
            _logger.error(ind_str + 'cannot construct sequence, element not found')

        self.update_optics_masks()
        

    def update_optics_masks(self):
        for element in self.sequence:
            print(element)

            get_transfer_function(element)

    def estimate_mesh(self):
        for element in self.sequence:
            element.mesh = 0

def get_transfer_function(element):
    """
    Matchs OpticalElement object to its mask. Several masks can correspond to one OpticalElement object 
    The function rewrites OpticalElement object parameters to Mask object parameters
    
    :param element: OpticsElement class
    """
    element.mask = Mask()

    if element.__class__ is None:
        raise ValueError('Optics element must belong to the OpticsElement class')

    elif element.__class__ is FreeSpace: #implementation of several masks
        if element.method in ['PropMask', 'PropMask_kf', 'PropMask_kt']:
            if element.mx == 1 and element.my == 1:
                element.mask = PropMask(element.l, element.method)                
            elif element.mx != 1 and element.my != 1:
                element.mask = Prop_mMask(element.l, element.mx, element.my, element.method)
            else:
                ValueError("mx and my must be positive non-zero values")
        elif element.type == 'Fraunhofer_Propagator':
            _logger.warn('Fraunhofer Propagation method has not implemented yet')
            pass    
        elif element.type == 'Fresnel_Propagator':
            _logger.warn('Fresnel Propagation method has not implemented yet')          
            pass                               
        else:            
            raise ValueError('Propagator method can be PropMask (see Tutorials), Fraunhofer_Propagator or Fresnel_Propagator')    
        
    elif element.__class__ is ThinLens: 
        element.mask = LensMask(element.fx, element.fy)
  
    elif element.__class__ is ApertureRect: 
        element.mask = ApertureRectMask(element.lx, element.ly, element.cx, element.cy)
        
    elif element.__class__ is ApertureEllips: #no checked
        element.mask = ApertureEllipsMask(element.ax, element.ay, element.cx, element.cy)
        
    elif element.__class__ is DispersiveSection: #no checked
        element.mask = PhaseDelayMask(element.coeff, element.E_ph0)
        
    elif element.__class__ is ImperfectMirrorSurface: #no checked
        element.mask = MirrorMask(element.height_profile, element.hrms, element.angle, element.plane, element.lx, element.ly)

    else:
        raise ValueError('Optics element must belong to one of the child OpticsElement classes')
