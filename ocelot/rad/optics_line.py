from copy import deepcopy
import logging
import numpy as np

from ocelot.common.globals import *
from ocelot.rad.transfer_function import *
from ocelot.rad.optics_elements import *


flatten = lambda *n: (e for a in n
                      for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))

class OpticsLine():
    def __init__(self, sequence, start=None, stop=None):
        self.sequence = list(flatten(sequence))
        self.stop = stop
        self.start = start
        self.update_optics_masks()

    def update_optics_masks(self):
        print("update mask")

        for element in self.sequence:
            print(element)

            get_transfer_function(element)

    def estimate_mesh(self):
        for element in self.sequence:
            element.mesh = 0
'''
def get_transfer_function(element):
    element.mask = Mask()

    if element.__class__ is None:
        raise ValueError('Optics element must belong to the OpticsElement class')
    
    elif element.__class__ is ApertureRect:
        mask = RectMask()
        mask.lx = element.lx
        mask.ly = element.ly
        element.mask = mask
        
    elif element.__class__ is ApertureEllips:
        mask = EllipsMask()
        mask.ax = element.ax
        mask.ay = element.ay
        element.mask = mask
        
    elif element.__class__ is FreeSpace:
        element.mask = DriftMask(element.l, element.mx, element.my)
    
    elif element.__class__ is ThinLens:
        element.mask = LensMask(element.fx, element.fy)
    
    elif element.__class__ is ImperfectMirror:
        mask = MirrorMask()
        mask.height_profile = element.height_profile
        mask.hrms = element.hrms
        mask.lx = element.lx
        mask.ly = element.ly
        element.mask = mask
'''