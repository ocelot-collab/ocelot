"""Lattice parser functions"""

import re
import sys

from ocelot import *
from ocelot.cpbd.io import *


flatten = lambda *n: (e for a in n for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))


class Parser():

    def __init__(self):
        pass
    

    def convert_elem_id(self, name):
    
        new_name = name
        new_name = new_name.replace('.', '_')
        new_name = new_name.replace(':', '_')
        new_name = new_name.replace('-', '_')

        return new_name


    def get_elements(self, cell):
        '''Create and return elements list from the cell'''

        sequence = list(flatten(cell))
        
        elements = {}
        for elem in sequence:
            elem_id = self.convert_elem_id(elem.id)
            elements[elem_id] = elem
        
        return elements


    def gui_lattice2input(self, gui_lattice, split=False):
        '''Prepare information about lattice to save into the python file'''

        lines = ['from ocelot import *\n']
        
        if gui_lattice.tws0 and isinstance(gui_lattice.tws0, Twiss):
            lines.append('\n#Initial Twiss parameters\n')
            lines.extend(twiss2input(gui_lattice.tws0))
        
        if gui_lattice.beam and isinstance(gui_lattice.beam, Beam):
            lines.append('\n#Beam parameters\n')
            lines.extend(beam2input(gui_lattice.beam))

        if gui_lattice.elements and isinstance(gui_lattice.elements, dict):
            lines.append('\n')
            lines.extend(elements2input(gui_lattice.lattice))

        if gui_lattice.cell and isinstance(gui_lattice.cell, (list, tuple, dict)):
            lines.append('\n# Lattice \n')
            lines.extend(cell2input(gui_lattice.lattice, True))

        lines.append('\n')

        return lines
