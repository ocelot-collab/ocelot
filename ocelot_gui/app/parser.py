"""Lattice parser functions"""

import re
import sys

from ocelot import *


flatten = lambda *n: (e for a in n
                      for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))


def convert_elem_id(name):
    
    new_name = name
    new_name = new_name.replace('.', '_')
    new_name = new_name.replace(':', '_')
    new_name = new_name.replace('-', '_')

    return new_name


class LatticeRead():

    def __init__(self):
        pass
    

    def parsing_elements(self, cell):
        
        sequence = list(flatten(cell))
        
        elements = {}
        for elem in sequence:
            elem_id = convert_elem_id(elem.id)
            elements[elem_id] = elem
        
        return elements


class LatticeSave():

    def __init__(self):
        self.param_list = {}
        self.param_list['Drift'] = ('l',)
        self.param_list['Quadrupole'] = ('l', 'k1')
        self.param_list['Bend'] = ('l', 'angle')
        self.param_list['SBend'] = ('l', 'angle')
        self.param_list['RBend'] = ('l', 'angle')
        self.param_list['Edge'] = ()
        self.param_list['Hcor'] = ('l', 'angle')
        self.param_list['Vcor'] = ('l', 'angle', 'tilt') # 'tilt' is used only for fix
        self.param_list['Marker'] = ()
        self.param_list['Monitor'] = ()
        self.param_list['Sextupole'] = ('l', 'k2')
        self.param_list['Octupole'] = ('l', 'k3')
        self.param_list['Undulator'] = ('lperiod', 'nperiods', 'Kx', 'Ky', 'l', 'ax') # 'l' and 'ax' are used only for fix
        self.param_list['Cavity'] = ('l', 'v', 'freq', 'phi')
        self.param_list['TDCavity'] = ('l', 'v', 'freq', 'phi')
        self.param_list['Solenoid'] = ('l', 'k')
        self.param_list['Multipole'] = ('kn',)
        self.param_list['Matrix'] = ('l',)
        self.param_list['UnknownElement'] = ('l',)


    def parsing(self, gui_lattice, split=False):

        self.gui_lattice = gui_lattice
        lines = 'from ocelot import *\n\n'

        if self.gui_lattice.tws0 and isinstance(self.gui_lattice.tws0, Twiss):
            lines += self.parsing_twiss(self.gui_lattice.tws0)

        if self.gui_lattice.beam and isinstance(self.gui_lattice.beam, Beam):
            lines += self.parsing_beam(self.gui_lattice.beam)

        if self.gui_lattice.elements and isinstance(self.gui_lattice.elements, dict):
            lines += self.parsing_elements(self.gui_lattice.elements)

        if self.gui_lattice.cell and isinstance(self.gui_lattice.cell, (list, tuple, dict)):
            lines += self.parsing_cell(self.gui_lattice.cell, split=split)

        return lines


    def parsing_elements(self, elements):

        elements_array = {}
        
        # sorting elements into array by type
        for key in elements:
            element_type = elements[key].__class__.__name__
            if element_type not in elements_array.keys():
                elements_array[element_type] = {}
            
            elements_array[element_type][key] = elements[key]

        # preparing elements list
        lines = ''
        for element_type in elements_array:
            
            for name in elements_array[element_type]:

                element = elements_array[element_type][name]

                lines += str(name) + ' = ' + str(element_type) + '('
                
                # prepare parameters from self.param_list
                for param in self.param_list[element_type]:
                    if element_type in ['Vcor'] and param == 'tilt':
                        pass
                    elif element_type in ['Undulator'] and (param == 'l' or param == 'ax'):
                        pass
                    else:
                        lines += str(param) + '=' + str(element.__dict__[param]) + ', '

                # prepare other parameters
                for param in element.__dict__:
                    if param not in self.param_list[element_type] and \
                       param != 'id' and \
                       isinstance(element.__dict__[param], (int, float)) and \
                       element.__dict__[param] != 0.0 and \
                       element.__dict__[param] != 0:
                        lines += str(param) + '=' + str(element.__dict__[param]) + ', '

                # fix for id - eid parameter
                lines += 'eid="' + str(element.__dict__['id']) + '", '

                lines = lines[0:-2]
                lines += ')\n'
        
            lines += '\n'

        return lines


    def parsing_twiss(self, tws):

        lines = 'tws0 = Twiss()\n'
        for param in tws.__dict__:
            if tws.__dict__[param] != 0.0 and \
                tws.__dict__[param] != 0 and \
                tws.__dict__[param] != '':

                lines += 'tws0.' + str(param) + ' = ' + str(tws.__dict__[param]) + '\n'

        lines += '\n'
        return lines


    def parsing_beam(self, beam):

        lines = 'beam = Beam()\n'
        for param in beam.__dict__:
            if beam.__dict__[param] != 0.0 and \
                beam.__dict__[param] != 0 and \
                beam.__dict__[param] != '' and \
                param != 'shape':

                lines += 'beam.' + str(param) + ' = ' + str(beam.__dict__[param]) + '\n'

        lines += '\n'
        return lines

    
    def parsing_cell(self, cell, split=False):
        lines = 'cell = ('
        
        sequence = list(flatten(cell))

        for i, elem in enumerate(sequence, start=1):
            
            elem_id = convert_elem_id(elem.id)
            
            lines += str(elem_id) + ', '
            
            if split and i % 10 == 0:
                lines += '\n'

        if len(sequence) > 0:
            lines = lines[0:-2]
        lines += ')\n'
        
        lines += '\n'
        return lines
    