"""Lattice parser functions"""

import re
import sys

from ocelot import *


class LatticeRead():
    
    def __init__(self):
        pass


    def parsing(self, lines):

        proc_lines = self.del_comments(lines)
        
        tws = self.parsing_twiss(proc_lines)
        beam = self.parsing_beam(proc_lines)
        elements = self.parsing_elements(proc_lines)
        cells, cells_order = self.parsing_cells(proc_lines, elements)

        return tws, beam, elements, cells, cells_order

    
    def del_comments(self, lines):

        # remove """ and ''' in one line comments
        lines = re.sub(r'\"{3}(.*)?\"{3}', '', lines)
        lines = re.sub(r'\'{3}(.*)?\'{3}', '', lines)
        
        proc_lines0 = []
        lines = lines.split("\n")
        for line in lines:
            proc_lines0.append(line.strip())

        # remove """ multiline comments
        proc_lines1 = []
        in_comment = False
        for line in proc_lines0:
            result = re.findall(r'(.*)\"{3}(.*)', line)
            if len(result) == 0:
                if not in_comment:
                    proc_lines1.append(line)
            else:
                if not in_comment:
                    proc_lines1.append(result[0][0])
                    in_comment = True
                else:
                    proc_lines1.append(result[0][1])
                    in_comment = False

        # remove ''' multiline comments
        proc_lines2 = []
        in_comment = False
        for line in proc_lines1:
            result = re.findall(r'(.*)\'{3}(.*)', line)
            if len(result) == 0:
                if not in_comment:
                    proc_lines2.append(line)
            else:
                if not in_comment:
                    proc_lines2.append(result[0][0])
                    in_comment = True
                else:
                    proc_lines2.append(result[0][1])
                    in_comment = False
        
        # remove # comments, empty lines, glue lines
        tmp = ''
        proc_lines3 = []
        for line in proc_lines2:

            # skip empty and comment line
            if line == '': continue
            if line[0] == '#': continue

            # glue lines
            if tmp != '':
                line = tmp + line
                tmp = ''
            if line[-1] == '\\':
                tmp = line[:-1]
                continue
            
            proc_lines3.append(line)
        
        return proc_lines3


    def parsing_twiss(self, lines):
        """parsing twiss parameters (only the first Twiss() defenition)"""

        tws = Twiss()
        tws_name = ''
        for line in lines:
            if tws_name == '':
                # find Twiss() variable name
                result = re.findall(r'(^[a-zA-Z]{1}[a-zA-Z0-9_]*)\s*=\s*Twiss\s*\(\s*\)', line)
                if len(result) != 0:
                    tws_name = result[0]
            else:
                # find twiss values
                result = re.findall(r'^' + str(tws_name) + r'\.([a-zA-Z]{1}[a-zA-Z0-9_]*)\s*=\s*([\+\-]?\d+\.?\d*)', line)
                if len(result) != 0 and hasattr(tws, result[0][0]):
                    tws.__dict__[result[0][0]] = float(result[0][1])
        
        return tws


    def parsing_beam(self, lines):
        """parsing beam parameters (only the first Beam() defenition)"""

        beam = Beam()
        beam_name = ''
        for line in lines:
            if beam_name == '':
                # find Beam() variable name
                result = re.findall(r'(^[a-zA-Z]{1}[a-zA-Z0-9_]*)\s*=\s*Beam\s*\(\s*\)', line)
                if len(result) != 0:
                    beam_name = result[0]
            else:
                # find beam values
                result = re.findall(r'^' + str(beam_name) + r'\.([a-zA-Z]{1}[a-zA-Z0-9_]*)\s*=\s*([\+\-]?\d+\.?\d*)', line)
                if len(result) != 0 and hasattr(beam, result[0][0]):
                    beam.__dict__[result[0][0]] = float(result[0][1])

        return beam


    def parsing_elements(self, lines):
        
        elements = {}
        for line in lines:
            result = re.findall(r'(^[a-zA-Z]{1}[a-zA-Z0-9_]*)\s*=\s*([a-zA-Z]+)\s*\((.*)\)', line)
            if len(result) != 0 and issubclass(getattr(sys.modules[__name__], result[0][1]), Element):
                elements[result[0][0]] = eval(result[0][1] + '(' + result[0][2] + ')')

        return elements
    

    def parsing_cells(self, lines, elements):
        """parsing sequences of elements"""

        cells = {}
        cells_order = []
        for line in lines:
            # find sequences of elements
            result = re.findall(r'(^[a-zA-Z]{1}[a-zA-Z0-9_]*)\s*=\s*\(([a-zA-Z0-9_\,\s]+)\)$', line)
            if len(result) != 0:

                cells[result[0][0]] = ()
                cells_order.append(result[0][0])

                sequence = result[0][1].split(",")
                for elem in sequence:
                    elem = elem.strip()
                    
                    if elem in elements.keys():
                        cells[result[0][0]] += (elem,)
                    
                    if elem in cells.keys():
                        cells[result[0][0]] += (elem,)
                continue

            #find reassignment of elements or sequences
            result = re.findall(r'(^[a-zA-Z]{1}[a-zA-Z0-9_]*)\s*=\s*[\(]?([a-zA-Z]{1}[a-zA-Z0-9_]*)[\)]?$', line)
            if len(result) != 0:
                
                if result[0][1] in elements.keys():
                    cells_order.append(result[0][0])
                    cells[result[0][0]] = (result[0][1],)
                    
                if result[0][1] in cells.keys():
                    cells_order.append(result[0][0])
                    cells[result[0][0]] = (result[0][1],)

        return (cells, cells_order)


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
        self.param_list['Vcor'] = ('l', 'angle')
        self.param_list['Marker'] = ()
        self.param_list['Monitor'] = ()
        self.param_list['Sextupole'] = ('l', 'k2')
        self.param_list['Octupole'] = ('l', 'k3')
        self.param_list['Undulator'] = ('lperiod', 'nperiods', 'Kx', 'Ky')
        self.param_list['Cavity'] = ('l', 'v', 'freq', 'phi')
        self.param_list['TDCavity'] = ('l', 'v', 'freq', 'phi')
        self.param_list['Solenoid'] = ('l', 'k0')
        self.param_list['Multipole'] = ('kn',)
        self.param_list['Matrix'] = ('l',)
        self.param_list['UnknownElement'] = ('l',)


    def parsing(self, gui_lattice):

        self.gui_lattice = gui_lattice
        lines = ''

        if self.gui_lattice.tws and isinstance(self.gui_lattice.tws, Twiss):
            lines += self.parsing_twiss(self.gui_lattice.tws)

        if self.gui_lattice.beam and isinstance(self.gui_lattice.beam, Beam):
            lines += self.parsing_beam(self.gui_lattice.beam)

        if self.gui_lattice.elements and isinstance(self.gui_lattice.elements, dict):
            lines += self.parsing_elements(self.gui_lattice.elements)

        if self.gui_lattice.cells and isinstance(self.gui_lattice.cells, dict):
            lines += self.parsing_cells(self.gui_lattice.cells, self.gui_lattice.cells_order)

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

        lines = 'tws = Twiss()\n'
        for param in tws.__dict__:
            if tws.__dict__[param] != 0.0 and \
               tws.__dict__[param] != 0 and \
               tws.__dict__[param] != '' and \
               tws.__dict__[param] != 'id':

                lines += 'tws.' + str(param) + '=' + str(tws.__dict__[param]) + '\n'

        lines += '\n'
        return lines


    def parsing_beam(self, beam):

        lines = 'beam = Beam()\n'
        for param in beam.__dict__:
            if beam.__dict__[param] != 0.0 and \
               beam.__dict__[param] != 0 and \
               beam.__dict__[param] != '':

                lines += 'beam.' + str(param) + '=' + str(beam.__dict__[param]) + '\n'

        lines += '\n'
        return lines


    def parsing_cells(self, cells, cells_order):
        lines = ''

        for cell in cells_order:
            lines += str(cell) + ' = ('

            for elem in cells[cell]:
                lines += str(elem) + ', '

            lines = lines[0:-2]
            lines += ')\n'

        lines += '\n'
        return lines
