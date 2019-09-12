"""Elegant <--> Ocelot lattice converter"""

import re
import sys
import operator
import numpy as np

from ocelot.cpbd.magnetic_lattice import *
from ocelot.cpbd.elements import *
from ocelot.cpbd.io import *


class ElegantLatticeConverter:
    """Main Elegant <--> Ocelot lattice converter class"""

    def __init__(self):
        pass
    
    
    def init_convert_matrix(self):
        """Init Elegant -> Ocelot convertion matrix"""

        self.elegant_matrix = {}
        self.elegant_matrix['DRIF'] = {'type': Drift, 'params':{'L':'l'}}
        self.elegant_matrix['DRIFT'] = {'type': Drift, 'params': {'L': 'l'}}
        self.elegant_matrix['LSCDRIFT'] = {'type': Drift, 'params': {'L': 'l'}}
        self.elegant_matrix['CSRDRIFT'] = {'type': Drift, 'params': {'L': 'l'}}
        self.elegant_matrix['SOLE'] = {'type': Solenoid, 'params': {'L': 'l'}}
        self.elegant_matrix['SBEN'] = {'type': SBend,
                                       'params':{'L': 'l', 'ANGLE': 'angle', "K1": "k1", "K2": "k2", "FINT": "fint",
                                                 'E1': 'e1', 'E2': 'e2', "HGAP": ["gap", "2"], 'TILT': 'tilt'}}
        self.elegant_matrix['SBEND'] = {'type': SBend,
                                        'params':{'L': 'l', 'ANGLE': 'angle', "K1": "k1", "K2": "k2", "FINT": "fint",
                                                 'E1': 'e1', 'E2': 'e2', "HGAP": ["gap", "2"], 'TILT': 'tilt'}}
        self.elegant_matrix['RBEN'] = {'type': RBend,
                                       'params':{'L': 'l', 'ANGLE': 'angle', "K1": "k1", "K2": "k2", "FINT": "fint",
                                                 'E1': 'e1', 'E2': 'e2', "HGAP": ["gap", "2"], 'TILT': 'tilt'}}
        self.elegant_matrix['QUAD'] = {'type': Quadrupole, 'params': {'L':'l', 'K1':'k1', "K2": "k2", 'TILT': 'tilt'}}
        self.elegant_matrix['SEXT'] = {'type': Sextupole, 'params': {'L':'l', 'K2':'k2'}}
        self.elegant_matrix['MONI'] = {'type': Monitor, 'params': {}}
        self.elegant_matrix['MARK'] = {'type': Marker, 'params': {}}
        self.elegant_matrix['WATCH'] = {'type': Marker, 'params': {}}
        self.elegant_matrix['RFCA'] = {'type': Cavity, 'params': {'L':'l', 'VOLT':['v','1.0e-9'], 'FREQ': 'freq', 'PHASE': 'phi'}}
        self.elegant_matrix['HKICK'] = {'type': Hcor, 'params': {'L':'l'}}
        self.elegant_matrix['VKICK'] = {'type': Vcor, 'params': {'L':'l'}}
        self.elegant_matrix['WIGGLER'] = {'type': Undulator, 'params': {'L': 'l', 'K': 'Kx', "POLES": ["nperiods", "0.5"]}}
        self.elegant_matrix['EMATRIX'] = {'type': Matrix, 'params': {'L':'l',
                                                'R11': 'rm11', 'R12': 'rm12', 'R13': 'rm13', 'R14': 'rm14', 'R15': 'rm15', 'R16': 'rm16',
                                                'R21': 'rm21', 'R22': 'rm22', 'R23': 'rm23', 'R24': 'rm24', 'R25': 'rm25', 'R26': 'rm26',
                                                'R31': 'rm31', 'R32': 'rm32', 'R33': 'rm33', 'R34': 'rm34', 'R35': 'rm35', 'R36': 'rm36',
                                                'R41': 'rm41', 'R42': 'rm42', 'R43': 'rm43', 'R44': 'rm44', 'R45': 'rm45', 'R46': 'rm46',
                                                'R51': 'rm51', 'R52': 'rm52', 'R53': 'rm53', 'R54': 'rm54', 'R55': 'rm55', 'R56': 'rm56',
                                                'R61': 'rm61', 'R62': 'rm62', 'R63': 'rm63', 'R64': 'rm64', 'R65': 'rm65', 'R66': 'rm66',}}
        self.elegant_matrix['CSRCSBEND'] = {'type': SBend, 'params': {'L':'l', 'ANGLE': 'angle', "K1": "k1", "K2": "k2", 'E1': 'e1', 'E2': 'e2', 'TILT': 'tilt'}}
        self.elegant_matrix['RFCW'] = {'type': Cavity, 'params': {'L': 'l', 'VOLT': ['v','1.0e-9'], 'FREQ': 'freq', 'PHASE': 'phi'}}


    def fix_convert_matrix(self):
        """Init and fix Ocelot -> Elegant convertion matrix"""

        self.init_convert_matrix()

        # Some elements are deleted from the matrix to exclude ambiguity in elements conversion
        del self.elegant_matrix['CSRCSBEND']
        del self.elegant_matrix['RFCW']


    def calc_rpn(self, expression, constants={}, info=''):
        """
        Calculation of the expression written in reversed polish notation
        """
        
        operators = {'+': operator.add, '-': operator.sub, '*': operator.mul, '/': operator.truediv}
        functions = {'SIN': np.sin, 'COS': np.cos, 'TAN': np.tan, 'ASIN': np.arcsin, 'ACOS': np.arccos, 'ATAN': np.arctan, 'SQRT': np.sqrt}
        
        stack = [0]
        for token in expression.split(" "):
            
            if token in operators:
                op2, op1 = stack.pop(), stack.pop()
                stack.append(operators[token](op1, op2))
                
            elif token in constants:
                stack.append(constants[token])
            
            elif token in functions:
                op = stack.pop()
                stack.append(functions[token](op))
                
            elif token:
                try:
                    stack.append(float(token))
                except ValueError:
                    print('********* ERROR! NO ' + token + ' in function conversion or constants list in calc_rpn function. ' + info)
                    sys.exit()
        
        return stack.pop()
    
    
    def convert_val(self, str, constants, info):
        """
        Convert string input value to float or change input string value by value from constants array
        """
        
        try:
            value = float(str)
        except ValueError:
            if str in constants:
                value = constants[str]
            else:
                value = self.calc_rpn(str, constants, 'In element '+info)
        
        return value

        
    def elegant2ocelot(self, file_name):
        """
        :filename - input Elegant lattice filename
        """
        
        # just in case double check
        # because in self.ocelot2elegant function these elements are deleted
        self.init_convert_matrix()

        with open(file_name) as file_link:
            data = file_link.read()

        # delete comments
        data = re.sub(r'![^\n]*\n', '\n', data)
        # merge splitted lines
        data = re.sub(r'&\s*\n', '', data)

        # element names correction
        # remove parentheses from the names
        bad_element_names = re.findall(r'.\[.*\].*:', data)
        for element in bad_element_names:
            element = element.replace(":", "")
            element = element.strip()
            element_new = element.replace("[", "")
            element_new = element_new.replace("]", "")
            data = data.replace(element, element_new)

        # element names correction
        bad_element_names = re.findall(r'"(.*)":', data)
        for element in bad_element_names:
            element_new = re.sub(r'[\W]+', '_', element)
            data = re.sub(r''+str(element), str(element_new), data)

        #data = re.sub(r'"', '', data)
        data = re.sub(r'\'', '"', data)
        lines = data.split('\n')

        # parsing lines
        elegant_elements_dict = {}
        used_cell_name = ''
        cell_flag = False
        constants = {}
        
        for line in lines:
            # print(line)
            string = line.upper().strip()

            # skip empty and comment lines
            if string == '' or string[0] == '!':
                continue
            
            # parsing constants
            if string[0] == '%' and string.find('STO') is not -1:
                expr = string[1:].split('STO')
                constants[expr[1].strip()] = self.calc_rpn(expr[0].strip(), constants, 'In string '+string)
                continue
            
            # delete spaces and quotes
            string0 = string
            string = ''
            del_space = True
            for s in string0:
                if s == '"':
                    del_space = not del_space
                    continue
                elif s == ' ' and del_space:
                    continue
                string += s

            # split line and add to the dictionary with element names and descriptions
            element_desc = string.split(':')
            element_params = ''

            for i in range(1, len(element_desc)):
                element_params += '_' + str(element_desc[i])
            element_params = element_params[1:]
            elegant_elements_dict.update({element_desc[0]:element_params})

            # finding used lattice cell
            if cell_flag is False and element_params != '':
                used_cell_name = element_desc[0]
            if cell_flag is False and element_desc[0][:3] == 'USE':
                used_cell_name = element_desc[0][4:]
                cell_flag = True

        # parsing found used cell
        if used_cell_name in elegant_elements_dict:
            elegant_cell = elegant_elements_dict[used_cell_name]
        else:
            print('********* ERROR! used unkown cell or element ' + used_cell_name + ' *********')
            sys.exit()
        # replace multiplications
        elegant_cell = self.replace_s_multiplications(elegant_cell)
        elegant_cell = self.replace_multiplications(elegant_cell)
        # replace cells by elements
        elegant_cell = self.replace_cells(elegant_cell, elegant_elements_dict)

        # create Ocelot cell
        ocelot_cell = ()
        elements_list = {}

        for elem in elegant_cell:

            if elem not in elements_list.keys():

                param = elegant_elements_dict[elem].split(',')
                param = [elem.strip() for elem in param]
                # convert element
                if param[0] in self.elegant_matrix.keys():
                    # create element
                    elements_list[elem] = self.elegant_matrix[param[0]]['type'](eid=elem)

                    # parse parameters
                    for data in param:
                        result = data.split('=')
                        if result[0] in self.elegant_matrix[param[0]]['params']:
                            tmp = self.elegant_matrix[param[0]]['params'][result[0]]
                            
                            val = self.convert_val(result[1], constants, elem)
                            
                            if tmp.__class__ == list:
                                elements_list[elem].__dict__[tmp[0]] = val * float(tmp[1])
                            else:
                                # fix for phi cavity
                                if elements_list[elem].__class__ == Cavity and tmp == 'phi':
                                    val = 90.0 - val

                                elements_list[elem].__dict__[tmp] = val
                    if elements_list[elem].__class__ == Undulator:
                        elements_list[elem].lperiod = elements_list[elem].l/elements_list[elem].nperiods
                # replace element by Drift (if it has L) or skip
                else:
                    elements_list[elem] = None
                    print_skiped = True
                    
                    for data in param:
                        if data[:2] == 'L=':
                            
                            val = self.convert_val(data[2:], constants, elem)
        
                            elements_list[elem] = Drift(eid=elem, l=val)
                            print_skiped = False
                            break

                    if print_skiped:
                        print('WARNING! Unknown element', elem, 'with type', param[0], 'was skiped')
                    else:
                        print('WARNING! Unknown element', elem, 'with type', param[0], 'was changed by Drift')
                        

            if elements_list[elem] is not None:
                ocelot_cell += (elements_list[elem],)

        return ocelot_cell


    def replace_s_multiplications(self, line):
        """Replace single multiplication of elements"""

        rep_data = re.findall(r'(\d+)\*([\w]+)', line)

        for arr in rep_data:
            old_element = arr[0]+'*'+arr[1]
            new_element = arr[1] + (',' + arr[1]) * (int(arr[0]) - 1)
            line = line.replace(old_element, new_element)

        return line


    def replace_multiplications(self, line):
        """Replace multiplication with brackets of elements"""

        rep_data = re.findall(r'(\d+)\*\(([\w\,]+)\)', line)

        while rep_data != []:

            for arr in rep_data:
                old_element = arr[0]+'*('+arr[1]+')'
                new_element = arr[1] + (',' + arr[1]) * (int(arr[0]) - 1)
                line = line.replace(old_element, new_element)

            rep_data = re.findall(r'(\d+)\*\(([\w\,]+)\)', line)

        return line


    def replace_cells(self, cell, elements_dict):
        """Replace cells by elements"""

        if cell[0:4] != 'LINE':
            return [cell]

        cell = cell[6:-1]

        elements_list = cell.split(',')
        elements_list = [elem.strip() for elem in elements_list]
        check_lines = True
        while check_lines:
            new_list = []
            check_lines = False

            for i in range(0, len(elements_list)):
                reverse = False
                if elements_list[i][0] == '-':
                    reverse = True
                    elements_list[i] = elements_list[i][1:]

                if elements_list[i] not in elements_dict:
                    print('********* ERROR! used undefined cell or element ' + elements_list[i] + ' *********')
                    sys.exit()

                if elements_dict[elements_list[i]][0:4] == 'LINE':
                    check_lines = True
                    tmp_list = elements_dict[elements_list[i]][6:-1].split(',')
                    if reverse:
                        start, stop, step = len(tmp_list)-1, -1, -1
                    else:
                        start, stop, step = 0, len(tmp_list), 1
                    for j in range(start, stop, step):
                        new_list.append(tmp_list[j])
                else:
                    new_list.append(elements_list[i])

            elements_list = new_list

        return elements_list


    def ocelot2elegant(self, lattice, file_name='lattice.lte'):
        """
        :lattice - Ocelot lattice objectl
        :filename - output Elegant lattice filename
        """

        # fix due to elements type doubling
        self.fix_convert_matrix()

        reverse_matrix = {}
        elements_arr = {}
        elegant_cell = []
        
        for elegant_type in self.elegant_matrix:
            reverse_matrix.update({self.elegant_matrix[elegant_type]['type'].__name__:elegant_type})
            elements_arr.update({self.elegant_matrix[elegant_type]['type'].__name__:[]})

        for elem in lattice.sequence:

            # add element if its class it the matrix
            if elem.__class__.__name__ in reverse_matrix.keys():

                if elem not in elements_arr[elem.__class__.__name__]:
                    elements_arr[elem.__class__.__name__].append(elem)

            # if element class not in the maxrix - replace it by Drift (L > 0.0) or skip it
            else:
                if elem.l > 0.0:
                    elements_arr['Drift'].append(Drift(l=elem.l, eid=elem.id))
                else:
                    continue

            # add element in the final cell
            elegant_cell.append(elem.id)

        #  save data to file
        lines = ''

        # save elements
        for element_class in elements_arr:
            for elem in elements_arr[element_class]:
                lines += elem.id + ': ' + reverse_matrix[element_class]

                # check parameter 'L' and print it first
                if 'L' in self.elegant_matrix[reverse_matrix[element_class]]['params']:
                    tmp = self.elegant_matrix[reverse_matrix[element_class]]['params']['L']
                    lines += ',L=' + str(elem.__dict__[tmp])

                # print other parameters
                for param in self.elegant_matrix[reverse_matrix[element_class]]['params']:
                    if param == 'L':
                        continue
                    tmp = self.elegant_matrix[reverse_matrix[element_class]]['params'][param]
                    if tmp.__class__ == list:
                        lines += ',' + param + '=' + str(elem.__dict__[tmp[0]]/float(tmp[1]))
                    else:
                        # fix for phi cavity
                        if elem.__class__ == Cavity and tmp == 'phi':
                            lines += ',' + param + '=' + str(90.0 - elem.__dict__[tmp])
                            lines += ',CHANGE_P0=1,END1_FOCUS=1,END2_FOCUS=1,BODY_FOCUS_MODEL="SRS"'
                        else:
                            lines += ',' + param + '=' + str(elem.__dict__[tmp])
                lines += '\n'

        # save cell
        lines += 'CELL: LINE=(' + ', '.join(elegant_cell) + ')\n'
        lines += 'USE,CELL\n'
        lines += 'RETURN\n'

        fff = open(file_name, 'w')
        fff.writelines(lines)
        fff.close()

        return 0



if __name__ == '__main__':
    pass
    
    """
    # example of Elegant - Ocelot convertion
    SC = ElegantLatticeConverter()
    read_cell = SC.elegant2ocelot('elbe.lte')
    lattice = MagneticLattice(read_cell)
    write_lattice(lattice, remove_rep_drifts=False)
    """
    
    """
    # example of Ocelot - Elegant convertion
    from lattice import *
    lattice = MagneticLattice(cell)
    SC = ElegantLatticeConverter()
    SC.ocelot2elegant(lattice)
    """
