'''
Elegant <--> Ocelot lattice converter
'''

import re
import sys

from ocelot.cpbd.magnetic_lattice import *
from ocelot.cpbd.elements import *
from ocelot.cpbd.io import *


class ElegantStructureConverter:
    '''
    Main Elegant <--> Ocelot lattice converter class
    '''

    def __init__(self):

        # elements for conversion matrix
        self.elegant_matrix = {}
        self.elegant_matrix['DRIF'] = {'type':Drift, 'params':{'L':'l'}}
        self.elegant_matrix['SBEN'] = {'type':SBend, 'params':{'L':'l', 'ANGLE':'angle', 'E1':'e1', 'E2':'e2'}}
        self.elegant_matrix['QUAD'] = {'type':Quadrupole, 'params':{'L':'l', 'K1':'k1', 'TILT':'tilt'}}
        self.elegant_matrix['SEXT'] = {'type':Sextupole, 'params':{'L':'l', 'K2':'k2'}}
        self.elegant_matrix['MONI'] = {'type':Marker, 'params':{}}


    def elegant2ocelot(self, file_name):
        '''
        filename - input Elegant lattice filename
        '''

        with open(file_name) as file_link:
            data = file_link.read()

        # merge splitted lines
        data = re.sub(r'&\s*\n', '', data)

        # element names correction
        bad_element_names = re.findall(r'"(.*)":', data)

        for element in bad_element_names:
            element_new = re.sub(r'[\W]+', '_', element)
            data = re.sub(r''+str(element), str(element_new), data)

        data = re.sub(r'"', '', data)
        data = re.sub(r'\'', '', data)

        lines = data.split('\n')

        # parsing lines
        elegant_elements_dict = {}
        used_cell_name = ''
        cell_flag = False
        for line in lines:
            string = line.upper().strip()
            string = re.sub(r'\s*', '', string)

            # skip comment line
            if string == '':
                continue
            if string[0] == '!':
                continue

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

                # convert element
                if param[0] in self.elegant_matrix.keys():
                    # creat element
                    elements_list[elem] = self.elegant_matrix[param[0]]['type'](eid=elem)

                    # parse parameters
                    for data in param:
                        result = data.split('=')
                        if result[0] in self.elegant_matrix[param[0]]['params']:
                            tmp = self.elegant_matrix[param[0]]['params'][result[0]]
                            elements_list[elem].__dict__[tmp] = float(result[1])

                # replace element by Drift (if it has L) or skip
                else:
                    elements_list[elem] = None
                    print_skiped = True
                    
                    for data in param:
                        if data[:2] == 'L=':
                            elements_list[elem] = Drift(eid=elem, l=float(data[2:]))
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
        '''
        Replace single multiplication of elements
        '''

        rep_data = re.findall(r'(\d+)\*([\w]+)', line)

        for arr in rep_data:
            old_element = arr[0]+'*'+arr[1]
            new_element = arr[1] + (',' + arr[1]) * (int(arr[0]) - 1)
            line = line.replace(old_element, new_element)

        return line


    def replace_multiplications(self, line):
        '''
        Replace multiplication with brackets of elements
        '''

        rep_data = re.findall(r'(\d+)\*\(([\w\,]+)\)', line)

        while rep_data != []:

            for arr in rep_data:
                old_element = arr[0]+'*('+arr[1]+')'
                new_element = arr[1] + (',' + arr[1]) * (int(arr[0]) - 1)
                line = line.replace(old_element, new_element)

            rep_data = re.findall(r'(\d+)\*\(([\w\,]+)\)', line)

        return line


    def replace_cells(self, cell, elements_dict):
        '''
        Replace cells by elements
        '''

        if cell[0:4] != 'LINE':
            return [cell]

        cell = cell[6:-1]

        elements_list = cell.split(',')

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
        '''
        lattice - Ocelot lattice objectl
        filename - output Elegant lattice filename
        '''

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
    '''
    # example of Elegant - Ocelot convertion
    SC = ElegantStructureConverter()
    read_cell = SC.elegant2ocelot('elegant_test.lte')

    lattice = MagneticLattice(read_cell)
    write_lattice(lattice, remove_rep_drifts=False)
    '''
    '''
    # example of Ocelot - Elegant convertion
    from lattice import *
    lattice = MagneticLattice(cell)

    SC = ElegantStructureConverter()
    SC.ocelot2elegant(lattice)
    '''
