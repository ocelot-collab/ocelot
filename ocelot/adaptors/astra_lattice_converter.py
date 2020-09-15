"""Astra <--> Ocelot lattice converter"""

import re

from ocelot.cpbd.magnetic_lattice import *
from ocelot.cpbd.elements import *
from ocelot.cpbd.io import *


class AstraLatticeConverter:

    def __init__(self):

        self.astra_matrix = {}

        # used elements for conversion
        self.astra_matrix['QUADRUPOLE'] = {'type':Quadrupole, 'eid_prefix':'Q_', 'position':'Q_POS', 'params':{'Q_K':'k1', 'Q_LENGTH':'l'}, 'default_vals':{'Q_BORE':0.0000001}}
        #self.astra_matrix['CAVITY'] = {'type':Cavity, 'position':'C_pos', 'params':{'Nue':['f', '1000000000'], 'Phi':'phi', 'MaxE':'v', 'FILE_EFieLD':None}}

        # used sections and flags in Astra lattice file
        self.astra_section_flags = []
        self.astra_section_flags.append({'NEWRUN':[]})
        #self.astra_section_flags.append({'CAVITY':['LEFieLD=.T']})
        self.astra_section_flags.append({'QUADRUPOLE':['LQUAD=.T', 'Loop=.F']})


    def astra2ocelot(self, filename):
        """
        :filename - input Astra lattice filename
        """

        with open(filename) as fp:
            lines = fp.read().split("\n")

        flag = False
        element_type = None
        element_name_array = []
        elements_list = []

        for line in lines:
            string = line.upper().strip()

            # skip comment line
            if string == '': continue
            if string[0] == '!': continue

            # parsing section name
            if string[0] == '&':
                if string[1:] in self.astra_matrix.keys():
                    flag = True
                    element_type = string[1:]
                else:
                    print('********* WARNING! section ' + string[1:] + ' is found but not parsed *********')

            # check the end of the section
            if flag and line == '/':
                flag = False
                element_type = None

            # in the section
            if flag:
                result = re.findall(r'(\w*)\((\d*)\)\s*=\s*([\+\-]?\d*[\.]?\d*)', string)

                if result == []: continue

                # parse element name
                if 'eid_prefix' in self.astra_matrix[element_type].keys() and (self.astra_matrix[element_type]['eid_prefix'] != '' or self.astra_matrix[element_type]['eid_prefix'] != None):
                    element_name = self.astra_matrix[element_type]['eid_prefix'] + result[0][1]
                else:
                    element_name = element_type + result[0][1]

                # create new element
                if element_name not in element_name_array:
                    element_name_array.append(element_name)
                    elements_list.append([self.astra_matrix[element_type]['type'](eid=element_name), None])

                element_index = element_name_array.index(element_name)

                # parsing element position
                if result[0][0] == self.astra_matrix[element_type]['position']:
                    elements_list[element_index][1] = float(result[0][2])

                # parsing element parameters
                if result[0][0] in self.astra_matrix[element_type]['params']:
                    param = self.astra_matrix[element_type]['params'][result[0][0]]
                    elements_list[element_index][0].__dict__[param] = float(result[0][2])

        elements_list = sorted(elements_list, key=lambda x: x[1])
        cell = lattice_format_converter(elements_list)

        return cell


    def ocelot2astra(self, lattice, pos_start=0.0, filename='lattice.in'):
        """
        :lattice - Ocelot lattice object
        :pos_start - longitudinal position of the beginning of the first element
        :filename - output Astra lattice filename
        """

        lines_arr = {}
        for key in self.astra_matrix:
            if self.astra_matrix[key] != {}:
                lines_arr.update({key:[]})

        elem_center = pos_start
        elem_end = pos_start

        # prepare data from the Ocelot lattice
        for i, elem in enumerate(lattice.sequence):

            elem_center = elem_end + 0.5 * elem.l
            elem_end += elem.l

            elem_class = elem.__class__.__name__.upper()
            if elem_class in self.astra_matrix:
                tmp_data = {}

                for j in self.astra_matrix[elem_class]['params']:
                    tmp_data[j] = elem.__dict__[self.astra_matrix[elem_class]['params'][j]]

                for j in self.astra_matrix[elem_class]['default_vals']:
                    tmp_data[j] = self.astra_matrix[elem_class]['default_vals'][j]

                tmp_data[self.astra_matrix[elem_class]['position']] = elem_center

                lines_arr[elem_class].append({'name':elem.id, 'data':tmp_data})


        # save data to file
        lines = ''
        for i in self.astra_section_flags:
            for section_name in i:
                lines += '&' + section_name + '\n'
                if i[section_name] != []:
                    for flags in i[section_name]:
                        lines += flags + '\n'

                if section_name in lines_arr:
                    element_number = 0

                    for j in lines_arr[section_name]:
                        lines += '\n! ' + section_name.lower() + ' ' + j['name'] + '\n'

                        element_number += 1
                        for k in j['data']:
                            lines += '  ' + k + '(' + str(element_number)  + ')=' + str(j['data'][k]) + '\n'

            lines += '/\n\n'

        f = open(filename, 'w')
        f.writelines(lines)
        f.close()

        return 0



if __name__ == '__main__':
    pass
    
    '''
    # example of Astra - Ocelot convertion
    SC = AstraLatticeConverter()
    cell = SC.astra2ocelot('astra_test.in')
    lattice = MagneticLattice(cell)
    lattice.write_lattice()
    '''
    
    '''
    # example of Ocelot - Astra convertion
    from lattice import *
    lattice = MagneticLattice(cell)
    SC = AstraLatticeConverter()
    SC.ocelot2astra(lattice)
    '''
