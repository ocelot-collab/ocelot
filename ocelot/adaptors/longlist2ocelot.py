
"""
athors: Ye. Fomin (NRC KI); S. Tomin (XFEL, NRC KI), 2016. 
"""

import xlrd
from ocelot.cpbd.magnetic_lattice import *
from ocelot.cpbd.elements import *
from ocelot.cpbd.io import *
from numpy import sin

class StructureConverter:
    
    def __init__(self):
        self.types = []
        #self.screens = ["OTRC", "OTRA", "OTRB"]

    def longlist_matrix_init(self):
        
        self.longlist_matrix = {}

        self.longlist_matrix['MAGNET'] = {}
        self.longlist_matrix['MAGNET']['SBEN'] = {'type': SBend, 'strength': 'angle', 'e1_lag': 'e1', 'e2_freq': 'e2'}
        self.longlist_matrix['MAGNET']['RBEN'] = {'type': RBend, 'strength': 'angle'}
        self.longlist_matrix['MAGNET']['QUAD'] = {'type': Quadrupole, 'strength': ['k1', '1./length']}
        self.longlist_matrix['MAGNET']['SEXT'] = {'type': Sextupole, 'strength': ['k2', '1./length']}
        self.longlist_matrix['MAGNET']['OCTU'] = {'type': Octupole, 'strength': ['k3', '1./length']}
        if "HKIC" in self.types:
            self.longlist_matrix['MAGNET']['HKIC'] = {'type': Hcor, 'strength': 'angle'}
        if "VKIC" in self.types:
            self.longlist_matrix['MAGNET']['VKIC'] = {'type': Vcor, 'strength': 'angle'}
        self.longlist_matrix['MAGNET']['SOLE'] = {'type': Solenoid, 'strength': 'k'}

        #self.longlist_matrix['FASTKICK'] = {}
        #self.longlist_matrix['FASTKICK']['HKIC'] = {'type': Hcor, 'strength': 'angle'}
        #self.longlist_matrix['FASTKICK']['VKIC'] = {'type': Vcor, 'strength': 'angle'}
        #self.longlist_matrix['FASTKICK']['RBEN'] = {'type': RBend, 'strength': 'angle'}
        
        #self.longlist_matrix['PMAGNET'] = {}
        #self.longlist_matrix['PMAGNET']['RBEN'] = {'type': RBend, 'strength': 'angle'}

        #self.longlist_matrix['MOVER'] = {}
        #self.longlist_matrix['MOVER']['HKIC'] = {'type': Hcor, 'strength': 'angle'}
        #self.longlist_matrix['MOVER']['VKIC'] = {'type': Vcor, 'strength': 'angle'}

        self.longlist_matrix['CAVITY'] = {}
        self.longlist_matrix['CAVITY']['LCAV'] = {'type': Cavity, 'strength': ['v', '1.e-3'], 'e1_lag': ['phi', '360.0'], 'e2_freq': ['f', '1000000']}

        self.longlist_matrix['DIAG'] = {}
        if "MONI" in self.types:
            self.longlist_matrix['DIAG']['MONI'] = {'type': Monitor}

        if "INSTR" in self.types:
            self.longlist_matrix['DIAG']['INSTR'] = {'type': Marker}

        if "MARK" in self.types:
            self.longlist_matrix['MARK'] = {}
            self.longlist_matrix['MARK']['MARK'] = {'type': Marker}

        self.longlist_matrix['UNDU'] = {}
        self.longlist_matrix['UNDU']['UNDULATOR'] = {'type': Undulator}

    def longlist_transform(self, row):
        
        # longlist parameters
        name_pos =     1+2
        ps_id_pos =    1+3
        group_pos =    1+4
        class_pos =    1+5
        length_pos =   1+7
        strength_pos = 1+8
        e1_lag_pos =   1+9
        e2_freq_pos =  1+10
        tilt_pos =     1+11
        s_pos =        1+12

        if row[0] == '': return None, 0.0

        element = None
        if row[group_pos] in self.longlist_matrix and row[class_pos] in self.longlist_matrix[row[group_pos]]:
                
            group_elem = self.longlist_matrix[row[group_pos]]
            class_elem = group_elem[row[class_pos]]
            
            element = class_elem['type'](eid=row[name_pos])
            element.ps_id = row[ps_id_pos]

            length = row[length_pos]
            if row[length_pos] != 0.0:
                element.l = length
            
            for key in class_elem.keys():
                if key != 'type':
                    pos = eval(key + '_pos')
                    if class_elem[key].__class__ == list:
                        element.__dict__[class_elem[key][0]] = row[pos] * eval(class_elem[key][1])
                    else:
                        element.__dict__[class_elem[key]] = row[pos]

            if row[tilt_pos] != 0:
                element.__dict__['tilt'] = row[tilt_pos]
                
        return element, row[s_pos]

    def sbend_l_correction(self, element):
        """
        correction SBEN length
        """

        if element.__class__ == SBend and element.angle != 0.0:
            if element.e1 != 0.0 and element.e2 == 0.0 or element.e1 == 0.0 and element.e2 != 0.0:
                element.l *= element.angle / sin(element.angle)
            else:
                element.l *= element.angle * 0.5 / sin(element.angle * 0.5)


    def Longlist2Ocelot(self, filename, sheet_name='LONGLIST', pos_start=3, pos_stop=None, sbend_l_corr=False):
        """
        pos_start - start line number into Excel file (1 is the first line, 3 is the first line with data)
        pos_stop - stop line number into Excel file
        """

        self.longlist_matrix_init()
        book = xlrd.open_workbook(filename)
        sheet = book.sheet_by_name(sheet_name)
        
        if pos_stop is None:
            pos_stop = sheet.nrows
        if pos_stop > sheet.nrows:
            pos_stop = sheet.nrows

        pos_start -= 1
        if pos_start < 0:
            pos_start = 0
        if pos_start > pos_stop:
            pos_start = pos_stop

        elements_list = []
        for row_index in range(pos_start, pos_stop):
            element, element_pos = self.longlist_transform(sheet.row_values(row_index))

            if sbend_l_corr: 
                self.sbend_l_correction(element)

            if element is not None:
                elements_list.append([element, element_pos])

        cell = lattice_format_converter(elements_list)

        return cell

      
if __name__ == '__main__':

    SC = StructureConverter()
    cell = SC.Longlist2Ocelot('component_list_8.5.1.xls', pos_stop=364, sbend_l_corr=True)

    lattice = MagneticLattice(cell)
    write_lattice(lattice)