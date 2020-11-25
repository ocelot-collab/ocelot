import sys

from ocelot.cpbd.elements import *
from ocelot.cpbd.beam import Beam, Twiss


class LatticeIO:
    """
    Utility class that holds all IO functions of MagneticLattice for storing it in a python file.
    """

    @staticmethod
    def save_lattice(lattice, tws0=None, file_name="lattice.py", remove_rep_drifts=True, power_supply=False):
        """
        saves lattice as python input file
        file_name - name of the file
        remove_rep_drifts - if True, remove the drifts with the same lengths from the lattice drifts definition
        """
        if remove_rep_drifts:
            lattice.rem_drifts()

        lines = LatticeIO.lat2input(lattice, tws0=tws0)

        if power_supply:
            lines = LatticeIO._write_power_supply_id(lattice, lines=lines)

        with open(file_name, 'w') as f:
            f.writelines(lines)

    @staticmethod
    def elements2input(lattice) -> str:
        """
        Generates a string, in a python readable format, that contains the elements in the lattice to store it in a python file.
        @param lattice: Input lattice
        @return: A string that contains the elements in the lattice in a python readable format
        """
        elements = LatticeIO._get_elements(lattice)
        elements_dict = LatticeIO._sort_elements(elements)
        lines = LatticeIO._print_elements(elements_dict)
        return lines

    @staticmethod
    def cell2input(lattice, split=False):
        """
        Generates a string, in a python readable format, that contains the cell of the lattice to store it in a python file.
        @param lattice: Input lattice
        @param split:
        @return: A string that contains the cell of the lattice in a python readable format
        """
        lines = []
        names = []
        for elem in lattice.sequence:
            # TODO: Edge and CouplerKick are not elements.
            if elem.__class__.__name__ not in ("Edge", "CouplerKick"):
                names.append(elem.name)

        new_names = []
        for i, name in enumerate(names):
            if split and i % 10 == 9:
                new_names.append('\n' + name)
            else:
                new_names.append(name)

        lines.append('cell = (' + ', '.join(new_names) + ')')

        return lines

    @staticmethod
    def lat2input(lattice, tws0=None):
        """
        Generates a string, in a python readable format, that contains the lattice to store it in a python file.
        @param lattice: Input lattice
        @return: A string that contains the lattice in a python readable format
        """
        lines = ['from ocelot import * \n']

        # prepare initial Twiss parameters
        if tws0 is not None and isinstance(tws0, Twiss):
            lines.append('\n#Initial Twiss parameters\n')
            lines.extend(LatticeIO.twiss2input(tws0))

        # prepare elements list
        lines.append('\n')
        lines.extend(LatticeIO.elements2input(lattice))

        # prepare cell list
        lines.append('\n# Lattice \n')
        lines.extend(LatticeIO.cell2input(lattice, True))

        lines.append('\n')

        return lines

    @staticmethod
    def twiss2input(twiss):
        """
        Generates a string, in a python readable format, that contains the Twiss parameter to store it in a python file.
        @param twiss: Input twiss
        @return: A string that contains Twiss parameter in a python readable format
        """
        lines = []
        tws_ref = Twiss()
        lines.append('tws0 = Twiss()\n')
        for param in twiss.__dict__:
            if twiss.__dict__[param] != tws_ref.__dict__[param]:
                lines.append('tws0.' + str(param) + ' = ' + str(twiss.__dict__[param]) + '\n')
        return lines

    @staticmethod
    def beam2input(beam):
        lines = []
        beam_ref = Beam()
        lines.append('beam = Beam()\n')
        for param in beam.__dict__:
            if beam.__dict__[param] != beam_ref.__dict__[param]:
                lines.append('beam.' + str(param) + ' = ' + str(beam.__dict__[param]) + '\n')

        return lines

    @staticmethod
    def _create_var_name(objects):
        alphabet = "abcdefgiklmn"
        ids = [obj.id for obj in objects]
        search_occur = lambda obj_list, name: [i for i, x in enumerate(obj_list) if x == name]
        for j, obj in enumerate(objects):
            inx = search_occur(ids, obj.id)
            if len(inx) > 1:
                for n, i in enumerate(inx):
                    name = ids[i]
                    name = name.replace('.', '_')
                    name = name.replace(':', '_')
                    name = name.replace('-', '_')
                    ids[i] = name + alphabet[n]
            else:
                name = ids[j]
                name = name.replace('.', '_')
                name = name.replace(':', '_')
                name = name.replace('-', '_')
                ids[j] = name
            obj.name = ids[j].lower()

        return objects

    @staticmethod
    def _get_elements(lattice):
        """
        Filters the elements in lattice and remove the fake elements (Edge, CouplerKick).
        @param lattice: input lattice
        @return: A list of elements
        """
        elements = []
        for element in lattice.sequence:
            element_type = element.__class__.__name__
            # TODO: Edge and CouplerKick are not elements.
            if element_type in ('Edge', "CouplerKick"):
                continue
            if element not in elements:
                elements.append(element)
        elements = LatticeIO._create_var_name(elements)
        return elements

    @staticmethod
    def _find_obj_and_create_name(lattice, types):
        objects = LatticeIO._find_objects(lattice, types=types)
        objects = LatticeIO._create_var_name(objects)
        return objects

    @staticmethod
    def _find_objects(lattice, types):
        """
        Function finds objects by types and adds it to list if object is unique.
        :param types: types of the Elements
        :return: list of elements
        """
        obj_id = []
        objs = []
        for elem in lattice.sequence:
            if elem.__class__ in types:
                if id(elem) not in obj_id:
                    objs.append(elem)
                    obj_id.append(id(elem))

        return objs

    @staticmethod
    def _write_power_supply_id(lattice, lines=[]):
        quads = LatticeIO._find_obj_and_create_name(lattice, types=["Quadrupole"])
        sexts = LatticeIO._find_obj_and_create_name(lattice, types=["Sextupole"])
        octs = LatticeIO._find_obj_and_create_name(lattice, types=["Octupole"])
        cavs = LatticeIO._find_obj_and_create_name(lattice, types=["Cavity"])
        bends = LatticeIO._find_obj_and_create_name(lattice, types=["Bend", "RBend", "SBend"])

        lines.append("\n# power supplies \n")
        for elem_group in [quads, sexts, octs, cavs, bends]:
            lines.append("\n#  \n")
            for elem in elem_group:
                if "ps_id" in dir(elem):
                    line = elem.name.lower() + ".ps_id = '" + elem.ps_id + "'\n"
                    lines.append(line)
        return lines

    @staticmethod
    def _sort_elements(elements):
        """
        Sort the elements by the element type.
        @param elements: A list of elements
        @return: A dict with the elements sorted by element type
        """
        elements_dict = {}
        for element in elements:
            element_type = element.__class__.__name__

            if element_type not in elements_dict:
                elements_dict[element_type] = []

            elements_dict[element_type].append(element)

        return elements_dict

    @staticmethod
    def _print_elements(elements_dict):
        """
        Creates a string, in a python readable format, of all elements in a Lattice sorted by the element types
        @param elements_dict: 
        @return: 
        """
        elements_order = []
        elements_order.append('Drift')
        elements_order.append('Quadrupole')
        elements_order.append('SBend')
        elements_order.append('RBend')
        elements_order.append('Bend')
        elements_order.append('Sextupole')
        elements_order.append('Octupole')
        elements_order.append('Multipole')
        elements_order.append('Hcor')
        elements_order.append('Vcor')
        elements_order.append('Undulator')
        elements_order.append('Cavity')
        elements_order.append('TDCavity')
        elements_order.append('Solenoid')
        elements_order.append('Monitor')
        elements_order.append('Marker')
        elements_order.append('Matrix')
        elements_order.append('Aperture')

        lines = []
        ordered_dict = {}
        unordered_dict = {}

        # sort on ordered and unordered elements dicts
        for type in elements_dict:
            if type in elements_order:
                ordered_dict[type] = elements_dict[type]
            else:
                unordered_dict[type] = elements_dict[type]

        # print ordered elements
        for type in elements_order:

            if type in ordered_dict:

                lines.append('\n# ' + type + 's\n')

                for element in ordered_dict[type]:
                    string = LatticeIO.element_def_string(element)
                    lines.append(string)

        # print remaining unordered elements
        for type in unordered_dict:

            lines.append('\n# ' + type + 's\n')

            for element in unordered_dict[type]:
                string = LatticeIO.element_def_string(element)
                lines.append(string)

        # delete new line symbol from the first line
        if lines != []:
            lines[0] = lines[0][1:]
        return lines

    @staticmethod
    def _matrix_def_string(element, params):
        """
        Creates a string, in a python readable format, for a matrix element to store it in a python file.
        This function is be used by element_def_string.
        @param element: input Element
        @return: A String that contains an matrix element in a python readable format
        """
        for key in element.__dict__:
            if isinstance(element.__dict__[key], np.ndarray):
                # r - elements
                if np.shape(element.__dict__[key]) == (6, 6):
                    for i in range(6):
                        for j in range(6):
                            val = element.__dict__[key][i, j]
                            if np.abs(val) > 1e-9:
                                params.append(key + str(i + 1) + str(j + 1) + '=' + str(val))
                # t - elements
                elif np.shape(element.__dict__[key]) == (6, 6, 6):
                    for i in range(6):
                        for j in range(6):
                            for k in range(6):
                                val = element.__dict__[key][i, j, k]
                                if np.abs(val) > 1e-9:
                                    params.append(key + str(i + 1) + str(j + 1) + str(k + 1) + '=' + str(val))
                # b - elements
                if np.shape(element.__dict__[key]) == (6, 1):
                    for i in range(6):
                        val = element.__dict__[key][i, 0]
                        if np.abs(val) > 1e-9:
                            params.append(key + str(i + 1) + '=' + str(val))
        return params

    @staticmethod
    def element_def_string(element) -> str:
        """
        Creates a string, in a python readable format, for an element to store it in a python file.
        @param element: input Element
        @return: A String that contains an element in a python readable format
        """
        params = []

        element_type = element.__class__.__name__
        element_ref = getattr(sys.modules[__name__], element_type)()
        params_order = element_ref.__init__.__code__.co_varnames
        argcount = element_ref.__init__.__code__.co_argcount

        for param in params_order[:argcount]:
            if param == 'self':
                continue

            # fix for parameter 'eid'
            if param == 'eid':
                params.append('eid=\'' + element.id + '\'')
                continue

            if isinstance(element.__dict__[param], np.ndarray):

                if not np.array_equal(element.__dict__[param], element_ref.__dict__[param]):
                    params.append(param + '=' + np.array2string(element.__dict__[param], separator=', '))
                continue

            if isinstance(element.__dict__[param], (int, float, complex)):

                # fix for parameters 'e1' and 'e2' in RBend element
                if element_type == 'RBend' and param in ('e1', 'e2'):
                    val = element.__dict__[param] - element.angle / 2.0
                    if val != 0.0:
                        params.append(param + '=' + str(val))
                    continue

                if element.__dict__[param] != element_ref.__dict__[param]:
                    params.append(param + '=' + str(element.__dict__[param]))
                continue

            if isinstance(element.__dict__[param], str):

                if element.__dict__[param] != element_ref.__dict__[param]:
                    params.append(param + '=\'' + element.__dict__[param] + '\'')
                continue

        if element.__class__.__name__ == "Matrix":
            params = LatticeIO._matrix_def_string(element, params)

        # join all parameters to element definition
        string = LatticeIO._pprinting(element, element_type, params)
        return string

    @staticmethod
    def _pprinting(element, element_type, params) -> str:
        string = element.name + ' = ' + element_type + '('
        n0 = len(string)
        n = n0
        for i, param in enumerate(params):
            n += len(params)
            if n > 250:
                string += "\n"
                string += " " * n0 + param + ", "
                n = n0 + len(param) + 2
            else:
                if i == len(params) - 1:
                    string += param
                else:
                    string += param + ", "
        string += ")\n"
        return string
