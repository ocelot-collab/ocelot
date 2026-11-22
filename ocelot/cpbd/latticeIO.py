import inspect
from numbers import Number


_ELEMENT_IMPORT_ORDER = (
    'UnknownElement',
    'Aperture',
    'Bend',
    'Cavity',
    'Drift',
    'Element',
    'Hcor',
    'Marker',
    'Matrix',
    'Monitor',
    'Multipole',
    'Octupole',
    'Pulse',
    'Quadrupole',
    'RBend',
    'SBend',
    'Sextupole',
    'Solenoid',
    'TDCavity',
    'TWCavity',
    'Undulator',
    'Vcor',
    'XYQuadrupole',
)


class LatticeIO:
    """
    Utility class that holds all IO functions of MagneticLattice for storing it in a python file.
    """

    @staticmethod
    def save_lattice(lattice, twiss0=None, file_name="lattice.py", remove_rep_drifts=True, power_supply=False,
                     **kwargs):
        """
        Save a lattice as a Python input file.

        Args:
            lattice (MagneticLattice): Input lattice.
            twiss0 (Twiss, optional): Initial Twiss parameters to write at the beginning of the file.
                The output variable is named ``twiss0``. The old keyword ``tws0`` is still accepted
                through ``kwargs`` for backward compatibility. Passing both names raises a `ValueError`.
            file_name (str): Name of the file to write.
            remove_rep_drifts (bool): If ``True``, remove repeated drifts from the lattice drift definitions.
            power_supply (bool): If ``True``, write power supply IDs into the file.
        """
        twiss0 = LatticeIO._resolve_twiss0_alias(twiss0, kwargs, "save_lattice")

        if remove_rep_drifts:
            lattice.rem_drifts()

        lines = LatticeIO.lat2input(lattice, twiss0=twiss0)

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
        if any(not hasattr(elem, "name") for elem in lattice.sequence):
            LatticeIO._create_var_name(LatticeIO._unique_sequence_elements(lattice.sequence))
        names = [elem.name for elem in lattice.sequence]

        new_names = []
        for i, name in enumerate(names):
            if split and i % 10 == 9:
                new_names.append('\n' + name)
            else:
                new_names.append(name)

        lines.append('cell = (' + ', '.join(new_names) + ')')

        return lines

    @staticmethod
    def lat2input(lattice, twiss0=None, **kwargs):
        """
        Generates a string, in a python readable format, that contains the lattice to store it in a python file.
        :param lattice: Input lattice
        :param twiss0: Initial Twiss parameters. The output variable is named ``twiss0``.
            The old keyword ``tws0`` is still accepted through ``kwargs`` for backward compatibility.
            Passing both names raises a `ValueError`.
        :return: A string that contains the lattice in a python readable format
        """
        twiss0 = LatticeIO._resolve_twiss0_alias(twiss0, kwargs, "lat2input")

        include_twiss = False
        if twiss0 is not None:
            from ocelot.cpbd.beam.core import Twiss
            include_twiss = isinstance(twiss0, Twiss)

        lines = LatticeIO._import_lines(lattice, include_twiss=include_twiss)

        # prepare initial Twiss parameters
        if include_twiss:
            lines.append('\n#Initial Twiss parameters\n')
            lines.extend(LatticeIO.twiss2input(twiss0))

        # prepare elements list
        lines.append('\n')
        lines.extend(LatticeIO.elements2input(lattice))

        # prepare cell list
        lines.append('\n# Lattice \n')
        lines.extend(LatticeIO.cell2input(lattice, True))

        lines.append('\n')

        return lines

    @staticmethod
    def _import_lines(lattice, include_twiss=False):
        """
        Return explicit imports needed by the generated lattice file.

        The generated file should not use ``from ocelot import *`` because that
        forces the full legacy facade to load. Known Ocelot element wrappers are
        imported from ``ocelot.cpbd.elements``; custom element classes are
        imported from their defining module when possible.
        """
        elements = LatticeIO._unique_sequence_elements(lattice.sequence)
        element_types = {element.__class__.__name__ for element in elements}

        lines = []
        known_elements = [name for name in _ELEMENT_IMPORT_ORDER if name in element_types]
        if known_elements:
            lines.append('from ocelot.cpbd.elements import ' + ', '.join(known_elements) + '\n')

        custom_imports = []
        known_element_set = set(_ELEMENT_IMPORT_ORDER)
        for element in elements:
            cls = element.__class__
            if cls.__name__ in known_element_set:
                continue
            module_name = cls.__module__
            if module_name == '__main__':
                continue
            custom_imports.append((module_name, cls.__name__))

        for module_name, class_name in sorted(set(custom_imports)):
            lines.append(f'from {module_name} import {class_name}\n')

        if include_twiss:
            lines.append('from ocelot.cpbd.beam import Twiss\n')

        return lines

    @staticmethod
    def _resolve_twiss0_alias(twiss0, kwargs, function_name):
        if "tws0" in kwargs:
            legacy_twiss0 = kwargs.pop("tws0")
            if twiss0 is not None:
                raise ValueError(f"{function_name} accepts either twiss0 or legacy tws0, not both.")
            twiss0 = legacy_twiss0

        if kwargs:
            name = next(iter(kwargs))
            raise TypeError(f"{function_name}() got an unexpected keyword argument '{name}'")

        return twiss0

    @staticmethod
    def twiss2input(twiss):
        """
        Generates Python input lines for Twiss parameters.

        The output variable is named ``twiss0``.

        :param twiss: Input twiss
        :return: A string that contains Twiss parameter in a python readable format
        """
        from ocelot.cpbd.beam.core import Twiss

        lines = []
        tws_ref = Twiss()
        lines.append('twiss0 = Twiss()\n')
        for param in twiss.__dict__:
            if twiss.__dict__[param] != tws_ref.__dict__[param]:
                lines.append('twiss0.' + str(param) + ' = ' + str(twiss.__dict__[param]) + '\n')
        return lines

    @staticmethod
    def beam2input(beam):
        from ocelot.cpbd.beam.core import Beam

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
        def search_occur(obj_list, name): return [i for i, x in enumerate(obj_list) if x == name]
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
        Collect unique lattice elements in sequence order.
        :param lattice: input lattice
        :return: A list of elements
        """
        return LatticeIO._create_var_name(LatticeIO._unique_sequence_elements(lattice.sequence))

    @staticmethod
    def _unique_sequence_elements(sequence):
        """Return unique sequence elements while preserving their first occurrence order."""
        return list(dict.fromkeys(sequence))

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
        from ocelot.cpbd.elements import Bend, Cavity, Octupole, Quadrupole, RBend, SBend, Sextupole

        quads = LatticeIO._find_obj_and_create_name(lattice, types=[Quadrupole])
        sexts = LatticeIO._find_obj_and_create_name(lattice, types=[Sextupole])
        octs = LatticeIO._find_obj_and_create_name(lattice, types=[Octupole])
        cavs = LatticeIO._find_obj_and_create_name(lattice, types=[Cavity])
        bends = LatticeIO._find_obj_and_create_name(lattice, types=[Bend, RBend, SBend])

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
        :param elements: A list of elements
        :return: A dict with the elements sorted by element type
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
        :param elements_dict:
        :return:
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
        :param element: input Element
        :return: A String that contains an matrix element in a python readable format
        """
        import numpy as np

        for key in ("r", "t", "b"):
            value = getattr(element, key)
            if np.shape(value) == (6, 6):
                for i in range(6):
                    for j in range(6):
                        val = value[i, j]
                        if np.abs(val) > 1e-9:
                            params.append(key + str(i + 1) + str(j + 1) + '=' + str(val))
            elif np.shape(value) == (6, 6, 6):
                for i in range(6):
                    for j in range(6):
                        for k in range(6):
                            val = value[i, j, k]
                            if np.abs(val) > 1e-9:
                                params.append(key + str(i + 1) + str(j + 1) + str(k + 1) + '=' + str(val))
            elif np.shape(value) == (6, 1):
                for i in range(6):
                    val = value[i, 0]
                    if np.abs(val) > 1e-9:
                        params.append(key + str(i + 1) + '=' + str(val))
        return params

    @staticmethod
    def _element_init_params(element):
        """Return public wrapper constructor parameters used for serialization."""
        parameters = inspect.signature(type(element).__init__).parameters.values()
        return [
            param.name
            for param in parameters
            if param.name not in ("self", "tm")
            and param.kind not in (inspect.Parameter.VAR_POSITIONAL, inspect.Parameter.VAR_KEYWORD)
        ]

    @staticmethod
    def element_def_string(element) -> str:
        """
        Creates a string, in a python readable format, for an element to store it in a python file.
        :param element: input Element
        :return: A String that contains an element in a python readable format
        """
        import numpy as np

        params = []

        element_type = element.__class__.__name__
        element_ref = type(element)()
        params_order = LatticeIO._element_init_params(element)
        # Alignment kwargs belong to the magnet atom and must survive a round trip.
        if "dtilt" in vars(element.element):
            params_order.extend(name for name in ("dx", "dy", "dtilt") if name not in params_order)

        for param in params_order:
            # fix for parameter 'eid'
            if param == 'eid':
                params.append('eid=\'' + element.id + '\'')
                continue

            value = getattr(element, param)
            value_ref = getattr(element_ref, param)

            if isinstance(value, np.ndarray):
                if not np.array_equal(value, value_ref):
                    params.append(param + '=' + np.array2string(value, separator=', '))
                continue

            if isinstance(value, Number):

                # fix for parameters 'e1' and 'e2' in RBend element
                if element_type == 'RBend' and param in ('e1', 'e2'):
                    val = value - element.angle / 2.0
                    if val != 0.0:
                        params.append(param + '=' + str(val))
                    continue

                if value != value_ref:
                    params.append(param + '=' + str(value))
                continue

            if isinstance(value, str):

                if value != value_ref:
                    params.append(param + '=\'' + value + '\'')
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
