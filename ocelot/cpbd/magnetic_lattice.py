from ocelot.cpbd.elements import *
from ocelot.cpbd.beam import Twiss
import logging
import numpy as np

_logger = logging.getLogger(__name__)


def lattice_format_converter(elements):
    """
    :param elements: lattice in the format: [[elem1, center_pos1], [elem2, center_pos2], [elem3, center_pos3], ... ]
    :return: lattice in the format: [elem1, drift1, elem2, drift2, elem3, drift3, ...]
    """
    cell = []
    drift_num = 0
    s_pos = 0.0
    for element in elements:
        if element[0].__class__ in (Edge, CouplerKick):
            continue
        element_start = element[1] - element[0].l / 2.0
        if element_start < s_pos - 1.0e-14:                 # 1.0e-14 is used as crutch for precision of float
            if element[0].l == 0.0:

                if s_pos - element_start > 1.0e-2:
                    print("************** WARNING! Element " + element[0].id + " was deleted")
                    continue
                element[1] = s_pos
                print("************** WARNING! Element " + element[0].id + " was moved from " + str(element_start) +
                      " to " + str(s_pos))
            else:
                dl = element[0].l / 2.0  - element[1] + s_pos
                if cell[-1].__class__ == Marker and cell[-2].__class__ == Drift and cell[-2].l > dl:
                    cell[-2].l -= dl
                    print("************** WARNING! Element " + cell[-1].id + " was deleted")
                    cell.pop()
                else:
                    print("************** ERROR! Element " + element[0].id + " has bad position (overlapping?)")
                    exit()

        if element_start > s_pos + 1.0e-12:                 # 1.0e-12 is used as crutch for precision of float
            drift_num += 1
            drift_l = round(element_start - s_pos, 10)      # round() is used as crutch for precision of float
            drift_eid = 'D_' + str(drift_num)
            cell.append(Drift(l=drift_l, eid=drift_eid))

        cell.append(element[0])
        s_pos = element[1] + element[0].l / 2.0
    return tuple(cell)


def merger(lat, remaining_types=[], remaining_elems=[], init_energy=0.):
    """
    Function to compress the lattice excluding elements by type or by individual elements

    :param lat: MagneticLattice
    :param remaining_types: list, the type of the elements which needed to be untouched
                            others will be "compress" to Matrix element
                            e.g. [Monitor, Quadrupole, Bend, Hcor]
    :param remaining_elems: list of elements (ids (str) or object)
    :param init_energy: initial energy
    :return: New MagneticLattice
    """
    _logger.debug("element numbers before: " + str(len(lat.sequence)))
    lattice_analysis = []
    merged_elems = []
    for elem in lat.sequence:
        if elem.__class__ in remaining_types or elem.id in remaining_elems or elem in remaining_elems:
            lattice_analysis.append(merged_elems)
            merged_elems = []
            lattice_analysis.append([elem])
        elif elem.__class__ == Edge and ((Bend in remaining_types) or (SBend in remaining_types) or (RBend in remaining_types)):
            continue
        else:
            merged_elems.append(elem)
    if len(merged_elems) != 0:
        lattice_analysis.append(merged_elems)

    seq = []
    E = init_energy
    for elem_list in lattice_analysis:
        if len(elem_list) == 1:
            E += elem_list[0].transfer_map.delta_e
            seq.append(elem_list[0])
        elif len(elem_list) == 0:
            continue
        else:
            delta_e = np.sum([elem.transfer_map.delta_e for elem in elem_list])
            lattice = MagneticLattice(elem_list, method=lat.method)
            R = lattice_transfer_map(lattice, energy=E)
            m = Matrix()
            m.r = lattice.R
            m.t = lattice.T
            m.b = lattice.B
            m.l = lattice.totalLen
            m.delta_e = delta_e
            E += delta_e
            seq.append(m)

    new_lat = MagneticLattice(seq, method=lat.method)
    _logger.debug("element numbers after: " + str(len(new_lat.sequence)))
    return new_lat


flatten = lambda *n: (e for a in n
                      for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))


class MagneticLattice:
    """
    sequence - list of the elements,
    start - first element of lattice. If None, then lattice starts from first element of sequence,
    stop - last element of lattice. If None, then lattice stops by last element of sequence,
    method - A dictionary that contains the method of the tracking. If nothing is set TransferMap will
    be used as the global method for all elements. Setting TransferMaps for specific elements is also possible.
    Notes: If the elements doesn't support the defined transfer map, the default transfer map will be used, which is
    defined in the specific element class.
    For example:
        {"global": SecondTM, "Undulator": UndulatorTestTM }
        Sets for all elements SecondTM as transfer map, expect for the Undulator elements.
    """
    def __init__(self, sequence, start=None, stop=None, method={'global': TransferMap}):
        self.sequence = list(flatten(sequence))
        self.method = method

        try:
            if start is not None:
                id1 = self.sequence.index(start)
            else:
                id1 = 0
            if stop is not None:
                id2 = self.sequence.index(stop) + 1
                self.sequence = self.sequence[id1:id2]
            else:
                self.sequence = self.sequence[id1:]
        except:
            print('cannot construct sequence, element not found')
            raise

        # create transfer map and calculate lattice length
        self.totalLen = 0
        if not EdgeUtil.check(self):
            EdgeUtil.add(self)

        if not CouplerKickUtil.check(self):
            CouplerKickUtil.add(self)

        self.update_transfer_maps()

        self.__hash__ = {}
        for e in self.sequence:
            self.__hash__[e] = e

    def __getitem__(self, el):
        try:
            return self.__hash__[el]
        except:
            return None

    def update_transfer_maps(self):
        self.totalLen = 0
        for i, element in enumerate(self.sequence):
            if element.__class__ == Undulator:
                if element.field_file is not None:
                    element.l = element.field_map.l * element.field_map.field_file_rep
                    if element.field_map.units == "mm":
                        element.l = element.l*0.001
            self.totalLen += element.l

            if element.__class__ == Edge:

                self.update_endings(lat_index=i, element=element, body_elements=(Bend, RBend, SBend), element_util=EdgeUtil)

            if element.__class__ == CouplerKick:
                self.update_endings(lat_index=i, element=element, body_elements=(Cavity, ), element_util=CouplerKickUtil)

            element.create_tm(self.method)

            _logger.debug("update: " + element.transfer_map.__class__.__name__)
            if 'pulse' in element.__dict__: element.transfer_map.pulse = element.pulse
        return self

    def update_endings(self, lat_index, element, body_elements, element_util):

        if element_util.suffix_1 in element.id:
            body = self.sequence[lat_index + 1]
            if body.__class__ not in body_elements:
                body = self.sequence[lat_index - 1]
                _logger.debug("Backtracking?")
            element_util.update_first(element, body)
        elif element_util.suffix_2 in element.id:
            body = self.sequence[lat_index - 1]
            if body.__class__ not in body_elements:
                body = self.sequence[lat_index + 1]
                _logger.debug("Backtracking?")
            element_util.update_last(element, body)
        else:
            _logger.error(element.__class__.__name__ + " is not updated. Use standard function to create and update MagneticLattice")

    def __str__(self):
        line = "LATTICE: length = " + str(self.totalLen) + " m \n"
        for e in self.sequence:
            line += "{0:15} length: {1:5.2f}      id: {2:10}\n".format(e.__class__.__name__, e.l, e.id)
        return line

    def find_indices(self, element):
        indx_elem = np.where([i.__class__ == element for i in self.sequence])[0]
        return indx_elem

    def find_drifts(self):
        drift_lengs = []
        drifts = []
        for elem in self.sequence:
            if elem.__class__ == Drift:
                elem_l = np.around(elem.l, decimals=6)
                if elem_l not in drift_lengs:
                    drifts.append(elem)
                    drift_lengs.append(elem_l)
        return drifts

    def _create_var_name(self, objects):
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

    def _get_elements(self):
        elements = []
        for element in self.sequence:
            element_type = element.__class__.__name__
            # TODO: Edge and CouplerKick are not elements.
            if element_type in ('Edge', "CouplerKick"):
                continue
            if element not in elements:
                elements.append(element)
        elements = self._create_var_name(elements)
        return elements

    def _find_obj_and_create_name(self, types):
        objects = self._find_objects(types=types)
        objects = self._create_var_name(objects)
        return objects

    def _find_objects(self, types):
        """
        Function finds objects by types and adds it to list if object is unique.
        :param types: types of the Elements
        :return: list of elements
        """
        obj_id = []
        objs = []
        for elem in self.sequence:
            if elem.__class__ in types:
                if id(elem) not in obj_id:
                    objs.append(elem)
                    obj_id.append(id(elem))

        return objs

    def _write_power_supply_id(self, lines=[]):
        quads = self._find_obj_and_create_name(types=[Quadrupole])
        sexts = self._find_obj_and_create_name(types=[Sextupole])
        octs = self._find_obj_and_create_name(types=[Octupole])
        cavs = self._find_obj_and_create_name(types=[Cavity])
        bends = self._find_obj_and_create_name(types=[Bend, RBend, SBend])

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

        elements_dict = {}
        for element in elements:
            element_type = element.__class__.__name__

            if element_type not in elements_dict:
                elements_dict[element_type] = []

            elements_dict[element_type].append(element)

        return elements_dict

    def _print_elements(self, elements_dict):

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
                    string = element.element_def_string()
                    lines.append(string)

        # print remaining unordered elements
        for type in unordered_dict:

            lines.append('\n# ' + type + 's\n')

            for element in unordered_dict[type]:
                string = element.element_def_string()
                lines.append(string)

        # delete new line symbol from the first line
        if lines != []:
            lines[0] = lines[0][1:]
        return lines

    def elements2input(self):
        elements = self._get_elements()
        elements_dict = self._sort_elements(elements)
        lines = self._print_elements(elements_dict)
        return lines

    def cell2input(self, split=False):

        lines = []
        names = []
        for elem in self.sequence:
            if elem.__class__ not in (Edge, CouplerKick):
                names.append(elem.name)

        new_names = []
        for i, name in enumerate(names):
            if split and i % 10 == 9:
                new_names.append('\n' + name)
            else:
                new_names.append(name)

        lines.append('cell = (' + ', '.join(new_names) + ')')

        return lines

    def lat2input(self, tws0=None):
        """
        returns python input string for the self in the lat object
        """

        lines = ['from ocelot import * \n']

        # prepare initial Twiss parameters
        if tws0 is not None and isinstance(tws0, Twiss):
            lines.append('\n#Initial Twiss parameters\n')
            lines.extend(tws0.twiss2input())

        # prepare elements list
        lines.append('\n')
        lines.extend(self.elements2input())

        # prepare cell list
        lines.append('\n# Lattice \n')
        lines.extend(self.cell2input(True))

        lines.append('\n')

        return lines

    def rem_drifts(self):
        drifts = {}
        for i, elem in enumerate(self.sequence):
            if elem.__class__ == Drift:
                if not (elem.l in drifts.keys()):
                    drifts[elem.l] = elem
                else:
                    self.sequence[i] = drifts[elem.l]

    def write_lattice(self, tws0=None, file_name="lattice.py", remove_rep_drifts=True, power_supply=False):
        """
        saves lattice as python input file
        file_name - name of the file
        remove_rep_drifts - if True, remove the drifts with the same lengths from the lattice drifts definition
        """
        if remove_rep_drifts:
            self.rem_drifts()

        lines = self.lat2input(tws0=tws0)

        if power_supply:
            lines = self._write_power_supply_id(lines=lines)

        f = open(file_name, 'w')
        f.writelines(lines)
        f.close()


class EndElements:
    suffix_1 = "_1"
    suffix_2 = "_2"

    @staticmethod
    def check(lattice):
        pass

    @staticmethod
    def add(lattice):
        pass

    @staticmethod
    def update_first(end_element, body):
        pass

    @staticmethod
    def update_last(end_element, body):
        pass


class EdgeUtil(EndElements):
    suffix_1 = "_e1"
    suffix_2 = "_e2"

    def __init__(self):
        super(EdgeUtil).__init__(self)

    @staticmethod
    def name_suffix_1():
        return "_1"

    @staticmethod
    def name_suffix_2():
        return "_2"

    @staticmethod
    def check(lattice):
        """
        if there are edges on the ends of dipoles return True, else False
        """
        if len(lattice.sequence) < 3:
            return False
        for i in range(len(lattice.sequence)-2):
            prob_edge1 = lattice.sequence[i]
            elem = lattice.sequence[i+1]
            prob_edge2 = lattice.sequence[i+2]
            if elem.__class__ in (SBend, RBend, Bend):  # , "hcor", "vcor"
                if prob_edge1.__class__ != Edge and prob_edge2.__class__ != Edge:
                    return False
        return True

    @staticmethod
    def add(lattice):
        n = 0
        for i in range(len(lattice.sequence)):
            elem = lattice.sequence[n]
            if elem.__class__ in (SBend, RBend, Bend) and elem.l != 0.:  # , "hcor", "vcor"

                e_name = elem.id

                if elem.id is None:
                    e_name = "b_" + str(i)

                e1 = Edge(l=elem.l, angle=elem.angle, k1=elem.k1, edge=elem.e1, tilt=elem.tilt,
                          dx=elem.dx, dy=elem.dy, h_pole=elem.h_pole1, gap=elem.gap, fint=elem.fint, pos=1,
                          eid=e_name + EdgeUtil.suffix_1)

                lattice.sequence.insert(n, e1)

                e2 = Edge(l=elem.l, angle=elem.angle, k1=elem.k1, edge=elem.e2, tilt=elem.tilt,
                          dx=elem.dx, dy=elem.dy, h_pole=elem.h_pole2, gap=elem.gap, fint=elem.fintx, pos=2,
                          eid=e_name + EdgeUtil.suffix_2)

                lattice.sequence.insert(n+2, e2)
                n += 2
            n += 1

    @staticmethod
    def update_first(edge, bend):
        if bend.l != 0.:
            edge.h = bend.angle/bend.l
        else:
            edge.h = 0
        edge.l = 0.
        edge.angle = bend.angle
        edge.k1 = bend.k1
        edge.edge = bend.e1
        edge.tilt = bend.tilt
        edge.dx = bend.dx
        edge.dy = bend.dy
        edge.h_pole = bend.h_pole1
        edge.gap = bend.gap
        edge.fint = bend.fint
        edge.pos = 1

    @staticmethod
    def update_last(edge, bend):
        if bend.l != 0.:
            edge.h = bend.angle/bend.l
        else:
            edge.h = 0
        edge.l = 0.
        edge.angle = bend.angle
        edge.k1 = bend.k1
        edge.edge = bend.e2
        edge.tilt = bend.tilt
        edge.dx = bend.dx
        edge.dy = bend.dy
        edge.h_pole = bend.h_pole2
        edge.gap = bend.gap
        edge.fint = bend.fintx
        edge.pos = 2


class CouplerKickUtil(EndElements):
    suffix_1 = "_ck1"
    suffix_2 = "_ck2"

    def __init__(self):
        super(CouplerKickUtil).__init__(self)

    @staticmethod
    def check(lattice):
        """
        if there are CouplerKicks on the ends of Cavities return True, else False
        """
        if len(lattice.sequence) < 3:
            return False
        for i in range(len(lattice.sequence)-2):
            prob_ck1 = lattice.sequence[i]
            elem = lattice.sequence[i+1]
            prob_ck2 = lattice.sequence[i+2]
            if elem.__class__ in (Cavity,):
                if prob_ck1.__class__ != CouplerKick and prob_ck2.__class__ != CouplerKick:
                    return False
        return True

    @staticmethod
    def add(lattice):
        n = 0
        for i in range(len(lattice.sequence)):
            elem = lattice.sequence[n]
            if elem.__class__ in (Cavity,) and elem.l != 0.:

                e_name = elem.id

                if elem.id is None:
                    e_name = "cav_" + str(i)

                e1 = CouplerKick(v=elem.v, phi=elem.phi, freq=elem.freq, vx=elem.vx_up, vy=elem.vy_up,
                                 vxx=elem.vxx_up, vxy=elem.vxy_up, eid=e_name + CouplerKickUtil.suffix_1)

                lattice.sequence.insert(n, e1)

                e2 = CouplerKick(v=elem.v, phi=elem.phi, freq=elem.freq, vx=elem.vx_down, vy=elem.vy_down,
                                 vxx=elem.vxx_down, vxy=elem.vxy_down, eid=e_name + CouplerKickUtil.suffix_2)

                lattice.sequence.insert(n+2, e2)
                n += 2
            n += 1

    @staticmethod
    def update_first(ckick, cavity):
        ckick.v = cavity.v
        ckick.phi = cavity.phi
        ckick.freq = cavity.freq
        ckick.vx = cavity.vx_up
        ckick.vy = cavity.vy_up
        ckick.vxx = cavity.vxx_up
        ckick.vxy = cavity.vxy_up

    @staticmethod
    def update_last(ckick, cavity):
        ckick.v = cavity.v
        ckick.phi = cavity.phi
        ckick.freq = cavity.freq
        ckick.vx = cavity.vx_down
        ckick.vy = cavity.vy_down
        ckick.vxx = cavity.vxx_down
        ckick.vxy = cavity.vxy_down


def merge_drifts(cell):
    """
    Merge neighboring Drifts in one Drift

    :param cell: list of element
    :return: new list of elements
    """
    new_cell = []
    L = 0.
    for elem in cell:

        if elem.__class__ in [Drift, UnknownElement]:
            L += elem.l
        else:
            if L != 0:
                new_elem = Drift(l=L)
                new_cell.append(new_elem)
            new_cell.append(elem)
            L = 0.
    if L != 0: new_cell.append(Drift(l=L))
    print("Merge drift -> Element numbers: before -> after: ", len(cell), "->", len(new_cell))
    return new_cell


def exclude_zero_length_element(cell, elem_type=[UnknownElement, Marker], except_elems=[]):
    """
    Exclude zero length elements some types in elem_type

    :param cell: list, sequence of elements
    :param elem_type: list, types of Elements which should be excluded
    :param except_elems: list, except elements
    :return: list, new sequence of elements
    """
    new_cell = []
    for elem in cell:
        if elem.__class__ in elem_type and elem.l == 0 and elem not in except_elems:
            continue
        new_cell.append(elem)
    print("Exclude elements -> Element numbers: before -> after: ", len(cell), "->", len(new_cell))
    return new_cell