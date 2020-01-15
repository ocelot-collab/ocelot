from ocelot.cpbd.optics import MethodTM, lattice_transfer_map
from ocelot.cpbd.elements import *
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
        if element[0].__class__ == Edge:
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
    method = MethodTM() - method of the tracking.
    """
    def __init__(self, sequence, start=None, stop=None, method=MethodTM()):
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
        if not self.check_edges():
            self.add_edges()
        self.update_transfer_maps()

        self.__hash__ = {}
        for e in self.sequence:
            self.__hash__[e] = e

    def __getitem__(self, el):
        try:
            return self.__hash__[el]
        except:
            return None

    def check_edges(self):
        """
        if there are edges on the ends of dipoles return True, else False
        """
        if len(self.sequence) < 3:
            return False
        for i in range(len(self.sequence)-2):
            prob_edge1 = self.sequence[i]
            elem = self.sequence[i+1]
            prob_edge2 = self.sequence[i+2]
            if elem.__class__ in (SBend, RBend, Bend):  # , "hcor", "vcor"
                if prob_edge1.__class__ != Edge and prob_edge2.__class__ != Edge:
                    return False
        return True

    def add_edges(self):
        n = 0
        for i in range(len(self.sequence)):
            elem = self.sequence[n]
            if elem.__class__ in (SBend, RBend, Bend) and elem.l != 0.:  # , "hcor", "vcor"

                e_name = elem.id

                if elem.id is None:
                    e_name = "b_" + str(i)

                e1 = Edge(l=elem.l, angle=elem.angle, k1=elem.k1, edge=elem.e1, tilt=elem.tilt, dtilt=elem.dtilt,
                          dx=elem.dx, dy=elem.dy, h_pole=elem.h_pole1, gap=elem.gap, fint=elem.fint, pos=1,
                          eid=e_name + "_e1")

                self.sequence.insert(n, e1)

                e2 = Edge(l=elem.l, angle=elem.angle, k1=elem.k1, edge=elem.e2, tilt=elem.tilt, dtilt=elem.dtilt,
                          dx=elem.dx, dy=elem.dy, h_pole=elem.h_pole2, gap=elem.gap, fint=elem.fintx, pos=2,
                          eid=e_name + "_e2")

                self.sequence.insert(n+2, e2)
                n += 2
            n += 1

    def update_edge_e1(self, edge, bend):
        if bend.l != 0.:
            edge.h = bend.angle/bend.l
        else:
            edge.h = 0
        edge.l = 0.
        edge.angle = bend.angle
        edge.k1 = bend.k1
        edge.edge = bend.e1
        edge.tilt = bend.tilt
        edge.dtilt = bend.dtilt
        edge.dx = bend.dx
        edge.dy = bend.dy
        edge.h_pole = bend.h_pole1
        edge.gap = bend.gap
        edge.fint = bend.fint
        edge.pos = 1

    def update_edge_e2(self, edge, bend):
        if bend.l != 0.:
            edge.h = bend.angle/bend.l
        else:
            edge.h = 0
        edge.l = 0.
        edge.angle = bend.angle
        edge.k1 = bend.k1
        edge.edge = bend.e2
        edge.tilt = bend.tilt
        edge.dtilt = bend.dtilt
        edge.dx = bend.dx
        edge.dy = bend.dy
        edge.h_pole = bend.h_pole2
        edge.gap = bend.gap
        edge.fint = bend.fintx
        edge.pos = 2

    def update_transfer_maps(self):
        self.totalLen = 0
        for i, element in enumerate(self.sequence):
            if element.__class__ == Undulator:
                if element.field_file != None:
                    element.l = element.field_map.l * element.field_map.field_file_rep
                    if element.field_map.units == "mm":
                        element.l = element.l*0.001
            self.totalLen += element.l

            if element.__class__ == Edge:

                if "_e1" in element.id:
                    bend = self.sequence[i+1]
                    if bend.__class__ not in (SBend, RBend, Bend):
                        bend = self.sequence[i - 1]
                        _logger.debug("Backtracking?")
                    self.update_edge_e1(element, bend)
                elif "_e2" in element.id:
                    bend = self.sequence[i-1]
                    if bend.__class__ not in (SBend, RBend, Bend):
                        bend = self.sequence[i + 1]
                        _logger.debug("Backtracking?")
                    self.update_edge_e2(element, bend)
                else:
                    _logger.error("EDGE is not updated. Use standard function to create and update MagneticLattice")
            element.transfer_map = self.method.create_tm(element)
            _logger.debug("update: " + element.transfer_map.__class__.__name__)
            if 'pulse' in element.__dict__: element.transfer_map.pulse = element.pulse
        return self


    def __str__(self):
        line = "LATTICE: length = " + str(self.totalLen) + " m \n"
        for e in self.sequence:
            line += "{0:15} length: {1:5.2f}      id: {2:10}\n".format(e.__class__.__name__, e.l, e.id)
        return line

    def find_indices(self, element):
        indx_elem = np.where([i.__class__ == element for i in self.sequence])[0]
        return indx_elem


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