from ocelot.cpbd.optics import MethodTM
from ocelot.cpbd.elements import *
from ocelot.common.logging import *
import numpy as np
from copy import deepcopy
logger = Logger()

def lattice_format_converter(elements):
    """
    :param elements: lattice in the format: [[elem1, center_pos1], [elem2, center_pos2], [elem3, center_pos3], ... ]
    :return: lattice in the format: [elem1, drift1, elem2, drift2, elem3, drift3, ...]
    """
    cell = []
    drift_num = 0
    s_pos = 0.0
    for element in elements:
        element_start = element[1] - element[0].l / 2.0
        if element_start < s_pos - 1.0e-14:                 # 1.0e-14 is used as crutch for precision of float
            if element[0].l == 0.0:
                if s_pos - element_start > 1.0e-2:
                    print("************** WARNING! Element " + element[0].id + " was deleted")
                    continue
                element[1] = s_pos
                print("************** WARNING! Element " + element[0].id + " was moved from " + str(element_start) + " to " + str(s_pos))
            else:
                dl = element[0].l / 2.0  - element[1] + s_pos
                if cell[-1].__class__ == Marker and cell[-2].__class__ == Drift and cell[-2].l > dl:
                    cell[-2].l -= dl
                    print("************** WARNING! Element " + cell[-1].id + " was deleted")
                    cell.pop()
                else:
                    print("************** ERROR! Element " + element[0].id + " has bed position")
                    exit()

        if element_start > s_pos + 1.0e-12:                 # 1.0e-12 is used as crutch for precision of float
            drift_num += 1
            drift_l = round(element_start - s_pos, 10)      # round() is used as crutch for precision of float
            drift_eid = 'D_' + str(drift_num)
            cell.append(Drift(l=drift_l, eid=drift_eid))

        cell.append(element[0])
        s_pos = element[1] + element[0].l / 2.0
    return tuple(cell)


def matrix2Matrix(elem, R):
    """
    Creation from R-matrix element Matrix
    :param elem: Matrix element
    :param R: R-matrix
    :return: Matrix Element with defined r-elements
    """
    elem.rm11, elem.rm21, elem.rm31 = R[0, 0], R[1, 0], R[2, 0]
    elem.rm12, elem.rm22, elem.rm32 = R[0, 1], R[1, 1], R[2, 1]
    elem.rm13, elem.rm23, elem.rm33 = R[0, 2], R[1, 2], R[2, 2]
    elem.rm14, elem.rm24, elem.rm34 = R[0, 3], R[1, 3], R[2, 3]
    elem.rm15, elem.rm25, elem.rm35 = R[0, 4], R[1, 4], R[2, 4]
    elem.rm16, elem.rm26, elem.rm36 = R[0, 5], R[1, 5], R[2, 5]
    elem.rm41, elem.rm51, elem.rm61 = R[3, 0], R[4, 0], R[5, 0]
    elem.rm42, elem.rm52, elem.rm62 = R[3, 1], R[4, 1], R[5, 1]
    elem.rm43, elem.rm53, elem.rm63 = R[3, 2], R[4, 2], R[5, 2]
    elem.rm44, elem.rm54, elem.rm64 = R[3, 3], R[4, 3], R[5, 3]
    elem.rm45, elem.rm55, elem.rm65 = R[3, 4], R[4, 4], R[5, 4]
    elem.rm46, elem.rm56, elem.rm66 = R[3, 5], R[4, 5], R[5, 5]


def shrinker(lat, remaining_types, init_energy=0.):
    """

    :param lat: MagneticLattice
    :param remaining_types: the type of the elements which needed to be untoched
                            others will be "compress" to Matrix element
                            e.g. [Monitor, Quadrupole, Bend, Hcor]
    :param init_energy: initial energy
    :return: New MagneticLattice
    """
    print("element numbers before: ", len(lat.sequence))
    R_seq = []
    new_elem = Matrix()
    _r = np.eye(6)
    _b = np.zeros(6)
    length = 0.
    count = 0
    E = init_energy
    for elem in lat.sequence:
        if elem.__class__ in remaining_types:
            if count != 0:
                new_elem.l = length
                matrix2Matrix(new_elem, _r)
                R_seq.append(new_elem)
                new_elem = Matrix()
                length = 0.
                _r = np.eye(6)
                _b = np.zeros(6)

                count = 0
            elem.transfer_map._r = elem.transfer_map.R(E)
            elem.transfer_map._b = elem.transfer_map.B(E)
            R_seq.append(elem)
        elif elem.__class__ == Edge and ((Bend in remaining_types) or (SBend in remaining_types) or (RBend in remaining_types)):
            continue
        else:
            _r = np.dot(elem.transfer_map.R(E), _r)
            _b =  np.dot(elem.transfer_map.R(E), _b) + elem.transfer_map.B(E)  #+dB #check
            length = elem.l + length
            new_elem.delta_e = elem.transfer_map.delta_e + new_elem.delta_e
            count += 1
        E += elem.transfer_map.delta_e
    if count != 0:
        matrix2Matrix(new_elem, _r)
        new_elem.l = length
        R_seq.append(new_elem)
    new_lat = MagneticLattice(R_seq)
    print("element numbers after: ", len(new_lat.sequence))
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
        #self.energy = energy
        self.sequence = list(flatten(sequence))
        self.method = method
        try:
            if start != None:
                id1 = self.sequence.index(start)
            else:
                id1 = 0
            if stop != None:
                id2 = self.sequence.index(stop) + 1
                self.sequence = self.sequence[id1:id2]
            else:
                self.sequence = self.sequence[id1:]
        except:
            print('cannot construct sequence, element not found')
            raise

        #self.transferMaps = {}
        # create transfer map and calculate lattice length
        self.totalLen = 0
        if not self.check_edges():
            self.add_edges()
        self.update_transfer_maps()

        self.__hash__ = {}
        #print 'creating hash'
        for e in self.sequence:
            #print e
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
        for i in range(len(self.sequence)-2):
            prob_edge1 = self.sequence[i]
            elem = self.sequence[i+1]
            prob_edge2 = self.sequence[i+2]
            if elem.__class__ in (SBend, RBend, Bend):  # , "hcor", "vcor"
                if prob_edge1.__class__ != Edge and prob_edge2.__class__ != Edge:
                    #print elem.type, prob_edge1.type, prob_edge2.type
                    return False
        return True

    def add_edges(self):
        n = 0
        for i in range(len(self.sequence)):
            elem = self.sequence[n]
            if elem.__class__ in (SBend, RBend, Bend) and elem.l != 0.:  # , "hcor", "vcor"

                e_name = elem.id

                if elem.id == None:
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
        edge.fint = bend.fintx
        edge.pos = 1

    def update_edge_e2(self, edge, bend):
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
        #E = self.energy
        self.totalLen = 0
        for i, element in enumerate(self.sequence):
            if element.__class__ == Undulator:
                if element.field_file != None:
                    element.l = element.field_map.l * element.field_map.field_file_rep
                    if element.field_map.units == "mm":
                        element.l = element.l*0.001
            self.totalLen += element.l
            #print(element.k1)
            if element.__class__ == Edge:

                if "_e1" in element.id:
                    bend = self.sequence[i+1]
                    if bend.__class__ not in (SBend, RBend, Bend):
                        bend = self.sequence[i - 1]
                        print("Backtracking?")
                    self.update_edge_e1(element, bend)
                elif "_e2" in element.id:
                    bend = self.sequence[i-1]
                    if bend.__class__ not in (SBend, RBend, Bend):
                        bend = self.sequence[i + 1]
                    self.update_edge_e2(element, bend)
                else:
                    print("EDGE is not updated. Use standard function to create and update MagneticLattice")
            element.transfer_map = self.method.create_tm(element)
            logger.debug("update: " + element.transfer_map.__class__.__name__)
            #print("update: ", element.transfer_map.__class__.__name__, element.l, element.id, element.transfer_map.R(0))
            if 'pulse' in element.__dict__: element.transfer_map.pulse = element.pulse
        return self

    def printElements(self):
        print('\nLattice\n')
        for e in self.sequence:
            print('-->',  e.id, '[', e.l, ']')

def merge_drifts(cell):
    """
    Merge neighboring Drifts in one Drift

    input: cell - list of element
    return: new_cell - new list of elements
    """

    new_cell = []
    L = 0.
    for elem in cell:

        if elem.__class__ in [Drift, UnknownElement]:
            L += elem.l
            #print(L)
        else:
            #print("new")
            if L != 0:
                new_elem = Drift(l=L)
                new_cell.append(new_elem)
            new_cell.append(elem)
            L = 0.
    if L != 0: new_cell.append(Drift(l=L))
    print("Merge drift -> Element numbers: before -> after: ", len(cell), "->", len(new_cell))
    return new_cell

def exclude_zero_length_element(cell, elem_type=[UnknownElement, Marker]):
    """
    Exclude zero length elements some types in elem_type
    input: cell
    return: new cell
    """
    new_cell = []
    for elem in cell:
        if elem.__class__ in elem_type and elem.l == 0:
            continue
        new_cell.append(elem)
    print("Exclude elements -> Element numbers: before -> after: ", len(cell), "->", len(new_cell))
    return new_cell