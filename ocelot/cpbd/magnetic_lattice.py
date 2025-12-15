from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.elements.marker import Marker
from ocelot.cpbd.elements.drift import Drift
from ocelot.cpbd.elements.monitor import Monitor
from ocelot.cpbd.elements.undulator import Undulator
from ocelot.cpbd.elements.unknown_element import UnknownElement
from ocelot.cpbd.elements.matrix import Matrix
from ocelot.cpbd.elements.bend import Bend
from ocelot.cpbd.elements.sbend import SBend
from ocelot.cpbd.elements.rbend import RBend
from ocelot.cpbd.latticeIO import LatticeIO
from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.optics import lattice_transfer_map, periodic_twiss
from ocelot.cpbd.tm_utils import transfer_maps_mult
from ocelot.common.globals import m_e_GeV
from ocelot.cpbd.beam import Twiss

import logging
import re
from collections import defaultdict
from typing import Mapping, Sequence, Tuple, Callable, Any, Generator, Iterator, Type, TypeVar, List
import numpy as np

_logger = logging.getLogger(__name__)

# Returned by insert_marker_by_type, insert_markers_by_name,
# insert_markers_by_predicate
MarkersInsertionReturnType = Mapping[Element, Sequence[Tuple[Marker, Marker]]]

# to solve typing issue in PyCharm
E = TypeVar('E', bound=OpticElement)


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
        if element_start < s_pos - 1.0e-14:  # 1.0e-14 is used as crutch for precision of float
            if element[0].l == 0.0:

                if s_pos - element_start > 1.0e-2:
                    print("************** WARNING! Element " + element[0].id + " was deleted")
                    continue
                element[1] = s_pos
                print("************** WARNING! Element " + element[0].id + " was moved from " + str(element_start) +
                      " to " + str(s_pos))
            else:
                dl = element[0].l / 2.0 - element[1] + s_pos
                if cell[-1].__class__ == Marker and cell[-2].__class__ == Drift and cell[-2].l > dl:
                    cell[-2].l -= dl
                    print("************** WARNING! Element " + cell[-1].id + " was deleted")
                    cell.pop()
                else:
                    print("************** ERROR! Element " + element[0].id + " has bad position (overlapping?)")
                    exit()

        if element_start > s_pos + 1.0e-12:  # 1.0e-12 is used as crutch for precision of float
            drift_num += 1
            drift_l = round(element_start - s_pos, 10)  # round() is used as crutch for precision of float
            drift_eid = 'D_' + str(drift_num)
            cell.append(Drift(l=drift_l, eid=drift_eid))

        cell.append(element[0])
        s_pos = element[1] + element[0].l / 2.0
    return tuple(cell)


def merger(lat, remaining_types=None, remaining_elems=None, init_energy=0.):
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
    if remaining_elems is None:
        remaining_elems = []
    if remaining_types is None:
        remaining_types = []
    _logger.debug("element numbers before: " + str(len(lat.sequence)))
    lattice_analysis = []
    merged_elems = []
    for elem in lat.sequence:
        if elem.__class__ in remaining_types or elem.id in remaining_elems or elem in remaining_elems:
            lattice_analysis.append(merged_elems)
            merged_elems = []
            lattice_analysis.append([elem])
        else:
            merged_elems.append(elem)
    if len(merged_elems) != 0:
        lattice_analysis.append(merged_elems)

    seq = []
    E = init_energy
    for elem_list in lattice_analysis:
        if len(elem_list) == 1:
            for tm in elem_list[0].tms:
                E += tm.get_delta_e()
            seq.append(elem_list[0])
        elif len(elem_list) == 0:
            continue
        else:
            magnetic_elems = [elem for elem in elem_list if elem.__class__ not in [Marker, Drift, Monitor]]
            if len(magnetic_elems) == 0:
                total_len = np.sum([elem.l for elem in elem_list])
                d = Drift(l=total_len)
                seq.append(d)
            else:
                delta_e = np.sum([tm.get_delta_e() for elem in elem_list for tm in elem.tms])
                lattice = MagneticLattice(elem_list, method=lat.method)
                m = Matrix()
                m.b, m.r, m.t = lattice.transfer_maps(energy=E)
                m.l = lattice.totalLen
                m.delta_e = delta_e
                E += delta_e
                seq.append(m)

    new_lat = MagneticLattice(seq, method=lat.method)
    _logger.debug("element numbers after: " + str(len(new_lat.sequence)))
    return new_lat


def flatten(iterable: Iterator[Any]) -> Generator[Any, None, None]:
    """Flatten arbitrarily nested iterable.
    Special case for strings that avoids infinite recursion.  Non iterables passed
    as arguments are yielded.

    :param iterable: Any iterable to be flattened.
    :raises: RecursionError

    """

    def _flatten(iterable):
        try:
            for item in iterable:
                try:
                    iter(item)
                except TypeError:
                    yield item
                else:
                    yield from _flatten(item)
        except TypeError:  # If iterable isn't actually iterable, then yield it.
            yield iterable

    try:
        yield from _flatten(iterable)
    except RecursionError:
        raise RecursionError("Maximum recursion reached.  Possibly trying"
                             " to flatten an infinitely nested iterable.")


class MagneticLattice:
    """
    Represents a magnetic lattice, which is a sequence of elements defining a beamline.

    Args:
        sequence (list): A list of elements forming the lattice.
        start (Element, optional): The first element of the lattice. If `None`, the lattice starts with the
            first element of the sequence. Defaults to `None`.
        stop (Element, optional): The last element of the lattice. If `None`, the lattice stops with the
            last element of the sequence. Defaults to `None`.
        method (dict, optional): A dictionary specifying the tracking method for the lattice. If no method is provided,
            `TransferMap` is used as the global default for all elements. Specific methods for individual elements
            can also be set. Defaults to `None`.

            Example:
                ```python
                from ocelot import *

                method = {"global": TransferMap} # default first order transfer map
                lat = MagneticLattice(cell, method=method)
                # or
                method = {"global": SecondTM, Octupole: KickTM, Undulator: RungeKuttaTM}
                lat = MagneticLattice(cell, method=method)
                ```

            In this example:
            - Sets `SecondTM` (second order transfer maps) as the global transfer map for all elements.
            - Assigns `KickTM` specifically for `Octupole` elements.
            - Assigns `RungeKuttaTM` specifically for `Undulator` elements.

    Notes:
        - If an element does not support the specified transfer map, the default transfer map defined in
          the element class is used.
        - For more details, refer to Section 7 of the tutorial:
          [Small Useful Features](https://nbviewer.org/github/ocelot-collab/ocelot/blob/dev/demos/ipython_tutorials/small_useful_features.ipynb).
    """

    def __init__(self, sequence, start: E = None, stop: E = None, method=None):
        if method is None:
            method = {'global': TransferMap}
        if isinstance(method, dict):
            self.method = method
        else:
            self.method = method.to_dict()

        self.sequence = list(flatten(sequence))
        self.sequence = self.get_sequence_part(start, stop)

        self.update_transfer_maps()

        self._elem_map = {e: e for e in self.sequence}

    def __getitem__(self, el):
        # Let KeyError propagate if not found
        return self._elem_map[el]

    def __iter__(self):
        return iter(self.sequence)

    @property
    def totalLen(self):
        return sum([e.l for e in self.sequence])

    def get_sequence_part(self, start: E | None, stop: E | None):
        """
        Return a sub-sequence of elements from `start` to `stop` (inclusive).
        Raises a ValueError if either element is not found, or if `stop` precedes `start`.

        Parameters
        ----------
        start : Element or None
            Element where the subsequence starts. If None, starts from the beginning.
        stop : Element or None
            Element where the subsequence ends. If None, continues to the end.

        Returns
        -------
        list
            Sub-list of elements from start to stop (inclusive).
        """
        seq = self.sequence

        try:
            id1 = seq.index(start) if start is not None else 0
        except ValueError:
            raise ValueError(f"Start element {getattr(start, 'id', start)} not found in lattice.")

        try:
            id2 = seq.index(stop) if stop is not None else len(seq) - 1
        except ValueError:
            raise ValueError(f"Stop element {getattr(stop, 'id', stop)} not found in lattice.")

        # --- new check: ensure stop is after start ---
        if id2 < id1:
            s_id = getattr(start, "id", start)
            e_id = getattr(stop, "id", stop)
            raise ValueError(
                f"Invalid range: stop element '{e_id}' appears before start element '{s_id}' "
                f"in the lattice sequence (indices {id1}>{id2})."
            )

        return seq[id1:id2 + 1]  # inclusive

    def update_transfer_maps(self):
        for i, element in enumerate(self.sequence):
            # TODO: This belongs to the Element Undulator
            if element.__class__ == Undulator:
                if element.field_file is not None:
                    element.l = element.field_map.l * element.field_map.field_file_rep
                    if element.field_map.units == "mm":
                        element.l = element.l * 0.001

            tm_class_type = self.method.get(element.__class__)
            if tm_class_type:
                element.set_tm(tm_class_type)
            else:
                tm_class_type = self.method.get('global')
                if tm_class_type:
                    element.set_tm(tm_class_type)

            _logger.debug(f"update: {','.join([tm.__class__.__name__ for tm in element.tms])}")
        return self

    def __str__(self):
        line = "LATTICE: length = " + str(self.totalLen) + " m \n"
        for e in self.sequence:
            line += "{0:15} length: {1:5.2f}      id: {2:10}\n".format(e.__class__.__name__, e.l, e.id)
        return line

    def find_indices(self, element):
        """Find index by class type: argument should be a class"""

        return self.find_indices_by_predicate(lambda elem: isinstance(elem, element))

    def find_indices_by_predicate(self, predicate: Callable[[Element], bool]) -> List[int]:
        """Get indices using some callable function.  Function should
        take one argument, an element of the sequence, and return
        either True or False.
        """
        return [index for (index, element) in enumerate(self.sequence) if predicate(element)]

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

    def rem_drifts(self):
        drifts = {}
        for i, elem in enumerate(self.sequence):
            if elem.__class__ == Drift:
                if not (elem.l in drifts.keys()):
                    drifts[elem.l] = elem
                else:
                    self.sequence[i] = drifts[elem.l]

    def save_as_py_file(self, file_name: str, tws0=None, remove_rep_drifts=True, power_supply=False):
        """
        Saves the lattice to a Python file.

        Args:
            file_name (str): The path and name of the Python file where the lattice will be stored.
            tws0 (Twiss, optional): A `Twiss` object. If provided, the Twiss parameters will be printed at the beginning
                of the lattice file. Defaults to `None`.
            remove_rep_drifts (bool, optional): If `True`, removes repeated drift elements from the lattice.
                Defaults to `True`.
            power_supply (bool, optional): If `True`, writes the power supply IDs into the file. Defaults to `False`.

        Returns:
            None
        """
        LatticeIO.save_lattice(self, tws0=tws0, file_name=file_name, remove_rep_drifts=remove_rep_drifts,
                               power_supply=power_supply)

    def transfer_maps(self, energy, output_at_each_step: bool = False, start: E = None,
                      stop: E = None):
        """
        Calculates the zero, first- and second-order transfer maps for the lattice.

        This method computes the full transfer maps (B, R, T) across the lattice or a specified segment.
        If `output_at_each_step` is set to `True`, it also returns the intermediate transfer maps after each element.

        :param start: (Element, optional): Element from which to start the transfer map calculation. Defaults to beginning.
        :param stop: Element at which to stop the transfer map calculation. Defaults to end.
        :param energy: the initial electron beam energy [GeV]
        :param output_at_each_step: [Bs], [Rs], [Ts], [S] return three list of matrices [Bs], [Rs], [Ts] after each element in the line
                                    and [S] position at the end of each transfer map.
        :return: B, R, T - matrices OR [Bs], [Rs], [Ts], [S] if output_at_each_step = True
        """
        sequence = self.get_sequence_part(start, stop)
        Ra = np.eye(6)
        Ta = np.zeros((6, 6, 6))
        Ba = np.zeros((6, 1))
        Bs, Rs, Ts, S = [], [], [], []
        s_pos = 0.
        E = energy
        for elem in sequence:
            for Rb, Bb, Tb, tm in zip(elem.R(E), elem.B(E), elem.T(E), elem.tms):
                Ba, Ra, Ta = transfer_maps_mult(Ba, Ra, Ta, Bb, Rb, Tb)
                E += tm.get_delta_e()
                Bs.append(Ba)
                Rs.append(Ra)
                Ts.append(Ta)
                s_pos += tm.length
                S.append(s_pos)
        if output_at_each_step:
            return Bs, Rs, Ts, S
        return Ba, Ra, Ta

    def survey(self, X0=0, Y0=0, Z0=0, theta0=0, phi0=0, psi0=0):
        """
        Calculates the 3D survey using the exact MAD-8 recursive vector/matrix method.

        :param X0, Y0, Z0: Initial global coordinates [m].
        :param theta0: Initial azimuth angle [rad].
        :param phi0: Initial elevation angle [rad].
        :param psi0: Initial roll angle [rad].

        :return: (mid_survey_data, end_survey_data)
                 Two lists of dictionaries containing survey data.
                 - mid_survey_data: Calculated at the geometric center of elements (Best for Tables).
                 - end_survey_data: Calculated at the exit of elements (Best for Plotting).
        """
        V = np.array([X0, Y0, Z0])
        S_pos = 0.0

        # Initial Matrices
        M_theta = np.array([[np.cos(theta0), 0, np.sin(theta0)], [0, 1, 0], [-np.sin(theta0), 0, np.cos(theta0)]])
        M_phi = np.array([[1, 0, 0], [0, np.cos(phi0), np.sin(phi0)], [0, -np.sin(phi0), np.cos(phi0)]])
        from ocelot.common.math_op import get_tilt_matrix
        M_psi = get_tilt_matrix(psi0)

        W = M_theta @ M_phi @ M_psi

        mid_survey_data = []
        end_survey_data = []

        # Helper: Now accepts v_start and v_end vectors
        def get_survey_data(w_mat, w_start, v_start, v_end, s_val, el):
            # 1. Scalars for Excel/Pandas
            xpd, ypd, zpd = w_mat[0, 2], w_mat[1, 2], w_mat[2, 2]
            val_phi = np.clip(ypd, -1.0, 1.0)
            phi = np.arcsin(val_phi)
            theta = np.arctan2(xpd, zpd)
            psi_angle = np.arctan2(w_mat[1, 0], w_mat[1, 1])

            length = getattr(el, 'l', 0.0) if el else 0.0
            tilt = getattr(el, 'tilt', 0.0) if el else 0.0

            return {
                # --- Excel/Table Data (Scalars) ---
                "LENGTH": length,
                "TILT": tilt,
                "S": s_val,
                "X": v_end[0], "Y": v_end[1], "Z": v_end[2],  # Current position
                "THETA": theta, "PHI": phi, "PSI": psi_angle,
                "XPD": xpd, "YPD": ypd, "ZPD": zpd,

                # --- Layout/Plotting Data (Vectors/Matrices) ---
                "W": w_mat.copy(),  # end or mid
                "W_start": w_start.copy(),
                "r_start": v_start.copy(),  # Vector [x,y,z] at element start
                "r_end": v_end.copy(),  # Vector [x,y,z] at element end

                "element": el
            }

        # Add Start Point
        # Start/End are same for the zero-length marker at 0.0
        start_point_data = get_survey_data(W, W, V, V, S_pos, None)
        mid_survey_data.append(start_point_data)
        end_survey_data.append(start_point_data)

        for elem in self.sequence:
            R_end, S_end, R_mid, S_mid = elem.get_transfer_geometry()

            W_start = W.copy()
            V_start = V.copy()

            # --- 1. MIDPOINT (For Tables) ---
            V_mid_global = V_start + W_start @ R_mid
            W_mid_global = W_start @ S_mid

            L = getattr(elem, 'l', 0.0)
            S_mid_val = S_pos + L / 2.0

            # For midpoint data, 'r_end' is the center, 'r_start' is entry
            mid_survey_data.append(get_survey_data(W_mid_global,W_start, V_start, V_mid_global, S_mid_val, elem))

            # --- 2. ENDPOINT (For Layouts/Connectivity) ---
            V_end_global = V_start + W_start @ R_end
            W_end_global = W_start @ S_end
            S_pos += L

            # For endpoint data, we have the full element segment
            end_survey_data.append(get_survey_data(W_end_global, W_start, V_start, V_end_global, S_pos, elem))

            # Update state
            V = V_end_global
            W = W_end_global

        return mid_survey_data, end_survey_data

    def survey_longlist(self, X0=0, Y0=0, Z0=0, theta0=0, phi0=0, psi0=0):
        """
        Calculates survey and converts angles/coordinates to the 'LongList' engineering convention.

        Mapping:
        - THETA (Azimuth)   -> Becomes PHI (Elevation)
        - PHI (Elevation)   -> Becomes -THETA (-Azimuth)
        - PSI (Roll)        -> Renamed to CHI
        - XPD               -> Becomes YPD (Consistent with angle swap)
        - YPD               -> Becomes -XPD
        """
        # 1. Get the raw physics data (Assuming this returns dicts with 'PSI', 'THETA', 'PHI')
        mid_phys, end_phys = self.survey(X0, Y0, Z0, theta0, phi0, psi0)

        # 2. Define the transformation logic
        def to_longlist_convention(data_list):
            new_list = []
            for item in data_list:
                # Shallow copy to protect original data
                new_item = item.copy()

                # --- 1. Rename PSI to CHI ---
                # .pop('PSI') gets the value AND removes 'PSI' from new_item
                # This ensures you don't have duplicate keys.
                new_item['CHI'] = new_item.pop('PSI', 0.0)

                # --- 2. Swap Angles ---
                mad_theta = item['THETA']
                mad_phi = item['PHI']

                new_item['THETA'] = mad_phi
                new_item['PHI'] = -mad_theta

                # --- 3. Swap Direction Cosines (Consistency) ---
                # If we swap angles, we must swap the unit vectors too.
                # If Theta -> Phi (X -> Y), then XPD -> YPD.
                # If Phi -> -Theta (Y -> -X), then YPD -> -XPD.
                mad_xpd = item['XPD']
                mad_ypd = item['YPD']

                new_item['XPD'] = mad_ypd
                new_item['YPD'] = -mad_xpd

                new_list.append(new_item)
            return new_list

        return to_longlist_convention(mid_phys), to_longlist_convention(end_phys)


    def print_sequence(self, start: E = None, stop: E = None):
        sequence = self.get_sequence_part(start, stop)
        lines = ["{:<17} {:<15} {:<10} {:<10}".format('id', 'length', 'start', 'end') ]
        s = 0.
        for elem in sequence:
            line = "{:<17} {:<15} {:<10} {:<10}".format(elem.id, np.round(elem.l, 4), np.round(s, 4), np.round(s + elem.l, 4))
            s += elem.l
            lines.append(line)
        return lines

    def periodic_twiss(self, tws=None):
        tws = Twiss(tws)
        R = self.transfer_maps(energy=tws.E)[1]
        tw_periodic = periodic_twiss(tws, R)
        return tw_periodic


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
    if L != 0:
        new_cell.append(Drift(l=L))
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


def insert_markers_by_name(sequence, string: str, regex=False,
                           before=True, after=True) -> MarkersInsertionReturnType:
    """Insert markers either side of elements in the magnetic lattice, selected
    based on the element name (either equality or with a regular expression).
    By default markers are placed either side of the matched elements.

    :param string: Element string or regex string (with regex=True) to select
    elements to wrap with markers.
    :param regex: Whether to interpret argument `string` as a regex or not.  If
    not, then `string` is checked for equality against the element name.
    :param: Place a marker before each matched element.
    :param: Place a marker after each matched element.
    :return: Dictionary of element instances to list of (marker_start,
    marker_end pairs, like
    {element_instance: [(start1, start2), ... (startn, endn)]}.

    """
    if regex:
        def fre(ele):
            return bool(re.match(string, ele.id))

        return insert_markers_by_predicate(sequence, fre)

    def f(ele):
        return ele.id == string

    return insert_markers_by_predicate(sequence, f, before=before, after=after)


def insert_markers_by_type(sequence, magnet_type: Element, before=True, after=True
                           ) -> MarkersInsertionReturnType:
    """Insert markers either side of elements in the magnetic lattice, selected
    based on the element
    By default markers are placed either side of the matched elements.

    :param string: Element string or regex string (with regex=True) to select
    elements to wrap with markers.
    :param: Place a marker before each matched element.
    :param: Place a marker after each matched element.
    :return: Dictionary of element instances to list of \
    (marker_start, marker_end pairs)
    """

    def f(ele):
        return isinstance(ele, magnet_type)

    return insert_markers_by_predicate(sequence, f, before=before, after=after)


def insert_markers_by_predicate(sequence, predicate: Callable[[Element], bool],
                                before=True, after=True,
                                before_suffix="_before",
                                after_suffix="_after",
                                matched_slice=None
                                ) -> MarkersInsertionReturnType:
    """Insert markers either side of elements in the magnetic lattice, selected
    based on the the provided predicate function.

    :param predicate: Function on each element to select elements to wrap with
    markers.
    :param: Place a marker before each matched element.
    :param: Place a marker after each matched element.
    :return: Dictionary of element instances to list of \
    (marker_start, marker_end pairs)
    """

    indices = [i for (i, ele) in enumerate(sequence) if predicate(ele)]
    out_dict = defaultdict(list)
    for i in reversed(indices):
        ele = sequence[i]
        marker_before = None
        marker_after = None
        if after:
            marker_id_after = f"{ele.id}{after_suffix}"
            marker_after = Marker(marker_id_after)
            sequence.insert(i + 1, marker_after)
        if before:
            marker_id_before = f"{ele.id}{before_suffix}"
            marker_before = Marker(marker_id_before)
            sequence.insert(i, marker_before)

        out_dict[ele].append((marker_before, marker_after))

    return dict(out_dict)
