"""
Utils for creating closed bumps in both planes and convert the correctors to dipoles
S.Tomin, DESY, 11.2020
"""
from typing import Union, List
import numpy as np
import copy
from ocelot.cpbd.elements import Element, Monitor, Marker, Hcor, Vcor, Bend, Edge
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.optics import lattice_transfer_map
from ocelot.cpbd.beam import Particle
from ocelot.cpbd.track import lattice_track


def bump_4cors(lat: MagneticLattice, cor_list: List[Element], marker: Union[Marker, Monitor],
                    x: float = 0.0001, xp: float = 0., energy: float = 14.) -> np.array:
    """
    Bump with 4 correctors.
    Function calculates correctors angles (strength) to make closed bump with (x, x') or (y, y') at the marker position.
    Calculated angles are applied to corresponding correctors (MagneticLattice is mutated).

    :param lat: MagneticLattice
    :param cor_list: list of 4 correctors, All for correctors must be in the same plane
    :param marker: Marker or Monitor. It must be between 2d and 3d corrector
    :param x: x or y beam coordinate at Marker position [m]
    :param xp: x' or y' beam coordinate at Marker position [rad]
    :param energy: Beam energy [GeV]
    :return: numpy array of corrector angles in [rad]
    """
    # lattice analysis, calculate R matrices between different devices
    if len(cor_list) != 4:
        raise TypeError("It must be 4 correctors")

    if len([cor for cor in cor_list if isinstance(cor, Hcor)]) == 4:
        x_idx, xp_idx = 0, 1
    elif len([cor for cor in cor_list if isinstance(cor, Vcor)]) == 4:
        x_idx, xp_idx = 2, 3
    else:
        raise TypeError("Not all correctors belong to Y- or X-plane")

    Rs = []
    sorted_list = [cor_list[0], cor_list[1], marker, cor_list[2], cor_list[3]]
    for i in range(4):
        sequence = copy.deepcopy(lat.get_sequence_part(start=sorted_list[i], stop=sorted_list[i + 1]))
        sequence[0].l /= 2.
        sequence[-1].l /= 2.
        lat1 = MagneticLattice(sequence)
        R1 = lattice_transfer_map(lat1, energy)
        Rs.append(R1)

    R21 = np.dot(Rs[1], Rs[0])

    M = np.array([[R21[x_idx, xp_idx], Rs[1][x_idx, xp_idx]],
                  [R21[xp_idx, xp_idx], Rs[1][xp_idx, xp_idx]]])

    Minv = np.linalg.inv(M)
    X = np.array([x, xp])
    cor_strengths = np.dot(Minv, X)
    angle_1, angle_2 = cor_strengths

    R43 = np.dot(Rs[3], Rs[2])

    Xf = np.array([0., 0.])

    M43 = np.array([[R43[x_idx, x_idx], R43[x_idx, xp_idx]],
                    [R43[xp_idx, x_idx], R43[xp_idx, xp_idx]]])
    X4 = np.dot(M43, X)
    angle_3 = (Xf[0] - X4[0]) / Rs[3][x_idx, xp_idx]
    angle_4 = -(X4[1] + Rs[3][xp_idx, xp_idx] * angle_3 - Xf[1])
    a = np.array([angle_1, angle_2, angle_3, angle_4])
    for i, cor in enumerate(cor_list):
        cor.angle = a[i]
    lat.update_transfer_maps()
    return a


def convert_cors2dipoles(lat: MagneticLattice, cor_list: List[Element], energy: float = 10.) -> MagneticLattice:
    """
    Function converts correctors with non zero angles to dipoles

    :param lat: MagneticLattice
    :param cor_list: list of correctors
    :param energy: beam energy
    :return: new MagneticLattice
    """

    # track a particle through lattice with correctors
    plist = lattice_track(lat, Particle(E=energy))
    splist = np.array([p.s for p in plist])

    # calculate edge positions for all elements
    L = 0.
    for elem in lat.sequence:
        elem.s0 = L
        L += elem.l
        elem.s1 = L

    # create dipoles
    dipoles = []
    for corr in cor_list:
        scor = [corr.s0, corr.s1]
        ins1 = [np.argmin(np.abs(splist - sbpm)) for sbpm in scor]
        if isinstance(corr, Hcor):
            p_str = "px"
            tilt = 0.
        else:
            p_str = "py"
            tilt = np.pi/2
        edges = []
        for i in range(2):
            edges.append(plist[ins1[i]].__dict__[p_str])
        dipoles.append(Bend(l=corr.l, angle=corr.angle, e1=edges[0], e2=edges[1], tilt=tilt, eid=corr.id))

    # create new sequence with substituted correctors by dipoles
    seq = []
    for i, elem in enumerate(lat.sequence):
        if elem in cor_list:
            for d in dipoles:
                if elem.id == d.id:
                    seq.append(d)
        else:
            if elem.__class__ is Edge:
                continue
            seq.append(elem)
    lat_new = MagneticLattice(seq)

    return lat_new