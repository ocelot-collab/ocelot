__author__ = 'Sergey Tomin'
"""
module contains lat2input function which creates python input string
(ocelot lattice) for a lattice object
author sergey.tomin
"""
from ocelot.cpbd.elements import *
import os, sys
from ocelot.adaptors.astra2ocelot import astraBeam2particleArray, particleArray2astraBeam
from ocelot.adaptors.csrtrack2ocelot import csrtrackBeam2particleArray, particleArray2csrtrackBeam
from ocelot.cpbd.beam import ParticleArray, Twiss, Beam


def save_particle_array2npz(filename, p_array):
    np.savez_compressed(filename, rparticles=p_array.rparticles,
                        q_array=p_array.q_array,
                        E=p_array.E, s=p_array.s)


def load_particle_array_from_npz(filename, print_params=False):
    """
    Load beam file in npz format and return ParticleArray

    :param filename:
    :return:
    """
    p_array = ParticleArray()
    with np.load(filename) as data:
        for key in data.keys():
            p_array.__dict__[key] = data[key]
    return p_array


def load_particle_array(filename, print_params=False):
    """
    Universal function to load beam file, *.ast (ASTRA), *.fmt1 (CSRTrack) or *.npz format

    Note that downloading ParticleArray from the astra file (.ast) and saving it back does not give the same distribution.
    The difference arises because the array of particles does not have a reference particle, and in this case
    the first particle is used as a reference.

    :param filename: path to file, filename.ast or filename.npz
    :return: ParticleArray
    """
    name, file_extension = os.path.splitext(filename)
    if file_extension == ".npz":
        parray = load_particle_array_from_npz(filename, print_params=False)
    elif file_extension in [".ast", ".001"]:
        parray = astraBeam2particleArray(filename, print_params=False)
    elif file_extension in [".fmt1"]:
        parray = csrtrackBeam2particleArray(filename)
    else:
        raise Exception("Unknown format of the beam file: " + file_extension + " but must be *.ast, *fmt1 or *.npz ")

    if print_params:
        print(parray)

    return parray


def save_particle_array(filename, p_array):
    """
    Universal function to save beam file, *.ast (ASTRA), *.fmt1 (CSRTrack) or *.npz format

    Note that downloading ParticleArray from the astra file (.ast) and saving it back does not give the same distribution.
    The difference arises because the array of particles does not have a reference particle, and in this case
    the first particle is used as a reference.

    :param filename: path to file, filename.ast or filename.npz
    :return: ParticleArray
    """
    name, file_extension = os.path.splitext(filename)
    if file_extension == ".npz":
        save_particle_array2npz(filename, p_array)
    elif file_extension == ".ast":
        particleArray2astraBeam(p_array, filename)
    elif file_extension == ".fmt1":
        particleArray2csrtrackBeam(p_array, filename)
    else:
        raise Exception("Unknown format of the beam file: " + file_extension + " but must be *.ast, *.fmt1 or *.npz")


def find_drifts(lat):
    drift_lengs = []
    drifts = []
    for elem in lat.sequence:
        if elem.__class__ == Drift:
            elem_l = np.around(elem.l, decimals=6)
            if elem_l not in drift_lengs:
                drifts.append(elem)
                drift_lengs.append(elem_l)
    return drifts


def find_objects(lat, types):
    """
    Function finds objects by types and adds it to list if object is unique.
    :param lat: MagneticLattice
    :param types: types of the Elements
    :return: list of elements
    """
    obj_id = []
    objs = []
    for elem in lat.sequence:
        if elem.__class__ in types:
            if id(elem) not in obj_id:
                objs.append(elem)
                obj_id.append(id(elem))

    return objs


def create_var_name(objects):
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


def find_obj_and_create_name(lat, types):
    objects = find_objects(lat, types=types)
    objects = create_var_name(objects)
    return objects


def lat2input(lattice, tws0=None):
    """
    returns python input string for the lattice in the lat object
    """

    lines = ['from ocelot import * \n']

    # prepare initial Twiss parameters
    if tws0 is not None and isinstance(tws0, Twiss):
        lines.append('\n#Initial Twiss parameters\n')
        lines.extend(twiss2input(tws0))

    # prepare elements list
    lines.append('\n')
    lines.extend(elements2input(lattice))

    # prepare cell list
    lines.append('\n# Lattice \n')
    lines.extend(cell2input(lattice, True))

    lines.append('\n')

    return lines


def get_elements(lattice):

    elements = []
        
    for element in lattice.sequence:
        element_type = element.__class__.__name__
        
        if element_type in ('Edge', "CouplerKick"):
            continue

        if element not in elements:
            elements.append(element)

    elements = create_var_name(elements)

    return elements


def matrix_def_string(element, params):
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

def element_def_string(element):

    #if element.__class__ == Matrix:
    #    return matrix_def_string(element)

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
                val = element.__dict__[param] - element.angle/2.0
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

    if element.__class__ is Matrix:
        params = matrix_def_string(element, params)

    # join all parameters to element definition
    string = pprinting(element, element_type, params)
    return string


def pprinting(element, element_type, params):
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


def print_elements(elements_dict):
    
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
                string = element_def_string(element)
                lines.append(string)

    # print remaining unordered elements
    for type in unordered_dict:
        
        lines.append('\n# ' + type + 's\n')

        for element in unordered_dict[type]:
            string = element_def_string(element)
            lines.append(string)

    # delete new line symbol from the first line
    if lines != []:
        lines[0] = lines[0][1:]
    return lines


def sort_elements(elements):

    elements_dict = {}
    for element in elements:
        element_type = element.__class__.__name__

        if element_type not in elements_dict:
            elements_dict[element_type] = []

        elements_dict[element_type].append(element)

    return elements_dict


def elements2input(lattice):

    elements = get_elements(lattice)
    elements_dict = sort_elements(elements)
    lines = print_elements(elements_dict)

    return lines


def cell2input(lattice, split=False):
    
    lines = []
    names = []
    for elem in lattice.sequence:
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


def twiss2input(tws):

    lines = []
    tws_ref = Twiss()
    lines.append('tws0 = Twiss()\n')
    for param in tws.__dict__:
        if tws.__dict__[param] != tws_ref.__dict__[param]:
            lines.append('tws0.' + str(param) + ' = ' + str(tws.__dict__[param]) + '\n')

    return lines


def beam2input(beam):
    
    lines = []
    beam_ref = Beam()
    lines.append('beam = Beam()\n')
    for param in beam.__dict__:
        if beam.__dict__[param] != beam_ref.__dict__[param]:
            lines.append('beam.' + str(param) + ' = ' + str(beam.__dict__[param]) + '\n')

    return lines


def rem_drifts(lat):
    drifts = {}
    for i, elem in enumerate(lat.sequence):
        if elem.__class__ == Drift:
            if not (elem.l in drifts.keys()):
                drifts[elem.l] = elem
            else:
                lat.sequence[i] = drifts[elem.l]

    return lat


def write_power_supply_id(lattice, lines=[]):
    quads = find_obj_and_create_name(lattice, types=[Quadrupole])
    sexts = find_obj_and_create_name(lattice, types=[Sextupole])
    octs = find_obj_and_create_name(lattice, types=[Octupole])
    cavs = find_obj_and_create_name(lattice, types=[Cavity])
    bends = find_obj_and_create_name(lattice, types=[Bend, RBend, SBend])

    lines.append("\n# power supplies \n")
    for elem_group in [quads, sexts, octs, cavs, bends]:
        lines.append("\n#  \n")
        for elem in elem_group:
            if "ps_id" in dir(elem):
                line = elem.name.lower() + ".ps_id = '" + elem.ps_id + "'\n"
                lines.append(line)
    return lines


def write_lattice(lattice, tws0=None, file_name="lattice.py", remove_rep_drifts=True, power_supply=False):
    """
    saves lattice as python imput file
    lattice - MagneticLattice
    file_name - name of the file
    remove_rep_drifts - if True, remove the drifts with the same lengths from the lattice drifts definition
    """
    if remove_rep_drifts:
        lattice = rem_drifts(lattice)

    lines = lat2input(lattice, tws0=tws0)

    if power_supply:
        lines = write_power_supply_id(lattice, lines=lines)

    f = open(file_name, 'w')
    f.writelines(lines)
    f.close()
