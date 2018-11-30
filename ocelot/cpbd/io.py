__author__ = 'Sergey Tomin'
"""
module contains lat2input function which creates python input string
(ocelot lattice) for a lattice object
author sergey.tomin
"""
from ocelot.cpbd.elements import *
import os
from ocelot.adaptors.astra2ocelot import astraBeam2particleArray, particleArray2astraBeam
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
    Universal function to load beam file, *.ast or *.npz format

    :param filename: path to file, filename.ast or filename.npz
    :return: ParticleArray
    """
    name, file_extension = os.path.splitext(filename)
    if file_extension == ".npz":
        return load_particle_array_from_npz(filename, print_params=print_params)
    elif file_extension in [".ast", ".001"]:
        return astraBeam2particleArray(filename, s_ref=-1, Eref=-1, print_params=print_params)
    else:
        raise Exception("Unknown format of the beam file: " + file_extension + " but must be *.ast or *.npz ")


def save_particle_array(filename, p_array, ref_index=0):
    """
    Universal function to save beam file, *.ast or *.npz format

    :param filename: path to file, filename.ast or filename.npz
    :param ref_index: index of ref particle
    :return: ParticleArray
    """
    name, file_extension = os.path.splitext(filename)
    if file_extension == ".npz":
        save_particle_array2npz(filename, p_array)
    elif file_extension == ".ast":
        particleArray2astraBeam(p_array, filename, ref_index)
    else:
        raise Exception("Unknown format of the beam file: " + file_extension + " but must be *.ast or *.npz")



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
                ids[i] = name + alphabet[n]
        else:
            name = ids[j]
            name = name.replace('.', '_')
            name = name.replace(':', '_')
            ids[j] = name
        obj.name = ids[j]

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


def elements2input(lattice):

    # find objects
    drifts = find_obj_and_create_name(lattice, types=[Drift])
    quads = find_obj_and_create_name(lattice, types=[Quadrupole])
    sexts = find_obj_and_create_name(lattice, types=[Sextupole])
    octs = find_obj_and_create_name(lattice, types=[Octupole])
    cavs = find_obj_and_create_name(lattice, types=[Cavity])
    sols = find_obj_and_create_name(lattice, types=[Solenoid])
    matrices = find_obj_and_create_name(lattice, types=[Matrix])
    marks = find_obj_and_create_name(lattice, types=[Marker])
    mons = find_obj_and_create_name(lattice, types=[Monitor])
    unds = find_obj_and_create_name(lattice, types=[Undulator])
    cors = find_obj_and_create_name(lattice, types=[Hcor, Vcor])
    bends = find_obj_and_create_name(lattice, types=[Bend, RBend, SBend])
    unkns = find_obj_and_create_name(lattice, types=[UnknownElement])
    tcavs = find_obj_and_create_name(lattice, types=[TDCavity])
    xyquads = find_obj_and_create_name(lattice, types=[XYQuadrupole])
    # prepare txt elements list
    lines = []

    if len(drifts) != 0:
        lines.append("\n# drifts \n")
    
    for drift in drifts:
        line = drift.name.lower() + " = Drift(l=" + str(drift.l) + ", eid='" + drift.id + "')\n"
        lines.append(line)

    if len(quads) != 0:
        lines.append("\n# quadrupoles \n")
    
    for quad in quads:
        line = quad.name.lower() + " = Quadrupole(l=" + str(quad.l) + ", k1=" + str(quad.k1) + ", tilt=" + str(
            quad.tilt) + ", eid='" + quad.id + "')\n"
        lines.append(line)

    if len(quads) != 0:
        lines.append("\n# quadrupoles with offsets \n")

    for quad in xyquads:
        line = quad.name.lower() + " = XYQuadrupole(l=" + str(quad.l) + ", k1=" + str(quad.k1) + \
               ", x_offs=" + str(quad.x_offs) + ", y_offs=" + str(quad.y_offs) + ", tilt=" + str(
            quad.tilt) + ", eid='" + quad.id + "')\n"
        lines.append(line)

    if len(bends) != 0:
        lines.append("\n# bending magnets \n")
    
    for bend in bends:
        if bend.__class__ == RBend:
            type = " = RBend(l="
        elif bend.__class__ == SBend:
            type = " = SBend(l="
        else:
            type = " = Bend(l="
        if bend.k1 == 0 or bend.k1 == None:
            k = ''
        else:
            k = ", k1 = " + str(bend.k1)

        line = bend.name.lower() + type + str(bend.l) + k + ", angle=" + str(bend.angle) + \
               ", e1=" + str(bend.e1) + ", e2=" + str(bend.e2) + ", gap=" + str(bend.gap) + \
               ", tilt=" + str(bend.tilt) + ", fint=" + str(bend.fint) + ", fintx=" + str(
            bend.fintx) + ", eid='" + bend.id + "')\n"
        lines.append(line)

    if len(cors) != 0:
        lines.append("\n# correctors \n")
    
    for cor in cors:
        if cor.__class__ == Hcor:
            type = " = Hcor(l="
        else:
            type = " = Vcor(l="

        line = cor.name.lower() + type + str(cor.l) + ", angle=" + str(cor.angle) + ", eid='" + cor.id + "')\n"
        lines.append(line)

    if len(marks) != 0:
        lines.append("\n# markers \n")
    
    for mark in marks:
        line = mark.name.lower() + " = Marker(eid='" + mark.id + "')\n"
        lines.append(line)

    if len(mons) != 0:
        lines.append("\n# monitors \n")
    
    for mon in mons:
        line = mon.name.lower() + " = Monitor(eid='" + mon.id + "')\n"
        lines.append(line)

    if len(sexts) != 0:
        lines.append("\n# sextupoles \n")
    
    for sext in sexts:
        line = sext.name.lower() + " = Sextupole(l=" + str(sext.l) + ", k2=" + str(sext.k2) + ", tilt=" + str(
            sext.tilt) + ", eid='" + sext.id + "')\n"
        lines.append(line)

    if len(octs) != 0:
        lines.append("\n# octupoles \n")
    
    for oct in octs:
        line = oct.name.lower() + " = Octupole(l=" + str(oct.l) + ", k3=" + str(oct.k3) + ", tilt=" + str(
            oct.tilt) + ", eid='" + oct.id + "')\n"
        lines.append(line)

    if len(unds) != 0:
        lines.append("\n# undulators \n")
    
    for und in unds:
        line = und.name.lower() + " = Undulator(lperiod=" + str(und.lperiod) + ", nperiods=" + str(
            und.nperiods) + ", Kx=" + str(und.Kx) + ", Ky=" + str(und.Ky) + ", eid='" + und.id + "')\n"
        lines.append(line)

    if len(cavs) != 0:
        lines.append("\n# cavities \n")
    
    for cav in cavs:
        line = cav.name.lower() + " = Cavity(l=" + str(cav.l) + ", v=" + str(cav.v) + \
               ", freq=" + str(cav.freq) + ", phi=" + str(cav.phi) + ", eid='" + cav.id + "')\n"
        lines.append(line)

    if len(tcavs) != 0:
        lines.append("\n# tdcavities \n")
    
    for tcav in tcavs:
        line = tcav.name.lower() + " = Cavity(l=" + str(tcav.l) + ", v=" + str(tcav.v) + \
               ", freq=" + str(tcav.freq) + ", phi=" + str(tcav.phi) + ", eid='" + tcav.id + "')\n"
        lines.append(line)

    if len(unkns) != 0:
        lines.append("\n# UnknowElements \n")
    
    for unkn in unkns:
        line = unkn.name.lower() + " = UnknownElement(l=" + str(unkn.l) + ", eid='" + unkn.id + "')\n"
        lines.append(line)

    if len(matrices) != 0:
        lines.append("\n# Matrices \n")
    
    for mat in matrices:
        line = mat.name.lower() + " = Matrix(l=" + str(mat.l) + \
               ", rm11=" + str(mat.rm11) + ", rm12=" + str(mat.rm12) + ", rm13=" + str(mat.rm13) + ", rm14=" + str(
            mat.rm14) + \
               ", rm21=" + str(mat.rm21) + ", rm22=" + str(mat.rm22) + ", rm23=" + str(mat.rm23) + ", rm24=" + str(
            mat.rm24) + \
               ", rm31=" + str(mat.rm31) + ", rm32=" + str(mat.rm32) + ", rm33=" + str(mat.rm33) + ", rm34=" + str(
            mat.rm34) + \
               ", rm41=" + str(mat.rm41) + ", rm42=" + str(mat.rm42) + ", rm43=" + str(mat.rm43) + ", rm44=" + str(
            mat.rm44) + \
               ", eid = '" + mat.id + "')\n"
        lines.append(line)

    if len(sols) != 0:
        lines.append("\n# Solenoids \n")

    for sol in sols:
        line = sol.name.lower() + " = Solenoid(l=" + str(sol.l) + ", k=" + str(sol.k) + ", eid='" + sol.id + "')\n"
        lines.append(line)

    if lines != []:
        lines[0] = lines[0][1:]
        
    return lines


def cell2input(lattice, split=False):
    
    lines = []
    names = []
    for elem in lattice.sequence:
        if elem.__class__ != Edge:
            names.append(elem.name.lower())

    new_names = []
    for i, name in enumerate(names):
        if split and i % 8 == 7:
            new_names.append("\n" + name)
        else:
            new_names.append(name)

    lines.append("cell = (" + ", ".join(new_names) + ")")
    
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
