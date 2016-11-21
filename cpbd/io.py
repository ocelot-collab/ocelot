__author__ = 'Sergey Tomin'
"""
module contains lat2input function which creates python input string
(ocelot natice) for a lattice object
author sergey.tomin
"""


from numpy import around, pi
from ocelot.cpbd.elements import *

def find_drifts(lat):
    drift_lengs = []
    drifts = []
    for elem in lat.sequence:
        if elem.__class__ == Drift:
            elem_l = around(elem.l, decimals = 6)
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
                name = name.replace('.','_')
                name = name.replace(':','_')
                ids[i] = name + alphabet[n]
        else:
            name = ids[j]
            name = name.replace('.','_')
            name = name.replace(':','_')
            ids[j] = name
        obj.name = ids[j]

    return objects


def find_obj_and_create_name(lat, types):
    objects = find_objects(lat, types=types)
    objects = create_var_name(objects)
    return objects

def lat2input(lat):
    """
    returns python input string for the lattice in the lat object
    """
    drifts = find_obj_and_create_name(lat, types = [Drift])
    quads = find_obj_and_create_name(lat, types = [Quadrupole])
    sexts = find_obj_and_create_name(lat, types = [Sextupole])
    octs = find_obj_and_create_name(lat, types = [Octupole])
    cavs = find_obj_and_create_name(lat, types = [Cavity])
    sols = find_obj_and_create_name(lat, types = [Solenoid])
    matrices = find_obj_and_create_name(lat, types = [Matrix])
    marks = find_obj_and_create_name(lat, types = [Marker])
    mons = find_obj_and_create_name(lat, types = [Monitor])
    unds = find_obj_and_create_name(lat, types = [Undulator])
    cors = find_obj_and_create_name(lat, types = [Hcor, Vcor])
    bends = find_obj_and_create_name(lat, types = [Bend, RBend, SBend])
    unkns = find_obj_and_create_name(lat, types = [UnknownElement])
    #end find objects

    lines = ["from ocelot import * \n"]
    lines.append("\n# drifts \n")
    for drift in drifts:
        line = drift.name.lower() + " = Drift(l=" + str(drift.l)+ ", eid='"+ drift.id+ "')\n"
        lines.append(line)

    lines.append("\n# quadrupoles \n")
    for quad in quads:
        line = quad.name.lower() + " = Quadrupole(l=" + str(quad.l) +", k1="+ str(quad.k1) +", tilt="+ str(quad.tilt)  +", eid='"+ quad.id+ "')\n"
        lines.append(line)

    lines.append("\n# bending magnets \n")
    for bend in bends:
        if bend.__class__ == RBend:
            type = " = RBend(l = "
        elif bend.__class__ == SBend:
            type = " = SBend(l = "
        else:
            type = " = Bend(l = "
        if bend.k1 == 0 or bend.k1 == None:
            k = ''
        else:
            k = ", k1 = "+ str(bend.k1)

        line = bend.name.lower() + type + str(bend.l) + k + ", angle=" + str(bend.angle)+\
               ", e1=" + str(bend.e1) + ", e2=" + str(bend.e2)+ ", gap=" + str(bend.gap)  + \
               ", tilt=" + str(bend.tilt) + ", fint=" + str(bend.fint)  + ", fintx=" + str(bend.fintx) +", eid='"+ bend.id+ "')\n"
        lines.append(line)

    lines.append("\n# correctors \n")
    for cor in cors:
        if cor.__class__ == Hcor:
            type = " = Hcor(l="
        else:
            type = " = Vcor(l="

        line = cor.name.lower() + type + str(cor.l) + ", angle="+ str(cor.angle) +", eid='"+ cor.id+ "')\n"
        lines.append(line)

    lines.append("\n# markers \n")
    for mark in marks:
        line = mark.name.lower() + " = Marker(eid='"+ mark.id+ "')\n"
        lines.append(line)

    lines.append("\n# monitor \n")
    for mon in mons:
        line = mon.name.lower() + " = Monitor(eid='"+ mon.id+ "')\n"
        lines.append(line)

    lines.append("\n# sextupoles \n")
    for sext in sexts:
        line = sext.name.lower() +" = Sextupole(l=" + str(sext.l) +", k2="+ str(sext.k2) +", tilt="+ str(sext.tilt) +", eid='"+ sext.id+ "')\n"
        lines.append(line)

    lines.append("\n# octupole \n")
    for oct in octs:
        line = oct.name.lower() + " = Octupole(l=" + str(oct.l) +", k3="+ str(oct.k3) +", tilt="+ str(oct.tilt) +", eid='"+ oct.id+ "')\n"
        lines.append(line)

    lines.append("\n# undulator \n")
    for und in unds:
        line = und.name.lower() +  " = Undulator(lperiod=" + str(und.lperiod) +", nperiods="+ str(und.nperiods) +", Kx="+ str(und.Kx) +", Ky="+ str(und.Ky) +", eid='"+ und.id+ "')\n"
        lines.append(line)

    lines.append("\n# cavity \n")
    for cav in cavs:
        line = cav.name.lower() + " = Cavity(l=" + str(cav.l) + ", v="+ str(cav.v) +\
               ", freq="+ str(cav.f) +", phi="+ str(cav.phi) + ", eid='"+ cav.id+ "')\n"
        lines.append(line)

    lines.append("\n# UnknowElement \n")
    for unkn in unkns:
        line = unkn.name.lower() +  " = UnknownElement(l=" + str(unkn.l) +  ", eid='"+ unkn.id+ "')\n"
        lines.append(line)
    #lines.append("\n# rfcavity \n")
    #for rfcav in rfcavs:
    #    line = rfcav.id.replace('.','_') + " = RFcavity(l = " + str(rfcav.l) +", volt = "+ str(rfcav.volt) +", lag = "+ str(rfcav.lag)+\
    #           ", harmon = "+ str(rfcav.harmon) +", eid = '"+ rfcav.id+ "')\n"
    #    lines.append(line)

    lines.append("\n# Matrices \n")

    for mat in matrices:
        line = mat.name.lower() + " = Matrix(l=" + str(mat.l) +\
               ", rm11="+ str(mat.rm11) + ", rm12="+ str(mat.rm12) + ", rm13="+ str(mat.rm13) + ", rm14="+ str(mat.rm14) +\
               ", rm21="+ str(mat.rm21) + ", rm22="+ str(mat.rm22) + ", rm23="+ str(mat.rm23) + ", rm24="+ str(mat.rm24) +\
               ", rm31="+ str(mat.rm31) + ", rm32="+ str(mat.rm32) + ", rm33="+ str(mat.rm33) + ", rm34="+ str(mat.rm34) +\
               ", rm41="+ str(mat.rm41) + ", rm42="+ str(mat.rm42) + ", rm43="+ str(mat.rm43) + ", rm44="+ str(mat.rm44) +\
               ", eid = '"+ mat.id+ "')\n"
        lines.append(line)

    lines.append("\n# Solenoids \n")

    for sol in sols:
        line = sol.name.lower() + " = Solenoid(l=" + str(sol.l) +", k="+ str(sol.k) +", eid='"+ sol.id+ "')\n"
        lines.append(line)

    lines.append("\n# lattice \n")
    names = []
    for elem in lat.sequence:
        if elem.__class__ != Edge:
            names.append(elem.name.lower())
    #names = map(lambda p: p.id, lat.sequence)
    new_names = []
    for i, name in enumerate(names):
        if i%8 == 7:
            new_names.append("\n"+name)
        else:
            new_names.append(name)
    line = "cell = (" + ", ".join(new_names) +")"
    lines.append(line)
    return lines


def rem_drifts(lat):
    drifts = {}
    for i, elem in enumerate(lat.sequence):
        if elem.__class__ == Drift:
            if not (elem.l in drifts.keys()):
                drifts[elem.l] = elem
            else:
                # print(cell[i],  drifts[elem.l])
                lat.sequence[i] = drifts[elem.l]

    return lat

def write_lattice(lattice, file_name="lattice.inp", remove_rep_drifts=True):
    """
    saves lattice as python imput file
    lattice - MagneticLattice
    file_name - name of the file
    remove_rep_drifts - if True, remove the drifts with the same lengths from the lattice drifts definition
    """
    lattice = rem_drifts(lattice)
    lines = lat2input(lattice)

    f = open(file_name, 'w')
    f.writelines(lines)
    f.close()
