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
    obj_id = []
    objs = []
    for elem in lat.sequence:
        if elem.__class__ in types:
            if elem.id not in obj_id:
                objs.append(elem)
                obj_id.append(elem.id)
    return objs



def lat2input(lat):
    """
    returns python input string for the lattice in the lat object
    """
    drifts = find_objects(lat, types = [Drift])
    quads = find_objects(lat, types = [Quadrupole])
    sexts = find_objects(lat, types = [Sextupole])
    octs = find_objects(lat, types = [Octupole])
    cavs = find_objects(lat, types = [Cavity])
    sols = find_objects(lat, types = [Solenoid])
    matrices = find_objects(lat, types = [Matrix])
    marks = find_objects(lat, types = [Marker])
    mons = find_objects(lat, types = [Monitor])
    unds = find_objects(lat, types = [Undulator])
    cors = find_objects(lat, types = [Hcor, Vcor])
    bends = find_objects(lat, types = [Bend, RBend, SBend])
    #end find objects

    lines = []
    lines.append("\n# drifts \n")
    for drift in drifts:
        line = drift.id.replace('.','_') + " = Drift(l = " + str(drift.l)+ ", eid = '"+ drift.id+ "')\n"
        lines.append(line)

    lines.append("\n# quadrupoles \n")
    for quad in quads:
        line = quad.id.replace('.','_') + " = Quadrupole(l = " + str(quad.l) +", k1 = "+ str(quad.k1) +", tilt = "+ str(quad.tilt)  +", eid = '"+ quad.id+ "')\n"
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
        fint = bend.fint1
        #if bend.fint1 == bend.fint2:
        #    fint = bend.fint1

        line = bend.id.replace('.','_') + type + str(bend.l) + k + ", angle = " + str(bend.angle)+ ", e1 = " + str(bend.e1) + ", e2 = " + str(bend.e2) + ", tilt = " + str(bend.tilt) + ", fint = " + str(fint) +", eid = '"+ bend.id+ "')\n"
        lines.append(line)

    lines.append("\n# correctors \n")
    for cor in cors:
        if cor.__class__ == Hcor:
            type = " = Hcor(l = "
        else:
            type = " = Vcor(l = "
        line = cor.id.replace('.','_') + type + str(cor.l) + ", angle = "+ str(cor.angle) +", eid = '"+ cor.id+ "')\n"
        lines.append(line)

    lines.append("\n# markers \n")
    for mark in marks:
        line = mark.id.replace('.','_') + " = Marker(eid = '"+ mark.id+ "')\n"
        lines.append(line)

    lines.append("\n# monitor \n")
    for mon in mons:
        line = mon.id.replace('.','_') + " = Monitor(eid = '"+ mon.id+ "')\n"
        lines.append(line)

    lines.append("\n# sextupoles \n")
    for sext in sexts:
        line = sext.id.replace('.','_') + " = Sextupole(l = " + str(sext.l) +", k2 = "+ str(sext.k2) +", tilt = "+ str(sext.tilt) +", eid = '"+ sext.id+ "')\n"
        lines.append(line)

    lines.append("\n# octupole \n")
    for oct in octs:
        line = oct.id.replace('.','_') + " = Octupole(l = " + str(oct.l) +", k3 = "+ str(oct.k3) +", tilt = "+ str(oct.tilt) +", eid = '"+ oct.id+ "')\n"
        lines.append(line)

    lines.append("\n# undulator \n")
    for und in unds:
        line = und.id.replace('.','_') + " = Undulator(lperiod = " + str(und.lperiod) +", nperiods = "+ str(und.nperiods) +", Kx = "+ str(und.Kx) +", Ky = "+ str(und.Ky) +", eid = '"+ und.id+ "')\n"
        lines.append(line)

    lines.append("\n# cavity \n")
    for cav in cavs:
        line = cav.id.replace('.','_') + " = Cavity(l = " + str(cav.l) + ", v = "+ str(cav.v) +\
               ", freq = "+ str(cav.f) +", phi = "+ str(cav.phi) + ", eid = '"+ cav.id+ "')\n"
        lines.append(line)

    #lines.append("\n# rfcavity \n")
    #for rfcav in rfcavs:
    #    line = rfcav.id.replace('.','_') + " = RFcavity(l = " + str(rfcav.l) +", volt = "+ str(rfcav.volt) +", lag = "+ str(rfcav.lag)+\
    #           ", harmon = "+ str(rfcav.harmon) +", eid = '"+ rfcav.id+ "')\n"
    #    lines.append(line)

    lines.append("\n# Matrices \n")

    for mat in matrices:
        line = mat.id.replace('.','_') + " = Matrix(l = " + str(mat.l) +\
               ", rm11 = "+ str(mat.rm11) + ", rm12 = "+ str(mat.rm12) + ", rm13 = "+ str(mat.rm13) + ", rm14 = "+ str(mat.rm14) +\
               ", rm21 = "+ str(mat.rm21) + ", rm22 = "+ str(mat.rm22) + ", rm23 = "+ str(mat.rm23) + ", rm24 = "+ str(mat.rm24) +\
               ", rm31 = "+ str(mat.rm31) + ", rm32 = "+ str(mat.rm32) + ", rm33 = "+ str(mat.rm33) + ", rm34 = "+ str(mat.rm34) +\
               ", rm41 = "+ str(mat.rm41) + ", rm42 = "+ str(mat.rm42) + ", rm43 = "+ str(mat.rm43) + ", rm44 = "+ str(mat.rm44) +\
               ", eid = '"+ mat.id+ "')\n"
        lines.append(line)

    lines.append("\n# Solenoids \n")

    for sol in sols:
        line = sol.id.replace('.','_') + " = Solenoid(l = " + str(sol.l) +", k = "+ str(sol.k) +", eid = '"+ sol.id+ "')\n"
        lines.append(line)

    lines.append("\n# lattice \n")
    names = []
    for elem in lat.sequence:
        #print elem.id, elem.type
        if elem.__class__ != Edge:
            names.append(elem.id.replace('.','_'))
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


def write_lattice(lattice, file_name="lattice.inp"):
    """
    saves lattice as python imput file
    """
    lines = lat2input(lattice)

    f = open(file_name, 'w')
    f.writelines(lines)
    f.close()
