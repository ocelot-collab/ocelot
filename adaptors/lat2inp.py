__author__ = 'Sergey Tomin'

from numpy import around

def find_drifts(lat):
    drift_lengs = []
    drifts = []
    for elem in lat.sequence:
        if elem.type == "drift":
            elem_l = around(elem.l, decimals = 6)
            if elem_l not in drift_lengs:
                drifts.append(elem)
                drift_lengs.append(elem_l)
    return drifts



def find_objects(lat, types):
    obj_id = []
    objs = []
    for elem in lat.sequence:
        if elem.type in types:
            if elem.id not in obj_id:
                objs.append(elem)
                obj_id.append(elem.id)
    return objs



def lat2input(lat):
    # start find objects
    drifts = find_objects(lat, types = ["drift"])
    quads = find_objects(lat, types = ["quadrupole"])
    sexts = find_objects(lat, types = ["sextupole"])
    octs = find_objects(lat, types = ["octupole"])
    cavs = find_objects(lat, types = ["cavity"])
    sols = find_objects(lat, types = ["solenoid"])
    matrices = find_objects(lat, types = ["matrix"])
    marks = find_objects(lat, types = ["marker"])
    bends = find_objects(lat, types = ["bend","rbend", "sbend"])
    #end find objects

    lines = []
    for drift in drifts:
        line = drift.id + " = Drift(l = " + str(drift.l)+ ", id = '"+ drift.id+ "')\n"
        lines.append(line)

    lines.append("\n# quadrupoles \n")
    for quad in quads:
        line = quad.id + " = Quadrupole(l = " + str(quad.l) +", k1 = "+ str(quad.k1) +", id = '"+ quad.id+ "')\n"
        lines.append(line)

    lines.append("\n# bending magnets \n")
    for bend in bends:
        if bend.type == "rbend":
            type = " = RBend(l = "
        elif bend.type == "sbend":
            type = " = SBend(l = "
        else:
            type = " = Bend(l = "
        if bend.k1 == 0 or bend.k1 == None:
            k = ''
        else:
            k = ", k1 = "+ str(bend.k1)
        line = bend.id + type + str(bend.l) + k + ", angle = "+ str(bend.angle)+ ", e1 = " + str(bend.e1) + ", e2 = " + str(bend.e2) + ", tilt = " + str(bend.tilt) +", id = '"+ bend.id+ "')\n"
        lines.append(line)

    lines.append("\n# markers \n")
    for mark in marks:
        line = mark.id + " = Marker(id = '"+ mark.id+ "')\n"
        lines.append(line)

    lines.append("\n# sextupoles \n")
    for sext in sexts:
        line = sext.id + " = Sextupole(l = " + str(sext.l) +", k2 = "+ str(sext.k2) +", tilt = "+ str(sext.tilt) +", id = '"+ sext.id+ "')\n"
        lines.append(line)

    lines.append("\n# octupole \n")
    for oct in octs:
        line = oct.id + " = Octupole(l = " + str(oct.l) +", k3 = "+ str(oct.k3) +", tilt = "+ str(oct.tilt) +", id = '"+ oct.id+ "')\n"
        lines.append(line)

    lines.append("\n# cavity \n")
    for cav in cavs:
        line = cav.id + " = Cavity(l = " + str(cav.l) +", volt = "+ str(cav.v) +", delta_e = "+ str(cav.delta_e)+\
               ", freq = "+ str(cav.f) +", volterr = "+ str(cav.volterr) +", id = '"+ cav.id+ "')\n"
        lines.append(line)

    lines.append("\n# Matrices \n")

    for mat in matrices:
        line = mat.id + " = Matrix(l = " + str(mat.l) +", rm11 = "+ str(mat.rm11) +", rm12 = "+ str(mat.rm12)+\
               ", rm13 = "+ str(mat.rm13) +", rm21 = "+ str(mat.rm21) + ", rm22 = "+ str(mat.rm22) +", rm33 = "+ str(mat.rm33) +\
               ", rm34 = "+ str(mat.rm34) +", rm43 = "+ str(mat.rm43) + ", rm44 = "+ str(mat.rm44) +\
               ", id = '"+ mat.id+ "')\n"
        lines.append(line)

    lines.append("\n# Solenoids \n")

    for sol in sols:
        line = sol.id + " = Solenoid(l = " + str(sol.l) +", k = "+ str(sol.k) +", id = '"+ sol.id+ "')\n"
        lines.append(line)

    lines.append("\n# lattice \n")
    names = []
    for elem in lat.sequence:
        if elem.type != "edge":
            names.append(elem.id)
    #names = map(lambda p: p.id, lat.sequence)
    new_names = []
    for i, name in enumerate(names):
        if i%8 == 7:
            new_names.append("\n"+name)
        else:
            new_names.append(name)
    line = "lattice = (" + ", ".join(new_names) +")"
    lines.append(line)
    return lines


def write_lattice(lattice, file_name = "lattice.inp"):
    lines = lat2input(lattice)

    f = open(file_name, 'w')
    f.writelines(lines)
    f.close()