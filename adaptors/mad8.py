__author__ = 'Sergey Tomin'

import numpy as np
def translate(lines):
    lines2 = []
    for line in lines:

        #print line

        line = line.replace('asin', "arcsin")
        line = line.replace('acos', "arccos")
        line = line.replace('phi0', "phi")
        #line = line.replace('deltae', "deltae")
        line = line.replace('^', "**")
        line = line.replace('solenoid', "Solenoid")
        line = line.replace('drift', "Drift")
        line = line.replace('quadrupole', "Quadrupole")
        line = line.replace('quad', "Quadrupole")
        line = line.replace('drif', "Drift")
        line = line.replace('csrcsbend', "SBend")
        line = line.replace('sbend', "Bend")
        line = line.replace('rbend', "RBend")
        line = line.replace('sben', "Bend")
        line = line.replace('bend', "Bend")

        line = line.replace('title', "#title")

        line = line.replace('octupole', "Octupole")
        line = line.replace('blmonitor', "Monitor")
        line = line.replace('imonitor', "Monitor")
        line = line.replace('monitor', "Monitor")
        #line = line.replace('moni', "Monitor")
        line = line.replace('wire', "Monitor")
        #line = line.replace('blmo', "Monitor")
        #line = line.replace('blmo', "Monitor")
        #line = line.replace('imon', "Monitor")
        line = line.replace('matrix', "Matrix")
        line = line.replace('lcavity', "Cavity")
        line = line.replace('lcav', "Cavity")

        line = line.replace('rfcavity', "Cavity")
        line = line.replace('sextupole', "Sextupole")
        line = line.replace('marker', "Marker")
        line = line.replace('mark', "Marker")
        #line = line.replace('rcol', "Marker")

        line = line.replace('profile', "Monitor")
        #line = line.replace('inst', "Monitor")
        line = line.replace('instrument', "UnknownElement")
        line = line.replace('rcollimator', "UnknownElement")
        line = line.replace('ecollimator', "UnknownElement")
        line = line.replace('ecol', "UnknownElement")
        line = line.replace('vkicker', "Vcor")
        line = line.replace('hkicker', "Hcor")
        line = line.replace('kicker', "UnknownElement")
        line = line.replace('hkic', "Hcor")
        line = line.replace('vkic', "Vcor")
        line = line.replace('hkick', "Hcor")
        line = line.replace('vkick', "Vcor")
        line = line.replace('sequence', "Sequence")
        line = line.replace('return', "#return")
        line = line.replace('->', ".")
        line = line.replace('centre', "'centre'")
        line = line.replace('value', "#")
        lines2.append(line)
        #print line
    return lines2


def find_objects(line, info):
    #info = {}
    if ":" in line and line[0] != "#":
        
        try:
            line_test = line.split("#")[0]
            line_test = line_test.split(":")[1]
        except:
            print 'crashed on', line
        #print line_test
        if line_test.find("line")>=0 or line_test.find("subroutine")>0:
            return line
        i = line.find("#")
        if i != -1: line = line[:i]
        line2 = line.replace(' ', '')
        words = line2.split(",")
        temp = words[0].split(":")
        if len(temp)<2:
            return line
        name = temp[0].replace(' ', '')
        type = temp[1].replace(' ', '')
        type = type.replace("\r", "")
        params = words[1:]
        params = [item for item in params if not ("type=" in item or "aper" in item)]
        params = [item.replace("hgap=", "gap=2*") for item in params]
        params = [item.replace("ks=", "k=") for item in params]

        if type in ["lcavity", "lcav"]:
            for i, p in enumerate(params):

                p = p.replace('deltae=', "v=1e-3*")
                p = p.replace('phi0=', "phi=360*")
                p = p.replace('freq=', "freq=1e3*")
                #if args[0] == "type": continue
                params[i] = p

        params = ', '.join(params)
        #print(params)
        if not("at" in params):
            info[name] = {"type": type, "params": params}

            if not (type in ["drift", "drif","sbend", "sben", "rbend", "bend", "rcol", "ecol",
                             "matrix", "quadrupole", "quad", "marker",
                             "hkicker", "vkicker", "hkic","vkic","hkick","vkick","monitor", "moni", "inst", "kicker", "mark", "prof",
                             "sextupole", "lcavity", 'lcav',"ecollimator",
                             "solenoid", 'octupole', "sequence",
                             "beta0", "wire", "blmo", "imon",
                             "instrument", "profile", "imonitor", "blmonitor", "rcollimator"]):
                #print "debug:",  name, ",", type, ",", params
                #print("info", info)
                type = info[type]["type"]

        if len(params )>0:
            line = name + " = "+ type +"(" + params + ", eid = '" + name + "')"
        else:
            line = name + " = "+ type +"(" + params + "eid = '" + name + "')"
        #print type
        if type == "matrix":
            #print line
            line = line.replace('rm(1, 1)', "rm11")
            line = line.replace('rm(1, 2)', "rm12")
            line = line.replace('rm(1, 3)', "rm13")
            line = line.replace('rm(2, 1)', "rm21")
            line = line.replace('rm(2, 2)', "rm22")
            line = line.replace('rm(3, 3)', "rm33")
            line = line.replace('rm(3, 4)', "rm34")
            line = line.replace('rm(4, 3)', "rm43")
            line = line.replace('rm(4, 4)', "rm44")


    return line

def find_subroutine(lines):
    new_lines = []
    n=0
    for line in lines:
        if ":" in line and line[0] != "#":
            if line.find("subroutine")>0:
                n+=1
                i = line.find("#")
                line2 = line[:i]
                line2 = line2.replace(' ', '')
                words = line2.split(",")
                temp = words[0].split(":")
                if len(temp)<2:
                    return line
                name = temp[0].replace(' ', '')
                type = temp[1].replace(' ', '')
                params = words[1:]

                line = "#"+name + " = "+ type +"(" + ', '.join(params) + '):'

                new_lines.append(line)
                continue
        if n>0:

            if line.find("endsubroutine")>0:
                n = 0
                line  = "#"+line
            line = " "*4 + "#" +line
        new_lines.append(line)
    return new_lines

def find_comments(lines):
    new_lines = []
    n=0
    for line in lines:
        if ":" in line and line[0] != "#":
            if line.find("subroutine")>0:
                n+=1
                i = line.find("#")
                line2 = line[:i]
                line2 = line2.replace(' ', '')
                words = line2.split(",")
                temp = words[0].split(":")
                if len(temp)<2:
                    return line
                name = temp[0].replace(' ', '')
                type = temp[1].replace(' ', '')
                params = words[1:]

                line = name + " = "+ type +"(" + ', '.join(params) + '):'
                #new_lines.append("pass")
                new_lines.append(line)
                continue
        if n>0:

            if line.find("endsubroutine")>0:
                n = 0
                line  = "#"+line
            line = " "*4+line
        new_lines.append(line)
    return new_lines

def find_functions(lines):
    new_lines = []
    #n=0
    for line in lines:
        if ":" in line and line[0] != "#":
            parts = line.split(":")
            if "(" in parts[0]:
                parts2 = parts[0].split("(")
                name = parts2[0]
                argums = parts2[1].replace(")", "")
                fpart = name + " = lambda "+ argums + ":"
                spart = parts[1].replace("line", "")
                spart = spart.replace("=", "")
                line = fpart + spart

        new_lines.append(line)
    return new_lines

def find_line(lines):
    new_lines = []
    #n=0
    for line in lines:
        if ":" in line and line[0] != "#":
            parts = line.split(":")
            if parts[1].find("line")>=0:
                line = parts[0] + parts[1].replace("line", "")
                #line = line_test.replace(":", "")
        new_lines.append(line)
    return new_lines


def xfel_line_transform(file):
    """
    replace ":=" by "="
    replace ": CONSTANT =" by " = "
    replace '!' by '#'
    if there is not ";" at the end line, collect the multiline
    all letters in lowercase
    """
    lines = []
    multiline = ''
    name_budget = []
    for line in file:
        #print line
        line = line.replace(" ", "")
        line = line.replace(':=', '=')
        line = line.replace('!', '#')
        ind = line.find("#")
        line = line[:ind]
        if len(line) == 0:
            continue
        if "CONSTANT" in line:
            line = line.replace(":CONSTANT=", "=")

        ind = line.find(":")
        if ind>0 and line[0] != "#":

            parts = line.split(":")
            part = parts[0].replace(" ", "")
            name_budget.append(part)
            part = part.replace(".", "_")
            line = part+line[ind:]

        names = np.array(name_budget)
        names = names[np.argsort(map(len, names))][::-1]
        for name in names:

            name_ = name.replace(".", "_")
            line = line.replace(name, name_)
        line = line.replace(';', '')
        line = line.lower()
        line = line.replace("[",".")
        line = line.replace(']', '')
        #print line
        lines.append(line)
        #print(line)
    #print name_budget
    return lines

def find_multiline(lines):
    multiline = ''
    new_lines = []
    n = 0
    for i, line in enumerate(lines):
        if line.find("&")>0 and line[0] != "#":
            n+=1
            line = line.replace("&", "")
            line = line.replace("\n", " ")
            multiline+= line
            continue
        if line.find("&")<0 and n > 0:
            multiline += line
            n = 0
            new_lines.append(multiline)
            multiline = ''
        else:
            new_lines.append(line)

    return new_lines


def lattice_str_from_mad8(name_file):
    f = open(name_file, "r")
    lines = xfel_line_transform(f)
    new_lines = find_multiline(lines)
    new_lines = find_subroutine(new_lines)
    lines = []
    info = {}
    for line in new_lines:
        line = find_objects(line, info)
        lines.append(line)
    lines = find_functions(lines)
    lines = find_line(lines)

    lines = translate(lines)
    f.close()
    return lines

def save_lattice_str(lines, name_file):
    f_new = open(name_file, "w")
    f_new.write("from ocelot import * \n")
    for line in lines:
        f_new.write(line+"\n")
    f_new.close()




def mad8saveline2lines(fname_saveline):
    f = open(fname_saveline, "r")
    lines = xfel_line_transform(f)
    lines2 = find_multiline(lines)
    lines = []
    info = {}
    for line in lines2:
        line = find_objects(line, info)
        lines.append(line)
    lines = find_line(lines)
    lines = translate(lines)
    f.close()
    return lines

def mad8saveline2ocelot(fname_saveline, fname_ocelot):
    lines = mad8saveline2lines(fname_saveline)
    save_lattice_str(lines, fname_ocelot)


if __name__ == "main":
    # TODO: check XFEL lattice

    from ocelot import *
    from ocelot.cpbd.magnetic_lattice import *
    from numpy import *
    from ocelot.gui.accelerator import *
    import re

    filename = "LCLS.saveline"

    mad8saveline2ocelot(fname_saveline=filename, fname_ocelot="test.py")

    from test import*

    method = MethodTM(params={Cavity: SlacCavityTM})

    lat = MagneticLattice(lsfel, method=method)


    cell = []
    for elem in lat.sequence:
        if elem.id == 'lh_und':
            elem = Undulator(nperiods=0.2531315/0.054, lperiod=0.054, Kx = 1.38523906872,  eid = 'lh_und')
        elif elem.id == 'us33':
            elem = Undulator(nperiods=1.592/0.03, lperiod=0.03, Kx=0, eid=elem.id)
        elif re.match("us+[0-3][0-9]", elem.id):
            elem = Undulator(lperiod=0.03, nperiods=55.5, Kx=3.5, Ky=0.0, eid=elem.id)

        cell.append(elem)

    cell = exclude_zero_length_element(cell, elem_type=[UnknownElement, Marker])
    cell = merge_drifts(cell)
    lat = MagneticLattice(cell, method=method)

    tws = Twiss()
    tws.E = 0.006
    tws.beta_y = 3.909300396E-01
    tws.alpha_y = 5.514326695E-03
    tws.beta_x = 1.557422201E+01
    tws.alpha_x = -3.081460532E+00

    write_lattice(lat, file_name="lcls_lattice.py")

    tws = twiss(lat, tws)
    plot_opt_func(lat, tws)
    plt.show()