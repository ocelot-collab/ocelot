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
        line = line.replace('monitor', "Monitor")
        line = line.replace('moni', "Monitor")
        line = line.replace('inst', "Monitor")
        line = line.replace('matrix', "Matrix")
        line = line.replace('lcavity', "Cavity")
        line = line.replace('lcav', "Cavity")

        line = line.replace('rfcavity', "Cavity")
        line = line.replace('sextupole', "Sextupole")
        line = line.replace('marker', "Marker")
        line = line.replace('mark', "Marker")
        line = line.replace('rcol', "Marker")
        line = line.replace('ecol', "Marker")
        line = line.replace('prof', "Monitor")
        line = line.replace('instrument', "UnknownElement")
        line = line.replace('rcollimator', "UnknownElement")
        line = line.replace('ecollimator', "UnknownElement")
        line = line.replace('ecol', "UnknownElement")
        line = line.replace('vkicker', "UnknownElement")
        line = line.replace('hkicker', "UnknownElement")
        line = line.replace('kicker', "UnknownElement")
        line = line.replace('hkic', "UnknownElement")
        line = line.replace('vkic', "UnknownElement")
        line = line.replace('hkick', "UnknownElement")
        line = line.replace('vkick', "UnknownElement")
        line = line.replace('sequence', "Sequence")
        line = line.replace('return', "#return")
        line = line.replace('->', ".")
        line = line.replace('centre', "'centre'")
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
        line2 = line[:i]
        line2 = line2.replace(' ', '')
        words = line2.split(",")
        temp = words[0].split(":")
        if len(temp)<2:
            return line
        name = temp[0].replace(' ', '')
        type = temp[1].replace(' ', '')
        #print name, type
        params = words[1:]
        if type == "lcavity":
            for i, p in enumerate(params):
                args = p.split("=")
                args[0] = args[0].replace('deltae', "delta_e")
                params[i] = "=".join(args)
            #print params
        params = ', '.join(params)
        if not("at" in params):
            info[name] = {"type": type, "params": params}

            if not (type in ["drift", "drif","sbend", "sben", "rbend", "bend", "rcol", "ecol",
                             "matrix", "quadrupole", "quad", "marker",
                             "hkicker", "vkicker", "hkic","vkic","hkick","vkick","monitor", "moni", "inst", "kicker", "mark", "prof",
                             "sextupole", "lcavity", 'lcav',"ecollimator",
                             "solenoid", 'octupole', "sequence"]):
                print "debug:",  name, ",", type, ",", params
                type = info[type]["type"]

        if len(params )>0:
            line = name + " = "+ type +"(" + params + ", id = '" + name + "')"
        else:
            line = name + " = "+ type +"(" + params + "id = '" + name + "')"
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
                #new_lines.append("    pass")
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
    replace '!' by '#'
    if there is not ";" at the end line, collect the multiline
    all letters in lowercase
    """
    lines = []
    multiline = ''
    name_budget = []
    for line in file:
        #print line
        line = line.replace(':=', '=')
        line = line.replace('!', '#')

        #line = line.replace('&', '\\')
        ind = line.find(":")
        if ind>0 and line[0] != "#":
            parts = line.split(":")
            #print parts
            #type = parts[1].split(",")[0]
            part = parts[0].replace(" ", "")
            name_budget.append(part)
            #print part, type
            part = part.replace(".", "_")
            line = part+line[ind:]
        #print line
        names = np.array(name_budget)
        names = names[np.argsort(map(len, names))][::-1]
        #if "i1.qih.20" in names:
        #    print names
        for name in names:

            name_ = name.replace(".", "_")
            line = line.replace(name, name_)
            #ind = line.find(name)
            #if ind>0:
            #    #print name, line, ind
            #    name_ = name.replace(".", "_")
            #    line = line[:ind]+name_+line[ind+len(name):]

        line = line.replace(';', '')
        line = line.lower()
        line = line.replace("[",".")
        line = line.replace(']', '')
        #print line
        lines.append(line)
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
    lines = []
    info = {}
    for line in new_lines:
        #print line
        line = find_objects(line, info)
        #print line
        lines.append(line)
    lines = find_functions(lines)
    lines = find_line(lines)
    lines = find_subroutine(lines)
    lines = translate(lines)
    f.close()
    return lines

def save_lattice_str(lines, name_file):
    #part_name = name_file.split(".")
    #part_name[0] += ".inp"
    f_new = open(name_file, "w")
    for line in lines:
        #print line
        f_new.write(line+"\n")
    f_new.close()