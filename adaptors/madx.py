__author__ = 'Sergey Tomin'

from ocelot.cpbd.elements import Drift

class MadObj:

    name = ""
    type = ""
    params = ""

    def parse_params(self):
        self.dic_p = {}
        for param in self.params:
            param = param.replace('=', ' ')
            param = param.split()
            self.dic_p[param[0]] = param[1]


def line_transform(file):
    """
    replace ":=" by "="
    replace '!' by '#'
    if there is not ";" at the end line, collect the multiline
    all letters in lowercase
    """
    lines = []
    multiline = ''
    for line in file:
        #print line
        line = line.replace(':=', '=')
        line = line.replace('!', '#')
        multiline += line
        if multiline.find(";")<0 and multiline[0] != "#":
            ind = line.find("#")
            multiline = multiline[:ind]
            continue
        else:
            line = multiline
            multiline = ''
        line = line.replace(';', '')
        line = line.lower()
        #print line
        lines.append(line)
    return lines

def find_objects(lines):
    """
    searching mad's objects. if there ara name and ":" it is object
    """
    mad_objs = []
    for line in lines:
        if ":" in line and line[0] != "#":
            madObj = MadObj()
            i = line.find("#")
            line2 = line[:i]
            words = line2.split(",")
            temp = words[0].split()
            madObj.type = temp[-1]
            madObj.name = temp[0].replace(":", "")
            madObj.params = words[1:]
            mad_objs.append(madObj)
    return mad_objs


def subs_objects(mad_objs):
    for i, mo in enumerate(mad_objs):
        for mo2 in mad_objs[i:]:
            if mo.name == mo2.type:
                mo2.type = mo.type
                params = mo.params + mo2.params
                mo2.params = params
    return mad_objs


def parse_obj(mad_objs):
    for mo in mad_objs:
        #print mo.name, mo.type, mo.params
        mo.parse_params()
        #print  mo.name, mo.type, mo.dic_p
    return mad_objs


def replace_objects(lines, mad_objs):
    new_lines = []
    multy_line = ''
    for line in lines:

        if ":" in line and line[0] != "#":
            i = line.find("#")
            line2 = line[:i]
            words = line2.split(",")
            temp = words[0].split()
            name = temp[0].replace(":", "")

            for mo in mad_objs:
                if name == mo.name:
                    line = ""
                    line = mo.name + " = " + mo.type + "("
                    for p in mo.dic_p:
                        line += p + " = " + mo.dic_p[p] +", "
                    line += ")"
        line = line.replace(', )', ')')
        new_lines.append(line)
    return new_lines

def translate(lines):
    lines2 = []
    for line in lines:
        line = line.replace('quadrupole', "Quadrupole")
        line = line.replace('sbend', "Bend")
        #line = line.replace('rbend', "RBend")
        #line = line.replace('bend', "Bend")
        line = line.replace('monitor', "Monitor")
        line = line.replace('matrix', "Matrix")
        line = line.replace('rfcavity', "RFcavity")
        line = line.replace('sextupole', "Sextupole")
        line = line.replace('marker', "Marker")
        line = line.replace('instrument', "Scavenger")
        line = line.replace('rcollimator', "Scavenger")
        line = line.replace('ecollimator', "Scavenger")
        line = line.replace('vkicker', "Scavenger")
        line = line.replace('hkicker', "Scavenger")
        line = line.replace('sequence', "Sequence")
        line = line.replace('return', "#return")
        line = line.replace('->', ".")
        line = line.replace('//', "!")
        line = line.replace('centre', "'centre'")
        line = line.replace('at =', "at=")
        lines2.append(line)
        #print line
    return lines2

def c2py(lines):
    lines2 = []
    for line in lines:
        #remove spases
        for i in range(8,2,-1):
            line = line.replace(' '*i, "\t")
        line = line.replace('\t', " "*4)

        # make operator body of "for" or "if" li
        line = line.replace('{', ": # start operator")
        line = line.replace('}', "#end operator")
        line = line.replace('print, text=', "print ")
        line = line.replace('print, text =', "print ")
        lines2.append(line)
    return lines2


def collect_sequence(lines):
    seq = []
    first = 1
    lines2 = []
    for line in lines:
        if line.find("at=")>0:
            if line[0] == "#":
                continue
            parts = line.split("at=")
            name = parts[0].replace(",", "")
            name = name.split()[0]
            #print name
            pos = parts[1].replace("\n", "")

            id = " '" + name + "' "
            line = ",".join([name,id, pos])

            ind = line.find("#")
            if ind>0:
                line = line[:ind]
                #print "ind = ", ind == True
            #print line
            if first:
                line = "lattice = [[" + line +"],"
                first = 0
            else:
                line = "[" + line + "],"
            #print line
            seq.append(line)
        line = line.replace("endSequence", "]")
        lines2.append(line)
    return lines2


def lattice_str_from_madx(name_file):
    f = open(name_file,"r")
    lines = line_transform(f)
    mad_objs = find_objects(lines)

    mad_objs = subs_objects(mad_objs)
    mo = parse_obj(mad_objs)
    new_lines = replace_objects(lines, mo)
    lines = translate(new_lines)
    lines = c2py(lines)
    lines = collect_sequence(lines)
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


def madx_seq2ocelot_seq(list_elem_pos, tot_length, exclude_elems = []):
    seq = []
    azimuth = 0.
    for i, term in enumerate(list_elem_pos):
        if term[1] in exclude_elems:
            continue
        element = term[0]
        #print element
        element.id = term[1]
        pos = term[2]
        drift = Drift(l = pos - element.l/2. - azimuth, id = "drift_" + str(i))
        azimuth = pos + element.l/2.
        seq.append(drift)
        seq.append(element)
        #print elem[0].l, elem[1], elem[2]
    len_last_drift = tot_length - list_elem_pos[-1][-1] - list_elem_pos[-1][0].l/2.
    drift = Drift(l = len_last_drift, id = "drift_last")
    seq.append(drift)
    return seq


if __name__ == "__main__":
    #name_file = "quadsex_20wig.dat"
    name_file = 'petra3.txt'
    f = open(name_file,"r")
    lines = line_transform(f)
    mad_objs = find_objects(lines)
    mad_objs = subs_objects(mad_objs)
    mo = parse_obj(mad_objs)
    new_lines = replace_objects(lines, mo)
    lines = translate(new_lines)
    lines = c2py(lines)
    lines = collect_sequence(lines)
    for line in lines:
        print line
    f.close()
    part_name = name_file.split(".")
    part_name[0] += ".py"
    f_new = open(part_name[0], "w")
    for line in lines:
        f_new.write(line+"\n")
    f_new.close()
























