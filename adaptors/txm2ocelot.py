__author__ = 'Sergey Tomin'


def translate(lines):
    lines2 = []
    for line in lines:
        line = line.replace('asin', "arcsin")
        line = line.replace('acos', "arccos")
        line = line.replace('phi0', "phi")
        line = line.replace('deltae', "delta_e")
        line = line.replace('^', "**")
        line = line.replace('solenoid', "Solenoid")
        line = line.replace('drift', "Drift")
        line = line.replace('quadrupole', "Quadrupole")
        line = line.replace('sbend', "Bend")
        line = line.replace('rbend', "RBend")
        line = line.replace('bend', "Bend")
        line = line.replace('ecol', "Scavenger")


        line = line.replace('octupole', "Scavenger")
        line = line.replace('monitor', "Monitor")
        line = line.replace('matrix', "Matrix")
        line = line.replace('lcavity', "Cavity")
        line = line.replace('rfcavity', "Cavity")
        line = line.replace('sextupole', "Sextupole")
        line = line.replace('marker', "Scavenger")
        line = line.replace('instrument', "Scavenger")
        line = line.replace('rcollimator', "Scavenger")
        line = line.replace('ecollimator', "Scavenger")
        line = line.replace('vkicker', "Scavenger")
        line = line.replace('hkicker', "Scavenger")
        line = line.replace('kicker', "Scavenger")
        line = line.replace('sequence', "Sequence")
        line = line.replace('return', "#return")
        line = line.replace('->', ".")
        line = line.replace('centre', "'centre'")
        lines2.append(line)
        #print line
    return lines2


def find_objects(line):
    if ":" in line and line[0] != "#":
        line_test = line.split("#")[0]
        line_test = line_test.split(":")[1]
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
        params = words[1:]
        if len(', '.join(params) )>0:
            line = name + " = "+ type +"(" + ', '.join(params) + ", id = '" + name + "')"
        else:
            line = name + " = "+ type +"(" + ', '.join(params) + "id = '" + name + "')"
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
            part = line[:ind]
            part = part.replace(" ", "")
            name_budget.append(part)
            #print part
            part = part.replace(".", "_")
            line = part+line[ind:]
        #print line
        for name in name_budget:
            ind = line.find(name)
            if ind>0:
                #print name, line, ind
                name_ = name.replace(".", "_")
                line = line[:ind]+name_+line[ind+len(name):]

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


# def xfel_find_objects(lines):
#     """
#     searching mad's objects. if there ara name and ":" it is object
#     """
#     mad_objs = []
#     for line in lines:
#         if ":" in line and line[0] != "#":
#             madObj = MadObj()
#             i = line.find("#")
#             line2 = line[:i]
#             words = line2.split(",")
#             temp = words[0].split(":")
#             name = temp[0].replace(' ', '')
#             type = temp[1].replace(' ', '')
#             madObj.type = type
#             madObj.name = name
#             madObj.params = words[1:]
#             mad_objs.append(madObj)
#     return mad_objs


def xfel2ocelot(name_file):
    f = open(name_file,"r")
    lines = xfel_line_transform(f)
    new_lines = find_multiline(lines)
    lines = []
    for line in new_lines:
        #print line
        line = find_objects(line)
        #print line
        lines.append(line)
    lines = find_functions(lines)
    lines = find_line(lines)
    lines = find_subroutine(lines)
    lines = translate(lines)
    #mad_objs = xfel_find_objects(lines)
    #mad_objs = subs_objects(mad_objs)
    #mo = parse_obj(mad_objs)
    #new_lines = replace_objects(lines, mo)
    #lines = translate(new_lines)
    #lines = c2py(lines)
    #lines = collect_sequence(lines)
    #for line in lines:
        #print line
    f.close()
    part_name = name_file.split(".")
    part_name[0] += ".inp"
    f_new = open(part_name[0], "w")
    for line in lines:
        #print line
        f_new.write(line+"\n")
    f_new.close()