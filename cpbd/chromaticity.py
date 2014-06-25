__author__ = 'Sergey Tomin'

from numpy import tan, linspace, array, pi, matrix
from optics import trace_z
from scipy.integrate import simps
from numpy.linalg import inv

def edge_chromaticity_old(lattice, tws_0):
    #tested !
    ksi_x_edge = 0.
    ksi_y_edge = 0.
    L = 0.
    for element in lattice.sequence:
        L += element.l
        if element.type == "rbend":
            r = element.l/element.angle
            tw_start = trace_z(lattice,tws_0,[(L - element.l)])
            tw_end = trace_z(lattice,tws_0,[L])
            #print "*************    ", tw_start, tw_end
            ksi_x_edge += (tw_start.beta_x + tw_end.beta_x)*tan(element.angle/2)/r
            ksi_y_edge += (tw_start.beta_y + tw_end.beta_y)*tan(-element.angle/2)/r
    return (ksi_x_edge, ksi_y_edge)

def edge_chromaticity(lattice, tws_0):
    # not tested!
    ksi_x_edge = 0.
    ksi_y_edge = 0.
    L = 0.
    for element in lattice.sequence:
        L += element.l
        if element.type == "edge":
            #r = element.l/element.angle
            tw_start = trace_z(lattice,tws_0,[(L - element.l)])
            #print element.id
            #tw_end = trace_z(lattice,tws_0,[L])
            #print "*************    ", tw_start, tw_end
            ksi_x_edge += tw_start[0].beta_x*tan(element.edge)*element.h
            ksi_y_edge += tw_start[0].beta_y*tan(-element.edge)*element.h
    return (ksi_x_edge, ksi_y_edge)

def chromaticity(lattice, tws_0, nsuperperiod = 1):
    edge_ksi_x, edge_ksi_y = edge_chromaticity(lattice, tws_0)
    tws_elem = tws_0
    #M = TransferMap()
    integr_x = 0.
    integr_y = 0.
    for elem in lattice.sequence:
        if elem.type == "rbend" or elem.type == "sbend" or elem.type == "bend" or elem.type == "quadrupole":
            bx = []
            by = []
            k = []
            h = []
            Z = []
            for z in linspace(0, elem.l,num = 20, endpoint=True):
                twiss_z = elem.transfer_map(z)*tws_elem
                bx.append(twiss_z.beta_x)
                by.append(twiss_z.beta_y)
                k.append(elem.k1)
                if  elem.type != "quadrupole" or elem.l == 0:
                    h.append(elem.angle/elem.l)
                else:
                    h.append(0.)
                Z.append(z)

            H2 = array(h)*array(h)
            X = array(bx)*(array(k)+ H2)
            Y = -array(by)*array(k)
            integr_x += simps(X, Z)
            integr_y += simps(Y, Z)

        tws_elem = elem.transfer_map*tws_elem
    ksi_x = -(integr_x - edge_ksi_x)/(4*pi)
    ksi_y = -(integr_y - edge_ksi_y)/(4*pi)
    return (ksi_x*nsuperperiod, ksi_y*nsuperperiod)


def sextupole_chromaticity(lattice, tws0, nsuperperiod = 1):
    #edge_ksi_x, edge_ksi_y = edge_chromaticity(lattice, tws_0)
    tws_elem = tws0

    integr_x = 0.
    integr_y = 0.
    for elem in lattice.sequence:
        if elem.type == "sextupole":
            bx = []
            by = []
            Dx = []
            Z = []

            for z in linspace(0, elem.l, num = int(elem.l/0.01)+1, endpoint=True):
                twiss_z = elem.transfer_map(z)*tws_elem
                bx.append(twiss_z.beta_x)
                by.append(twiss_z.beta_y)
                Dx.append(twiss_z.Dx)

                Z.append(z)

            X = array(bx)*array(Dx)
            Y = array(by)*array(Dx)
            integr_x += simps(X, Z)*elem.k2
            integr_y += simps(Y, Z)*elem.k2

        tws_elem = elem.transfer_map*tws_elem
    chrom_sex_x = (integr_x)/(4*pi)
    chrom_sex_y = -(integr_y)/(4*pi)
    return (chrom_sex_x*nsuperperiod, chrom_sex_y*nsuperperiod)



def sextupole_id(lattice):
    sex = {}
    for element in lattice.sequence:
        if element.type == "sextupole":
            sex[element.id] = element.ms
    sex_name = sex.keys()
    if len(sex_name)>2:
        exit("more than two families of sextupole")
    return sex

def calculate_sex_strength(lattice, tws_0, ksi, ksi_comp, nsuperperiod):
    # we suppose that sextupole is the thin element l = 0
    ksi_x, ksi_y = ksi
    ksi_x_comp, ksi_y_comp = ksi_comp
    sex = sextupole_id(lattice)
    sex_name = sex.keys()
    print "Chromatism compensate: Before: ", sex
    m1x = 0.
    m1y = 0.
    m2x = 0.
    m2y = 0.
    #L = 0.
    tws_elem = tws_0
    sex_dict_stg = {}
    for element in lattice.sequence:
        if element.type == "sextupole":
            if element.id == sex_name[0]:
                m1x += tws_elem.Dx*tws_elem.beta_x/(4*pi)*nsuperperiod
                m1y -= tws_elem.Dx*tws_elem.beta_y/(4*pi)*nsuperperiod
            elif element.id == sex_name[1]:
                m2x += tws_elem.Dx*tws_elem.beta_x/(4*pi)*nsuperperiod
                m2y -= tws_elem.Dx*tws_elem.beta_y/(4*pi)*nsuperperiod
        tws_elem = element.transfer_map*tws_elem
    M = matrix([ [m1x, m2x], [m1y, m2y] ])
    ksi = matrix([ [ksi_x_comp - ksi_x], [ksi_y_comp - ksi_y] ])
    KSI = inv(M)*ksi
    sex_dict_stg[sex_name[0]] = KSI[0,0]
    sex_dict_stg[sex_name[1]] = KSI[1,0]
    return sex_dict_stg


def compensate_chromaticity(lattice, tws_0, ksi_x_comp = 0, ksi_y_comp = 0,  nsuperperiod = 1):
    ksi_comp = (ksi_x_comp, ksi_y_comp)
    ksi = chromaticity(lattice, tws_0, nsuperperiod)
    sex_dict_stg = calculate_sex_strength(lattice, tws_0, ksi, ksi_comp, nsuperperiod)
    print "Chromatism compensate: After:  ", sex_dict_stg
    #print sex_dict_stg
    sex_name = sex_dict_stg.keys()
    for element in lattice.sequence:
        if element.type == "sextupole":
            if element.id == sex_name[0]:
                element.ms = sex_dict_stg[element.id]
            elif element.id == sex_name[1]:
                element.ms = sex_dict_stg[element.id]
    lattice.update_transfer_maps()