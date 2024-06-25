
__author__ = 'Sergey Tomin'

from numpy.linalg import eig
from scipy.integrate import simpson

from ocelot.cpbd.beam import *
from ocelot.cpbd.elements import *
from ocelot.cpbd.optics import *
from ocelot.cpbd.transformations.transformation import TMTypes


def natural_chromaticity(lattice, tws_0, nsuperperiod=1):
    edge_ksi_x, edge_ksi_y = 0, 0
    tws_elem = tws_0
    #M = TransferMap()
    integr_x = 0.
    integr_y = 0.
    for elem in lattice.sequence:
        for tm in elem.first_order_tms:
            if elem.__class__ in [SBend, RBend, Bend, Quadrupole] and tm.tm_type == TMTypes.MAIN:
                bx = []
                by = []
                k = []
                h = []
                Z = []
                for z in np.linspace(0, elem.l, num=5, endpoint=True):
                    section_tms = elem.get_section_tms(delta_l=z, ignore_edges=True)
                    twiss_z = tws_elem
                    for s_tm in section_tms:
                        twiss_z = s_tm * twiss_z
                    bx.append(twiss_z.beta_x)
                    by.append(twiss_z.beta_y)
                    k.append(elem.k1)
                    # print elem.l
                    if elem.__class__ != Quadrupole and elem.l != 0:
                        h.append(elem.angle/elem.l)
                    else:
                        h.append(0.)
                    Z.append(z)

                H2 = np.array(h)*np.array(h)
                X = np.array(bx)*(np.array(k) + H2)
                Y = -np.array(by)*np.array(k)
                integr_x += simpson(X, x=Z)
                integr_y += simpson(Y, x=Z)
            elif elem.__class__ == Multipole:
                twiss_z = tm * tws_elem
                integr_x += twiss_z.beta_x*elem.kn[1]
                integr_y -= twiss_z.beta_y*elem.kn[1]
            tws_elem = tm * tws_elem
    ksi_x = -(integr_x - edge_ksi_x)/(4*pi)
    ksi_y = -(integr_y - edge_ksi_y)/(4*pi)
    return np.array([ksi_x*nsuperperiod, ksi_y*nsuperperiod])


def sextupole_chromaticity(lattice, tws0, nsuperperiod=1):
    tws_elem = tws0

    integr_x = 0.
    integr_y = 0.
    for elem in lattice.sequence:
        if elem.__class__ == Sextupole:
            bx = []
            by = []
            Dx = []
            Z = []

            for z in np.linspace(0, elem.l, num=5, endpoint=True):
                delta_tms = elem.get_section_tms(start_l=0.0, delta_l=z, first_order_only=True)
                twiss_z = tws_elem
                for tm in delta_tms:
                    twiss_z = tm*twiss_z
                    bx.append(twiss_z.beta_x)
                    by.append(twiss_z.beta_y)
                    Dx.append(twiss_z.Dx)

                    Z.append(z)

            X = np.array(bx)*np.array(Dx)
            Y = np.array(by)*np.array(Dx)
            integr_x += simpson(X, x=Z)*elem.k2
            integr_y += simpson(Y, x=Z)*elem.k2

        for tm in elem.first_order_tms:
            tws_elem = tm*tws_elem
    chrom_sex_x = (integr_x)/(4*pi)
    chrom_sex_y = -(integr_y)/(4*pi)
    return np.array([chrom_sex_x*nsuperperiod, chrom_sex_y*nsuperperiod])


def chromaticity(lattice, tws_0, nsuperperiod=1):
    return natural_chromaticity(lattice, tws_0, nsuperperiod) + sextupole_chromaticity(lattice, tws_0, nsuperperiod)


def sextupole_id(lattice):
    sex = {}
    for element in lattice.sequence:
        if element.__class__ == Sextupole:
            if element.id == None:
                print("Please add ID to sextupole")
            sex[element.id] = element.k2*element.l
        if element.__class__ == Multipole and element.n == 3:
            if element.id == None:
                print("Please add ID to sextupole")
            sex[element.id] = element.kn[2]
    sex_name = sex

    if len(sex_name) < 2:
        exit("Only " + str(len(sex_name))+" families of sextupole are found")
    if len(sex_name) > 2:
        exit("more than two families of sextupole")
    return sex


def calculate_sex_strength(lattice, tws_0, ksi, ksi_comp, nsuperperiod):
    '''
    compensation with two families
    '''
    ksi_x, ksi_y = ksi
    ksi_x_comp, ksi_y_comp = ksi_comp
    sex = sextupole_id(lattice)
    sex_name = list(sex)
    print("Chromaticity compensation: Before: ", sex)
    m1x = 0.
    m1y = 0.
    m2x = 0.
    m2y = 0.
    #L = 0.
    tws_elem = tws_0
    sex_dict_stg = {}
    for element in lattice.sequence:
        for tm in element.first_order_tms:
            # if element.type == "sextupole":
            if element.id == sex_name[0]:
                m1x += tws_elem.Dx*tws_elem.beta_x/(4*pi)*nsuperperiod
                m1y -= tws_elem.Dx*tws_elem.beta_y/(4*pi)*nsuperperiod
            elif element.id == sex_name[1]:
                m2x += tws_elem.Dx*tws_elem.beta_x/(4*pi)*nsuperperiod
                m2y -= tws_elem.Dx*tws_elem.beta_y/(4*pi)*nsuperperiod
            tws_elem = tm * tws_elem
    M = np.array([[m1x, m2x], [m1y, m2y]])
    ksi = np.array([[ksi_x_comp - ksi_x], [ksi_y_comp - ksi_y]])
    KSI = np.dot(inv(M), ksi)
    sex_dict_stg[sex_name[0]] = KSI[0, 0]
    sex_dict_stg[sex_name[1]] = KSI[1, 0]
    return sex_dict_stg


def compensate_chromaticity(lattice,  ksi_x_comp=0, ksi_y_comp=0,  nsuperperiod=1):
    '''
    old chromaticity compensation with 2 sextupole families
    '''
    tws0 = Twiss()
    tws = twiss(lattice, tws0)
    tws_0 = tws[0]
    ksi_comp = (ksi_x_comp, ksi_y_comp)
    ksi = natural_chromaticity(lattice, tws_0, nsuperperiod)
    print("ksi_x = ", ksi[0])
    print("ksi_y = ", ksi[1])
    sex_dict_stg = calculate_sex_strength(lattice, tws_0, ksi, ksi_comp, nsuperperiod)
    print("Chromaticity compensatation: After:  ", sex_dict_stg)
    # print sex_dict_stg
    sex_name = list(sex_dict_stg)
    for element in lattice.sequence:
        if element.__class__ == Sextupole:
            if element.id == sex_name[0]:
                element.k2 = sex_dict_stg[element.id] / element.l
                #if element.l != 0: element.k2 = element.ms / element.l
            elif element.id == sex_name[1]:
                element.k2 = sex_dict_stg[element.id] / element.l
                #if element.l != 0: element.k2 = element.ms / element.l
        if element.__class__ == Multipole and element.n == 3:
            if element.id == sex_name[0]:
                element.kn[2] = sex_dict_stg[element.id]
            elif element.id == sex_name[1]:
                element.kn[2] = sex_dict_stg[element.id]

    lattice.update_transfer_maps()


# TODO: finish DR, DZ guys and cromaticity from T matrices


def DZ(lattice, energy):
    tws0 = Twiss()
    tws = twiss(lattice, tws0)
    R = lattice_transfer_map(lattice, energy)
    # print np.array(R[:4, :4])- np.eye(4),R[:4, 5]

    x = np.dot(R, [1, 1, 1, 1, 0, 0])
    x2 = np.dot(R, [1, 1, 1, 1, 0, 0.001])
    #print (x2 - x)/0.001
    # print R

    w, v = eig(R)
    # print R
    v1 = np.dot(R, v[0].real)
    # print "v0 = ", v[0].real
    # print "v*R = ", v1
    DZ = np.dot(inv(-R[:4, :4] + np.eye(4)), R[:4, 5])
    DZ = np.append(DZ, [0, 1])
    print("DZ0 = ", DZ)
    DR_new = np.zeros(6)
    DZ0 = deepcopy(DZ)
    DQ1n = 0.
    DQ2n = 0.
    R_lat = np.eye(6)
    for elem in lattice.sequence:
        R = elem.transfer_map.R(energy)
        R_lat = np.dot(R, R_lat)
        #cosmx = (R[0, 0] + R[1, 1])/2.
        #cosmy = (R[2, 2] + R[3, 3])/2.

        # if abs(cosmx) >= 1 or abs(cosmy) >= 1:
        #    logger.warn("************ periodic solution does not exist. return None ***********")
        #    #print("************ periodic solution does not exist. return None ***********")
        #    return None
        # nsinmx = np.sign(R[0, 1])*sqrt(R[0,1]*R[1, 0] - (R[0, 0]*R[1,1])/4.)
        # nsinmy = np.sign(R[2, 3])*sqrt(R[2,3]*R[3, 2] - (R[2, 2]*R[3,3])/4.)

        DZ = np.dot(R_lat, DZ0)
        DR = np.zeros((6, 6))
        T = elem.transfer_map.t_mat_z_e(elem.l, energy)
        for i in range(6):
            for j in range(6):
                t = 0.
                for m in range(6):
                    t += T[i, j, m]*DZ[m]
                DR[i, j] = 2*t  # np.dot(lattice.T[i, j, :], DZ)
        print(elem.l)
        print(DR)
        DR_new = DR + DR_new
        #DQ1n += (DR[0, 0] + DR[1, 1])/(2*sin(tws[-1].mux))/2/pi
        #DQ2n += (DR[2, 2] + DR[3, 3])/(2*sin(tws[-1].muy))/2/pi
    print("DZ = ",  DZ)
    # print "DZ2 = ", dot(R, DZ)
    print(DQ1n*6, DQ2n*6)
    DR = np.zeros((6, 6))

    for i in range(6):
        for j in range(6):
            t = 0.
            for m in range(6):
                t += lattice.T[i, j, m]*DZ[m]
            DR[i, j] = 2*t  # np.dot(lattice.T[i, j, :], DZ)
    # print DR
    DR = DR_new

    DQ1 = (DR[0, 0] + DR[1, 1])/(2*sin(tws[-1].mux))/2/pi
    DQ2 = (DR[2, 2] + DR[3, 3])/(2*sin(tws[-1].muy))/2/pi
    print(DQ1*6, DQ2*6)
    print(tws[-1].mux/2/pi*6)
    print(tws[-1].muy/2/pi*6)
