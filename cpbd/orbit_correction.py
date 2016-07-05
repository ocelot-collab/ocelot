__author__ = 'Sergey Tomin'


import pickle
from time import sleep

import matplotlib.pyplot as plt
from numpy import diag, shape
from numpy.linalg import svd
from scipy.interpolate import splrep, splev

from ocelot.cpbd.match import closed_orbit
from ocelot.cpbd.track import *
from ocelot.gui.accelerator import *
import copy

def show_currents( elems, alpha):
    print( "******* displaying currents - START ********")
    for elem in elems:
        if elem.dI == 0:
            continue
        n = len(elem.id)
        n2 = len(str(elem.I + elem.dI))
        n3 = len(str(elem.I))
        print( elem.id, " "*(10-n) + "<-- ", elem.I + elem.dI,  " "*(18-n2)+ " was = ", elem.I, " "*(18-n3) + " dI = ", elem.dI, "x", alpha)
    print ("******* displaying currents - END ********")
"""
class BPM(object):
    def __init__(self, id = None):
        self.id = id
        self.type = "monitor"
        self.__x = 0
        self.__y = 0
        self.dx = 0
        self.dy = 0
        self.s = 0
        self.beta_x = 0
        self.beta_y = 0
        self.phi_x = 0
        self.phi_y = 0

    @property
    def x(self):
        #print 'getter x'
        return self.__x

    @x.setter
    def x(self, value):
        a = 1.
        b = c = 0.
        self.__x = a*value + b + c

    @property
    def y(self):
        #print 'getter x'
        return self.__y

    @y.setter
    def y(self, value):
        a = 1.
        b = c = 0.
        self.__y = a*value + b + c

    def save(self):
        self.x0 = self.x
        self.y0 = self.y
"""

class Response_matrix:
    def __init__(self):
        self.cor_names = []
        self.bpm_names = []
        self.matrix = []
        self.mode = "radian"  # or "ampere"

    def save(self, filename):
        dict_rmatrix = {}
        dict_rmatrix["cor_names"] = self.cor_names
        dict_rmatrix["bpm_names"] = self.bpm_names
        dict_rmatrix["matrix"] = self.matrix
        dict_rmatrix["mode"] = self.mode
        pickle.dump(dict_rmatrix, open(filename, "wb"))

    def load(self, filename):
        dict_rmatrix = pickle.load(open(filename, "rb"))
        self.cor_names = dict_rmatrix["cor_names"]
        self.bpm_names = dict_rmatrix["bpm_names"]
        self.matrix = dict_rmatrix["matrix"]
        self.mode = dict_rmatrix["mode"]
        #print (self.mode)
        return 1

    def extract(self, cor_list, bpm_list):
        cor_list = np.array(cor_list)
        bpm_list = np.array(bpm_list)
        cors1 = np.array(self.cor_names)
        cors2 = cor_list
        bpms1 = np.array(self.bpm_names)
        bpms2 = bpm_list
        nb1 = len(bpms1)

        c_names = cors1[np.in1d(cors1, cors2)]

        c_i1 = np.where(np.in1d(cors1, cors2))[0]
        c_i2 = np.where(np.in1d(cors2, c_names))[0]
        #print bpms1, np.in1d(bpms1, bpms2)
        b_names = bpms1[np.in1d(bpms1, bpms2)]
        #print b_names
        b_i1 = np.where(np.in1d(bpms1, bpms2))[0]
        b_i2 = np.where(np.in1d(bpms2, b_names))[0]

        if not np.array_equal(c_names, cor_list):
            print (" in origin response matrix does not exist correctors:")
            #print c_names
            print (cors2[np.in1d(cors2, c_names, invert=True)])
        if not np.array_equal(b_names, bpm_list):
            print (" in origin response matrix does not exist BPMs:")
            print (bpm_list[b_i2[:]])

        extr_matrix = np.zeros((len(b_names)*2, len(c_names)))
        #plane = ["X", "Y"]
        for n in range(2):
            #print "****************   ", plane[n], "   ****************"
            for i, c in enumerate(c_names):
                for j, b in enumerate(b_names):
                    #print b_i1[j],  nb1*n, c_i1[i]
                    x1 = self.matrix[b_i1[j] + nb1*n, c_i1[i]]
                    extr_matrix[j + n*len(b_names), i] = x1

        rmatrix = Response_matrix()
        rmatrix.cor_names = c_names
        rmatrix.bpm_names = b_names
        rmatrix.matrix = extr_matrix
        rmatrix.mode = self.mode
        return rmatrix


    #def measure(self, mi, dp):
    #    I0 = mi.init_corrector_vals(self.cor_names)
    #    for cor in self.cor_names:
    #        I0 = mi.init_corrector_vals([cor])[0]
    #        print cor, " I0 = ", I0
    #        I0array.append(I0)


    def compare(self, rmatrix, absolut = 0.001, relative = 0.1):
        cors1 = np.array(self.cor_names)
        cors2 = np.array(rmatrix.cor_names)
        bpms1 = np.array(self.bpm_names)
        bpms2 = np.array(rmatrix.bpm_names)
        nb1 = len(bpms1)
        nb2 = len(bpms2)
        #c_names = np.intersect1d(cors1, cors2)
        c_names = cors1[np.in1d(cors1, cors2)]
        c_i1 = np.where(np.in1d(cors1, cors2))[0]
        c_i2 = np.where(np.in1d(cors2, cors1))[0]
        #b_names = np.intersect1d(bpms1, bpms2)
        b_names = bpms1[np.in1d(bpms1, bpms2)]
        b_i1 = np.where(np.in1d(bpms1, bpms2))[0]
        b_i2 = np.where(np.in1d(bpms2, bpms1))[0]
        plane = ["X", "Y"]
        for n in range(2):
            print ("****************   ", plane[n], "   ****************")
            counter = 0
            for i, c in enumerate(c_names):
                for j, b in enumerate(b_names):
                    #print b_i1[j],  nb1*n, c_i1[i]
                    x1 = self.matrix[b_i1[j] + nb1*n, c_i1[i]]
                    x2 = rmatrix.matrix[b_i2[j] + nb2*n, c_i2[i]]
                    if abs(x1 - x2) <absolut:
                        continue
                    if abs(x1 - x2)/max(np.abs([x1, x2])) < relative:
                        continue
                    l_x1 = len(str(x1))
                    print (plane[n], c, " "*(10 - len(c)), b, " "*(10 - len(b)), "r1: ", x1," "*(18 - l_x1),"r2: ", x2)
                    counter += 1
            print("shown", counter, "elements of", len(c_names)*len(b_names))

    def show(self, list_cor=None, list_bpm=None):
        print (" "*10,)
        for bpm in self.bpm_names:
            print (bpm,)
        print()
        for i in range(shape(self.matrix)[1]):
            print (self.cor_names[i] + " "*(10 - len(self.cor_names[i])),)
            #print np.array_str(self.matrix[:, i], precision=2, suppress_small=True)
            for j in range(shape(self.matrix)[0]):
                print ("%.2f" % self.matrix[j, i],)
            print()



class Orbit:
    def __init__(self, lattice):
        #self.monitors = []
        #self.correctors = []
        self.lat = lattice
        #self.h_resp = []
        #self.v_resp = []
        self.bpms = []
        self.hcors = []
        self.vcors = []
        self.nu_x = 0.
        self.nu_y = 0.
        self.resp = []
        self.mode = "radian" # or "ampere"
        #if lattice != None:
        self.create_BPM()
        self.create_COR()

    def create_BPM(self, bpm_list=None):
        """
        Search bpm in the lattice and create list of bpms
        :param lattice: class MagneticLattice
        :return: self.bpms - list of BPMs (class BPM)
        """
        self.bpms = []
        L = 0.
        for elem in self.lat.sequence:
            if elem.__class__ == Monitor:
                if bpm_list is None or elem.id in bpm_list:
                    try:
                        elem.weight
                    except:
                        elem.weight = 1.
                    elem.s = L+elem.l/2.
                    elem.x_ref = 0.
                    elem.y_ref = 0.
                    self.bpms.append(elem)
            L += elem.l
        if len(self.bpms) == 0:
            print("there are not monitors")
        return self.bpms

    def set_ref_pos(self):
        for bpm in self.bpms:
            bpm.x_ref = bpm.x
            bpm.y_ref = bpm.y

    def minus_reference(self):
        for bpm in self.bpms:
            #print bpm.x, bpm.x_ref
            bpm.x = bpm.x - bpm.x_ref
            bpm.y = bpm.y - bpm.y_ref

    def set_bpm_signal(self, x_bpm, y_bpm):
        for i, bpm in enumerate(self.bpms):
            #print bpm.x, bpm.x_ref
            bpm.x = x_bpm[i]
            bpm.y = y_bpm[i]

    def create_COR(self, cor_list=None):
        """
        Search correctors (horizontal and vertical) in the lattice and create list of hcors and list of vcors
        :param lattice: class MagneticLattice
        :return:
        """
        self.hcors = []
        self.vcors = []
        L = 0.
        for elem in self.lat.sequence:
            if elem.__class__ == Vcor:
                if cor_list is None or elem.id in cor_list:
                    elem.s = L+elem.l/2.
                    self.vcors.append(elem)
            elif elem.__class__ == Hcor:
                if cor_list is None or elem.id in cor_list:
                    elem.s = L+elem.l/2.
                    self.hcors.append(elem)
            L += elem.l
        if len(self.hcors) == 0:
            print("there are not horizontal correctors")
        if len(self.vcors) == 0:
            print("there are not vertical correctors")


    def create_types(self, types, remove_elems=[]):
        self.htypes = []
        self.vtypes = []
        L = 0.
        for elem in self.lat.sequence:
            L += elem.l
            if elem.__class__ in types:
                if "_U" in elem.id:
                    continue
                if elem.id in remove_elems:
                    continue
                elem.s = L - elem.l/2.
                self.htypes.append(elem)
                self.vtypes.append(elem)

    def export_response_matrix(self, r_matrix):
        self.create_BPM(bpm_list=r_matrix.bpm_names)
        self.create_COR(cor_list=r_matrix.cor_names)
        self.resp = r_matrix.matrix
        self.mode = r_matrix.mode

    def read_virtual_orbit(self, p_init=None):
        """
        searching closed orbit by function closed_orbit(lattice) and searching coordinates of beam at the bpm possitions
        :param lattice: class MagneticLattice
        :return: orbit.bpms
        """
        X = []
        Y = []
        if p_init == None:
            self.particle0 = closed_orbit(self.lat)
        else:
            self.particle0 = p_init
        #print "particle2 = ", self.particle0.s, self.particle0.x
        p = copy.copy(self.particle0)
        navi = Navigator()
        L = 0.
        for bpm in self.bpms:
            #print("energy = ", p.E)
            dz = bpm.s - L
            tracking_step(self.lat, [p], dz, navi)
            bpm.x = p.x
            bpm.y = p.y
            bpm.E = p.E
            L = bpm.s
            X.append(p.x)
            Y.append(p.y)
        #print("energy = ", p.E)
        return array(X), array(Y)

    def response_matrix(self, mi, dp, timeout=0.5, delta_i=0.01):
        resp = np.zeros((len(self.bpms)*2, len(self.hcors)+len(self.vcors)))
        plane = ["X", "Y"]
        bpm_names = [b.id for b in self.bpms]
        for n, cor_list in enumerate([self.hcors, self.vcors]):
            #print "cor_list = ", cor_list
            for cor in cor_list:
                #print "cor = ", cor
                i = mi.init_corrector_vals([cor.id])
                cor.I = i[0]
                try:
                    cor.dI
                except:
                    #print cor.id, " delta_i = ", delta_i
                    cor.dI = delta_i
                print ("X:  ", cor.id, "I = ", cor.I)

            show_currents(cor_list, alpha=1.)
            inp = raw_input("Start measurement of response matrix for " + plane[n]+":? ")
            if inp == "yes":
                #resp = np.zeros(len(self.bpms)*2, len(cor_list))
                X0, Y0 = mi.get_bpms_xy(bpm_names)
                XY0 = np.append(X0, Y0)
                for i, cor in enumerate(cor_list):
                    print (i, "/", len(cor_list), cor.id)
                    mi.set_value(cor.id, cor.I + cor.dI)
                    sleep(timeout)
                    X, Y = mi.get_bpms_xy(bpm_names)
                    XY = np.append(X, Y)
                    #print "XY = ", XY, XY0
                    #print (XY - XY0)/cor.dI
                    resp[:, i + n*len(self.hcors)] = (XY - XY0)/cor.dI
                    mi.set_value(cor.id, cor.I - cor.dI)
                    sleep(timeout)
                    X0, Y0 = mi.get_bpms_xy(bpm_names)
                    XY0 = np.append(X0, Y0)
        rmatrix = Response_matrix()
        rmatrix.bpm_names = [b.id for b in self.bpms]
        rmatrix.cor_names = np.append(np.array([c.id for c in self.hcors]), np.array([c.id for c in self.vcors]))
        rmatrix.matrix = resp
        rmatrix.mode = "ampere"
        return rmatrix

    def optical_func_params(self, tw_init=None):
        """
        Optical function parameters for correctors and bpms. It is needed for calculation of ideal response matrix:
        defining beta functions on the azimuth of correctors and bpms: beta_x, beta_y;
        defining phase shift between origin of lattice and element: mu_x, mu_y;
        defining tunes of beta functions of whole lattice: nu_x, nu_y = mu(totalLen)/(2*pi)
        :param lattice: class MagneticLattice
        :return:
        """
        if tw_init == None:
            tw_init = Twiss()

        tws = twiss(self.lat, tw_init, nPoints=int(self.lat.totalLen/0.05))
        s = array([tw.s for tw in tws])
        tck_mux = splrep(s, array([tw.mux for tw in tws]))
        tck_muy = splrep(s, array([tw.muy for tw in tws]))

        beta_x = array([tw.beta_x for tw in tws])
        beta_y = array([tw.beta_y for tw in tws])
        self.nu_x = tws[-1].mux/2./pi
        self.nu_y = tws[-1].muy/2./pi

        tck_bx = splrep(s, beta_x)
        tck_by = splrep(s, beta_y)
        energy_s = [tw.E for tw in tws]
        tck_E = splrep(s, energy_s)
        #s_bpm = [bpm.s for bpm in self.bpms]

        for bpm in self.bpms:
            bpm.phi_x = splev([0, bpm.s], tck_mux)[1]
            bpm.phi_y = splev([0, bpm.s], tck_muy)[1]

            bpm.beta_x = splev([0, bpm.s], tck_bx)[1]
            bpm.beta_y = splev([0, bpm.s], tck_by)[1]
            bpm.E = splev([0, bpm.s], tck_E)[1]

        for hcor in self.hcors:
            hcor.phi_x = splev([0, hcor.s], tck_mux)[1]
            hcor.phi_y = splev([0, hcor.s], tck_muy)[1]

            hcor.beta_x = splev([0, hcor.s], tck_bx)[1]
            hcor.beta_y = splev([0, hcor.s], tck_by)[1]
            hcor.E = splev([0, hcor.s], tck_E)[1]

        for vcor in self.vcors:
            vcor.phi_x = splev([0, vcor.s], tck_mux)[1]
            vcor.phi_y = splev([0, vcor.s], tck_muy)[1]

            vcor.beta_x = splev([0, vcor.s], tck_bx)[1]
            vcor.beta_y = splev([0, vcor.s], tck_by)[1]
            vcor.E = splev([0, vcor.s], tck_E)[1]

    def measure_response_matrix(self, p_init=None, match_ic=False, order=1):
        """
        :param lattice:
        :param p_init:
        :param match_ic: matching initial coordinates of particles
        :return:
        """
        print ("measure = ", p_init.x, p_init.y, p_init.px, p_init.py, p_init.E)
        shift = 0.0001
        m = len(self.bpms)
        nx = len(self.hcors)
        ny = len(self.vcors)
        if match_ic:
            real_resp = zeros((m*2, nx + ny+4))
        else:
            real_resp = zeros((m*2, nx + ny))
        self.read_virtual_orbit(p_init=copy.deepcopy(p_init))
        bpms = copy.deepcopy(self.bpms)

        for ix, hcor in enumerate(self.hcors):
            print("measure X - ", ix, "/", nx)
            hcor.angle = shift
            self.lat.update_transfer_maps()
            self.read_virtual_orbit(p_init=copy.deepcopy(p_init))

            for j, bpm in enumerate(self.bpms):
                real_resp[j, ix] = (bpm.x - bpms[j].x)/shift
                real_resp[j+m, ix] =(bpm.y - bpms[j].y)/shift
            hcor.angle = 0
        self.lat.update_transfer_maps()

        for iy, vcor in enumerate(self.vcors):
            print("measure Y - ", iy,"/",ny)
            vcor.angle = shift
            self.lat.update_transfer_maps()
            self.read_virtual_orbit(p_init=copy.deepcopy(p_init))

            for j, bpm in enumerate(self.bpms):
                real_resp[j, iy+nx] = (bpm.x - bpms[j].x)/shift
                real_resp[j+m, iy+nx] = (bpm.y - bpms[j].y)/shift
            vcor.angle = 0
        self.lat.update_transfer_maps()
        self.read_virtual_orbit(p_init=copy.deepcopy(p_init))
        if match_ic:
            for i, par in enumerate(["x", "px", "y", "py"]):
                print(i)
                p_i = Particle(E = p_init.E)
                p_i.__dict__[par] = 0.0001
                p2 = copy.deepcopy(p_i)
                print ("measure = ", p2.x, p2.y, p2.px, p2.py, p2.E)
                self.read_virtual_orbit(p_init=p2, order=order)
                for j, bpm in enumerate(self.bpms):
                    real_resp[j, nx + ny + i] = (bpm.x - bpms[j].x)/0.0001
                    real_resp[j+m, nx + ny + i] = (bpm.y - bpms[j].y)/0.0001
        self.resp = real_resp
        rmatrix = Response_matrix()
        rmatrix.bpm_names = [b.id for b in self.bpms]
        rmatrix.cor_names = np.append(np.array([c.id for c in self.hcors]), np.array([c.id for c in self.vcors]))
        rmatrix.matrix = self.resp
        rmatrix.mode = "radian"
        return rmatrix

    def ring_response_matrix(self, tw_init=None):
        """
        calculation of ideal response matrix
        :param lattice: class MagneticLattice
        :param tw_init: if tw_init == None, function tries to find periodical solution
        :return: orbit.resp
        """
        self.optical_func_params(tw_init=tw_init)

        m = len(self.bpms)
        nx = len(self.hcors)
        ny = len(self.vcors)
        h_resp = zeros((m, nx))
        v_resp = zeros((m, ny))
        sin_pnu_x = sin(pi*self.nu_x)
        sin_pnu_y = sin(pi*self.nu_y)
        for i, bpm in enumerate(self.bpms):
            kx = sqrt(bpm.beta_x)/(2.*sin_pnu_x)
            ky = sqrt(bpm.beta_y)/(2.*sin_pnu_y)
            for j, hcor in enumerate(self.hcors):
                mu_x = abs(bpm.phi_x - hcor.phi_x)
                h_resp[i,j] = kx*sqrt(hcor.beta_x)*cos(mu_x - pi*self.nu_x)
            for n, vcor in enumerate(self.vcors):
                mu_y = abs(bpm.phi_y - vcor.phi_y)
                v_resp[i,n] = ky*sqrt(vcor.beta_y)*cos(mu_y - pi*self.nu_y)

        m = len(self.bpms)
        kx = len(self.hcors)
        ky = len(self.vcors)
        self.resp = zeros((2*m, kx + ky))
        self.resp[:m,:kx] = h_resp[:,:]
        self.resp[m:,kx:] = v_resp[:,:]
        #print "shape = ", shape(self.resp)
        return self.resp

    def linac_response_matrix(self, tw_init=None):
        """
        calculation of ideal response matrix
        :param lattice: class MagneticLattice
        :param tw_init: if tw_init == None, function tries to find periodical solution
        :return: orbit.resp
        """
        self.optical_func_params(tw_init=tw_init)

        m = len(self.bpms)
        nx = len(self.hcors)
        ny = len(self.vcors)
        h_resp = zeros((m, nx))
        v_resp = zeros((m, ny))

        for i, bpm in enumerate(self.bpms):
            kx = sqrt(bpm.beta_x)#/(2.*sin_pnu_x)
            ky = sqrt(bpm.beta_y)#/(2.*sin_pnu_y)
            for j, hcor in enumerate(self.hcors):
                if hcor.s < bpm.s:
                    mu_x = (bpm.phi_x - hcor.phi_x)
                    h_resp[i, j] = kx*sqrt(hcor.beta_x)*sin(mu_x)*sqrt(hcor.E/bpm.E)

            for n, vcor in enumerate(self.vcors):
                if vcor.s < bpm.s:
                    mu_y = (bpm.phi_y - vcor.phi_y)
                    v_resp[i, n] = ky*sqrt(vcor.beta_y)*sin(mu_y)*sqrt(vcor.E/bpm.E)

        m = len(self.bpms)
        kx = len(self.hcors)
        ky = len(self.vcors)

        self.resp = zeros((2*m, kx + ky))
        self.resp[:m,:kx] = h_resp[:, :]
        self.resp[m:,kx:] = v_resp[:, :]

        rmatrix = Response_matrix()
        rmatrix.bpm_names = [b.id for b in self.bpms]
        rmatrix.cor_names = np.append(np.array([c.id for c in self.hcors]), np.array([c.id for c in self.vcors]))
        rmatrix.matrix = self.resp
        rmatrix.mode = "radian"
        return rmatrix

    def apply_svd(self, resp_matrix, misallign, weight=None):
        #print resp_matrix
        if weight is None:
            weight = eye(len(misallign))
        resp_matrix_w = dot(weight, resp_matrix)
        misallign_w = dot(weight, misallign)
        U, s, V = svd(resp_matrix_w)
        #print (s)
        s_inv = zeros(len(s))
        for i in range(len(s)):
            #if s[i]<1./max(s):
            if s[i] < 1.e-4:
                s_inv[i] = 0.
            else:
                s_inv[i] = 1./s[i]
        Sinv = zeros((shape(U)[0], shape(V)[0]))
        Sinv[:len(s), :len(s)] = diag(s_inv)
        Sinv = transpose(Sinv)
        A = dot(transpose(V), dot(Sinv, transpose(U)))
        angle = dot(A, misallign_w)
        return angle

    def correction(self, p_init=None):
        m = len(self.bpms)
        monitors = zeros(2*m)
        weights = eye(len(monitors))
        for i, bpm in enumerate(self.bpms):
            monitors[i] = bpm.x
            monitors[i+m] = bpm.y
            weights[i, i] = bpm.weight
            weights[i+m, i+m] = bpm.weight

        start = time()
        angle = self.apply_svd(self.resp, monitors, weight=weights)

        print("correction = ", time() - start)
        for i, cor in enumerate(np.append(self.hcors, self.vcors)):
            if self.mode == "ampere":
                #print "ampere"
                cor.dI = -angle[i]
            else:
                #print len(np.append(self.hcors, self.vcors)), i, len(angle)
                cor.angle -= angle[i]
                #print(cor.angle)
        """
        for i, vcor in enumerate(self.vcors):
            vcor.angle -= angle[i+len(self.hcors)]
        """
        """
        ix = 0
        iy = 0
        for elem in lattice.sequence:
            if ix<len(self.hcors) and elem.id == self.hcors[ix].id:
                elem.angle -= angle[ix]
                #print "x:", elem.angle
                ix += 1

            if iy<len(self.vcors) and elem.id == self.vcors[iy].id:
                elem.angle -= angle[iy+len(self.hcors)]
                #print "y:", elem.angle
                iy += 1
        """
        #print "ix = ", ix, "iy =", iy, len(angle)
        self.lat.update_transfer_maps()
        if p_init is not None:
            p_init.x = -angle[-4]
            p_init.px = -angle[-3]
            p_init.y  = -angle[-2]
            p_init.py = -angle[-1]
        return 0

    def elem_correction(self, elem_response, elem_types,  remove_elems=[]):
        m = len(self.bpms)
        monitors = zeros(2*m)
        self.create_types( elem_types, remove_elems=remove_elems)
        for i, bpm in enumerate(self.bpms):
            monitors[i] = bpm.x
            monitors[i+m] = bpm.y
        start = time()

        poss = self.apply_svd(elem_response, monitors)
        print("correction = ", time() - start)
        ix = 0
        iy = 0
        for elem in self.lat.sequence:
            if ix<len(self.htypes) and elem.id == self.htypes[ix].id:
                print ("quad, ", elem.dx, poss[ix])
                elem.dx += poss[ix]
                #self.hquads[ix].dx -= poss[ix]
                ix += 1

            if iy<len(self.vtypes) and elem.id == self.vtypes[iy].id:
                elem.dy += poss[iy+len(self.htypes)]
                #self.vquads[iy].dy -= poss[iy+len(self.hquads)]
                iy += 1
        p = Particle(x=poss[-4], px=poss[-3], y=poss[-2], py=poss[-1])
        #print poss[-5:]
        self.lat.update_transfer_maps()
        return p
    """
    def save_rmatrix(self, filename):
        dict_rmatrix = {}
        cors = np.append(self.hcors, self.vcors)
        cor_names = [cor.id for cor in cors]
        bpm_names = [bpm.id for bpm in self.bpms]
        dict_rmatrix["cor_names"] = cor_names
        dict_rmatrix["bpm_names"] = bpm_names
        dict_rmatrix["matrix"] = self.resp
        pickle.dump(dict_rmatrix, open(filename, "wb"))

    def read_rmatrix(self, filename):
        dict_rmatrix = pickle.load(open(filename, "rb"))
        cor_names = dict_rmatrix["cor_names"]
        bpm_names = dict_rmatrix["bpm_names"]
        rmatrix = dict_rmatrix["matrix"]
        orbit = Orbit(self.lat)
        orbit.create_COR(cor_list=cor_names)
        orbit.create_BPM(bpm_list=bpm_names)
        orbit.resp = rmatrix
        return orbit

    """

    def show_orbit(self, title=" "):

        s_bpm = np.array([p.s for p in self.bpms])
        x_bpm = np.array([p.x for p in self.bpms])
        y_bpm = np.array([p.y for p in self.bpms])

        ax = plot_API(self.lat)
        ax.set_title(title)
        ax.plot(s_bpm, x_bpm*1000.,  "ro-", label="X/mm")
        ax.plot(s_bpm, y_bpm*1000.,   "bo-", label="Y/mm")
        ax.legend()
        plt.draw()

    def calc_track(self,lattice):
        part_list = trace_obj(lattice, self.particle0, nPoints = None)
        self.x_track = map(lambda p: p.x, part_list)
        self.y_track = map(lambda p: p.y, part_list)
        self.s_track = map(lambda p: p.s, part_list)
"""
def draw_orbit(orbit, lattice, traject=True):
    if traject:
        bpm_draw = "o"
        part_list = trace_obj(lattice, orbit.particle0, nPoints = None)
        #plt.plot(map(lambda p: p.s, part_list), map(lambda p: p.y, part_list), "-")
    else:
        bpm_draw = "o-"
    #plt.plot(map(lambda p: p.s, orbit.bpms), map(lambda p: p.y, orbit.bpms), bpm_draw)
    #plt.grid(True)
    return part_list
"""
"""
def other_method(lat_err, p_list_1):
    navi = Navigator()
    p = p_list_1[0]
    p_list = [p]

    for i in range(len(p_list_1)-1):
        dz = p_list_1[i+1].s - p_list_1[i].s
        p = single_track(lat_err, p, dz, navi)
        p_list.append(p)
    return p_list


def real_response_matirx():
    pass
"""


def change_corrector(corrector, lattice):
    for elem in lattice.sequence:
        if elem.id == corrector.id:
            elem.angle += corrector.dI
            #print "change ", elem.angle
            elem.transfer_map = create_transfer_map(elem)
            #print elem.transfer_map.b(1)
    return lattice#.update_transfer_maps()

def restore_corrector(corrector, lattice):
    for elem in lattice.sequence:
        if elem.id == corrector.id:
            elem.angle -= corrector.dI
            elem.transfer_map = create_transfer_map(elem)
    return lattice#.update_transfer_maps()

def change_quad_position(quad, lattice, dx=0., dy=0.):
    for elem in lattice.sequence:
        if elem.id == quad.id:
            elem.dx += dx
            elem.dy += dy
            elem.transfer_map = create_transfer_map(elem)
    return lattice.update_transfer_maps()


def measure_response_matrix(orbit, lattice):

    m = len(orbit.bpms)
    real_resp = zeros((m*2, len(orbit.hcors)+len(orbit.vcors)))
    orbit.read_virtual_orbit( lattice)
    bpms = copy.deepcopy(orbit.bpms)
    for ix, hcor in enumerate(orbit.hcors):
        print("measure X - ", ix,"/",len(orbit.hcors))
        lattice = change_corrector(hcor, lattice)

        orbit.read_virtual_orbit(lattice)

        for j, bpm in enumerate(orbit.bpms):

            real_resp[j, ix] = (bpm.x - bpms[j].x)/hcor.dI
            real_resp[j+m, ix] = (bpm.y - bpms[j].y)/hcor.dI
        lattice = restore_corrector(hcor, lattice)

    for iy, vcor in enumerate(orbit.vcors):

        lattice = change_corrector(vcor, lattice)

        orbit.read_virtual_orbit(lattice)

        for j, bpm in enumerate(orbit.bpms):
            real_resp[j, iy+len(orbit.hcors)] = (bpm.x - bpms[j].x)/vcor.dI
            real_resp[j+m, iy+len(orbit.hcors)] = (bpm.y - bpms[j].y)/vcor.dI
        lattice = restore_corrector(vcor, lattice)
    return real_resp

def quad_response_matrix(orbit, lattice):

    m = len(orbit.bpms)
    nx = len(orbit.hquads)
    ny = len(orbit.vquads)
    print(nx, ny, m)
    real_resp = zeros((m*2, nx + ny))
    orbit.read_virtual_orbit(lattice)
    bpms = copy.deepcopy(orbit.bpms)
    for ix, hquad in enumerate(orbit.hquads):
        print("measure X - ", ix,"/",nx)
        lattice = change_quad_position(hquad, lattice, dx = 0.001, dy = 0)
        orbit.read_virtual_orbit(lattice)

        for j, bpm in enumerate(orbit.bpms):
            real_resp[j, ix] = (bpm.x - bpms[j].x)/0.001
            real_resp[j+m, ix] = (bpm.y - bpms[j].y)/0.001

            #if real_resp[j, ix] == 0 or real_resp[j+m, ix] == 0:

                #print bpm.x ,bpm.y, j, j+m, ix
        lattice = change_quad_position(hquad, lattice, dx = -0.001, dy = 0)

    for iy, vquad in enumerate(orbit.vquads):
        lattice = change_quad_position(vquad, lattice, dx = 0., dy = 0.001)
        orbit.read_virtual_orbit(lattice)

        for j, bpm in enumerate(orbit.bpms):
            real_resp[j, iy+nx] = (bpm.x - bpms[j].x)/0.001
            real_resp[j+m, iy+nx] = (bpm.y - bpms[j].y)/0.001

            #if real_resp[j, iy+nx] == 0 or real_resp[j+m, iy+nx] == 0:
            #print bpm.x ,bpm.y, j, j+m, iy+nx
        lattice = change_quad_position(vquad, lattice, dx = 0., dy = -0.001)
    return real_resp

def elem_response_matrix(orbit, lattice, p_init, elem_types, remove_elem):
    shift = 0.001
    m = len(orbit.bpms)
    orbit.create_types(lattice, elem_types, remove_elem)
    nx = len(orbit.htypes)
    ny = len(orbit.vtypes)
    print(nx, ny, m)
    real_resp = zeros((m*2, nx + ny +4))
    orbit.read_virtual_orbit(lattice, p_init=copy.deepcopy(p_init))
    bpms = copy.deepcopy(orbit.bpms)
    for ix, hquad in enumerate(orbit.htypes):
        print("measure X - ", ix, "/", nx)
        hquad.dx += shift
        lattice.update_transfer_maps()
        orbit.read_virtual_orbit(lattice, p_init=copy.deepcopy(p_init))

        for j, bpm in enumerate(orbit.bpms):
            real_resp[j, ix] = (bpm.x - bpms[j].x)/shift
            real_resp[j+m, ix] = (bpm.y - bpms[j].y)/shift

        hquad.dx -= shift
        lattice.update_transfer_maps()

    for iy, vquad in enumerate(orbit.vtypes):
        print("measure Y - ", iy,"/",ny)
        vquad.dy += shift
        lattice.update_transfer_maps()
        orbit.read_virtual_orbit(lattice, p_init=copy.deepcopy(p_init))
        #plt.plot([bpm.s for bpm in orbit.bpms], [bpm.x for bpm in orbit.bpms], "r")
        #plt.plot([bpm.s for bpm in orbit.bpms], [bpm.y for bpm in orbit.bpms], "b")
        #plt.show()
        for j, bpm in enumerate(orbit.bpms):
            real_resp[j, iy+nx] = (bpm.x - bpms[j].x)/shift
            real_resp[j+m, iy+nx] = (bpm.y - bpms[j].y)/shift
        vquad.dy -= shift
        lattice.update_transfer_maps()

    for i, par in enumerate(["x", "px", "y", "py"]):
        print(i)
        p_i = Particle(E = p_init.E)
        p_i.__dict__[par] = 0.0001
        #print p_i.x, p_i.px, p_i.y, p_i.py, p_i.E
        p2 = copy.deepcopy(p_i)
        orbit.read_virtual_orbit(lattice, p_init=p2)
        #print ("energy = ", p2.E)
        #plt.plot([bpm.s for bpm in orbit.bpms], [bpm.x for bpm in orbit.bpms], "r")
        #plt.plot([bpm.s for bpm in orbit.bpms], [bpm.y for bpm in orbit.bpms], "b")
        #plt.show()
        for j, bpm in enumerate(orbit.bpms):
            real_resp[j, nx + ny + i] = (bpm.x - bpms[j].x)/0.0001
            real_resp[j+m, nx + ny + i] = (bpm.y - bpms[j].y)/0.0001
            #print j+m, nx + ny + i, (bpm.x - bpms[j].x)/0.00001
    #print real_resp[:,-5:]
    return real_resp


def test(lattice, errors):
    lat_errors = errors_seed(lattice, errors, nsuperperiods = 6)

    orbit = Orbit()
    orbit.lattice_analysis(lat_errors)
    orbit.ring_response_matrix(lat_errors)
    #print orbit.resp
    #real_resp = measure_response_matrix(orbit, lat_errors)
    #print real_resp
    #orbit.resp = real_resp
    orbit.read_virtual_orbit(lat_errors)
    draw(orbit, lat_errors, traject = True)

    lat_corrected = orbit.correction( lat_errors)
    orbit.read_virtual_orbit(lat_corrected)
    draw(orbit,lat_corrected, traject = True)


    lat_corrected2 = orbit.correction(lat_errors)
    for hcor in orbit.hcors:
        print (hcor.id, hcor.angle)

    orbit.read_virtual_orbit(lat_corrected2)
    draw(orbit, lat_corrected2, traject = True)
    plt.show()




if __name__ == "__main__":
    from ocelot.cpbd.elements import *
    exec( open("../../repository/siberia2/sibir2_correct.inp" ))
    err_list = {"quadrupole": {"offset": 0.02e-3, "dtilt": 0.001},
            "sbend":{"dtilt": 0.001}, "rbend":{"dtilt": 0.001}, "bend":{"dtilt": 0.001}}
    lat = MagneticLattice(superperiod)
    #print lat.totalLen
    #tw0 = Twiss(beam)
    #tws = twiss(lat, tw0, nPoints=1000)
    #plt.plot(map(lambda p: p.s, tws), map(lambda tws: tws.beta_x, tws), 'b-')
    #plt.plot(map(lambda p: p.s, tws), map(lambda tws: tws.beta_y, tws), 'b-')
    #plt.grid()
    #plt.show()

    test(lat, err_list)