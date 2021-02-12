from ocelot.cpbd.optics import *
from ocelot.cpbd.match import closed_orbit
from ocelot.cpbd.track import tracking_step
import copy
import os
import numpy as np
from scipy.interpolate import splrep, splev
import json
import time
from threading import Thread
import pandas as pd

import logging
logger = logging.getLogger(__name__)


class MeasureResponseMatrix:
    def __init__(self, lattice, hcors, vcors, bpms):
        self.lat = lattice
        self.hcors = hcors
        self.vcors = vcors
        self.bpms = bpms

    def calculate(self):
        pass

    def read_virtual_orbit(self, p_init=None, write2bpms=True):
        """
        searching closed orbit by function closed_orbit(lattice) and searching coordinates of beam at the bpm positions

        :param lattice: class MagneticLattice
        :return: orbit.bpms
        """
        particles = []
        X = []
        Y = []
        if p_init == None:
            self.particle0 = closed_orbit(self.lat)
        else:
            self.particle0 = p_init
        p = copy.copy(self.particle0)
        navi = Navigator(self.lat)
        L = 0.
        for bpm in self.bpms:
            #print("energy = ", p.E)
            dz = bpm.s - L
            tracking_step(self.lat, [p], dz, navi)
            if write2bpms:
                bpm.x = p.x
                bpm.y = p.y
                bpm.E = p.E
                bpm.p = p.p
            L = bpm.s
            X.append(p.x)
            Y.append(p.y)

            particles.append(copy.copy(p))
        #print("energy = ", p.E)
        return np.array(X), np.array(Y),

    def read_virtual_dispersion(self, E0):
        X0, Y0 = self.read_virtual_orbit(p_init=Particle(p=0.000, E=E0))
        X1, Y1 = self.read_virtual_orbit(p_init=Particle(p=0.01, E=E0), write2bpms=False)
        #p = np.array([bpm.p for bpm in self.bpms])
        Dx0 = (X1 - X0) / 0.01
        Dy0 = (Y1 - Y0) / 0.01
        for i, bpm in enumerate(self.bpms):
            bpm.Dx = Dx0[i]
            bpm.Dy = Dy0[i]
        return Dx0, Dy0

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

        tws = twiss(self.lat, tw_init, nPoints=int(self.lat.totalLen / 0.05))
        s = np.array([tw.s for tw in tws])
        tck_mux = splrep(s, np.array([tw.mux for tw in tws]))
        tck_muy = splrep(s, np.array([tw.muy for tw in tws]))

        beta_x = np.array([tw.beta_x for tw in tws])
        beta_y = np.array([tw.beta_y for tw in tws])
        self.nu_x = tws[-1].mux / 2. / np.pi
        self.nu_y = tws[-1].muy / 2. / np.pi

        tck_bx = splrep(s, beta_x)
        tck_by = splrep(s, beta_y)
        energy_s = [tw.E for tw in tws]
        tck_E = splrep(s, energy_s)
        # s_bpm = [bpm.s for bpm in self.bpms]

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


class RingRM(MeasureResponseMatrix):

    def __init__(self, lattice, hcors, vcors, bpms):
        super(RingRM, self).__init__(lattice, hcors, vcors, bpms)

    def calculate(self, tw_init=None):
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
        h_resp = np.zeros((m, nx))
        v_resp = np.zeros((m, ny))
        sin_pnu_x = np.sin(np.pi * self.nu_x)
        sin_pnu_y = np.sin(np.pi * self.nu_y)
        for i, bpm in enumerate(self.bpms):
            kx = np.sqrt(bpm.beta_x) / (2. * sin_pnu_x)
            ky = np.sqrt(bpm.beta_y) / (2. * sin_pnu_y)
            for j, hcor in enumerate(self.hcors):
                mu_x = abs(bpm.phi_x - hcor.phi_x)
                h_resp[i, j] = kx * np.sqrt(hcor.beta_x) * np.cos(mu_x - np.pi * self.nu_x)
            for n, vcor in enumerate(self.vcors):
                mu_y = abs(bpm.phi_y - vcor.phi_y)
                v_resp[i, n] = ky * np.sqrt(vcor.beta_y) * np.cos(mu_y - np.pi * self.nu_y)
        m = len(self.bpms)
        kx = len(self.hcors)
        ky = len(self.vcors)
        self.resp = np.zeros((2 * m, kx + ky))
        self.resp[:m, :kx] = h_resp[:, :]
        self.resp[m:, kx:] = v_resp[:, :]
        # print "shape = ", shape(self.resp)
        return self.resp


class LinacOpticalRM(MeasureResponseMatrix):

    def __init__(self, lattice, hcors, vcors, bpms):
        super(LinacOpticalRM, self).__init__(lattice, hcors, vcors, bpms)

    def calculate(self, tw_init=None):
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
        h_resp = np.zeros((m, nx))
        v_resp = np.zeros((m, ny))

        for i, bpm in enumerate(self.bpms):
            kx = np.sqrt(bpm.beta_x)#/(2.*sin_pnu_x)
            ky = np.sqrt(bpm.beta_y)#/(2.*sin_pnu_y)
            for j, hcor in enumerate(self.hcors):
                if hcor.s < bpm.s:
                    mu_x = (bpm.phi_x - hcor.phi_x)
                    h_resp[i, j] = kx*np.sqrt(hcor.beta_x)*np.sin(mu_x)*np.sqrt(hcor.E/bpm.E)


            for n, vcor in enumerate(self.vcors):
                if vcor.s < bpm.s:
                    mu_y = (bpm.phi_y - vcor.phi_y)
                    v_resp[i, n] = ky*np.sqrt(vcor.beta_y)*np.sin(mu_y)*np.sqrt(vcor.E/bpm.E)

        m = len(self.bpms)
        kx = len(self.hcors)
        ky = len(self.vcors)

        self.resp = np.zeros((2*m, kx + ky))
        self.resp[:m,:kx] = h_resp[:, :]
        self.resp[m:,kx:] = v_resp[:, :]

        #rmatrix = ResponseMatrix()
        #rmatrix.bpm_names = [b.id for b in self.bpms]
        #rmatrix.cor_names = np.append(np.array([c.id for c in self.hcors]), np.array([c.id for c in self.vcors]))
        #rmatrix.matrix = self.resp
        #rmatrix.mode = "radian"
        return self.resp


class LinacSimRM(MeasureResponseMatrix):

    def __init__(self, lattice, hcors, vcors, bpms):
        super(LinacSimRM, self).__init__(lattice, hcors, vcors, bpms)

    def calculate(self, tw_init=None):
        """
        calculation of ideal response matrix

        :param lattice: class MagneticLattice
        :param tw_init: if tw_init == None, function tries to find periodical solution
        :return: orbit.resp
        """
        #self.optical_func_params(tw_init=tw_init)
        match_ic = False  # for future, fitting the initial conditions
        m = len(self.bpms)
        nx = len(self.hcors)
        ny = len(self.vcors)

        add_i = 0
        if match_ic:
            add_i = 4
        self.resp = np.zeros((2 * m, nx + ny + add_i))
        s = [bpm.s for bpm in self.bpms]
        X0, Y0 = self.read_virtual_orbit(p_init=Particle( E=tw_init.E))

        Or0 = np.append(X0, Y0)
        for j, cor in enumerate([item for sublist in [self.hcors, self.vcors] for item in sublist]):
            print(j,"/", nx+ny, cor.id)
            cor.angle = 0.0001
            self.lat.update_transfer_maps()
            #start = time()
            X1, Y1 = self.read_virtual_orbit(p_init=Particle(E=tw_init.E))
            #print(time() - start)

            Or1 = np.append(X1, Y1)
            self.resp[:, j] = (Or1 - Or0)/cor.angle
            cor.angle = 0.00
        self.lat.update_transfer_maps()
        if match_ic:
            for i, par in enumerate(["x", "px", "y", "py"]):
                #print(i)
                p_i = Particle(E = tw_init.E)
                p_i.__dict__[par] = 0.0001
                #print ("measure = ", p_i.x, p_i.y, p_i.px, p_i.py, p_i.E)
                X1, Y1 = self.read_virtual_orbit(p_init=p_i)
                Or1 = np.append(X1, Y1)
                self.resp[:, nx + ny + i] = (Or1 - Or0) /0.0001
        X1, Y1 = self.read_virtual_orbit(p_init=Particle(E=tw_init.E))
        #rmatrix = ResponseMatrix()
        #rmatrix.bpm_names = [b.id for b in self.bpms]
        #rmatrix.cor_names = np.append(np.array([c.id for c in self.hcors]), np.array([c.id for c in self.vcors]))
        #rmatrix.matrix = self.resp
        #rmatrix.mode = "radian"
        return self.resp


class LinacRmatrixRM(MeasureResponseMatrix):

    def __init__(self, lattice, hcors, vcors, bpms):
        super(LinacRmatrixRM, self).__init__(lattice, hcors, vcors, bpms)

    def calculate(self, tw_init=None):
        """
        calculation of ideal response matrix

        :param lattice: class MagneticLattice
        :param tw_init: if tw_init == None, initial beam energy is ZERO
        :return: orbit.resp
        """
        if tw_init is None:
            logger.warning("tw_init is None. Initial beam energy is assuemed ZERO")
            Einit = 0
        else:
            Einit = tw_init.E

        m = len(self.bpms)
        nx = len(self.hcors)
        ny = len(self.vcors)
        self.resp = np.zeros((2 * m, nx + ny))

        for j, cor in enumerate([item for sublist in [self.hcors, self.vcors] for item in sublist]):
            print(j, "/", nx + ny, cor.id)
            Ra = np.eye(6)
            E = Einit
            for i, elem in enumerate(self.lat.sequence):
                if i < cor.lat_inx:
                    E += elem.transfer_map.delta_e
                    continue

                Rb = elem.transfer_map.R(E)
                Ra = np.dot(Rb, Ra)
                E += elem.transfer_map.delta_e
                if elem in self.bpms:

                    n = self.bpms.index(elem)

                    if cor.__class__ == Hcor:
                        self.resp[n, j] = Ra[0, 1]

                    else:
                        self.resp[n + m, j] = Ra[2, 3]
        return self.resp


class LinacDisperseSimRM(MeasureResponseMatrix):

    def __init__(self, lattice, hcors, vcors, bpms):
        super(LinacDisperseSimRM, self).__init__(lattice, hcors, vcors, bpms)

    def calculate(self, tw_init=None):
        """
        calculation of ideal dispersive response matrix

        :param lattice: class MagneticLattice
        :param tw_init: if tw_init == None, function tries to find periodical solution
        :return: orbit.resp
        """
        if tw_init == None:
            print("ADD INITIAL TWISS TO LinacDisperseSimRM.calculate(tw_init)")

        m = len(self.bpms)
        nx = len(self.hcors)
        ny = len(self.vcors)

        self.resp = np.zeros((2 * m, nx + ny))
        s = [bpm.s for bpm in self.bpms]
        Dx0, Dy0 = self.read_virtual_dispersion(E0=tw_init.E)
        D0 = np.append(Dx0, Dy0)
        for j, cor in enumerate([item for sublist in [self.hcors, self.vcors] for item in sublist]):
            print(j, "/", nx + ny, cor.id)
            cor.angle = 0.0005
            self.lat.update_transfer_maps()
            cor.transfer_map = self.lat.method.create_tm(cor)
            start = time.time()
            Dx1, Dy1 = self.read_virtual_dispersion(E0=tw_init.E)
            #if np.max(Dx1)>1e+100 or np.max(Dy1) > 1e+100:
            print(time.time() - start)

            D1 = np.append(Dx1, Dy1)
            self.resp[:, j] = (D1 - D0) / cor.angle
            cor.angle = 0.00
            cor.transfer_map = self.lat.method.create_tm(cor)
        #self.lat.update_transfer_maps()
        return self.resp


class LinacDisperseTmatrixRM(MeasureResponseMatrix):

    def __init__(self, lattice, hcors, vcors, bpms):
        super(LinacDisperseTmatrixRM, self).__init__(lattice, hcors, vcors, bpms)

    def calculate(self, tw_init=None):
        """
        calculation of ideal response matrix

        :param lattice: class MagneticLattice
        :param tw_init: if tw_init == None, function tries to find periodical solution
        :return: orbit.resp
        """
        if tw_init == None:
            print("ADD INITIAL TWISS TO LinacDisperseTmatrixRM.calculate(tw_init)")
        m = len(self.bpms)
        nx = len(self.hcors)
        ny = len(self.vcors)
        self.resp = np.zeros((2 * m, nx + ny))

        for j, cor in enumerate([item for sublist in [self.hcors, self.vcors] for item in sublist]):
            print(j, "/", nx + ny, cor.id)
            Ra = np.eye(6)
            Ta = np.zeros((6, 6, 6))
            E = tw_init.E
            for i, elem in enumerate(self.lat.sequence):
                if i < cor.lat_inx:
                    E += elem.transfer_map.delta_e
                    continue
                #Tc = np.zeros((6, 6, 6))
                Rb = elem.transfer_map.R(E)
                Tb = deepcopy(elem.transfer_map.t_mat_z_e(elem.l, E))
                #Ra = dot(Rb, Ra)
                Ra, Ta = transfer_maps_mult(Ra, Ta, Rb, Tb)
                E += elem.transfer_map.delta_e
                if elem in self.bpms:

                    n = self.bpms.index(elem)

                    if cor.__class__ == Hcor:
                        self.resp[n, j] = Ta[0, 1, 5]

                    else:
                        self.resp[n + m, j] = Ra[2, 3, 5]
        return self.resp


class ResponseMatrixJSON:
    def __init__(self, method=None):
        self.cor_names = []
        self.bpm_names = []
        self.matrix = []
        self.method = method
        self.mode = "radian"  # or "ampere"

        self.tw_init = None   # for self.run()
        self.filename = None  # for self.run()

    def calculate(self, tw_init=None):
        """
        rewrites cor_name, bpm_name and matrix

        :param method:
        :return:
        """
        if self.method != None:
            hcors = self.method.hcors
            vcors = self.method.vcors
            bpms = self.method.bpms
            self.cor_names = np.append([cor.id for cor in hcors], [cor.id for cor in vcors])
            self.bpm_names = [bpm.id for bpm in bpms]
            self.matrix = self.method.calculate(tw_init=tw_init)
        else:
            print("ResponseMatrix.method = None, Add the method, e.g. MeasureResponseMatrix")

    def get_matrix(self):
        return self.matrix

    def extract(self, cor_list, bpm_list):
        cor_list = np.array(cor_list)
        bpm_list = np.array(bpm_list)
        #print("EXTRACT given: ", cor_list)
        #print("EXTRACT given: ", bpm_list)
        cors1 = np.array(self.cor_names)
        cors2 = cor_list
        bpms1 = np.array(self.bpm_names)
        bpms2 = bpm_list
        nb1 = len(bpms1)
        #print("EXTRACT self: ", cors1)
        #print("EXTRACT self: ", bpms1)
        c_names = cors1[np.in1d(cors1, cors2)]

        c_i1 = np.where(np.in1d(cors1, cors2))[0]
        c_i2 = np.where(np.in1d(cors2, c_names))[0]
        b_names = bpms1[np.in1d(bpms1, bpms2)]
        b_i1 = np.where(np.in1d(bpms1, bpms2))[0]
        b_i2 = np.where(np.in1d(bpms2, b_names))[0]

        if not np.array_equal(c_names, cor_list):
            print (" Origin response matrix has no correctors:")
            print (cors2[np.in1d(cors2, c_names, invert=True)])
        if not np.array_equal(b_names, bpm_list):
            print (" Origin response matrix has no BPMs:")
            print (bpm_list[b_i2[:]])

        extr_matrix = np.zeros((len(b_names)*2, len(c_names)))
        for n in range(2):
            for i, c in enumerate(c_names):
                for j, b in enumerate(b_names):
                    x1 = self.matrix[b_i1[j] + nb1*n, c_i1[i]]
                    extr_matrix[j + n*len(b_names), i] = x1

        return extr_matrix

    def inject(self, cor_list, bpm_list, inj_matrix):
        """
        Update some elements of the response matrix

        :param cor_list:
        :param bpm_list:
        :param inj_matrix:
        :return:
        """
        cor_list = np.array(cor_list)
        bpm_list = np.array(bpm_list)
        #print("EXTRACT given: ", cor_list)
        #print("EXTRACT given: ", bpm_list)
        cors1 = np.array(self.cor_names)
        cors2 = cor_list
        bpms1 = np.array(self.bpm_names)
        bpms2 = bpm_list
        nb1 = len(bpms1)
        #print("EXTRACT self: ", cors1)
        #print("EXTRACT self: ", bpms1)
        c_names = cors1[np.in1d(cors1, cors2)]

        c_i1 = np.where(np.in1d(cors1, cors2))[0]
        c_i2 = np.where(np.in1d(cors2, c_names))[0]
        b_names = bpms1[np.in1d(bpms1, bpms2)]
        b_i1 = np.where(np.in1d(bpms1, bpms2))[0]
        b_i2 = np.where(np.in1d(bpms2, b_names))[0]

        if not np.array_equal(c_names, cor_list):
            logger.warning(" ResponseMatrix.inject: Origin response matrix has no correctors:")
            print (cors2[np.in1d(cors2, c_names, invert=True)])
        if not np.array_equal(b_names, bpm_list):
            logger.warning(" ResponseMatrix.inject: Origin response matrix has no BPMs:")
            print (bpm_list[b_i2[:]])

        #extr_matrix = np.zeros((len(b_names)*2, len(c_names)))

        for n in range(2):
            for i, c in enumerate(c_names):
                for j, b in enumerate(b_names):
                    #print("old matrix elem = ",b_i1[j] + nb1*n,  c_i1[i], self.matrix[b_i1[j] + nb1*n, c_i1[i]],
                    #      "inj_mat_elem =", j + n*len(b_names), i, inj_matrix[j + n*len(b_names), i])
                    self.matrix[b_i1[j] + nb1*n, c_i1[i]] = inj_matrix[j + n*len(b_names), i]
                    #extr_matrix[j + n*len(b_names), i] = x1

        return self.matrix

    def dump(self, filename):
        dict_rmatrix = {}
        dict_rmatrix["cor_names"] = list(self.cor_names)
        dict_rmatrix["bpm_names"] = list(self.bpm_names)
        dict_rmatrix["matrix"] = list(self.matrix.flatten())

        dict_rmatrix["method_name"] = self.method.__class__.__name__ if self.method != None else "None"
        dict_rmatrix["mode"] = self.mode

        directory = os.path.dirname(filename)
        if not os.path.exists(directory):
            os.makedirs(directory)

        with open(filename, 'w+') as f:
            json.dump(dict_rmatrix, f)

    def load(self, filename):
        with open(filename, 'r') as f:
            dict_rmatrix = json.load(f)
        self.cor_names = dict_rmatrix["cor_names"]
        self.bpm_names = dict_rmatrix["bpm_names"]
        self.method_name = dict_rmatrix["method_name"]
        r_matrix = np.array(dict_rmatrix["matrix"])
        self.matrix = r_matrix.reshape(2*len(self.bpm_names),len(self.cor_names))
        self.mode = dict_rmatrix["mode"]
        return 1

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
        for i in range(np.shape(self.matrix)[1]):
            print (self.cor_names[i] + " "*(10 - len(self.cor_names[i])),)
            #print np.array_str(self.matrix[:, i], precision=2, suppress_small=True)
            for j in range(np.shape(self.matrix)[0]):
                print ("%.2f" % self.matrix[j, i],)
            print()


class ResponseMatrix:
    def __init__(self, method=None):
        self.cor_names = []
        self.bpm_names = []
        self.matrix = []
        self.method = method
        self.mode = "radian"  # or "ampere"
        self.df = None
        self.tw_init = None   # for self.run()
        self.filename = None  # for self.run()

    def bpm2x_name(self, bpm_id):
        """
        Transform bpm id to a name how it use in a control system to get horizontal beam position
        :param bpm_id:
        :return: channel for X position
        """
        return bpm_id + ".X"

    def bpm2y_name(self, bpm_id):
        """
        Transform bpm id to a name how it use in a control system to get vertical beam position

        :param bpm_id:
        :return: channel for Y position
        """
        return bpm_id + ".Y"

    def xy_names2bpm_id(self, xy_names):
        """
        transform BPM channels to bpm ids

        :param xy_names:
        :return:
        """
        bpm_ids = [bpm.replace(".X", "") for bpm in xy_names if ".X" in bpm]
        return bpm_ids

    def calculate(self, tw_init=None):
        """
        rewrites cor_name, bpm_name and matrix

        :param method:
        :return:
        """
        if self.method != None:
            hcors = self.method.hcors
            vcors = self.method.vcors
            bpms = self.method.bpms
            self.cor_names = np.append([cor.id for cor in hcors], [cor.id for cor in vcors])
            self.bpm_names = [bpm.id for bpm in bpms]
            self.matrix = self.method.calculate(tw_init=tw_init)
            self.df = self.data2df(matrix=self.matrix, bpm_names=self.bpm_names , cor_names=self.cor_names)
        else:
            print("ResponseMatrix.method = None, Add the method, e.g. MeasureResponseMatrix")

    def get_matrix(self):
        return self.matrix

    def extract_df_slice(self, cor_list, bpm_list):
        bpm_x = [self.bpm2x_name(bpm) for bpm in bpm_list]
        bpm_y = [self.bpm2y_name(bpm) for bpm in bpm_list]
        rows = bpm_x + bpm_y
        cols = list(cor_list)
        cor_list_exist = np.array([item in self.df.columns for item in cor_list])
        bpm_list_exist = np.array([item in self.df.index for item in rows])

        if all(cor_list_exist) and all(bpm_list_exist):
            df_slice = self.df.loc[rows, cols]
            return df_slice
        else:
            print("correctors are not in the RM")
            print(np.array(cor_list)[~cor_list_exist])
            print()
            print("BPMs are not in the RM")
            print(np.array(rows)[~bpm_list_exist])
            return None

    def extract(self, cor_list, bpm_list):
        df_slice = self.extract_df_slice(cor_list, bpm_list)
        return df_slice.values

    def retrieve_from_scan(self, df_scan):
        from sklearn.linear_model import LinearRegression

        bpm_x = [self.bpm2x_name(bpm) for bpm in self.bpm_names]
        bpm_y = [self.bpm2y_name(bpm) for bpm in self.bpm_names]
        bpm_names_xy = bpm_x + bpm_y

        x = df_scan.loc[:, self.cor_names].values
        y = df_scan.loc[:, bpm_names_xy].values

        reg = LinearRegression().fit(x, y)
        #x_test = np.eye(np.shape(x)[1])
        rm = reg.coef_
        #df_rm = pd.DataFrame(rm.T, columns=self.cor_names, index=bpm_x+bpm_y)
        self.df = self.data2df(matrix=rm, bpm_names=self.bpm_names, cor_names=self.cor_names)
        return self.df

    def clean_rm(self, coupling=True):
        if self.method != None:
            hcors = self.method.hcors
            vcors = self.method.vcors
            bpms = self.method.bpms
            for hcor in hcors:
                for bpm in bpms:
                    if bpm.s < hcor.s:

                        self.df.loc[self.bpm2x_name(bpm.id), hcor.id] = 0
                    if not coupling:
                        self.df.loc[self.bpm2y_name(bpm.id), hcor.id] = 0
            for vcor in vcors:
                for bpm in bpms:
                    if bpm.s < vcor.s:
                        self.df.loc[self.bpm2y_name(bpm.id), vcor.id] = 0
                    if not coupling:
                        self.df.loc[self.bpm2x_name(bpm.id), vcor.id] = 0

        else:
            print("ResponseMatrix.method = None, Add the method, e.g. MeasureResponseMatrix")

    def inject(self, cor_list, bpm_list, inj_matrix):
        """
        Update some elements of the response matrix

        :param cor_list:
        :param bpm_list:
        :param inj_matrix:
        :return:
        """
        df_slice = self.data2df(matrix=inj_matrix, bpm_names=bpm_list, cor_names=cor_list)
        self.df.update(df_slice)
        self.matrix = self.df.values
        return self.matrix

    def data2df(self, matrix, bpm_names, cor_names):
        bpm_x = [self.bpm2x_name(bpm) for bpm in bpm_names]
        bpm_y = [self.bpm2y_name(bpm) for bpm in bpm_names]
        df = pd.DataFrame(matrix, columns=cor_names, index=bpm_x + bpm_y)
        return df

    def df2data(self):
        self.cor_names = list(self.df.columns.values)
        bpms_all = list(self.df.index.values)
        self.bpm_names = self.xy_names2bpm_id(bpms_all)
        self.matrix = self.df.values

    def dump(self, filename):
        df = self.data2df(matrix=self.matrix, bpm_names=self.bpm_names, cor_names=self.cor_names)
        directory = os.path.dirname(filename)
        if directory != "" and not os.path.exists(directory):
            os.makedirs(directory)
        df.to_pickle(filename)

    def load(self, filename):
        self.df = pd.read_pickle(filename)
        self.df2data()
        return 1

    def compare(self, rmatrix, absolut=0.001, relative=0.1):
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
                    if abs(x1 - x2) < absolut:
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
        for i in range(np.shape(self.matrix)[1]):
            print (self.cor_names[i] + " "*(10 - len(self.cor_names[i])),)
            #print np.array_str(self.matrix[:, i], precision=2, suppress_small=True)
            for j in range(np.shape(self.matrix)[0]):
                print ("%.2f" % self.matrix[j, i],)
            print()