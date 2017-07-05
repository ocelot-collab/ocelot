__author__ = 'Sergey Tomin'


#import pickle
from time import sleep

#import matplotlib.pyplot as plt
#from ocelot.gui.accelerator import *
from numpy import diag, shape
from numpy.linalg import svd
from scipy.interpolate import splrep, splev
from scipy.optimize import linprog
from ocelot.cpbd.match import closed_orbit
from ocelot.cpbd.track import *

from ocelot.cpbd.response_matrix import *
import copy
import json



class OrbitSVD:
    def __init__(self, resp_matrix, orbit, weights=None, epsilon_x=0.001, epsilon_y=0.001):
        self.resp_matrix = resp_matrix
        self.orbit = orbit
        self.weights = weights
        self.epsilon_x = epsilon_x
        self.epsilon_y = epsilon_y

    def apply(self):
        # print resp_matrix
        if self.weights is None:
            self.weights = eye(len(self.orbit))
        #print(np.shape(self.weights), np.shape(self.resp_matrix))
        resp_matrix_w = dot(self.weights, self.resp_matrix)
        misallign_w = dot(self.weights, self.orbit)
        #U, s, V = svd(resp_matrix_w)

        U, s, V = svd(self.resp_matrix)
        # print (s)
        s_inv = zeros(len(s))
        s_max = max(s)
        for i in range(len(s)):
            #print("S[",i,"]=", s[i], "s max = ", s_max)
            if i < int(len(s)/2.):
                epsilon = self.epsilon_x
            else:
                epsilon = self.epsilon_y
            if s[i] < s_max * epsilon:
                s_inv[i] = 0.
            else:
                s_inv[i] = 1. / s[i]
        #print(s_inv)
        Sinv = zeros((shape(U)[0], shape(V)[0]))
        Sinv[:len(s), :len(s)] = diag(s_inv)
        Sinv = transpose(Sinv)
        A = dot(transpose(V), dot(Sinv, transpose(U)))
        #angle = dot(A, misallign_w)
        angle = dot(A, self.orbit)
        #print(A)
        #print(angle)
        return angle

class LInfinityNorm(OrbitSVD):
    def __init__(self, resp_matrix, orbit, weights=None, epsilon_x=0.001, epsilon_y=0.001):
        OrbitSVD.__init__(self, resp_matrix=resp_matrix, orbit=orbit, weights=weights,
                          epsilon_x=epsilon_x, epsilon_y=epsilon_y)

    def apply(self):
        m, n = np.shape(self.resp_matrix)
        f = np.zeros(n + 1)
        f[-1] = 1
        Ane = np.vstack((np.hstack((self.resp_matrix, -np.ones((m, 1)))), np.hstack((-self.resp_matrix, -np.ones((m, 1))))))
        bne = np.vstack((+self.orbit, -self.orbit))
        res = linprog(f, A_ub=Ane, b_ub=bne)
        x = res["x"][:-1]
        return x


class NewOrbit:
    def __init__(self, lattice, rm_method=None, disp_rm_method=None, empty=False):
        self.lat = lattice
        self.bpms = []
        self.hcors = []
        self.vcors = []
        self.nu_x = 0.
        self.nu_y = 0.
        self.rm_method = rm_method
        self.disp_rm_method = disp_rm_method
        self.response_matrix = None
        self.disp_response_matrix = None
        self.mode = "radian" # or "ampere"

        if not empty:
            self.create_bpms()
            self.create_correctors()

        if rm_method != None and (not empty):
            self.setup_response_matrix()

        if disp_rm_method != None and (not empty):
            self.setup_disp_response_matrix()

    def setup_response_matrix(self):
        method = self.rm_method(lattice=self.lat, hcors=self.hcors, vcors=self.vcors, bpms=self.bpms)
        self.response_matrix = ResponseMatrix(method=method)

    def setup_disp_response_matrix(self):
        method = self.disp_rm_method(lattice=self.lat, hcors=self.hcors, vcors=self.vcors, bpms=self.bpms)
        self.disp_response_matrix = ResponseMatrix(method=method)

    #def update_devices_in_RMs(self):
    #
    #    if self.response_matrix != None:
    #        self.response_matrix.hcors = self.hcors
    #        self.response_matrix.vcors = self.vcors
    #        self.response_matrix.bpms = self.bpms
    #
    #    if self.disp_response_matrix != None:
    #        self.disp_response_matrix.hcors = self.hcors
    #        self.disp_response_matrix.vcors = self.vcors
    #        self.disp_response_matrix.bpms = self.bpms

    def create_bpms(self, bpm_list=None):
        """
        Search bpm in the lattice and create list of bpms
        :param lattice: class MagneticLattice
        :return: self.bpms - list of BPMs (class BPM)
        """
        self.bpms = []
        L = 0.
        for i, elem in enumerate(self.lat.sequence):
            if elem.__class__ == Monitor:
                if bpm_list is None or elem.id in bpm_list:
                    try:
                        elem.weight
                    except:
                        elem.weight = 1.
                    elem.s = L + elem.l / 2.
                    elem.x_ref = 0.
                    elem.y_ref = 0.
                    elem.Dx = 0.
                    elem.Dy = 0.
                    elem.Dx_des = 0.
                    elem.Dy_des = 0.
                    elem.lat_inx = i
                    self.bpms.append(elem)
            L += elem.l
        if len(self.bpms) == 0:
            print("there are not monitors")
        return self.bpms

    def create_correctors(self, cor_list=None):
        """
        Search correctors (horizontal and vertical) in the lattice and create list of hcors and list of vcors
        :param lattice: class MagneticLattice
        :return:
        """
        self.hcors = []
        self.vcors = []
        L = 0.
        for i, elem in enumerate(self.lat.sequence):
            if elem.__class__ == Vcor:
                if cor_list is None or elem.id in cor_list:
                    elem.s = L+elem.l/2.
                    elem.lat_inx = i
                    self.vcors.append(elem)
            elif elem.__class__ == Hcor:
                if cor_list is None or elem.id in cor_list:
                    elem.s = L+elem.l/2.
                    elem.lat_inx = i
                    self.hcors.append(elem)
            L += elem.l
        if len(self.hcors) == 0:
            print("there are not horizontal correctors")
        if len(self.vcors) == 0:
            print("there are not vertical correctors")

    def get_ref_orbit(self):
        for bpm in self.bpms:
            bpm.x_ref = 0.
            bpm.y_ref = 0.

    def get_orbit(self):

        #self.get_ref_orbit()

        m = len(self.bpms)
        orbit = zeros(2 * m)
        for i, bpm in enumerate(self.bpms):
            #print("get_orbit = ",bpm.id, bpm.x,  bpm.x_ref)
            orbit[i] = bpm.x - bpm.x_ref
            orbit[i+m] = bpm.y - bpm.y_ref
        return orbit


    def get_dispersion(self):
        m = len(self.bpms)
        disp = zeros(2 * m)
        for i, bpm in enumerate(self.bpms):
            disp[i] = bpm.Dx - bpm.Dx_des
            disp[i + m] = bpm.Dy - bpm.Dy_des
        return disp

    def combine_matrices(self, mat1, mat2):
        """

        :param mat1:
        :param mat2:
        :return:
        """
        n1, m1 = np.shape(mat1)
        n2, m2 = np.shape(mat2)
        rm = np.zeros((n1+n2, m1+m2))
        rm[:n1, :m1] = mat1[:, :]
        rm[n1:, m1:] = mat2[:, :]
        return rm

    def correction(self, alpha=0,  epsilon_x=0.001, epsilon_y=0.001, p_init=None, print_log=True):
        #TODO: initial condition for particle was removed. Add it again
        cor_list = [cor.id for cor in np.append(self.hcors, self.vcors)]
        bpm_list = [bpm.id for bpm in self.bpms]
        orbit = (1 - alpha) * self.get_orbit()

        RM = (1 - alpha) * self.response_matrix.extract(cor_list=cor_list, bpm_list=bpm_list)
        #print("RM = ", np.shape(RM))
        if alpha != 0:
            disp = alpha * self.get_dispersion()
        else:
            disp = alpha * np.array(orbit)

        if self.disp_response_matrix != None:
            DRM = alpha * self.disp_response_matrix.extract(cor_list=cor_list, bpm_list=bpm_list)
        else:
            DRM = np.zeros_like(RM)
        #print("DRM = ", np.shape(DRM))

        rmatrix = self.combine_matrices(RM, DRM)


        # trying to minimize strength of the correctors. does not work actually
        #rmatrix = self.combine_matrices(rmatrix, 10*np.eye(np.shape(DRM)[0]))

        #print("rmatrix = ", rmatrix)
        orbit = np.append(orbit, disp)

        # trying to minimize strength of the correctors. does not work actually
        #orbit = np.append(orbit, np.zeros(np.shape(DRM)[0]))

        # bpm weights
        #bpm_weights = np.eye(len(orbit))
        bpm_weights = np.array([bpm.weight for bpm in self.bpms])
        bpm_weights_diag = np.diag( np.append(bpm_weights, [bpm_weights, bpm_weights, bpm_weights]))
        #print("bpm_weights = ", np.shape(bpm_weights), len(orbit))

        self.orbit_svd = OrbitSVD(resp_matrix=rmatrix, orbit=orbit, weights=bpm_weights_diag, epsilon_x=epsilon_x, epsilon_y=epsilon_x)

        #self.orbit_svd = LInfinityNorm(resp_matrix=rmatrix, orbit=orbit, weights=bpm_weights_diag, epsilon_x=epsilon_x,
        #                          epsilon_y=epsilon_x)
        angle = self.orbit_svd.apply()
        ncor = len(cor_list)
        for i, cor in enumerate(np.append(self.hcors, self.vcors)):
            if print_log:
                print("correction:", cor.id," angle before: ", cor.angle*1000, "  after:", angle[i]*1000,angle[ncor+i]*1000)
            cor.angle -= ((1 - alpha) * angle[i] + alpha* angle[ncor + i])

        self.lat.update_transfer_maps()
        if p_init is not None:
            p_init.x = -angle[-4]
            p_init.px = -angle[-3]
            p_init.y  = -angle[-2]
            p_init.py = -angle[-1]
        return 0



class Orbit:
    def __init__(self, lattice, empty=False):
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
        if not empty:
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
        for i, elem in enumerate(self.lat.sequence):
            if elem.__class__ == Monitor:
                if bpm_list is None or elem.id in bpm_list:
                    try:
                        elem.weight
                    except:
                        elem.weight = 1.
                    elem.s = L+elem.l/2.
                    elem.x_ref = 0.
                    elem.y_ref = 0.
                    elem.lat_inx = i
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
        for i, elem in enumerate(self.lat.sequence):
            if elem.__class__ == Vcor:
                if cor_list is None or elem.id in cor_list:
                    elem.s = L+elem.l/2.
                    elem.lat_inx = i
                    self.vcors.append(elem)
            elif elem.__class__ == Hcor:
                if cor_list is None or elem.id in cor_list:
                    elem.s = L+elem.l/2.
                    elem.lat_inx = i
                    self.hcors.append(elem)
            L += elem.l
        if len(self.hcors) == 0:
            print("there are not horizontal correctors")
        if len(self.vcors) == 0:
            print("there are not vertical correctors")


    def create_types(self, types, remove_elems):
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

            particles.append(copy.copy(p))
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
        rmatrix = ResponseMatrix()
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
        #print ("measure = ", p_init.x, p_init.y, p_init.px, p_init.py, p_init.E)
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
        rmatrix = ResponseMatrix()
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

        rmatrix = ResponseMatrix()
        rmatrix.bpm_names = [b.id for b in self.bpms]
        rmatrix.cor_names = np.append(np.array([c.id for c in self.hcors]), np.array([c.id for c in self.vcors]))
        rmatrix.matrix = self.resp
        rmatrix.mode = "radian"
        return rmatrix

    def linac_response_matrix_r(self, tw_init=None):
        """
        calculation of ideal response matrix

        :param lattice: class MagneticLattice
        :param tw_init: if tw_init == None, function tries to find periodical solution
        :return: orbit.resp
        """
        m = len(self.bpms)
        nx = len(self.hcors)
        ny = len(self.vcors)
        self.resp = zeros((2 * m, nx + ny))


        for j, cor in enumerate([item for sublist in [self.hcors, self.vcors] for item in sublist]):
            print(j,"/", nx+ny, cor.id)
            Ra = np.eye(6)
            E = tw_init.E
            for i, elem in enumerate(self.lat.sequence):
                if i<cor.lat_inx:
                    #print(i, cor.lat_inx, E)
                    E += elem.transfer_map.delta_e
                    continue

                Rb = elem.transfer_map.R(E)
                #print("Rb = ", Rb)
                Ra = dot(Rb, Ra)
                E += elem.transfer_map.delta_e
                if elem in self.bpms:
                    #print("R12 = ", cor.lat_inx, cor.id, cor.s, Ra[0,1])
                    #exit(0)
                    n = self.bpms.index(elem)
                    #v_resp[n, j] = Ra[0,1]
                    #v_resp[n, j] = Ra[0, 1]
                    if cor.__class__ == Hcor:
                        self.resp[n, j] = Ra[0,1]
                        #print(Ra[0,1])
                    else:
                        self.resp[n+m, j] = Ra[2, 3]

        rmatrix = ResponseMatrix()
        rmatrix.bpm_names = [b.id for b in self.bpms]
        rmatrix.cor_names = np.append(np.array([c.id for c in self.hcors]), np.array([c.id for c in self.vcors]))
        rmatrix.matrix = self.resp
        rmatrix.mode = "radian"
        return rmatrix

    def linac_response_matrix_meas(self, tw_init=None):
        """
        calculation of ideal response matrix

        :param lattice: class MagneticLattice
        :param tw_init: if tw_init == None, function tries to find periodical solution
        :return: orbit.resp
        """
        #self.optical_func_params(tw_init=tw_init)
        match_ic = False
        m = len(self.bpms)
        nx = len(self.hcors)
        ny = len(self.vcors)
        #h_resp = zeros((m, nx))
        #v_resp = zeros((m, ny))
        add_i = 0
        if match_ic:
            add_i = 4
        self.resp = zeros((2 * m, nx + ny + add_i))
        s = [bpm.s for bpm in self.bpms]
        X0, Y0 = self.read_virtual_orbit(p_init=Particle( E=tw_init.E))
        #plt.plot(s, X0,"ro-",  s, Y0,"bo-")
        #plt.show()
        Or0 = np.append(X0, Y0)
        for j, cor in enumerate([item for sublist in [self.hcors, self.vcors] for item in sublist]):
            print(j,"/", nx+ny, cor.id)
            cor.angle = 0.0001
            self.lat.update_transfer_maps()
            start = time()
            X1, Y1 = self.read_virtual_orbit(p_init=Particle(E=tw_init.E))
            print(time() - start)

            Or1 = np.append(X1, Y1)
            self.resp[:, j] = (Or1 - Or0)/cor.angle
            cor.angle = 0.00
        self.lat.update_transfer_maps()
        if match_ic:
            for i, par in enumerate(["x", "px", "y", "py"]):
                print(i)
                p_i = Particle(E = tw_init.E)
                p_i.__dict__[par] = 0.0001
                #p2 = copy.deepcopy(p_i)
                print ("measure = ", p_i.x, p_i.y, p_i.px, p_i.py, p_i.E)
                X1, Y1 = self.read_virtual_orbit(p_init=p_i)
                #plt.plot(s, X1, "ro-", s, Y1, "bo-")
                #plt.show()
                Or1 = np.append(X1, Y1)
                self.resp[:, nx + ny + i] = (Or1 - Or0) /0.0001
        X1, Y1 = self.read_virtual_orbit(p_init=Particle(E=tw_init.E))
        rmatrix = ResponseMatrix()
        rmatrix.bpm_names = [b.id for b in self.bpms]
        rmatrix.cor_names = np.append(np.array([c.id for c in self.hcors]), np.array([c.id for c in self.vcors]))
        rmatrix.matrix = self.resp
        rmatrix.mode = "radian"
        return rmatrix

    def disp_measur(self, E0):
        X0, Y0 = self.read_virtual_orbit(p_init=Particle(p=0.0, E=E0))
        X1, Y1 = self.read_virtual_orbit(p_init=Particle(p=0.01, E=E0))
        E = np.array([bpm.E for bpm in self.bpms])
        Dx0 = (X1 - X0) * E / E0 / 0.01
        Dy0 = (Y1 - Y0) * E / E0 / 0.01
        return Dx0, Dy0

    # def disp_measur2(self, E0):
    #     p_list = np.array(lattice_track(self.lat, Particle(p=0.0, E=E0)))
    #     p_list2 = np.array(lattice_track(self.lat, Particle(p=0.001, E=E0)))
    #     s = np.array([bpm.s for bpm in self.bpms])
    #     p_s = np.array([p.s for p in p_list])
    #     indces = np.where(np.in1d(p_s, s))[0][::2]
    #     #print(len(s), len(indces))
    #     #print(indces)
    #     E = np.array([p.E for p in p_list[indces]])
    #     X0 = np.array([p.x for p in p_list[indces]])
    #     X1 = np.array([p.x for p in p_list2[indces]])
    #     Y0 = np.array([p.y for p in p_list[indces]])
    #     Y1 = np.array([p.y for p in p_list2[indces]])
    #     Dx0 = (X1 - X0) * E / E0 / 0.01
    #     Dy0 = (Y1 - Y0) * E / E0 / 0.01
    #     return Dx0, Dy0

    def linac_disp_response_matrix(self, tw_init):
        """
        calculation of ideal response matrix

        :param lattice: class MagneticLattice
        :param tw_init: if tw_init == None, function tries to find periodical solution
        :return: orbit.resp
        """
        #self.optical_func_params(tw_init=tw_init)

        m = len(self.bpms)
        nx = len(self.hcors)
        ny = len(self.vcors)
        #h_resp = zeros((m, nx))
        #v_resp = zeros((m, ny))

        self.resp = zeros((2 * m, nx + ny))
        s = [bpm.s for bpm in self.bpms]
        Dx0, Dy0 = self.disp_measur(E0=tw_init.E)
        #plt.plot(s, Dx0,"ro-",  s, Dy0,"bo-")
        #plt.show()
        D0 = np.append(Dx0, Dy0)
        for j, cor in enumerate([item for sublist in [self.hcors, self.vcors] for item in sublist]):
            print(j,"/", nx+ny, cor.id)
            cor.angle = 0.001
            self.lat.update_transfer_maps()
            # cor.transfer_map = self.lat.method.create_tm(cor)
            start = time()
            Dx1, Dy1 = self.disp_measur(E0=tw_init.E)
            print(time() - start)
            #plt.plot(s, Dx1, "ro-", s, Dy1, "bo-")
            #plt.show()
            D1 = np.append(Dx1, Dy1)
            self.resp[:, j] = (D1 - D0)/cor.angle
            cor.angle = 0.00
            # cor.transfer_map = self.lat.method.create_tm(cor)

        rmatrix = ResponseMatrix()
        rmatrix.bpm_names = [b.id for b in self.bpms]
        rmatrix.cor_names = np.append(np.array([c.id for c in self.hcors]), np.array([c.id for c in self.vcors]))
        rmatrix.matrix = self.resp
        rmatrix.mode = "radian"
        return rmatrix


    def apply_svd(self, resp_matrix, misallign, weight=None, alpha=1.e-4):
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
            if s[i] < alpha:
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
        print(monitors)
        start = time()
        angle = self.apply_svd(self.resp, monitors, weight=weights)
        #print(len(angle))
        print("correction = ", time() - start)
        for i, cor in enumerate(np.append(self.hcors, self.vcors)):
            if self.mode == "ampere":
                #print "ampere"
                cor.dI = -angle[i]
            else:
                print("correction", cor.angle*1000, angle[i]*1000, cor.id)
                cor.angle -= angle[i]
                #print(cor)
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
        #for i, cor in enumerate(np.append(self.hcors, self.vcors)):
        #    print(cor.id, cor.angle, cor)
        if p_init is not None:
            p_init.x = -angle[-4]
            p_init.px = -angle[-3]
            p_init.y  = -angle[-2]
            p_init.py = -angle[-1]
        return 0

    def disp_correction(self, monitors, p_init=None, alpha=1e-3):
        m = len(self.bpms)
        #monitors = zeros(2*m)
        weights = eye(len(monitors))
        for i, bpm in enumerate(self.bpms):
            #monitors[i] = bpm.x
            #monitors[i+m] = bpm.y
            weights[i, i] = bpm.weight
            weights[i+m, i+m] = bpm.weight

        start = time()
        angle = self.apply_svd(self.resp, monitors, weight=weights, alpha=alpha)

        print("correction = ", time() - start)
        for i, cor in enumerate(np.append(self.hcors, self.vcors)):
            if self.mode == "ampere":
                #print "ampere"
                cor.dI = -angle[i]
            else:
                #print len(np.append(self.hcors, self.vcors)), i, len(angle)
                cor.angle -= angle[i]
                print(cor.id, ".angle=", cor.angle)

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
        print(monitors)
        print(elem_response)
        poss = self.apply_svd(elem_response, monitors)
        print("correction = ", time() - start)
        ix = 0
        iy = 0
        for elem in self.lat.sequence:
            if ix<len(self.htypes) and elem.id == self.htypes[ix].id:
                print ("X: ", elem.__class__.__name__, elem.id, elem.dx, poss[ix])
                elem.dx = poss[ix]
                #self.hquads[ix].dx -= poss[ix]
                ix += 1

            if iy<len(self.vtypes) and elem.id == self.vtypes[iy].id:
                print("Y: ", elem.__class__.__name__, elem.id, elem.dy, poss[iy])
                elem.dy = poss[iy+len(self.htypes)]
                #self.vquads[iy].dy -= poss[iy+len(self.hquads)]
                iy += 1
        p = Particle(x=poss[-4], px=poss[-3], y=poss[-2], py=poss[-1])
        #print poss[-5:]
        self.lat.update_transfer_maps()
        return p

    def misalignment_rm(self, p_init, elem_types, remove_elem):
        shift = 0.001
        m = len(self.bpms)
        self.create_types(elem_types, remove_elem)
        nx = len(self.htypes)
        ny = len(self.vtypes)
        print(nx, ny, m)
        real_resp = zeros((m * 2, nx + ny + 4))
        X0, Y0 = self.read_virtual_orbit(p_init=copy.deepcopy(p_init))
        Or0 = np.append(X0, Y0)
        #bpms = copy.deepcopy(self.bpms)
        for ix, hquad in enumerate(self.htypes):
            print("measure X - ", ix, "/", nx)
            hquad.dx += shift
            self.lat.update_transfer_maps()
            X1, Y1 = self.read_virtual_orbit(p_init=copy.deepcopy(p_init))
            Or1 = np.append(X1, Y1)
            real_resp[:, ix] = (Or1 - Or0) / shift

            hquad.dx -= shift
            self.lat.update_transfer_maps()

        for iy, vquad in enumerate(self.vtypes):
            print("measure Y - ", iy, "/", ny)
            vquad.dy += shift
            self.lat.update_transfer_maps()
            X1, Y1 = self.read_virtual_orbit(p_init=copy.deepcopy(p_init))
            Or1 = np.append(X1, Y1)
            real_resp[:, iy] = (Or1 - Or0) / shift
            #plt.plot([bpm.s for bpm in self.bpms], [bpm.x for bpm in self.bpms], "r")
            #plt.plot([bpm.s for bpm in self.bpms], [bpm.y for bpm in self.bpms], "b")
            #plt.show()
            vquad.dy -= shift
            self.lat.update_transfer_maps()

        X1, Y1 = self.read_virtual_orbit(p_init=copy.deepcopy(p_init))

        for i, par in enumerate(["x", "px", "y", "py"]):
            print(i)
            p_i = Particle(E=p_init.E)
            p_i.__dict__[par] = 0.0001
            # print p_i.x, p_i.px, p_i.y, p_i.py, p_i.E
            p2 = copy.deepcopy(p_i)
            X1, Y1 = self.read_virtual_orbit(p_init=p2)
            Or1 = np.append(X1, Y1)
            # print ("energy = ", p2.E)
            # plt.plot([bpm.s for bpm in orbit.bpms], [bpm.x for bpm in orbit.bpms], "r")
            # plt.plot([bpm.s for bpm in orbit.bpms], [bpm.y for bpm in orbit.bpms], "b")
            # plt.show()
            real_resp[:, nx + ny + i] = (Or1 - Or0) / 0.0001
            #for j, bpm in enumerate(self.bpms):
            #    real_resp[:, nx + ny + i] = (Or1 - Or0) / 0.0001
                #real_resp[j + m, nx + ny + i] = (bpm.y - bpms[j].y) / 0.0001
                # print j+m, nx + ny + i, (bpm.x - bpms[j].x)/0.00001
        # print real_resp[:,-5:]
        X1, Y1 = self.read_virtual_orbit(p_init=copy.deepcopy(p_init))
        return real_resp

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

#def show_currents( elems, alpha):
#    print( "******* displaying currents - START ********")
#    for elem in elems:
#        if elem.dI == 0:
#            continue
#        n = len(elem.id)
#        n2 = len(str(elem.I + elem.dI))
#        n3 = len(str(elem.I))
#        print( elem.id, " "*(10-n) + "<-- ", elem.I + elem.dI,  " "*(18-n2)+ " was = ", elem.I, " "*(18-n3) + " dI = ", elem.dI, "x", alpha)
#    print ("******* displaying currents - END ********")




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