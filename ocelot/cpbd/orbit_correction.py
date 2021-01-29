__author__ = 'Sergey Tomin'

import numpy as np
from numpy.linalg import svd
from scipy.interpolate import splrep, splev
from scipy.optimize import linprog
from ocelot.cpbd.match import closed_orbit
from ocelot.cpbd.track import *
from ocelot.cpbd.response_matrix import *
from scipy.linalg import block_diag, lstsq
import copy
import json
from time import sleep, time
import logging

logger = logging.getLogger(__name__)


class OrbitSVD:
    def __init__(self, epsilon_x=0.001, epsilon_y=0.001):
        self.epsilon_x = epsilon_x
        self.epsilon_y = epsilon_y

    def apply(self, resp_matrix, orbit, weights=None):
        if weights is None:
            weights = np.eye(len(orbit))
        resp_matrix = np.dot(weights, resp_matrix)
        misallign = np.dot(weights, orbit)

        U, s, V = svd(resp_matrix)
        s_inv = np.zeros(len(s))
        s_max = max(s)
        for i in range(len(s)):
            #print("S[",i,"]=", s[i], "s max = ", s_max)
            if i < int(len(s)/2.):
                epsilon = self.epsilon_x
            else:
                epsilon = self.epsilon_y
            if s[i] <= s_max * epsilon:
                s_inv[i] = 0.
            else:
                s_inv[i] = 1. / s[i]
        Sinv = np.zeros((np.shape(U)[0], np.shape(V)[0]))
        Sinv[:len(s), :len(s)] = np.diag(s_inv)
        Sinv = np.transpose(Sinv)
        A = np.dot(np.transpose(V), np.dot(Sinv, np.transpose(U)))
        angle = np.dot(A, misallign)
        #a, b, _, _ = np.linalg.lstsq(resp_matrix, misallign, rcond=self.epsilon_x)
        #print(np.shape(a), np.shape(angle))
        #print(f"angle = {angle}, a = {a}")
        logger.debug("max(abs(angle)) = " + str(np.max(np.abs(angle))) + " min(abs(angle)) = " + str(np.min(np.abs(angle))))
        return angle


class MICADO(OrbitSVD):
    """
    Iterative method based on CERN-ISR-MA/73-17, 1973.
    Strategy:
    1. each corrector is tested singly. For each corrector we have residual vector: r = b + A * x.
        Corrector (x_j) and column A_j (response matrix) which correspond to minimum residual are exchanged with first
        position.
    2. Two correctors are tested. the first one which is defined from step 1 and second one iteratively is tested
        j = 2, ... n
    The iteration is stopped when the peak-to-peak amplitude of the residual vector is smaller than a fixed value given
    in advance.
    """
    def __init__(self, epsilon_x=0.001, epsilon_y=0.001, epsilon_ksi=1e-5):
        super(MICADO, self).__init__(epsilon_x=epsilon_x, epsilon_y=epsilon_y)
        self.epsilon_ksi = epsilon_ksi

    @staticmethod
    def swap_columns(m, col1, col2):
        """
        method to swap columns in matrix m
        :param m: matrix
        :param col1: index
        :param col2: index
        :return:
        """
        c1 = np.copy(m[:, col1])
        c2 = np.copy(m[:, col2])
        m[:, col1] = c2[:]
        m[:, col2] = c1[:]

    def solver_svd(self, resp_matrix, orbit, weights):
        svd_method = OrbitSVD(self.epsilon_x, self.epsilon_y)
        angles = svd_method.apply(resp_matrix, orbit, weights)
        return angles

    def solver_lstsq(self, resp_matrix, orbit, weights):
        if weights is not None:
            # weights = np.eye(len(orbit))
            resp_matrix = np.dot(weights, resp_matrix)
            orbit = np.copy(np.dot(weights, orbit))
            
        #misallign = np.dot(weights, orbit)
        if np.shape(resp_matrix)[1] == 1:
            A = np.dot(resp_matrix.T, resp_matrix)
            Ainv = 1./A if A[0, 0] != 0 else A
            a = np.dot(np.dot(Ainv, resp_matrix.T), orbit)
        else:
            a, _, _, _ = np.linalg.lstsq(resp_matrix, orbit, rcond=self.epsilon_x)
            #a, _, _, _ = lstsq(resp_matrix, misallign, cond=self.epsilon_x, lapack_driver="gelsy")
        return a

    def apply(self, resp_matrix, orbit, weights=None):
        weights = None
        bpm_num = np.shape(resp_matrix)[0]
        mask = np.arange(np.shape(resp_matrix)[1])
        self.orb_residual = []
        angle = np.zeros(np.shape(resp_matrix)[1])
        angles_part = []
        resp_matrix_part = np.empty((bpm_num, 0))
        resp_vector = np.empty((bpm_num, 0))
        start = time()
        for n in range(np.shape(resp_matrix)[1]):

            resp_matrix_part = np.hstack((resp_matrix_part, resp_vector))
            residual = []
            for i in range(n, np.shape(resp_matrix)[1]):
                resp_vector = resp_matrix[:, i].reshape(bpm_num, -1)
                resp_matrix_tmp = np.hstack((resp_matrix_part, resp_vector))

                #angles_part = self.solver_svd(resp_matrix_tmp, orbit, weights)
                angles_part = self.solver_lstsq(resp_matrix_tmp, orbit, weights)

                new_orbit = np.dot(resp_matrix_tmp, angles_part)
                res = np.sqrt(np.sum((new_orbit - orbit)**2))
                residual.append(res)

            index = np.argmin(residual) + n
            self.orb_residual.append(min(residual))

            resp_vector = np.copy(resp_matrix[:, index].reshape(bpm_num, -1))

            self.swap_columns(resp_matrix, n, index)
            mask[[n, index]] = mask[[index, n]]
            if len(self.orb_residual) > 1 and self.orb_residual[-2] - self.orb_residual[-1] < self.epsilon_ksi:
                logger.info(" MICADO: number of correctors " + str(n))
                angle[:len(angles_part)] = angles_part[:]
                break
        print( time() - start, " sec")
        return angle[np.argsort(mask)]


class LInfinityNorm(OrbitSVD):
    def __init__(self, epsilon_x=0.001, epsilon_y=0.001):
        OrbitSVD.__init__(self, epsilon_x=epsilon_x, epsilon_y=epsilon_y)

    def apply(self, resp_matrix, orbit, weights=None):
        m, n = np.shape(resp_matrix)
        f = np.zeros(n + 1)
        f[-1] = 1
        Ane = np.vstack((np.hstack((resp_matrix, -np.ones((m, 1)))), np.hstack((-resp_matrix, -np.ones((m, 1))))))
        bne = np.vstack((+orbit, -orbit))
        res = linprog(f, A_ub=Ane, b_ub=bne)
        x = res["x"][:-1]
        return x


class Orbit(object):
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
        self.orbit_solver = OrbitSVD()

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
            logger.warning(" create_correctors: there are not horizontal correctors")
        if len(self.vcors) == 0:
            logging.warning(" create_correctors: there are not vertical correctors")

    def get_ref_orbit(self):
        for bpm in self.bpms:
            bpm.x_ref = 0.
            bpm.y_ref = 0.

    def get_orbit(self):

        #self.get_ref_orbit()

        m = len(self.bpms)
        orbit = np.zeros(2 * m)
        for i, bpm in enumerate(self.bpms):
            #print("get_orbit = ",bpm.id, bpm.x,  bpm.x_ref)
            orbit[i] = bpm.x - bpm.x_ref
            orbit[i+m] = bpm.y - bpm.y_ref
        return orbit


    def get_dispersion(self):
        m = len(self.bpms)
        disp = np.zeros(2 * m)
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
        #m = block_diag(mat1, mat2)
        return rm

    def correction(self, alpha=0, beta=0, p_init=None, print_log=True):
        """
        Method to find corrector kicks using SVD. bpm weights are ignored for a moment but everything ready to immplement.

        :param alpha: 0 - 1, trade off between orbit and dispersion correction, 0 - only orbit, 1 - only dispersion
        :param epsilon_x: cut s-matrix diag for x-plane, if s[i] < s_max * epsilon: s_inv[i] = 0. else s_inv[i] = 1/s[i]
        :param epsilon_y: cut s-matrix diag for y-plane, if s[i] < s_max * epsilon: s_inv[i] = 0. else s_inv[i] = 1/s[i]
        :param beta: weight for suppress large kicks
        :param p_init: particle initial conditions. Removed in that version.
        :param print_log:
        :return:
        """
        #TODO: initial condition for particle was removed. Add it again
        cor_list = [cor.id for cor in np.append(self.hcors, self.vcors)]
        bpm_list = [bpm.id for bpm in self.bpms]
        orbit = (1 - alpha) * self.get_orbit()

        RM = (1 - alpha) * self.response_matrix.extract(cor_list=cor_list, bpm_list=bpm_list)
        logger.debug(" shape(RM) = " + str(np.shape(RM)))
        # dispersion
        if alpha != 0:
            disp = alpha * self.get_dispersion()
        else:
            disp = alpha * np.array(orbit)

        if self.disp_response_matrix is not None:
            DRM = alpha * self.disp_response_matrix.extract(cor_list=cor_list, bpm_list=bpm_list)
        else:
            DRM = np.zeros_like(RM)
        logger.debug(" shape(DRM) = " + str(np.shape(DRM)))

        # Add beta coefficient for minimizing strength of the correctors.
        b_mat = beta * np.eye(np.shape(DRM)[0])
        b_orbit = np.zeros(np.shape(DRM)[0])

        rmatrix = block_diag(RM, DRM, b_mat)
        orbit = np.append(orbit, [disp, b_orbit])
        logger.debug(" Combine: shape(RM + DRM + beta) = " + str(np.shape(rmatrix)) +
                     " shape(orbit) = " + str(np.shape(orbit)))


        # add bpm weights
        bpm_weights = np.array([bpm.weight for bpm in self.bpms])
        bpm_weights_diag = np.diag(np.append(bpm_weights, [bpm_weights, bpm_weights, bpm_weights]))
        logger.debug(" shape(bpm weight) = " + str(np.shape(bpm_weights_diag)))

        #if beta > 0:
        bpm_weights_diag = self.combine_matrices(bpm_weights_diag, np.diag(np.append(bpm_weights, [bpm_weights])))
        logger.debug(" beta > 0: shape(bpm weight) = " + str(np.shape(bpm_weights_diag)))

        #self.orbit_correction_method = self.get_correction_solver(resp_matrix=rmatrix, orbit=orbit,
        #                                                      weights=bpm_weights_diag, epsilon_x=epsilon_x,
        #                                                      epsilon_y=epsilon_y)

        #self.orbit_svd = LInfinityNorm(resp_matrix=rmatrix, orbit=orbit, weights=bpm_weights_diag, epsilon_x=epsilon_x,
        #                          epsilon_y=epsilon_x)
        angle = self.orbit_solver.apply(resp_matrix=rmatrix, orbit=orbit, weights=bpm_weights_diag)
        ncor = len(cor_list)
        for i, cor in enumerate(np.append(self.hcors, self.vcors)):
            if print_log:
                print("correction:", cor.id," angle before: ", cor.angle*1000, "  after:", angle[i]*1000,angle[ncor+i]*1000)
            cor.angle -= ((1 - alpha) * angle[i] + alpha * angle[ncor + i])

        self.lat.update_transfer_maps()
        if p_init is not None:
            p_init.x = -angle[-4]
            p_init.px = -angle[-3]
            p_init.y  = -angle[-2]
            p_init.py = -angle[-1]
        return 0


class NewOrbit(Orbit):
    def __init__(self, lattice, rm_method=None, disp_rm_method=None, empty=False):
        super(NewOrbit, self).__init__(lattice, rm_method=rm_method, disp_rm_method=disp_rm_method, empty=empty)
