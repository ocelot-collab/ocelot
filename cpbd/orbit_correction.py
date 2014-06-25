__author__ = 'Sergey Tomin'

from ocelot.cpbd.beam import Particle, Twiss
from ocelot.cpbd.errors import *
from ocelot.cpbd.matching import closed_orbit
from ocelot.cpbd.optics import *
import copy
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev
from scipy.integrate import simps, quad
from time import time
from numpy.linalg import svd
from numpy import diag, transpose, linspace, shape


class BPM(object):
    def __init__(self, id = None):
        self.id = id
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


class Orbit:
    def __init__(self, lattice = None):
        #self.monitors = []
        #self.correctors = []
        self.lattice = 0
        self.h_resp = []
        self.v_resp = []
        self.bpms = []
        self.hcors = []
        self.vcors = []
        self.nu_x = 0.
        self.nu_y = 0.
        self.resp = []
        if lattice != None:
            self.lattice_analysis(lattice)

    def lattice_analysis(self, lattice):
        """
        Analysis of lattice and creation of bpms and correctors lists
        self.create_BPM(lattice)
        self.create_COR(lattice)
        :param lattice: class MagneticLattice
        :return: orbit
        """
        self.create_BPM(lattice)
        self.create_COR(lattice)
        return self

    def create_BPM(self, lattice):
        """
        Search bpm in the lattice and create list of bpms
        :param lattice: class MagneticLattice
        :return: self.bpms - list of BPMs (class BPM)
        """
        self.bpms = []
        L = 0.
        for elem in lattice.sequence:
            if elem.type == "monitor":
                bpm = BPM(id = elem.id)
                bpm.s = L+elem.l/2.
                self.bpms.append(bpm)
            L += elem.l
        if len(self.bpms) == 0:
            print "there is not monitors"
        return self.bpms


    def create_COR(self, lattice):
        """
        Search correctors (horizontal and vertical) in the lattice and create list of hcors and list of vcors
        :param lattice: class MagneticLattice
        :return:
        """
        self.hcors = []
        self.vcors = []
        L = 0.
        for elem in lattice.sequence:
            if elem.type == "vcor":
                vcor = copy.copy(elem)
                vcor.s = L+elem.l/2.
                vcor.dI = 0.0001
                self.vcors.append(vcor)
            elif elem.type == "hcor":
                hcor = copy.copy(elem)
                hcor.s = L+elem.l/2.
                hcor.dI = 0.0001
                self.hcors.append(hcor)
            L += elem.l
        if len(self.hcors) == 0:
            print "there is not horizontal corrector"
        if len(self.vcors) == 0:
            print "there is not vertical corrector"


    def read_virtual_orbit(self, lattice):
        """
        searching closed orbit by function closed_orbit(lattice) and searching coordinates of beam at the bpm possitions
        :param lattice: class MagneticLattice
        :return: orbit.bpms
        """
        X = []
        Y = []
        self.particle0 = closed_orbit(lattice)
        p = self.particle0
        navi = Navigator(lattice)
        L = 0.
        for bpm in self.bpms:
            dz = bpm.s - L
            p = single_track(lattice, p, dz, navi)
            bpm.x = p.x
            bpm.y = p.y
            L = bpm.s
            X.append(p.x)
            Y.append(p.y)
        return array(X), array(Y)


    def optical_func_params(self, lattice):
        """
        Optical function parameters for correctors and bpms. It is needed for calculation of ideal response matrix:
        defining beta functions on the azimuth of correctors and bpms: beta_x, beta_y;
        defining phase shift between origin of lattice and element: mu_x, mu_y;
        defining tunes of beta functions of whole lattice: nu_x, nu_y = mu(totalLen)/(2*pi)
        :param lattice: class MagneticLattice
        :return:
        """
        tw0 = Twiss()
        tws = twiss(lattice, tw0, nPoints=int(lattice.totalLen/0.05))

        s = array(map(lambda tw:tw.s, tws))
        beta_x = array(map(lambda tw:tw.beta_x, tws))
        beta_y = array(map(lambda tw:tw.beta_y, tws))
        self.nu_x = simps(1./beta_x, s)/2./pi
        self.nu_y = simps(1./beta_y, s)/2./pi

        tck_x = splrep(s, beta_x)
        tck_y = splrep(s, beta_y)

        for bpm in self.bpms:
            z = linspace(0., bpm.s, num = int(bpm.s*100.))
            f_beta_x = splev(z, tck_x)
            f_beta_y = splev(z, tck_y)
            fx = 1./f_beta_x
            fy = 1./f_beta_y
            bpm.phi_x = simps(fx, z) #+ phi_x
            bpm.phi_y = simps(fy, z)# + phi_y
            bpm.beta_x = f_beta_x[-1]
            bpm.beta_y = f_beta_y[-1]

        for hcor in self.hcors:
            z = linspace(0., hcor.s, num = int(hcor.s*100.))
            f_beta_x = splev(z, tck_x)
            fx = 1./f_beta_x
            hcor.phi_x = simps(fx, z)
            hcor.beta_x = f_beta_x[-1]

        for vcor in self.vcors:
            z = linspace(0., vcor.s, num = int(vcor.s*100.))
            f_beta_y = splev(z, tck_y)
            fy = 1./f_beta_y
            vcor.phi_y = simps(fy, z)#
            vcor.beta_y = f_beta_y[-1]

    def read_response_matrix(self, dictionary):
        Energy = dictionary["energy"]
        print Energy
        m = len(self.bpms)
        nx = len(self.hcors)
        ny = len(self.vcors)
        real_resp = zeros((m*2, nx + ny))
        #self.read_virtual_orbit( lattice)
        #bpms = copy.deepcopy(orbit.bpms)
        for ix, hcor in enumerate(self.hcors):

            response = dictionary[hcor.id]

            for j, bpm in enumerate(self.bpms):

                real_resp[j, ix] = response[bpm.id][0]*(Energy/(hcor.angle_coef*hcor.length_iron*30.))
                real_resp[j+m, ix] = response[bpm.id][1]*(Energy/(hcor.angle_coef*hcor.length_iron*30.))

        for iy, vcor in enumerate(self.vcors):

            response = dictionary[vcor.id]

            for j, bpm in enumerate(self.bpms):
                real_resp[j, iy+nx] = response[bpm.id][0]*(Energy/(vcor.angle_coef*vcor.length_iron*30.))
                real_resp[j+m, iy+nx] = response[bpm.id][1]*(Energy/(vcor.angle_coef*vcor.length_iron*30.))
        self.resp = real_resp
        return self.resp

    def ideal_response_matrix(self, lattice):
        """
        calculation of ideal response matrix
        :param lattice: class MagneticLattice
        :return: orbit.resp
        """
        self.optical_func_params(lattice)

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
                #if mu_x<0.:
                #    print "mu_x = ", mu_x
                    #mu_x = 2*pi*self.nu_x + mu_x
                #mu_x = abs(mu_x)
                h_resp[i,j] = kx*sqrt(hcor.beta_x)*cos(mu_x - pi*self.nu_x)
            for n, vcor in enumerate(self.vcors):
                mu_y = abs(bpm.phi_y - vcor.phi_y)
                #if mu_y<0.:
                #    mu_y = 2*pi*self.nu_y + mu_y
                v_resp[i,n] = ky*sqrt(vcor.beta_y)*cos(mu_y - pi*self.nu_y)

        #self.h_resp = h_resp
        #self.v_resp = v_resp
        m = len(self.bpms)
        kx = len(self.hcors)
        ky = len(self.vcors)
        self.resp = zeros((2*m, kx + ky))
        self.resp[:m,:kx] = h_resp[:,:]
        self.resp[m:,kx:] = v_resp[:,:]
        #print "shape = ", shape(self.resp)
        return self.resp

    def apply_svd(self, resp_matrix, misallign):
        #print resp_matrix
        U, s,V = svd(resp_matrix)
        s_inv = 1./s
        for i in xrange(len(s)):
            if s_inv[i]>1.e+5:
                s_inv[i] = 0.

        Sinv = zeros((shape(U)[0],shape(V)[0]))
        Sinv[:len(s), :len(s)] = diag(s_inv)
        Sinv = transpose(Sinv)
        A = dot(transpose(V),dot(Sinv,transpose(U)))
        angle = dot(A, misallign)
        return angle

    def correction(self, lattice):
        m = len(self.bpms)
        monitors = zeros(2*m)
        for i, bpm in enumerate(self.bpms):
            monitors[i] = bpm.x
            monitors[i+m] = bpm.y
        start = time()
        angle = self.apply_svd(self.resp, monitors)
        print "correction = ", time() - start
        ix = 0
        iy = 0
        for elem in lattice.sequence:
            if ix<len(self.hcors) and elem.id == self.hcors[ix].id:
                elem.angle -= angle[ix]
                self.hcors[ix].angle -= angle[ix]
                ix += 1

            if iy<len(self.vcors) and elem.id == self.vcors[iy].id:
                elem.angle -= angle[iy+len(self.hcors)]
                self.vcors[iy].angle -= angle[iy+len(self.hcors)]
                iy += 1

        return lattice.update_transfer_maps()

    def quad_correction(self, lattice, quad_response):
        m = len(self.bpms)
        monitors = zeros(2*m)
        for i, bpm in enumerate(self.bpms):
            monitors[i] = bpm.x
            monitors[i+m] = bpm.y
        start = time()
        poss = self.apply_svd(quad_response, monitors)
        print "correction = ", time() - start
        ix = 0
        iy = 0
        for elem in lattice.sequence:
            if ix<len(self.hquads) and elem.id == self.hquads[ix].id:
                elem.dx -= poss[ix]
                #self.hquads[ix].dx -= poss[ix]
                ix += 1

            if iy<len(self.vquads) and elem.id == self.vquads[iy].id:
                elem.dy -= poss[iy+len(self.hquads)]
                #self.vquads[iy].dy -= poss[iy+len(self.hquads)]
                iy += 1

        return lattice.update_transfer_maps()


    def calc_track(self,lattice):
        part_list = trace_particle(lattice, self.particle0, nPoints = None)
        self.x_track = map(lambda p: p.x, part_list)
        self.y_track = map(lambda p: p.y, part_list)
        self.s_track = map(lambda p: p.s, part_list)

def draw_orbit(orbit,lattice, traject = True):
    if traject:
        bpm_draw = "o"
        part_list = trace_particle(lattice, orbit.particle0, nPoints = None)
        #plt.plot(map(lambda p: p.s, part_list), map(lambda p: p.y, part_list), "-")
    else:
        bpm_draw = "o-"
    #plt.plot(map(lambda p: p.s, orbit.bpms), map(lambda p: p.y, orbit.bpms), bpm_draw)
    #plt.grid(True)
    return  part_list

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
            elem.transfer_map = create_transfer_map(elem, energy = lattice.energy)
            #print elem.transfer_map.b(1)
    return lattice#.update_transfer_maps()

def restore_corrector(corrector, lattice):
    for elem in lattice.sequence:
        if elem.id == corrector.id:
            elem.angle -= corrector.dI
            elem.transfer_map = create_transfer_map(elem, energy = lattice.energy)
    return lattice#.update_transfer_maps()

def change_quad_position(quad, lattice, dx = 0., dy = 0.):
    for elem in lattice.sequence:
        if elem.id == quad.id:
            elem.dx += dx
            elem.dy += dy
            elem.transfer_map = create_transfer_map(elem, energy = lattice.energy)
    return lattice.update_transfer_maps()


def measure_response_matrix(orbit, lattice):

    m = len(orbit.bpms)
    real_resp = zeros((m*2, len(orbit.hcors)+len(orbit.vcors)))
    orbit.read_virtual_orbit( lattice)
    bpms = copy.deepcopy(orbit.bpms)
    for ix, hcor in enumerate(orbit.hcors):
        print "measure X - ", ix,"/",len(orbit.hcors)
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
    print nx, ny, m
    real_resp = zeros((m*2, nx + ny))
    orbit.read_virtual_orbit(lattice)
    bpms = copy.deepcopy(orbit.bpms)
    for ix, hquad in enumerate(orbit.hquads):
        print "measure X - ", ix,"/",nx
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


def test(lattice, errors):
    lat_errors = errors_seed(lattice, errors, nsuperperiods = 6)

    orbit = Orbit()
    orbit.lattice_analysis(lat_errors)
    orbit.ideal_response_matrix(lat_errors)
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
        print hcor.id, hcor.angle

    orbit.read_virtual_orbit(lat_corrected2)
    draw(orbit, lat_corrected2, traject = True)
    plt.show()




if __name__ == "__main__":
    from xframework.cpbd.elements import *
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