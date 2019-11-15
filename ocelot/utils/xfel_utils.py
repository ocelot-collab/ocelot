'''
functions common to fel decks
'''

# import scipy.special as sf
# import scipy.integrate as integrate
# from numpy.polynomial.chebyshev import *
import os
import time
import numpy as np
from numpy import inf, complex128, complex64
from copy import copy, deepcopy
import ocelot
from ocelot import ocelog 
from ocelot.common.globals import *  #import of constants like "h_eV_s" and "speed_of_light" 
from ocelot.common.py_func import * 
from ocelot.common.math_op import *
from ocelot.rad.undulator_params import *

from ocelot.gui import *

# from ocelot.optics.utils import *
# from ocelot.rad.undulator_params import *
# from ocelot.rad.fel import *
from ocelot.adaptors.genesis import *
# from ocelot.adaptors.genesis4 import *

import multiprocessing
nthread = multiprocessing.cpu_count()

from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.elements import *
from ocelot.cpbd.optics import *

_logger = logging.getLogger(__name__)

'''
SELF-SEEDING - relevant
'''





def dfl_st_cpl(dfl, theta_b, inp_axis='y', s_start=None):
    _logger.info('introducing spatio-temporal coupling')
    start = time.time()
    direction = 1
    if s_start == None:
        s_start = n_moment(dfl.scale_z(), dfl.int_z(), 0, 1)

    # dfl = deepcopy(dfl)
    shift_z_scale = dfl.scale_z() - s_start
    shift_m_scale = shift_z_scale / np.tan(theta_b)
    shift_m_scale[np.where(shift_m_scale > 0)] = 0
    shift_pix_scale = np.floor(shift_m_scale / dfl.dx).astype(int)

    # pix_start=np.where(shift_m_scale==0)[0][0]

    fld = dfl.fld
    if inp_axis == 'y':
        for i in np.where(shift_m_scale != 0)[0]:
            fld[i, :, :] = np.roll(fld[i, :, :], -direction * shift_pix_scale[i], axis=0)
            # if direction==1:
                # fld[i,:abs(shift_pix_scale[i]),:]=0

    elif inp_axis == 'x':
        for i in np.where(shift_m_scale != 0)[0]:
            fld[i, :, :] = np.roll(fld[i, :, :], -direction * shift_pix_scale[i], axis=1)
            # if direction==1:
                # fld[i,:,:abs(shift_pix_scale[i])]=0

    dfl.fld = fld
    t_func = time.time() - start
    _logger.debug(ind_str + 'done in {:.2f} sec'.format(t_func))
    # print('      done in %.2f ' % t_func + 'sec')
    return dfl


def dfl_hxrss_filt(dfl, trf, s_delay, st_cpl=1, enforce_padn=None, res_per_fwhm=6, fft_method='mp', dump_proj=0, debug=1):
    # needs optimizing?
    # import matplotlib.pyplot as plt

    nthread = multiprocessing.cpu_count()
    if nthread > 8:
        nthread = int(nthread * 0.9)  # not to occupy all CPUs on login server
    # print('  HXRSS dfl filtering')
    _logger.info('HXRSS dfl filtering')
    start = time.time()
    # klpos, krpos, cwidth = FWHM(trf.k, 1.0-np.abs(trf.tr))
    
    trf_m = np.abs(trf.tr) # temp array to measure fwhm
    trf_m /= np.amax(trf_m)
    trf_m = 1 - trf_m
    idx = fwhm3(trf_m)[2]
    cwidth = abs(trf.k[idx[1]] - trf.k[idx[0]]) # width the crystal transfer function hole (amplitude)
    # cwidth = fwhm(trf.k, trf_m)
    if hasattr(trf,'compound'):
        if trf.compound:
            cwidth = trf.dk
    
    dfl.to_domain('st')
    dk_old = 2 * pi / dfl.Lz()
    dk = cwidth / res_per_fwhm
    padn = np.ceil(dk_old / dk).astype(int)
    if np.mod(padn, 2) == 0 and padn != 0:  # check for odd
        padn = int(padn + 1)
    
    _logger.debug(ind_str + 'calculated padn = {}'.format(padn))
    if enforce_padn!=None:
        padn=enforce_padn
        _logger.debug(ind_str + 'enforced padn = {}'.format(padn))
        
    
    dfl_z = dfl.scale_z()
    dfl_z_mesh = dfl_z[-1]-dfl_z[0]
    if s_delay > dfl_z_mesh * padn:
        raise Exception('s_delay %.2e is larger that the padded dfl %.2e Consider using enforce_padn > %.2f' %(s_delay, dfl_z_mesh * padn, s_delay / dfl_z_mesh))
        
    if dump_proj:
        _logger.debug(ind_str + 'shape_before '+str(dfl.shape()))
        dfl_pad_z(dfl, padn)
        _logger.debug(ind_str + 'shape_after '+str(dfl.shape()))
        t1 = time.time()
        t_l_scale = dfl.scale_z()
        # t_l_int_b=dfl.int_z()
        t2 = time.time()

        dfl.fft_z(method=fft_method, nthread=multiprocessing.cpu_count())

        t3 = time.time()
        f_l_scale = dfl.scale_z()  # frequency_large_scale (wavelength in m)
        # f_l_int_b=dfl.int_z()
        t4 = time.time()

        f_l_filt = dfl_trf(dfl, trf, mode='tr', dump_proj=dump_proj)

        t5 = time.time()
        f_l_int_a = dfl.int_z()
        t6 = time.time()

        dfl.fft_z(method=fft_method, nthread=multiprocessing.cpu_count())

        t7 = time.time()
        t_l_int_a = dfl.int_z()
        t_l_pha_a = dfl.ang_z_onaxis()
        t8 = time.time()

        if st_cpl:
            dfl_st_cpl(dfl, trf.thetaB)

        dfl_shift_z(dfl, s_delay, set_zeros=0)
        _logger.debug(ind_str + 'shape_before '+str(dfl.shape()))
        dfl_pad_z(dfl, -padn)
        _logger.debug(ind_str + 'shape_after '+str(dfl.shape()))

        t_func = time.time() - start
        t_proj = t2 + t4 + t6 + t8 - (t1 + t3 + t5 + t7)
        _logger.debug(ind_str + 'done in {:.2f} sec, (inkl. {:.2f} sec for proj calc)'.format(t_func, t_proj))
        return ((t_l_scale, None, t_l_int_a, t_l_pha_a), (f_l_scale, f_l_filt, None, f_l_int_a))  # f_l_int_b,t_l_int_b,

    else:

        dfl_pad_z(dfl, padn)
        dfl.fft_z(method=fft_method, nthread=multiprocessing.cpu_count())
        dfl_trf(dfl, trf, mode='tr', dump_proj=dump_proj)
        dfl.fft_z(method=fft_method, nthread=multiprocessing.cpu_count())
        if st_cpl:
            dfl_st_cpl(dfl, trf.thetaB)
        dfl_shift_z(dfl, s_delay, set_zeros=0)
        dfl_pad_z(dfl, -padn)

        t_func = time.time() - start
        _logger.debug(ind_str + 'done in {:.2f} sec'.format(t_func))
        return ()


def save_xhrss_dump_proj(dump_proj, filePath):
    # saves the dfl_hxrss_filt radiation projections dump to text files

    (t_l_scale, _, t_l_int_a, t_l_pha_a), (f_l_scale, f_l_filt, _, f_l_int_a) = dump_proj

    f = open(filePath + '.t.txt', 'wb')
    header = 'Distance Power Phase'
    np.savetxt(f, np.c_[t_l_scale, t_l_int_a, t_l_pha_a], header=header, fmt="%e", newline='\n', comments='')
    f.close()

    f = open(filePath + '.f.txt', 'wb')
    header = 'Wavelength Spectrum Filter_Abs Filter_Ang'
    np.savetxt(f, np.c_[f_l_scale, f_l_int_a, np.abs(f_l_filt), np.unwrap(np.angle(f_l_filt))], header=header, fmt="%.8e", newline='\n', comments='')
    f.close()
    
    
def tap_exp(n, n0, a0, a1, a2):
    '''
    exponential tapering
    '''
    if n <= n0:
        return a0
    
    return a0 * (  1 + a1 * (n - n0)**a2 )

def tap_pol_old(n, n0, a0, a1, a2):
    '''
    piecewise-quadratic tapering function
    '''
    # if n0.__class__ is int:
        # n0 = [0, n0]
    # if a0.__class__ is int:
        # a0 = [a0]
    # if a1.__class__ is int:
        # a1 = [a1]
    # if a2.__class__ is int:
        # a2 = [a2]
        
    for i in range(1,len(n0)):
        if n < n0[i]:
            return a0 + (n-n0[i-1])*a1[i-1] + (n-n0[i-1])**2 * a2[i-1]
        a0 += (n0[i]-n0[i-1])*a1[i-1] + (n0[i]-n0[i-1])**2 * a2[i-1]
    
    return 1.0

def tap_pol(n, n0, a0, a1, a2):
    '''
    piecewise-quadratic tapering function
    '''
    if n < n0:
        return a0
    else:
        return a0 + (n-n0)*a1 + (n-n0)**2 * a2


def create_fel_lattice(und_N = 35,
                    und_L = 5,
                    und_l = 0.04,
                    und_Kx = 0,
                    und_Ky = 0,
                    inters_L = 1.08,
                    inters_K = 'K_und',
                    inters_phi=0,
                    quad_L = 0.4,
                    quad_K = 0,
                    phs_L = 0.0,
                    quad_start = 'd',
                    **kwargs):
    if quad_L > inters_L:
        _logger.warning('Quarrupole cannot be longer than intersection')

    # und_n = np.floor(und_L/und_l).astype(int)
    und_n = und_L/und_l

    und= Undulator(nperiods=und_n, lperiod=und_l, Kx=und_Kx, Ky=und_Ky, eid = "und")
    qf = Quadrupole (l=quad_L, eid = "qf")
    qd = Quadrupole (l=quad_L, eid = "qd")
    qfh = Quadrupole (l=qf.l / 2.)
    qdh = Quadrupole (l=qd.l / 2.)


    phs = UnknownElement(l=0) #phase shifter (defines expected retardation of electrons in intersections between undulators)
    phs.phi = inters_phi
    phs.K = inters_K #overrides phi, would be K of free space, identical to rms K_und if "K_und"

    d1 = Drift(l=(inters_L - quad_L) / 2, eid = "d1") #drift
    d2 = Drift(l=(inters_L - quad_L) / 2, eid = "d2")
    
    if und_N < 2:
        cell_N = 0
        cell_N_last = 0
    else:
        cell_N = np.floor((und_N - 1)/2).astype(int)
        cell_N_last = int((und_N - 1)/2%1)

    if quad_start == 'd':
        cell = (und, d1, qf, phs, d2, und, d1, qd, phs, d2) 
        extra_fodo = (und, d2, qdh)
        lat = (und, d1, qd, phs, d2) + cell_N * cell + cell_N_last * (und,)
    elif quad_start == 'f':
        cell = (und, d1, qd, phs, d2, und, d1, qf, phs, d2) 
        extra_fodo = (und, d2, qfh)
        lat = (und, d1, qf, phs, d2) + cell_N * cell + cell_N_last * (und,)

    return (MagneticLattice(lat), extra_fodo, cell)

def create_fel_lattice_tmp(und_N = 34,
                    und_L = 5,
                    und_l = 0.04,
                    und_Kx = 0,
                    und_Ky = 0,
                    inters_L = 1.08,
                    inters_K = 'K_und',
                    inters_phi=0,
                    quad_L = 0.1,
                    quad_K = 0,
                    phs_L = 0.0,
                    quad_start = 'd',
                    **kwargs):
    if quad_L > inters_L:
        _logger.warning('Quadrupole cannot be longer than intersection')

    # und_n = np.floor(und_L/und_l).astype(int)
    und_n = und_L/und_l

    und= Undulator(nperiods=und_n, lperiod=und_l, Kx=und_Kx, Ky=und_Ky, eid = "und")
    qf = Quadrupole (l=quad_L / 2, eid = "qf", k1=1)
    qd = Quadrupole (l=quad_L / 2, eid = "qd", k1=-1)
    qfh = Quadrupole (l=qf.l / 2., k1=1)
    qdh = Quadrupole (l=qd.l / 2., k1=-1)


    phs = UnknownElement(l=0) #phase shifter (defines expected retardation of electrons in intersections between undulators)
    phs.phi = inters_phi
    phs.K = inters_K #overrides phi, would be K of free space, identical to rms K_und if "K_und"

    d1 = Drift(l=(inters_L - quad_L) / 2, eid = "d1") #drift
    d2 = Drift(l=(inters_L - quad_L) / 2, eid = "d2")
    
    if und_N < 2:
        cell_N = 0
        cell_N_last = 0
    else:
        cell_N = np.floor((und_N)/2).astype(int)
        cell_N_last = int((und_N)/2%1)

    if quad_start == 'd':
        cell = (qd, phs, d2, und, d1, qf, qf, phs, d2, und, d1, qd) 
        # extra_fodo = (und, d2, qfh)
        lat = cell_N * cell + cell_N_last * (und,)
    elif quad_start == 'f':
        cell = (qf, phs, d2, und, d1, qd, qd, phs, d2, und, d1, qf) 
        # extra_fodo = (und, d2, qdh)
        lat = cell_N * cell + cell_N_last * (und,)

    return (MagneticLattice(lat), None, cell)


def create_fel_beamline(beamline = 'sase1', inters_phi=0, inters_K = "K_und"):
    if beamline in ['sase1', 1, 'sase2', 2, 'EuXFEL_SASE1', 'EuXFEL_SASE2']:
        return create_fel_lattice(und_N = 35,
                        und_L = 5-0.04,
                        und_l = 0.04,
                        und_Kx = 0,
                        und_Ky = 0,
                        inters_L = 1.08+0.04,
                        inters_K = inters_K,
                        inters_phi=inters_phi,
                        quad_L = 0.4,
                        quad_K = 0,
                        phs_L = 0.0,
                        quad_start = 'd',
                            )
    elif beamline in ['sase3', 3, 'EuXFEL_SASE3']:
        return create_fel_lattice(und_N = 21,
                        und_L = 5.032,#5
                        und_l = 0.068,
                        und_Kx = 0,
                        und_Ky = 0,
                        inters_L = 1.088,#1.08,
                        inters_K = inters_K,
                        inters_phi=inters_phi,
                        quad_L = 0.408,
                        quad_K = 0,
                        phs_L = 0.0,
                        quad_start = 'd',
                            )
    else:
        raise ValueError('Unknown beamline')

def create_exfel_lattice(*agrs, **kwargs):
    _logger.warning('will be deprecated, use create_fel_beamline instead')
    return create_fel_beamline(*agrs, **kwargs)

def prepare_el_optics(beam, lat_pkg, E_photon=None, beta_av=None, s=None):
    from ocelot.rad.undulator_params import Ephoton2K
    # if s is None:
        # jj = beam.I / (beam.beta_x * beam.beta_y * beam.emit_x * beam.emit_y)
        # jj = beam.I
        # s = beam.s[jj.argmax()]
    
    if s is None:
        beam_match = get_beam_peak(beam)
    else:
        beam_match = beam.get_s(s)
    
    
    lat, extra_fodo, cell = lat_pkg
    indx_und = np.where([i.__class__ == Undulator for i in lat.sequence])[0]
    und = lat.sequence[indx_und[0]]
    l_fodo= MagneticLattice(cell).totalLen / 2
    
    if E_photon is not None:
        und.Kx = Ephoton2K(E_photon, und.lperiod, beam_match.E)
        if np.isnan(und.Kx):
            gamma = beam_match.E/m_e_GeV
            Emax=2*gamma**2/und.lperiod*h_eV_s*speed_of_light
            _logger.error('requested photon energy {}eV is beyond reach for given lperiod {}m and beam energy {}GeV\nmax energy for given parameters is {}'.format(E_photon, und.lperiod, beam_match.E, Emax))
            raise ValueError('requested photon energy is beyond reach for given parameters')
    
    if beta_av is None:
        l_period = und.lperiod
        K_peak = np.max([und.Kx, und.Ky])
        if und.Kx != und.Ky:
            iwityp = 0 # planar undulator
        elif und.Kx == und.Ky:
            iwityp = 1 # helical undulator
        else:
            raise ValueError('unknown undulator: neither planar nor helical, estimation method not applicable')
        fel = beam2fel(beam_match, l_period, K_peak, iwityp=iwityp)
        beta_av = fel.beta_opt(apply=0, method=fel.method)
        if beta_av <= l_fodo * 2:
            beta_av = l_fodo * 2 * 1.001
        raise_min_beta = True
        _logger.info('optimal minimal average beta function was estimated and set to {:.3f}m'.format(beta_av))
    else:
        raise_min_beta = False
        

    rematch_beam_lat(beam_match, lat_pkg, beta_av, raise_min_beta=raise_min_beta)
    transform_beam_twiss(beam, Twiss(beam_match), s=s)



'''
legacy
'''


def detune_k(lat, sig):
    lat2 = deepcopy(lat)
    n = 0
    for i in range(len(lat2.sequence)):
        # print (lat2.sequence[i].__class__)
        if lat2.sequence[i].__class__ == Undulator:
            lat2.sequence[i] = deepcopy(lat.sequence[i])
            lat2.sequence[i].Kx = lat2.sequence[i].Kx * (1 + np.random.randn() * sig)
            n += 1

    return lat2


def detune_E(inp, beam, sig):
    # energy modulation
    inp.gamma0 = (beam.E + np.random.randn() * sig) / 0.000510998
    beam.emit_x = beam.emit_xn / inp.gamma0
    beam.emit_y = beam.emit_yn / inp.gamma0
    inp.rxbeam = np.sqrt(beam.emit_x * beam.beta_x)
    inp.rybeam = np.sqrt(beam.emit_y * beam.beta_y)


def taper(lat, k):
    lat2 = deepcopy(lat)
    n = 0
    for i in range(len(lat2.sequence)):
        if lat2.sequence[i].__class__ == Undulator:
            # print lat2.sequence[i].id, lat2.sequence[i].Kx
            lat2.sequence[i] = deepcopy(lat.sequence[i])
            # MOD BY GG. #lat2.sequence[i].Kx = lat2.sequence[i].Kx * k(n+1)
            lat2.sequence[i].Kx = k(n + 1)  # /np.sqrt(0.5) ##MOD BY GG.
            n += 1

    return lat2


def update_beam(beam_new, g, n_interp):
    '''
    check and rewrite!
    '''
    beam = deepcopy(beam_new)
    # g0 = np.array(map(lambda x : g.sliceValues[x]['energy'][-1], range(1,g.nSlices+1)) )
    # dg = np.array(map(lambda x : g.sliceValues[x]['e-spread'][-1], range(1,g.nSlices+1)) )
    g0 = g.el_energy[:, -1]  # * (0.511e-3)
    dg = g.el_e_spread[:, -1]

    print (len(g0))
    print (g.nSlices)

    print (len(beam_new.z))

    I = np.array(g.I)

    if n_interp == 0:
        n_interp = g.nSlices

    beam_new.z = np.linspace(beam.z[0], beam.z[-1], n_interp)
    z2 = np.linspace(beam.z[0], beam.z[-1], g.nSlices)
    beam_new.I = np.interp(beam_new.z, beam.z, beam.I)

    zmax, Imax = peaks(beam_new.z, beam_new.I, n=1)
    beam_new.idx_max = np.where(beam_new.z == zmax)[0][0]

    beam_new.ex = np.interp(beam_new.z, beam.z, beam.ex)
    beam_new.ey = np.interp(beam_new.z, beam.z, beam.ey)
    beam_new.zsep = beam.zsep * len(beam.z) / len(beam_new.z)
    #beam_new.g0 = np.interp(beam_new.z, beam.z, beam.g0)
    # print ("_______________________________")
    # print (g0)
    # print(beam.E)
    # print(beam.E/(0.511e-3))
    # print ("_______________________________")
    beam_new.g0 = g0 + beam.E / (0.511e-3)  # potential problem here, no beam.gamma_rel
    print (len(beam_new.z))
    print (len(beam_new.g0))
    print (len(beam.z))
    beam_new.g0 = np.interp(beam_new.z, z2, beam_new.g0)
    beam_new.dg = dg
    beam_new.dg = np.interp(beam_new.z, z2, beam_new.dg)

    beam_new.eloss = np.interp(beam_new.z, beam.z, beam.eloss)

    beam_new.betax = np.interp(beam_new.z, beam.z, beam.betax)
    beam_new.betay = np.interp(beam_new.z, beam.z, beam.betay)
    beam_new.alphax = np.interp(beam_new.z, beam.z, beam.alphax)
    beam_new.alphay = np.interp(beam_new.z, beam.z, beam.alphay)

    beam_new.x = np.interp(beam_new.z, beam.z, beam.x)
    beam_new.px = np.interp(beam_new.z, beam.z, beam.px)
    beam_new.y = np.interp(beam_new.z, beam.z, beam.y)
    beam_new.py = np.interp(beam_new.z, beam.z, beam.py)


#old one
def rematch(beta_mean, l_fodo, qdh, lat, extra_fodo, beam, qf, qd):
    '''
    requires l_fodo to be defined in the lattice
    '''

    k, betaMin, betaMax, __ = fodo_parameters(betaXmean=beta_mean, L=l_fodo, verbose=True)

    k1 = k[0] / qdh.l

    tw0 = Twiss(beam)
    
    _logger.debug('before rematching k=%f %f   beta=%f %f alpha=%f %f' % (qf.k1, qd.k1, tw0.beta_x, tw0.beta_y, tw0.alpha_x, tw0.alpha_y))
    # print('before rematching k=%f %f   beta=%f %f alpha=%f %f' % (qf.k1, qd.k1, tw0.beta_x, tw0.beta_y, tw0.alpha_x, tw0.alpha_y))

    extra = MagneticLattice(extra_fodo)
    tws = twiss(extra, tw0)
    tw2 = tws[-1]

    tw2m = Twiss(tw2)
    tw2m.beta_x = betaMin[0]
    tw2m.beta_y = betaMax[0]
    tw2m.alpha_x = 0.0
    tw2m.alpha_y = 0.0
    tw2m.gamma_x = (1 + tw2m.alpha_x * tw2m.alpha_x) / tw2m.beta_x
    tw2m.gamma_y = (1 + tw2m.alpha_y * tw2m.alpha_y) / tw2m.beta_y

    #k1 += 0.5

    qf.k1 = k1
    qd.k1 = -k1
    qdh.k1 = -k1

    lat.update_transfer_maps()
    extra.update_transfer_maps()

    R1 = lattice_transfer_map(extra, beam.E)
    Rinv = np.linalg.inv(R1)

    m1 = TransferMap()

    m1.R = lambda e: Rinv

    tw0m = m1.map_x_twiss(tw2m)
    _logger.debug('after rematching k=%f %f   beta=%f %f alpha=%f %f' % (qf.k1, qd.k1, tw0m.beta_x, tw0m.beta_y, tw0m.alpha_x, tw0m.alpha_y))
    # print ('after rematching k=%f %f   beta=%f %f alpha=%f %f' % (qf.k1, qd.k1, tw0m.beta_x, tw0m.beta_y, tw0m.alpha_x, tw0m.alpha_y))

    beam.beta_x, beam.alpha_x = tw0m.beta_x, tw0m.alpha_x
    beam.beta_y, beam.alpha_y = tw0m.beta_y, tw0m.alpha_y


def rematch_beam_lat(beam, lat_pkg, beta_mean, raise_min_beta=False):
    
    lat, extra_fodo, cell = lat_pkg
    l_fodo= MagneticLattice(cell).totalLen / 2
    
    if beta_mean <= 2 * l_fodo:
        if raise_min_beta:
            _logger.warning('Desired average beta function {}m was not larger that cell length and was set to that length {}m'.format(beta_mean, 2*l_fodo*1.001))
            beta_mean = 2 * l_fodo * 1.001
        else:
            _logger.error('Desired average beta function {}m should be larger that cell length {}m'.format(beta_mean, 2*l_fodo))
            raise ValueError('Desired average beta function {}m should be larger that cell length {}m'.format(beta_mean, 2*l_fodo))

    indx_q = np.where([i.__class__ == Quadrupole for i in lat.sequence])[0]

    qd = lat.sequence[indx_q[0]]
    qf = lat.sequence[indx_q[1]]
    
    indx_qh = np.where([i.__class__ == Quadrupole for i in extra_fodo])[0]
    qdh = extra_fodo[indx_qh[-1]]
    '''
    requires l_fodo to be defined in the lattice
    '''

    k, betaMin, betaMax, __ = fodo_parameters(betaXmean=beta_mean, L=l_fodo, verbose=False)

    k1 = k[0] / qdh.l

    tw0 = Twiss(beam)

    _logger.debug('before rematching k=%f %f   beta=%f %f alpha=%f %f' % (qf.k1, qd.k1, tw0.beta_x, tw0.beta_y, tw0.alpha_x, tw0.alpha_y))
    # print('before rematching k=%f %f   beta=%f %f alpha=%f %f' % (qf.k1, qd.k1, tw0.beta_x, tw0.beta_y, tw0.alpha_x, tw0.alpha_y))

    extra = MagneticLattice(extra_fodo)
    tws = twiss(extra, tw0)
    tw2 = tws[-1]

    tw2m = Twiss(tw2)
    tw2m.beta_x = betaMin[0]
    tw2m.beta_y = betaMax[0]
    tw2m.alpha_x = 0.0
    tw2m.alpha_y = 0.0
    tw2m.gamma_x = (1 + tw2m.alpha_x * tw2m.alpha_x) / tw2m.beta_x
    tw2m.gamma_y = (1 + tw2m.alpha_y * tw2m.alpha_y) / tw2m.beta_y

    #k1 += 0.5

    qf.k1 = k1
    qd.k1 = -k1
    qdh.k1 = -k1

    lat.update_transfer_maps()
    extra.update_transfer_maps()

    R1 = lattice_transfer_map(extra, beam.E)
    Rinv = np.linalg.inv(R1)

    m1 = TransferMap()

    m1.R = lambda e: Rinv

    tw0m = m1.map_x_twiss(tw2m)
    _logger.debug('after rematching k=%f %f   beta=%f %f alpha=%f %f' % (qf.k1, qd.k1, tw0m.beta_x, tw0m.beta_y, tw0m.alpha_x, tw0m.alpha_y))
    # print ('after rematching k=%f %f   beta=%f %f alpha=%f %f' % (qf.k1, qd.k1, tw0m.beta_x, tw0m.beta_y, tw0m.alpha_x, tw0m.alpha_y))

    beam.beta_x, beam.alpha_x = tw0m.beta_x, tw0m.alpha_x
    beam.beta_y, beam.alpha_y = tw0m.beta_y, tw0m.alpha_y


'''
CHEDULED FOR REMOVAL
'''


def get_data_dir():
    host = socket.gethostname()

    if host.startswith('it-hpc'):
        return '/data/netapp/xfel/iagapov/xcode_data/'
    return '/tmp/'


def checkout_run(run_dir, run_id, prefix1, prefix2, save=False, debug=1):
    print ('    checking out run from ' + prefix1 + '.gout to ' + prefix2 + '.gout')
    old_file = run_dir + 'run.' + str(run_id) + prefix1 + '.gout'
    new_file = run_dir + 'run.' + str(run_id) + prefix2 + '.gout'

    if save:
        os.system('cp ' + old_file + ' ' + new_file)
        os.system('cp ' + old_file + '.dfl ' + new_file + '.dfl 2>/dev/null')  # 2>/dev/null to supress error messages if no such file
        os.system('cp ' + old_file + '.dpa ' + new_file + '.dpa 2>/dev/null')
        os.system('cp ' + old_file + '.beam ' + new_file + '.beam 2>/dev/null')
        os.system('cp ' + run_dir + 'tmp.gen' + ' ' + run_dir + 'geninp.' + str(run_id) + prefix2 + '.inp 2>/dev/null')
        os.system('cp ' + run_dir + 'lattice.inp' + ' ' + run_dir + 'lattice.' + str(run_id) + prefix2 + '.inp 2>/dev/null')
    else:
        if debug > 0:
            print ('      moving *.out file')
        os.system('mv ' + old_file + ' ' + new_file)
        if debug > 0:
            print ('      moving *.dfl file')
        os.system('mv ' + old_file + '.dfl ' + new_file + '.dfl 2>/dev/null')  # 2>/dev/null to supress error messages if no such file
        if debug > 0:
            print ('      moving *.dpa file')
        os.system('mv ' + old_file + '.dpa ' + new_file + '.dpa 2>/dev/null')
        if debug > 0:
            print ('      moving *.beam file')
        os.system('mv ' + old_file + '.beam ' + new_file + '.beam 2>/dev/null')
        if debug > 0:
            print ('      moving input files')
        os.system('mv ' + run_dir + 'tmp.gen' + ' ' + run_dir + 'geninp.' + str(run_id) + prefix2 + '.inp 2>/dev/null')
        os.system('mv ' + run_dir + 'lattice.inp' + ' ' + run_dir + 'lattice.' + str(run_id) + prefix2 + '.inp 2>/dev/null')
    print ('      done')
    # os.system('rm ' + run_dir + 'run.' +str(run_id) + '.gout*')


class FelSimulator(object):
    '''
    configurable to e.g. semi-empirical models
    '''

    def __init__(self):
        self.engine = 'genesis'

    def run(self):
        if self.engine == 'test_1d':
            w1 = read_signal(file_name=self.input, npad=self.npad, E_ref=self.E_ev)
            return w1, None
        if self.engine == 'test_3d':
            ''' produced  sliced field '''
            w1 = read_signal(file_name=self.input, npad=self.npad, E_ref=self.E_ev)
            s3d = Signal3D()
            s3d.fs = [w1, deepcopy(w1), deepcopy(w1), deepcopy(w1)]
            s3d.mesh_size = (2, 2)
            return s3d, None
        if self.engine == 'test_genesis':
            ''' read test sliced field '''
            g = read_out_file(self.input)
            print ('read sliced field ', g('ncar'), g.nSlices)
            slices = readRadiationFile(fileName=self.input + '.dfl', npoints=g('ncar'))
            s3d = Signal3D()
            s3d.slices = slices
            s3d.mesh_size = (int(g('ncar')), int(g('ncar')))
            s3d.g = g
            return s3d, None
