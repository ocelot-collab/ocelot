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
from ocelot import *
from ocelot.common.math_op import *

# from ocelot.optics.utils import *
# from ocelot.rad.undulator_params import *
# from ocelot.rad.fel import *
from ocelot.adaptors.genesis import *

import multiprocessing
import pyfftw
nthread = multiprocessing.cpu_count()


def background(command):
    '''
    start command as background process
    the argument shohuld preferably be a string in triple quotes '
    '''
    import subprocess
    imports = 'from ocelot.adaptors.genesis import *; from ocelot.gui.genesis_plot import *; from ocelot.utils.xfel_utils import *; '
    subprocess.Popen(["python3", "-c", imports + command])


def copy_this_script(scriptName, scriptPath, folderPath):
    cmd = 'cp ' + scriptPath + ' ' + folderPath + 'exec_' + scriptName + ' '
    os.system(cmd)

'''
SELF-SEEDING - relevant
'''


def dfl_prop(dfl, z, fine=1, debug=1):
    '''
    Fourier propagator for fieldfile

    can handle wide spectrum
      (every slice in freq.domain is propagated 
       according to its frequency)
    no kx**2+ky**2<<k0**2 limitation

    dfl is the RadiationField() object
    z is the propagation distance in [m] 
    fine==0 is a flag for ~2x faster propagation. 
        no Fourier transform to frequency domain is done
        assumes no angular dispersion (true for plain FEL radiation)
        assumes narrow spectrum at center of xlamds (true for plain FEL radiation)

    returns RadiationField() object

    z>0 ==> forward
    '''
    if debug > 0:
        print('    propagating dfl file by %.2f meters' % (z))

    start = time.time()

    dfl_out = deepcopy(dfl)
    domain_xy = dfl.domain_xy
    domain_z = dfl.domain_z

    # switch to inv-space/freq domain
    if dfl_out.domain_xy == 's':
        dfl_out = dfl_fft_xy(dfl_out, debug=debug)
    if dfl_out.domain_z == 't' and fine:
        dfl_out = dfl_fft_z(dfl_out, debug=debug)

    if fine:
        k_x, k_y = np.meshgrid(dfl_out.scale_kx(), dfl_out.scale_ky())
        for i in range(dfl_out.Nz()):
            k = dfl_out.scale_kz()[i]
            H = exp(1j * z * (sqrt(k**2 - k_x**2 - k_y**2) - k))
            dfl_out.fld[i, :, :] *= H
    else:
        k_x, k_y = np.meshgrid(dfl_out.scale_kx(), dfl_out.scale_ky())
        k = 2 * np.pi / dfl_out.xlamds
        H = exp(1j * z * (sqrt(k**2 - k_x**2 - k_y**2) - k))
        for i in range(dfl_out.Nz()):
            dfl_out.fld[i, :, :] *= H

    # switch to original domain
    if domain_xy == 's':
        dfl_out = dfl_fft_xy(dfl_out, debug=debug)
    if domain_z == 't' and fine:
        dfl_out = dfl_fft_z(dfl_out, debug=debug)

    t_func = time.time() - start
    if debug > 0:
        print('      done in %.2f ' % t_func + 'sec')

    return dfl_out


def dfl_waistscan(dfl, z_pos, projection=0, debug=1):
    '''
    propagates the RadaitionField object dfl 
    through the sequence of positions z_pos
    and calculates transverse distribution parameters
    such as peak photon density and sizes in both dimentions

    if projection==1, then size of projection is calculated
        otherwise - size across the central line passing through the mesh center
    '''
    if debug > 0:
        print('    scanning dfl waist in range %s meters' % (z_pos))
    start = time.time()

    sc_res = WaistScanResults()
    sc_res.xlamds = dfl.xlamds
    sc_res.filePath = dfl.filePath

    for z in z_pos:

        if debug > 0:
            print('      scanning at z = %.2f m' % (z))

        I_xy = dfl_prop(dfl, z, fine=0, debug=0).int_xy()  # integrated xy intensity

        scale_x = dfl.scale_x()
        scale_y = dfl.scale_y()
        center_x = (shape(I_xy)[0] + 1) / 2
        center_y = (shape(I_xy)[1] + 1) / 2

        if projection:
            I_x = np.sum(I_xy, axis=1)
            I_y = np.sum(I_xy, axis=0)
        else:
            I_x = I_xy[:, center_y]
            I_y = I_xy[center_x, :]

        sc_res.z_pos = np.append(sc_res.z_pos, z)
        sc_res.phdens_max = np.append(sc_res.phdens_max, np.amax(I_xy))
        sc_res.phdens_onaxis = np.append(sc_res.phdens_onaxis, I_xy[center_x, center_y])
        sc_res.fwhm_x = np.append(sc_res.fwhm_x, fwhm(scale_x, I_x))
        sc_res.fwhm_y = np.append(sc_res.fwhm_y, fwhm(scale_y, I_y))
        sc_res.std_x = np.append(sc_res.std_x, std_moment(scale_x, I_x))
        sc_res.std_y = np.append(sc_res.std_y, std_moment(scale_y, I_y))

        sc_res.z_max_phdens = sc_res.z_pos[np.argmax(sc_res.phdens_max)]

    t_func = time.time() - start
    if debug > 0:
        print('      done in %.2f ' % t_func + 'sec')

    return sc_res


def dfl_interp(dfl, interpN=(1, 1), interpL=(1, 1), newN=(None, None), newL=(None, None), method='cubic', debug=1):
    ''' 
    2d interpolation of the coherent radiation distribution 
    interpN and interpL define the desired interpolation coefficients for  
    transverse point __density__ and transverse mesh __size__ correspondingly 
    newN and newL define the final desire number of points and size of the mesh 
    when newN and newL are not None interpN and interpL values are ignored 
    coordinate convention is (x,y) 
    '''
    from scipy.interpolate import interp2d

    if debug > 0:
        print ('    interpolating radiation file')
    start_time = time.time()

    # in case if interpolation is the same in toth dimentions
    if np.size(interpN) == 1:
        interpN = (interpN, interpN)
    if np.size(interpL) == 1:
        interpL = (interpL, interpL)
    if np.size(newN) == 1:
        newN = (newN, newN)
    if np.size(newL) == 1:
        newL = (newL, newL)

    if debug > 1:
        print('newL=', newL)
        print('newN=', newN)

    if (interpN == (1, 1) and interpL == (1, 1) and newN == (None, None) and newL == (None, None)) or \
       (interpN == (1, 1) and interpL == (1, 1) and newN == (dfl.Nx(), dfl.Ny()) and newL == (dfl.Lx(), dfl.Ly())):
        print('      skip (no interpolation required)')
        return dfl

    # calculate new mesh parameters only if not defined explicvitly
    if newN == (None, None) and newL == (None, None):
        interpNx = interpN[0]
        interpNy = interpN[1]
        interpLx = interpL[0]
        interpLy = interpL[1]

        if interpNx == 0 or interpLx == 0 or interpNy == 0 or interpLy == 0:
            print('interpolation values cannot be 0')
            return None
            # place exception
        elif interpNx == 1 and interpLx == 1 and interpNy == 1 and interpLy == 1:
            return dfl
            print('no interpolation required, returning original')
        else:
            Nx2 = int(dfl.Nx() * interpNx * interpLx)
            if Nx2 % 2 == 0 and Nx2 > dfl.Nx():
                Nx2 -= 1
            if Nx2 % 2 == 0 and Nx2 < dfl.Nx():
                Nx2 += 1
            Ny2 = int(dfl.Ny() * interpNy * interpLy)
            if Ny2 % 2 == 0 and Ny2 > dfl.Ny():
                Ny2 -= 1
            if Ny2 % 2 == 0 and Ny2 < dfl.Ny():
                Ny2 += 1

            Lx2 = dfl.Lx() * interpLx
            Ly2 = dfl.Ly() * interpLy

    else:
        # redo to maintain mesh density
        if newN[0] != None:
            Nx2 = newN[0]
        else:
            Nx2 = dfl.Nx()

        if newN[1] != None:
            Ny2 = newN[1]
        else:
            Ny2 = dfl.Ny()

        if newL[0] != None:
            Lx2 = newL[0]
        else:
            Lx2 = dfl.Lx()

        if newL[1] != None:
            Ly2 = newL[1]
        else:
            Ly2 = dfl.Ly()

    xscale1 = np.linspace(-dfl.Lx() / 2, dfl.Lx() / 2, dfl.Nx())
    yscale1 = np.linspace(-dfl.Ly() / 2, dfl.Ly() / 2, dfl.Ny())
    xscale2 = np.linspace(-Lx2 / 2, Lx2 / 2, Nx2)
    yscale2 = np.linspace(-Ly2 / 2, Ly2 / 2, Ny2)

    ix_min = np.where(xscale1 >= xscale2[0])[0][0]
    ix_max = np.where(xscale1 <= xscale2[-1])[-1][-1]
    iy_min = np.where(yscale1 >= yscale2[0])[0][0]
    iy_max = np.where(yscale1 <= yscale2[-1])[-1][-1]
    if debug > 1:
        print('      energy before interpolation ' + str(dfl.E()))
    #interp_func = rgi((zscale1,yscale1,xscale1), dfl.fld, fill_value=0, bounds_error=False, method='nearest')
    fld2 = []
    for fslice in dfl.fld:
        re_func = interp2d(xscale1, yscale1, np.real(fslice), fill_value=0, bounds_error=False, kind=method)
        im_func = interp2d(xscale1, yscale1, np.imag(fslice), fill_value=0, bounds_error=False, kind=method)
        fslice2 = re_func(xscale2, yscale2) + 1j * im_func(xscale2, yscale2)
        P1 = sum(abs(fslice[iy_min:iy_max, ix_min:ix_max])**2)
        P2 = sum(abs(fslice2)**2)
        fslice2 = fslice2 * sqrt(P1 / P2)
        fld2.append(fslice2)

    dfl2 = deepcopy(dfl)
    # dfl2=RadiationField()
    dfl2.fld = np.array(fld2)
    dfl2.dx = Lx2 / dfl2.Nx()
    dfl2.dy = Ly2 / dfl2.Ny()
    # dfl2.fileName=dfl.fileName+'i'
    # dfl2.filePath=dfl.filePath+'i'
    if debug > 1:
        print('      energy after interpolation ' + str(F2.E()))
    if debug > 0:
        print('      done in %.2f sec' % (time.time() - start_time))

    return dfl2


def dfl_shift_z(dfl, s, set_zeros=1):
    '''
    shift the radiation within the window in time domain
    '''
    # set_zeros - to set the values from out of initial time window to zeros
    assert dfl.domain_z == 't', 'dfl_shift_z works only in time domain!'
    shift_n = int(s / dfl.dz)
    print('    shifting dfl by %.2f um (%.0f slices)' % (s * 1e6, shift_n))
    start = time.time()
    dfl.fld = np.roll(dfl.fld, shift_n, axis=0)
    if set_zeros:
        if shift_n > 0:
            dfl.fld[:shift_n, :, :] = 0
        if shift_n < 0:
            dfl.fld[shift_n:, :, :] = 0

    t_func = time.time() - start
    print('      done in %.2f ' % t_func + 'sec')
    return dfl


def dfl_pad_z(dfl, padn):
    assert np.mod(padn, 1) == 0, 'pad should be integer'
    start = time.time()

    if padn > 1:
        print('    padding dfl by ' + str(padn))
        padn_n = (padn - 1) / 2 * dfl.Nz()  # number of slices to add before and after
        dfl_pad = RadiationField((dfl.Nz() + 2 * padn_n, dfl.Ny(), dfl.Nx()))
        dfl_pad.copy_param(dfl)
        dfl_pad.fld[padn_n:-padn_n, :, :] = dfl.fld
    elif padn < -1:
        padn = abs(padn)
        print('    de-padding dfl by ' + str(padn))
        padn_n = dfl.Nz() / padn * ((padn - 1) / 2)
        dfl_pad = RadiationField()
        dfl_pad.copy_param(dfl)
        dfl_pad.fld = dfl.fld[padn_n:-padn_n, :, :]
    else:
        print('    padding dfl by ' + str(padn))
        print('      pass')
        return dfl

    t_func = time.time() - start
    if t_func < 60:
        print('      done in %.2f ' % t_func + 'sec')
    else:
        print('      done in %.2f ' % t_func / 60 + 'min')
    return dfl_pad


def dfl_fft_z(dfl, method='mp', nthread=multiprocessing.cpu_count(), debug=1):  # move to somewhere else
    if debug > 0:
        print('      calculating fft_z from ' + dfl.domain_z + ' domain with ' + method)
    start = time.time()
    dfl_fft = RadiationField(dfl.shape())
    dfl_fft.copy_param(dfl)

    if nthread < 2:
        method = 'np'

    if dfl.domain_z == 't':
        if method == 'np':
            dfl_fft.fld = np.fft.fft(dfl.fld, axis=0)
        elif method == 'mp':
            fft = pyfftw.builders.fft(dfl.fld, axis=0, overwrite_input=False, planner_effort='FFTW_ESTIMATE', threads=nthread, auto_align_input=False, auto_contiguous=False, avoid_copy=True)
            dfl_fft.fld = fft()
        else:
            raise ValueError('fft method should be "np" or "mp"')
        dfl_fft.fld = np.fft.ifftshift(dfl_fft.fld, 0)
        dfl_fft.fld /= sqrt(dfl_fft.Nz())
        dfl_fft.domain_z = 'f'
    elif dfl.domain_z == 'f':
        dfl_fft.fld = np.fft.fftshift(dfl.fld, 0)
        if method == 'np':
            dfl_fft.fld = np.fft.ifft(dfl_fft.fld, axis=0)
        elif method == 'mp':
            fft = pyfftw.builders.ifft(dfl_fft.fld, axis=0, overwrite_input=False, planner_effort='FFTW_ESTIMATE', threads=nthread, auto_align_input=False, auto_contiguous=False, avoid_copy=True)
            dfl_fft.fld = fft()
        else:
            raise ValueError("fft method should be 'np' or 'mp'")
        dfl_fft.fld *= sqrt(dfl_fft.Nz())
        dfl_fft.domain_z = 't'
    else:
        raise ValueError("domain_z value should be 't' or 'f'")

    if debug > 0:
        t_func = time.time() - start
        if t_func < 60:
            print('        done in %.2f sec' %(t_func))
        else:
            print('        done in %.2f min' %(t_func / 60))
    return dfl_fft


def dfl_fft_xy(dfl, method='mp', nthread=multiprocessing.cpu_count(), debug=1):  # move to somewhere else
    if debug > 0:
        print('      calculating fft_xy from ' + dfl.domain_xy + ' domain with ' + method)
    start = time.time()
    dfl_fft = RadiationField(dfl.shape())
    dfl_fft.copy_param(dfl)

    if nthread < 2:
        method = 'np'

    if dfl.domain_xy == 's':
        if method == 'np':
            dfl_fft.fld = np.fft.fft2(dfl.fld, axes=(1, 2))
        elif method == 'mp':
            fft = pyfftw.builders.fft2(dfl.fld, axes=(1, 2), overwrite_input=False, planner_effort='FFTW_ESTIMATE', threads=nthread, auto_align_input=False, auto_contiguous=False, avoid_copy=True)
            dfl_fft.fld = fft()
        else:
            raise ValueError("fft method should be 'np' or 'mp'")
        dfl_fft.fld = np.fft.fftshift(dfl_fft.fld, axes=(1, 2))
        dfl_fft.fld /= sqrt(dfl_fft.Nx() * dfl_fft.Ny())
        dfl_fft.domain_xy = 'k'
    elif dfl.domain_xy == 'k':
        dfl_fft.fld = np.fft.ifftshift(dfl.fld, axes=(1, 2))
        if method == 'np':
            dfl_fft.fld = np.fft.ifft2(dfl_fft.fld, axes=(1, 2))
        elif method == 'mp':
            fft = pyfftw.builders.ifft2(dfl_fft.fld, axes=(1, 2), overwrite_input=False, planner_effort='FFTW_ESTIMATE', threads=nthread, auto_align_input=False, auto_contiguous=False, avoid_copy=True)
            dfl_fft.fld = fft()
        else:
            raise ValueError("fft method should be 'np' or 'mp'")
        dfl_fft.fld *= sqrt(dfl_fft.Nx() * dfl_fft.Ny())
        dfl_fft.domain_xy = 's'

    else:
        raise ValueError("domain_xy value should be 's' or 'k'")

    if debug > 0:
        t_func = time.time() - start
        if t_func < 60:
            print('        done in %.2f ' % t_func + 'sec')
        else:
            print('        done in %.2f ' % t_func / 60 + 'min')
    return dfl_fft


def dfl_trf(dfl, trf, mode):
    dfl_out = deepcopy(dfl)
    assert dfl_out.domain_z == 'f', 'dfl_trf works only in frequency domain!'
    print('    multiplying dfl by trf')
    start = time.time()
    # assert trf.__class__==TransferFunction,'Wrong TransferFunction class'
    assert dfl_out.domain_z == 'f', 'wrong dfl domain (must be frequency)!'
    if mode == 'tr':
        filt = trf.tr
    elif mode == 'ref':
        filt = trf.ref
    else:
        raise AttributeError('Wrong z_domain attribute')
    filt_lamdscale = 2 * pi / trf.k
    if min(dfl_out.scale_z()) > max(filt_lamdscale) or max(dfl_out.scale_z()) < min(filt_lamdscale):
        raise ValueError('frequency scales of dfl and transfer function do not overlap')

    filt_interp_re = np.flipud(np.interp(np.flipud(dfl_out.scale_z()), np.flipud(filt_lamdscale), np.flipud(np.real(filt))))
    filt_interp_im = np.flipud(np.interp(np.flipud(dfl_out.scale_z()), np.flipud(filt_lamdscale), np.flipud(np.imag(filt))))
    filt_interp = filt_interp_re - 1j * filt_interp_im
    del filt_interp_re, filt_interp_im
    dfl_out.fld = dfl_out.fld * filt_interp[:, np.newaxis, np.newaxis]

    t_func = time.time() - start
    print('      done in %.2f ' % t_func + 'sec')
    return dfl_out, filt_interp


def dfl_st_cpl(dfl, theta_b, inp_axis='y', s_start=None):

    print('    introducing spatio-temporal coupling')
    start = time.time()
    direction = 1
    if s_start == None:
        s_start = n_moment(dfl.scale_z(), dfl.int_z(), 0, 1)

    dfl2 = deepcopy(dfl)
    shift_z_scale = dfl2.scale_z() - s_start
    shift_m_scale = shift_z_scale / tan(theta_b)
    shift_m_scale[np.where(shift_m_scale > 0)] = 0
    shift_pix_scale = np.floor(shift_m_scale / dfl2.dx).astype(int)

    # pix_start=np.where(shift_m_scale==0)[0][0]

    fld = dfl2.fld
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

    dfl2.fld = fld
    t_func = time.time() - start
    print('      done in %.2f ' % t_func + 'sec')
    return dfl2


def dfl_hxrss_filt(dfl, trf, ev_seed, s_delay, st_cpl=1, enforce_padn=None, res_per_fwhm=6, fft_method='mp', dump_proj=0, debug=1):
    # needs optimizing?
    # tmp
    import matplotlib.pyplot as plt

    nthread = multiprocessing.cpu_count()
    if nthread > 8:
        nthread = int(nthread * 0.9)  # not to occupy all CPUs on login server
    print('  HXRSS dfl filtering')
    start = time.time()
    # klpos, krpos, cwidth = FWHM(trf.k, 1.0-np.abs(trf.tr))
    cwidth = fwhm(trf.k, 1.0 - np.abs(trf.tr))
    dk_old = 2 * pi / dfl.Lz()
    dk = cwidth / res_per_fwhm
    padn = np.ceil(dk_old / dk).astype(int)
    if np.mod(padn, 2) == 0 and padn != 0:  # check for odd
        padn = int(padn + 1)
    
    if enforce_padn!=None:
        padn=enforce_padn
        
        
    if dump_proj:

        dfl = dfl_pad_z(dfl, padn)

        t1 = time.time()
        t_l_scale = dfl.scale_z()
        # t_l_int_b=dfl.int_z()
        t2 = time.time()

        dfl = dfl_fft_z(dfl, method=fft_method, nthread=multiprocessing.cpu_count())

        t3 = time.time()
        f_l_scale = dfl.scale_z()  # frequency_large_scale (wavelength in m)
        # f_l_int_b=dfl.int_z()
        t4 = time.time()

        dfl, f_l_filt = dfl_trf(dfl, trf, mode='tr')

        t5 = time.time()
        f_l_int_a = dfl.int_z()
        t6 = time.time()

        dfl = dfl_fft_z(dfl, method=fft_method, nthread=multiprocessing.cpu_count())

        t7 = time.time()
        t_l_int_a = dfl.int_z()
        t_l_pha_a = dfl.ang_z_onaxis()
        t8 = time.time()

        # if debug>2:###
            # plt.figure('before st-c')
            # plt.plot(dfl.scale_z(),dfl.ang_z_onaxis())
            # plt.show()

        if st_cpl:
            dfl = dfl_st_cpl(dfl, trf.thetaB)

        # if debug>2:###
           # plt.figure('after st-c')
           # plt.plot(dfl.scale_z(),dfl.ang_z_onaxis())
           # plt.show()

        dfl = dfl_shift_z(dfl, s_delay, set_zeros=0)

        dfl = dfl_pad_z(dfl, -padn)

        t_func = time.time() - start
        t_proj = t2 + t4 + t6 + t8 - (t1 + t3 + t5 + t7)
        print('    done in %.2f sec, (inkl. %.2f sec for proj calc)' % (t_func, t_proj))
        return dfl, ((t_l_scale, None, t_l_int_a, t_l_pha_a), (f_l_scale, f_l_filt, None, f_l_int_a))  # f_l_int_b,t_l_int_b,

    else:

        dfl = dfl_pad_z(dfl, padn)
        dfl = dfl_fft_z(dfl, method=fft_method, nthread=multiprocessing.cpu_count())
        dfl, _ = dfl_trf(dfl, trf, mode='tr')
        dfl = dfl_fft_z(dfl, method=fft_method, nthread=multiprocessing.cpu_count())
        if st_cpl:
            dfl = dfl_st_cpl(dfl, trf.thetaB)
        dfl = dfl_shift_z(dfl, s_delay, set_zeros=0)
        dfl = dfl_pad_z(dfl, -padn)

        t_func = time.time() - start
        print('    done in %.2f ' % t_func + 'sec')
        return dfl, ()


def save_xhrss_dump_proj(dump_proj, filePath):
    # saves the dfl_hxrss_filt radiation projections dump to text files

    (t_l_scale, _, t_l_int_a, t_l_pha_a), (f_l_scale, f_l_filt, _, f_l_int_a) = dump_proj

    f = open(filePath + '.t.txt', 'wb')
    header = 'Distance Power Phase'
    np.savetxt(f, np.c_[t_l_scale, t_l_int_a, t_l_pha_a], header=header, fmt="%e", newline='\n', comments='')
    f.close()

    f = open(filePath + '.f.txt', 'wb')
    header = 'Wavelength Spectrum Filter_Abs Filter_Ang'
    np.savetxt(f, np.c_[f_l_scale, f_l_int_a, np.abs(f_l_filt), np.angle(f_l_filt)], header=header, fmt="%e", newline='\n', comments='')
    f.close()


def save_trf(trf, attr, flePath):
    if hasattr(trf, attr):
        filt = getattr(trf, attr)
    else:
        raise ValueError('no attribute', attr, 'in fransfer function')

    f = open(flePath, 'wb')
    header = 'Energy[eV] Filter_Abs Filter_Ang'
    np.savetxt(f, np.c_[trf.ev(), np.abs(trf.tr), np.angle(trf.tr)], header=header, fmt="%e", newline='\n', comments='')
    f.close()

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


def rematch(beta_mean, l_fodo, qdh, lat, extra_fodo, beam, qf, qd):
    '''
    requires l_fodo to be defined in the lattice
    '''

    k, betaMin, betaMax, __ = fodo_parameters(betaXmean=beta_mean, L=l_fodo, verbose=True)

    k1 = k[0] / qdh.l

    tw0 = Twiss(beam)

    print('before rematching k=%f %f   beta=%f %f alpha=%f %f' % (qf.k1, qd.k1, tw0.beta_x, tw0.beta_y, tw0.alpha_x, tw0.alpha_y))

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
    print ('after rematching k=%f %f   beta=%f %f alpha=%f %f' % (qf.k1, qd.k1, tw0m.beta_x, tw0m.beta_y, tw0m.alpha_x, tw0m.alpha_y))

    beam.beta_x, beam.alpha_x = tw0m.beta_x, tw0m.alpha_x
    beam.beta_y, beam.alpha_y = tw0m.beta_y, tw0m.alpha_y


def rematch_beam_lat(beam, lat, extra_fodo, l_fodo, beta_mean):

    isquad = find([i.__class__ == Quadrupole for i in lat.sequence])
    qd = lat.sequence[isquad[0]]
    qf = lat.sequence[isquad[1]]
    qdh = deepcopy(qd)
    qdh.l /= 2
    '''
    requires l_fodo to be defined in the lattice
    '''

    k, betaMin, betaMax, __ = fodo_parameters(betaXmean=beta_mean, L=l_fodo, verbose=False)

    k1 = k[0] / qdh.l

    tw0 = Twiss(beam)

    print('before rematching k=%f %f   beta=%f %f alpha=%f %f' % (qf.k1, qd.k1, tw0.beta_x, tw0.beta_y, tw0.alpha_x, tw0.alpha_y))

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
    print ('after rematching k=%f %f   beta=%f %f alpha=%f %f' % (qf.k1, qd.k1, tw0m.beta_x, tw0m.beta_y, tw0m.alpha_x, tw0m.alpha_y))

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
    old_file = run_dir + '/run.' + str(run_id) + prefix1 + '.gout'
    new_file = run_dir + '/run.' + str(run_id) + prefix2 + '.gout'

    if save:
        os.system('cp ' + old_file + ' ' + new_file)
        os.system('cp ' + old_file + '.dfl ' + new_file + '.dfl 2>/dev/null')  # 2>/dev/null to supress error messages if no such file
        os.system('cp ' + old_file + '.dpa ' + new_file + '.dpa 2>/dev/null')
        os.system('cp ' + old_file + '.beam ' + new_file + '.beam 2>/dev/null')
        os.system('cp ' + run_dir + '/tmp.gen' + ' ' + run_dir + '/geninp.' + str(run_id) + prefix2 + '.inp 2>/dev/null')
        os.system('cp ' + run_dir + '/lattice.inp' + ' ' + run_dir + '/lattice.' + str(run_id) + prefix2 + '.inp 2>/dev/null')
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
        os.system('mv ' + run_dir + '/tmp.gen' + ' ' + run_dir + '/geninp.' + str(run_id) + prefix2 + '.inp 2>/dev/null')
        os.system('mv ' + run_dir + '/lattice.inp' + ' ' + run_dir + '/lattice.' + str(run_id) + prefix2 + '.inp 2>/dev/null')
    print ('      done')
    # os.system('rm ' + run_dir + '/run.' +str(run_id) + '.gout*')


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
