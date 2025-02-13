from __future__ import annotations

__author__ = 'Sergey Tomin'

import os
from dataclasses import dataclass, astuple
import copy
import logging
from time import time
import multiprocessing as mp
from typing import Union, List, Tuple, Optional, Any, Iterable

from scipy.stats import truncnorm
import pandas as pd
import numpy as np

from ocelot.cpbd.transformations.transformation import TMTypes
from ocelot.cpbd.optics import *
from ocelot.cpbd.beam import *
from ocelot.cpbd.errors import *
from ocelot.cpbd.elements import *
from ocelot.cpbd.io import is_an_mpi_process, ParameterScanFile
from ocelot.cpbd.physics_proc import CopyBeam
from ocelot.cpbd.navi import Navigator

_logger = logging.getLogger(__name__)

try:
    from scipy.signal import argrelextrema
    extrema_chk = 1
except:
    extrema_chk = 0

try:
    from mpi4py import MPI
except ImportError:
    IS_MPI = False
else:
    if is_an_mpi_process():
        IS_MPI = True
        COMM = MPI.COMM_WORLD
        N_CORES = COMM.Get_size()
        RANK = COMM.Get_rank()
    else:
        IS_MPI = False

def aperture_limit(lat, xlim = 1, ylim = 1):
    tws=twiss(lat, Twiss(), nPoints=1000)
    bxmax = max([tw.beta_x for tw in tws])
    bymax = max([tw.beta_y for tw in tws])
    bx0 = tws[0].beta_x
    by0 = tws[0].beta_y
    px_lim = float(xlim)/np.sqrt(bxmax*bx0)
    py_lim = float(ylim)/np.sqrt(bymax*by0)
    xlim = float(xlim)*np.sqrt(bx0/bxmax)
    ylim = float(ylim)*np.sqrt(by0/bymax)

    return xlim, ylim, px_lim, py_lim


def arg_peaks(data, extrema_chk = extrema_chk):
    """
    the function search peaks of spectrum and return positions of all peaks
    if extrema_chk == 1 uses numpy module
    if extrema_chk == 0 uses independent code (see below)
    """
    if extrema_chk == 1:
        return argrelextrema(data, np.greater)[0]
    else:
        diff_y = np.diff(data)
        extrm_y = np.diff(np.sign(diff_y))
        return np.where(extrm_y<0)[0]+1


def spectrum(data1D):
    """
    input: 1D sample data
    output: frequency and fourier transform
    """
    len_data1D = len(data1D)
    ft = np.abs(np.fft.fft(data1D))
    ft_shift = np.fft.fftshift(ft)
    freq = np.fft.fftshift(np.fft.fftfreq(len_data1D))
    return freq, ft_shift


def find_nearest(positions, value):
    """
    input: 1D array and value
    the function searches nearest value in the array to the given value
    """
    idx = (np.abs(positions-value)).argmin()
    return positions[idx]


def find_highest(sorted_posns, value, diap):
    """
    input: 1D array and value
    the function searches highest value in the array to the given value
    """
    poss = []
    for pos in sorted_posns:
        if value-diap <= pos <= value+diap:
            poss.append(pos)
    return poss[-1]


def nearest_particle(track_list, xi,yi):
    #x_array = np.unique(np.sort(map(lambda pxy: pxy.x, pxy_list)))
    y_array = np.unique(np.sort([pxy.y for pxy in track_list]))
    yi = find_nearest(y_array, yi)
    x_array_i = []
    for pxy in track_list:
        if pxy.y == yi:
            x_array_i.append(pxy.x)
    xi = find_nearest(np.array(x_array_i), xi)

    #print "inside nearest_particle, xi, yi : ", xi, yi
    for pxy in track_list:
        #print "inside nearest_particle: ", pxy.x, pxy.y
        if pxy.x == xi and pxy.y == yi:
            return pxy


def harmonic_position(data1D, nu = None, diap = 0.1, nearest = False):
    """
    function searches three highest harmonics and return:
    a. the highest if nu == None
    b. the nearest harmonics to the nu (if nu != None)
    """
    freq, ft_shift = spectrum(data1D)
    ft_maxi = arg_peaks(ft_shift)
    if len(ft_maxi) == 0:
        return -0.001
    freq_peaks = freq[ft_maxi][int(len(ft_maxi)/2):]
    peaks = ft_shift[ft_maxi][int(len(ft_maxi)/2):]

    main_3 = freq_peaks[np.argsort(peaks)]
    if nearest:
        return find_nearest(main_3, nu)

    if nu is None:
        return main_3[-1]
    if diap is None:
        main_3 = main_3[-5:]
        nearest_nu = find_nearest(main_3, nu)
    else:
        nearest_nu = find_highest(main_3, nu, diap)
    return nearest_nu


def freq_analysis(track_list, lat, nturns, harm=True, diap=0.10, nearest=False, nsuperperiods=1):

    def beta_freq(lat):

        tws = twiss(lat, Twiss())
        nux = tws[-1].mux/2./pi*nsuperperiods
        nuy = tws[-1].muy/2./pi*nsuperperiods
        print ("freq. analysis: Qx = ", nux, " Qy = ", nuy)
        nux = abs(int(nux+0.5) - nux)
        nuy = abs(int(nuy+0.5) - nuy)
        print("freq. analysis: nux = ", nux)
        print("freq. analysis: nuy = ", nuy)
        return nux, nuy

    nux, nuy = None, None
    if harm is True:
        nux, nuy = beta_freq(lat)
    #fma(pxy_list, nux = nux, nuy = nuy)
    for n, pxy in enumerate(track_list):
        if pxy.turn == nturns-1:
            if len(pxy.p_list) == 1:
                #print len(pxy.p_list)
                print ("For frequency analysis coordinates are needed for each turns. Check tracking option 'save_track' must be True ")
                return track_list
            x = [p[0] for p in pxy.p_list]
            y = [p[2] for p in pxy.p_list]
            pxy.mux = harmonic_position(x, nux, diap, nearest)
            pxy.muy = harmonic_position(y, nuy, diap, nearest)

    return track_list


class Track_info:
    def __init__(self, particle, x=0., y=0.):
        self.particle = particle
        self.turn = 0
        self.x = particle.x #initail coordinate
        self.y = particle.y #initail coordinate
        #self.x_array = [p.x]
        #self.y_array = [p.y]
        self.mux = -0.001
        self.muy = -0.001
        self.p_list = [[particle.x, particle.px, particle.y, particle.py, particle.tau, particle.p]]

    def get_x(self):
        return np.array([p[0] for p in self.p_list])

    def get_xp(self):
        return np.array([p[1] for p in self.p_list])

    def get_y(self):
        return np.array([p[2] for p in self.p_list])

    def get_yp(self):
        return np.array([p[3] for p in self.p_list])


def contour_da(track_list, nturns, lvl = 0.9):
    """
    the function defines contour of DA. If particle "lived" > lvl*nturns then we set up nturns
    if particle "lived" < lvl*nturns then we set up 0
    """
    if lvl>1:
        lvl = 1
    elif lvl<=0:
        lvl = 1

    ctr_da = []
    for pxy in track_list:
        if pxy.turn >= lvl*(nturns-1):
            ctr_da.append(nturns)
        else:
            ctr_da.append(0)
    return np.array(ctr_da)


def stable_particles(track_list, nturns):
    pxy_list_sbl = []
    for pxy in track_list:
        if pxy.turn >= nturns-1:
            pxy_list_sbl.append(pxy)

    return np.array(pxy_list_sbl)


def phase_space_transform(x,y, tws):
    """
    curved line of second order
    a11*x**2 + a22*y**2 + 2*a12*x*y + 2*a13*x + 2*a23*y + a33 = 0
    gamma*x**2 + 2*alpha*x*x' + beta*x'**2 = const
    """

    angle = np.arctan(2*tws.alpha_x/(tws.gamma_x-tws.beta_x))/2.
    x = x*np.cos(angle) - y*np.sin(angle)
    y = x*np.sin(angle) + y*np.cos(angle)
    return x,y


def create_track_list(x_array, y_array, p_array, energy=0.):
    """
    the function create list of Pxy
    """
    track_list = []
    for p in p_array:
        for y in (y_array):
            for x in (x_array):
                particle = Particle(x=x, y=y, p=p, E=energy)
                pxy = Track_info(particle, x, y)
                track_list.append(pxy)

    return track_list


def ellipse_track_list(beam, n_t_sigma = 3, num = 1000, type = "contour"):
    beam.sizes()
    #sigma_x = sqrt((sigma_e*tws0.Dx)**2 + emit*tws0.beta_x)
    #sigma_xp = sqrt((sigma_e*tws0.Dxp)**2 + emit*tws0.gamma_x)
    if type == "contour":
        t = np.linspace(0,2*pi, num)
        x = n_t_sigma*beam.sigma_x*np.cos(t)
        y = n_t_sigma*beam.sigma_xp*np.sin(t)
    else:
        x = truncnorm( -n_t_sigma,  n_t_sigma, loc=0, scale=beam.sigma_x).rvs(num)
        y = truncnorm( -n_t_sigma,  n_t_sigma, loc=0, scale=beam.sigma_xp).rvs(num)
    tws0 = Twiss(beam)
    x_array, xp_array = phase_space_transform(x,y, tws0)
    track_list = []
    for x,y in zip(x_array + beam.x, xp_array + beam.xp):
        p = Particle(x = x, px = y, p=-0.0)
        pxy = Track_info(p, x, y)
        track_list.append(pxy)

    return track_list


def track_nturns(lat, nturns, track_list, nsuperperiods=1, save_track=True, print_progress=True):
    xlim, ylim, px_lim, py_lim = aperture_limit(lat, xlim = 1, ylim = 1)
    navi = Navigator(lat)

    t_maps = navi.get_map(lat.totalLen)
    track_list_const = copy.copy(track_list)
    p_array = ParticleArray()
    p_list = [p.particle for p in track_list]
    p_array.list2array(p_list)

    for i in range(nturns):
        if print_progress: print(i)
        for n in range(nsuperperiods):
            for tm in t_maps:
                tm.apply(p_array)
            p_indx = p_array.rm_tails(xlim, ylim, px_lim, py_lim)

            track_list = np.delete(track_list, p_indx)
        for n, pxy in enumerate(track_list):
            pxy.turn = i
            if save_track:
                pxy.p_list.append(p_array.rparticles[:, n].tolist())
    return np.array(track_list_const)


def track_nturns_mpi(mpi_comm, lat, nturns, track_list, errors=None, nsuperperiods=1, save_track=True):
    size = mpi_comm.Get_size()
    rank = mpi_comm.Get_rank()
    lat_copy = create_copy(lat, nsuperperiods = nsuperperiods)
    nsuperperiods = 1
    if errors != None:
        if rank == 0:
            lat, errors = errors_seed(lat_copy, errors)
        else:
            errors = None

        errors = mpi_comm.bcast(errors, root=0)
        for i, elem in enumerate(lat_copy.sequence):
            elem.dx = errors[0][i]
            elem.dy = errors[1][i]
            elem.dtilt = errors[2][i]

    lat = MagneticLattice(lat_copy.sequence, method=lat_copy.method)

    if size == 1:
        # it is made to prevent memory crash in mpi_comm.gather() for one-tread case and for case of big pxy_list
        # (for instance, number of pxy in the list - nx*ny = 120*60 and nturns = 1000 (nturns means length of pxy class))
        # for instance, for case nturns = 500 is all ok
        # but for nturns = 1000 program crashes with error in mpi_comm.gather()
        # the same situation if treads not so much - solution increase number of treads.
        print("nsuperperiods = ", nsuperperiods)
        track_list = track_nturns(lat, nturns, track_list, nsuperperiods, save_track=save_track)
        return track_list

    if rank == 0:
        # dividing data into chunks
        chunks_track_list = [[] for _ in range(size)]
        N = len(track_list)
        for i, x in enumerate(track_list):
            chunks_track_list[int(size*i/N)].append(x)
    else:
        track_list = None
        chunks_track_list = None

    start = time()
    track_list = mpi_comm.scatter(chunks_track_list, root=0)
    print(" scatter time = ", time() - start, " sec, rank = ", rank, "  len(pxy_list) = ", len(track_list) )
    start = time()
    track_list = track_nturns(lat, nturns, track_list, nsuperperiods, save_track =save_track)
    print( " scanning time = ", time() - start, " sec, rank = ", rank)
    start = time()
    out_track_list = mpi_comm.gather(track_list, root=0)
    print(" gather time = ", time() - start, " sec, rank = ", rank)

    if rank == 0:
        start = time()
        track_list = []
        for i, chank in enumerate(out_track_list):
            for pxy in chank:
                track_list.append(pxy)
        print(" time exec = ", time() - start)
        return track_list


def fma(lat, nturns, x_array, y_array, nsuperperiods = 1):
    from mpi4py import MPI
    mpi_comm = MPI.COMM_WORLD
    rank = mpi_comm.Get_rank()
    track_list = create_track_list(x_array, y_array, p_array=[0])
    track_list = track_nturns_mpi(mpi_comm, lat, nturns, track_list, nsuperperiods=nsuperperiods)
    if rank == 0:
        nx = len(x_array)
        ny = len(y_array)
        ctr_da = contour_da(track_list, nturns)
        #ctr_da = tra.countour_da()
        track_list = freq_analysis(track_list, lat, nturns, harm = True)
        da_mux = np.array(map(lambda pxy: pxy.mux, track_list))
        da_muy = np.array(map(lambda pxy: pxy.muy, track_list))
        return ctr_da.reshape(ny,nx), da_mux.reshape(ny,nx), da_muy.reshape(ny,nx)


def da_mpi(lat, nturns, x_array, y_array, errors=None, nsuperperiods=1):
    from mpi4py import MPI
    mpi_comm = MPI.COMM_WORLD
    rank = mpi_comm.Get_rank()

    track_list = create_track_list(x_array, y_array, p_array=[0])
    track_list = track_nturns_mpi(mpi_comm, lat, nturns, track_list, errors=errors, nsuperperiods=nsuperperiods, save_track=False)

    if rank == 0:
        da = np.array(map(lambda track: track.turn, track_list))#.reshape((len(y_array), len(x_array)))
        nx = len(x_array)
        ny = len(y_array)
        return da.reshape(ny, nx)


def tracking_step(lat, particle_list, dz, navi):
    """
    tracking for a fixed step dz
    :param lat: Magnetic Lattice
    :param particle_list: ParticleArray or Particle list
    :param dz: step in [m]
    :param navi: Navigator
    :return: None
    """
    if lat is not navi.lat:
        _logger.error("MagneticLattice in the Navigator and tracking_step() is not the same")
        raise Exception("MagneticLattice in the Navigator and tracking_step() is not the same")

    if navi.z0 + dz > navi.lat.totalLen:
        dz = navi.lat.totalLen - navi.z0

    t_maps = navi.get_map(dz)
    for tm in t_maps:
        start = time()
        tm.apply(particle_list)
        _logger.debug(" tracking_step -> tm.class: " + tm.__class__.__name__  + "  l= "+ str(tm.length))
        _logger.debug(" tracking_step -> tm.apply: time exec = " + str(time() - start) + "  sec")
    return


def track(
        lattice,
        p_array,
        navi=None,
        print_progress=True,
        calc_tws=True,
        bounds=None,
        return_df=False,
        overwrite_progress=True,
        slice=None
        ) -> Tuple[Union[List[Twiss], pd.DataFrame], ParticleArray]:

    """
    tracking through the lattice

    :param lattice: Magnetic Lattice
    :param p_array: ParticleArray
    :param navi: Navigator, if None default Navigator wo any PhysProc
    :param print_progress: True, print tracking progress
    :param calc_tws: True, during the tracking twiss parameters are calculated from the beam distribution
    :param bounds: None, optional, [left_bound, right_bound] - bounds in units of std(p_array.tau())
    :return: twiss_list, ParticleArray. In case calc_tws=False, twiss_list is list of empty Twiss classes.
    """
    if navi is None:
        navi = Navigator(lattice)
    if navi.lat is not lattice:
        _logger.warning("MagneticLattice is not the same in lattice argument and in Navigator")
    tw0 = get_envelope(p_array, bounds=bounds, slice=slice) if calc_tws else Twiss()
    tws_track = [tw0]
    L = 0.

    for t_maps, dz, proc_list, phys_steps in navi.get_next_step():
        for tm in t_maps:
            tm.apply(p_array)
            _logger.debug("tracking_step -> tm.class: %s  l = %s", tm.__class__.__name__, tm.length)

        #part = p_array[0]
        for p, z_step in zip(proc_list, phys_steps):
            p.z0 = navi.z0
            p.apply(p_array, z_step)
        #p_array[0] = part
        if p_array.n == 0:
            _logger.debug(" Tracking stop: p_array.n = 0")
            return tws_track, p_array
        tw = get_envelope(p_array, bounds=bounds, slice=slice) if calc_tws else Twiss()
        L += dz
        tw.s += L
        tws_track.append(tw)

        if print_progress:
            names = [type(p).__name__ for p in proc_list] # process names
            names = ', '.join(names)
            msg = f"z = {navi.z0} / {lattice.totalLen}. Applied: {names}"
            end = "\n"
            if overwrite_progress:
                msg = f"\r{msg}"
                end = ""
            print(msg, end=end)

    # finalize PhysProcesses
    for p in navi.get_phys_procs():
        p.finalize()

    if return_df:
        return twiss_iterable_to_df(tws_track), p_array

    return tws_track, p_array


def lattice_track(lat, p):
    """
    Function tracks Particle through lattice and save Particles after each element

    :param lat: MagneticLattice
    :param p: Particle
    :return: list of Particles along the lattice
    """
    plist = [copy.copy(p)]
    for elem in lat.sequence:
        for tm in elem.tms:
            tm.apply([p])
        plist.append(copy.copy(p))
    return plist


def merge_drifts(lat):
    print( "before merging: len(sequence) = ", len(lat.sequence) )
    L = 0.
    seq = []
    new_elem = None
    for elem in lat.sequence:
        #next_elem = lat.sequence[i+1]
        if elem.__class__ == Drift:
            L += elem.l
            new_elem = Drift(l=L, eid=elem.id)
        else:
            if new_elem != None:
                seq.append(new_elem)
                L = 0.
                new_elem = None
            seq.append(elem)
    if new_elem != None:
        seq.append(new_elem)
    print( "after merging: len(sequence) = ", len(seq) )
    return MagneticLattice(sequence=seq)


def update_effective_beta(beam, lat):
    tws0 = Twiss()
    beta_x_eff = []
    beta_y_eff = []
    for beam_sl in beam:
        tws0.beta_x = beam_sl.beta_x
        tws0.beta_y = beam_sl.beta_y
        tws0.alpha_x = beam_sl.alpha_x
        tws0.alpha_y = beam_sl.alpha_y

        tws = twiss(lat, tws0)
        bx = [tw.beta_x for tw in tws]
        by = [tw.beta_y for tw in tws]

        beta_x_eff.append(np.mean(bx))
        beta_y_eff.append(np.mean(by))

    beam.beta_x_eff = np.array(beta_x_eff)
    beam.beta_y_eff = np.array(beta_y_eff)


class ParameterScanner:
    """A class for performing parameter scans (multiple calls to track with some
    different input ParticleArray or other variation (e.g. a dipole angle, a
    PhysProc hyperparameter---anything modifiable from a Navigator instance).

    :param navigator: The base Navigator instance to be used for the tracking.  Any
    modifications can be done in prepare_navigator.
    :type Navigator:
    :param parameter_values: The list of values to be passed to prepare_navigator
    and also written to the output file.
    :param parray0: ParticleArray or list of ParticleArray instances to be used as
    input for the tracking.  If an iterable of insstances is provided, it must be
    equal in length to parameter_values.
    :param parameter_name: Optional metadata to be written to the outputfile
    providing the name of the written parameter_values.
    :param markers: List of Marker instances.

    """
    # For now user should ensure nprocesses is not great than number of jobs in mpi
    # case.  otherwise empty output groups are made.
    def __init__(self,
                 navigator: Navigator,
                 parameter_values: Iterable[Any],
                 parray0: Union[ParticleArray, Iterable[ParticleArray]],
                 parameter_name: str = "parameter",
                 markers: Optional[Iterable[Marker]] = None,
                 ):
        self.navigator = navigator
        self.parameter_values = parameter_values
        self.parray0 = parray0
        self.parameter_name = parameter_name
        self.markers = markers
        if self.markers is None:
            self.markers = []

        # If an iterable of input is provided.
        if (not isinstance(self.parray0, ParticleArray)
                and len(self.parray0) != len(self.parameter_values)):
            raise ValueError(
                "parameter_values length doesn't match number of parrays"
            )

    def prepare_navigator(self,
                          _value: Any = None,
                          _parray0: Optional[ParticleArray] = None,
                          _job_index: Optional[int] = None) -> Navigator:
        return copy.deepcopy(self.navigator)

    def generate_unique_navigators_with_parray0s(self):
        # Instantiate a navigator for each parameter.
        parray0s = self._prepare_parray0s()
        navigators = []
        for job_index, v in enumerate(self.parameter_values):
            parray0 = parray0s[job_index]
            navi = self.prepare_navigator(v, parray0, job_index)
            navigators.append(navi)
        self._attach_dump_processes_to_markers(navigators)
        return navigators, parray0s

    def _attach_dump_processes_to_markers(self, navigators: Iterable[Navigator]):
        sequence = self.navigator.lat.sequence
        marker_locations = [sequence.index(m) for m in self.markers]
        for navi in navigators:
            sequence = navi.lat.sequence
            for i in marker_locations:
                marker = sequence[i]
                process = CopyBeam(marker.id)
                navi.add_physics_proc(process, marker, marker)

    def prepare_track_payloads(self, navigators, parray0s):
        payloads = []
        for i in range(self.njobs):
            navigator = navigators[i]
            parray0 = parray0s[i]
            # Prepare arguments for the processes.  Each one gets a differently
            # prepared navigator and a copy of the parray.
            payload = TrackPayload(navigator.lat, parray0.copy(), navigator,
                                      job_index=i)

            payloads.append(payload)
        return payloads

    @property
    def njobs(self) -> int:
        return len(self.parameter_values)

    def _prepare_parray0s(self) -> ParticleArray:
        # Either given a list of ParticleArray instances or a single one.  If a
        # list it has to have the same length as parameter_values.
        if isinstance(self.parray0, ParticleArray):
            return len(self.parameter_values) * [self.parray0]
        if len(self.parameter_values) != len(self.parray0):
            raise ValueError(f"Parameter values len != [parray0]")
        return self.parray0 # Then it is an iterable of ParticleArray instances

    def scan(self, filename: str, nproc: int = 1):
        """Run the parameter scan, possibly across multiple cores.  MPI is
        automatically detected if used, where nproc will then have no effect.


        :param nproc: Number of processes to spawn.  If set to -1, then spawn a
        number of cores equal to the number of CPUs on this computer.

        """

        # Map of scanned parameter values to corresponding track arg pack to be
        # run.
        navigators, parray0s = self.generate_unique_navigators_with_parray0s()
        args = self.prepare_track_payloads(navigators, parray0s)
        with ParameterScanFile(filename, "w") as psf:
            psf.init_from_parameter_scanner(self)
            if is_an_mpi_process():
                self._run_mpi(args, psf)
            else:
                self._run_pool(args, psf, nproc)

    def _run_pool(self, track_payloads: Iterable[TrackPayload], psf, nproc: int) -> None:
        # TODO Better:
        # https://stackoverflow.com/questions/15704010/write-data-to-hdf-file-using-multiprocessing
        if nproc == -1: # Spawn as many jobs as possibly or necessary.
            nproc = min(self.njobs, mp.cpu_count())

        # outqueue = mp.Queue()
        # inqueue = mp.Queue()
        # # Process dedicated to writing output
        # output_proc = mp.Process(target=_handle_process_output, args=(outqueue, ))

        with mp.Pool(nproc) as p:
            results = p.map(TrackPayload.as_function, track_payloads)

        parray0s = [a.parray0 for a in track_payloads]
        parray1s = [x[0] for x in results]
        all_jobs_dumps = [x[1] for x in results]

        # Write input and output arrays, maybe also output array for no physics.
        for job_index, parray0 in enumerate(parray0s):
            psf.write_parray0(job_index, parray0)
            psf.write_parray1(job_index, parray1s[job_index])

        # Write marker parrays.
        for job_index, one_jobs_dumps in enumerate(all_jobs_dumps):
            for dump in one_jobs_dumps:
                psf.write_parray_marker(job_index, dump.name, dump.parray)

    def _run_mpi(self, payloads: Iterable[TrackPayload], psf) -> None:
        job_indices = self._get_job_indices_for_this_mpi_core()
        this_cores_track_args = [payloads[i] for i in job_indices]
        if not job_indices: # If given no jobs to run on this core, do nothing
            return

        for payload in payloads:
            job_index = payload.job_index
            parray0 = payload.parray0.copy()
            parray1, dumps = payload.run_and_get_dumps()

            psf.write_parray0(job_index, parray0)
            psf.write_parray1(job_index, parray1)

            for dump in dumps:
                psf.write_parray_marker(job_index, dump.name, dump.parray)

            _logger.info(
                f"Written marker distributions for %s at %s",
                job_index, RANK
            )

        _logger.info("Finished all parameters to be scanned"
                     f" at RANK={RANK} and results written to {psf.filename}")

    def _get_job_indices_for_this_mpi_core(self) -> List[int]:
        all_job_indices = range(self.njobs)
        this_cores_job_indices = np.array_split(all_job_indices, N_CORES)[RANK]
        # Converting to a list actually matters because bool(np.array([0])) is
        # False whereas bool([0]) is True...!
        this_cores_job_indices = list(this_cores_job_indices)
        return this_cores_job_indices


@dataclass
class TrackPayload:
    lattice: MagneticLattice
    parray0: ParticleArray
    navigator: Navigator
    calc_tws: bool = False
    bounds: Optional[Tuple[float, float]] = None
    return_df: bool = True
    print_progress: bool = True
    overwrite_progress: bool = False
    job_index: int = 0
    kwargs: Optional[Mapping[str, Any]] = None

    def run_and_get_dumps(self):
        _logger.info("Starting tracking for job number %s.", self.job_index)

        _, parray1 = track(lattice=self.lattice,
                           p_array=self.parray0,
                           navi=self.navigator,
                           print_progress=self.print_progress,
                           calc_tws=self.calc_tws,
                           bounds=self.bounds,
                           return_df=self.return_df,
                           overwrite_progress=self.overwrite_progress)

        marker_dump_processes = []
        for process in self.navigator.inactive_processes:
            if isinstance(process, CopyBeam):
                marker_dump_processes.append(process)

        _logger.info("Finished tracking for job number %s.", self.job_index)

        return parray1, marker_dump_processes

    @staticmethod
    def as_function(payload):
        """For use with pool.map"""
        return payload.run_and_get_dumps()


class UnitStepScanner(ParameterScanner):
    """Simple ParameterScanner subclass for scanning the Navigator unit step."""
    def prepare_navigator(self, unit_step: float, _parray0, _job_index) -> Navigator:
        navi = super().prepare_navigator(unit_step, _parray0, _job_index)
        navi.unit_step = unit_step
        return navi
