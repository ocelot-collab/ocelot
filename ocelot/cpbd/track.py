__author__ = 'Sergey Tomin'

from ocelot.cpbd.optics import *
from ocelot.cpbd.beam import *
from ocelot.cpbd.errors import *
from ocelot.cpbd.elements import *
from time import time
from scipy.stats import truncnorm
import copy
import sys
import logging

_logger = logging.getLogger(__name__)

try:
    from scipy.signal import argrelextrema
    extrema_chk = 1
except:
    extrema_chk = 0


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
        if value-diap<=pos<=value+diap:
            poss.append(pos)
    #print poss
    return poss[-1]
    #idx = (np.abs(sorted_posns-value)).argmin()
    #return sorted_posns[idx]


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

    main_3 =  freq_peaks[np.argsort(peaks)]

    if nearest:
        return find_nearest(main_3, nu)

    if nu == None:
        return main_3[-1]
    if diap == None:
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
    if harm == True:
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

    t_maps = get_map(lat, lat.totalLen, navi)
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
                pxy.p_list.append(p_array.rparticles[:, n])
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
    if navi.z0 + dz > lat.totalLen:
        dz = lat.totalLen - navi.z0

    t_maps = get_map(lat, dz, navi)
    for tm in t_maps:
        start = time()
        tm.apply(particle_list)
        _logger.debug(" tracking_step -> tm.class: " + tm.__class__.__name__  + "  l= "+ str(tm.length))
        _logger.debug(" tracking_step -> tm.apply: time exec = " + str(time() - start) + "  sec")
    return


def track(lattice, p_array, navi, print_progress=True, calc_tws=True, bounds=None):
    """
    tracking through the lattice

    :param lattice: Magnetic Lattice
    :param p_array: ParticleArray
    :param navi: Navigator
    :param print_progress: True, print tracking progress
    :param calc_tws: True, during the tracking twiss parameters are calculated from the beam distribution
    :param bounds: None, optional, [left_bound, right_bound] - bounds in units of std(p_array.tau())
    :return: twiss_list, ParticleArray. In case calc_tws=False, twiss_list is list of empty Twiss classes.
    """

    tw0 = get_envelope(p_array, bounds=bounds) if calc_tws else Twiss()
    tws_track = [tw0]
    L = 0.

    while np.abs(navi.z0 - lattice.totalLen) > 1e-10:
        if navi.kill_process:
            _logger.info("Killing tracking ... ")
            return tws_track, p_array

        dz, proc_list, phys_steps = navi.get_next()
        tracking_step(lat=lattice, particle_list=p_array, dz=dz, navi=navi)
        #part = p_array[0]
        for p, z_step in zip(proc_list, phys_steps):
            p.z0 = navi.z0
            p.apply(p_array, z_step)
        #p_array[0] = part
        if p_array.n == 0:
            _logger.debug(" Tracking stop: p_array.n = 0")
            return tws_track, p_array
        tw = get_envelope(p_array, bounds=bounds) if calc_tws else Twiss()
        L += dz
        tw.s += L
        tws_track.append(tw)

        if print_progress:
            poc_names = [p.__class__.__name__ for p in proc_list]
            sys.stdout.write( "\r" + "z = " + str(navi.z0)+" / "+str(lattice.totalLen) + " : applied: " + ", ".join(poc_names)  )
            sys.stdout.flush()

    # finalize PhysProcesses
    for p in navi.get_phys_procs():
        p.finalize()

    return tws_track, p_array


def lattice_track(lat, p):
    plist = [copy.copy(p)]

    for elem in lat.sequence:
        elem.transfer_map.apply([p])
        #print(p)
        if not (elem.__class__ in [Bend, RBend, SBend] and elem.l != 0.): #, "hcor", "vcor"
            if elem.__class__ == Edge:
                #print elem.pos
                if elem.pos == 1:
                    continue
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

