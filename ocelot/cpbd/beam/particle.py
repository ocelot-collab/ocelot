import numpy as np
from copy import deepcopy
from typing import TypeVar

import ocelot.common.globals as glb
from ocelot.common.ocelog import *
from . import analysis
from . import noise
from . import core

_logger = logging.getLogger(__name__)

TypeParticleArray = TypeVar("TypeParticleArray", bound="ParticleArray")


class ParticleArray:
    """
    array of particles of fixed size; for optimized performance
    (x, x' = px/p0),(y, y' = py/p0),(ds = c*tau, p = dE/(p0*c))
    p0 - ntum
    """

    class LostParticleRecorder:
        """
        Stores information about particles that are getting lost.

        Attributes:
            lost_particles: List of indices of deleted particle. Notes: The indices of the initial ParticleArray are used.
            lp_to_pos_hist: Histogram of number of lost particles to position s. [(pos, num of lost particles), ...]
        """

        def __init__(self, n):
            self.lost_particles = []
            self.lp_to_pos_hist = []
            self._current_particle = np.arange(n)

        def add(self, inds, position):
            self.lost_particles += self._current_particle[inds].tolist()
            self.lp_to_pos_hist.append((position, len(inds)))
            self._current_particle = np.delete(self._current_particle, inds)

        def initial_idx_2_p_idx(self, idx):
            found_idx = np.where(self._current_particle == idx)[0]
            if found_idx.size > 0:
                return found_idx[0]
            return None

    @classmethod
    def random(cls, n, sigma_x=0.000121407185261, sigma_px=1.80989470506e-05, sigma_y=0.000165584800564,
               sigma_py=4.00994225888e-05):
        # generate beam file
        x = np.random.randn(n) * sigma_x
        px = np.random.randn(n) * sigma_px
        y = np.random.randn(n) * sigma_y
        py = np.random.randn(n) * sigma_py

        # covariance matrix for [tau, p] for beam compression in BC
        cov_t_p = [[1.30190131e-06, 2.00819771e-05],
                   [2.00819771e-05, 3.09815718e-04]]
        long_dist = np.random.multivariate_normal((0, 0), cov_t_p, n)
        tau = long_dist[:, 0]
        dp = long_dist[:, 1]

        p_array = cls(n=n)
        p_array.E = 0.130  # GeV
        p_array.rparticles[0] = x
        p_array.rparticles[1] = px
        p_array.rparticles[2] = y
        p_array.rparticles[3] = py
        p_array.rparticles[4] = tau
        p_array.rparticles[5] = dp

        Q = 5e-9

        p_array.q_array = np.ones(n) * Q / n
        return p_array

    def __init__(self, n=0):
        self.rparticles = np.zeros((6, n))
        self.q_array = np.zeros(n)  # charge
        self.s = 0.0
        self.E = 0.0
        self.lost_particle_recorder = self.LostParticleRecorder(n)

    def rm_tails(self, xlim, ylim, px_lim, py_lim):
        """
        Method removes particles outside range [-xlim, +xlim], [-px_lim, +px_lim] ...

        """
        x = abs(self.x())
        px = abs(self.px())
        y = abs(self.y())
        py = abs(self.py())
        ind_angles = np.append(np.argwhere(px > px_lim), np.argwhere(py > py_lim))
        p_idxs = np.unique(np.append(np.argwhere(x > xlim), np.append(np.argwhere(y > ylim),
                                                                      np.append(np.argwhere(x != x),
                                                                                np.append(np.argwhere(y != y),
                                                                                          ind_angles)))))
        self.delete_particles(p_idxs)
        return p_idxs

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            result = ParticleArray()
            result.rparticles = self.rparticles[..., idx]
            result.q_array = self.q_array[idx]
            return result

        return Particle(x=self.rparticles[0, idx], px=self.rparticles[1, idx],
                        y=self.rparticles[2, idx], py=self.rparticles[3, idx],
                        tau=self.rparticles[4, idx], p=self.rparticles[5, idx],
                        s=self.s)

    def __setitem__(self, idx, p):
        self.rparticles[0, idx] = p.x
        self.rparticles[1, idx] = p.px
        self.rparticles[2, idx] = p.y
        self.rparticles[3, idx] = p.py
        self.rparticles[4, idx] = p.tau
        self.rparticles[5, idx] = p.p
        self.s = p.s

    def list2array(self, p_list):
        self.rparticles = np.zeros((6, len(p_list)))
        self.q_array = np.zeros(len(p_list))
        for i, p in enumerate(p_list):
            self[i] = p
        self.s = p_list[0].s
        self.E = p_list[0].E
        self.lost_particle_recorder = self.LostParticleRecorder(len(p_list))

    def array2list(self):
        p_list = []
        for i in range(self.size()):
            p_list.append(self[i])
        return p_list

    def array2ex_list(self, p_list):

        for i, p in enumerate(p_list):
            p.x = self.rparticles[0, i]
            p.px = self.rparticles[1, i]
            p.y = self.rparticles[2, i]
            p.py = self.rparticles[3, i]
            p.tau = self.rparticles[4, i]
            p.p = self.rparticles[5, i]
            p.E = self.E
            p.s = self.s
        return p_list

    def size(self):
        return int(self.rparticles.size / 6)

    def x(self):
        return self.rparticles[0]

    def px(self):
        return self.rparticles[1]  # xp

    def y(self):
        return self.rparticles[2]

    def py(self):
        return self.rparticles[3]  # yp

    def tau(self):
        return self.rparticles[4]

    def p(self):
        return self.rparticles[5]

    @property
    def t(self):
        return self.rparticles[4]

    @t.setter
    def t(self, value):
        self.rparticles[4] = value

    @property
    def n(self):
        return np.shape(self.rparticles)[1]

    @property
    def pz(self) -> float:
        """pz/p0 - the z-components of the macroparticle momenta normalised with
        respect to the reference momentum p0."""
        return np.sqrt((self.momenta / self.p0c) ** 2 - self.px() ** 2 - self.py() ** 2)

    @property
    def p0c(self) -> float:
        """Get macroparticle reference momentum * speed of light in GeV."""
        return np.sqrt(self.E ** 2 - glb.m_e_GeV ** 2)

    @property
    def energies(self) -> float:
        """Get all macroparticle energies in GeV."""
        return self.p() * self.p0c + self.E

    @property
    def momenta(self) -> float:
        """Get all macroparticle momenta in GeV/c."""
        return np.sqrt(self.energies ** 2 - glb.m_e_GeV ** 2)

    @property
    def gamma(self) -> float:
        """Get all macroparticle relativistic gamma factors."""
        return self.energies / glb.m_e_GeV

    @property
    def beta(self) -> float:
        """Get all macroparticle relativistic betas (v/c)."""
        return np.sqrt(1 - self.gamma ** -2)

    @property
    def total_charge(self) -> float:
        """Get total charge of the ParticleArray"""
        return sum(self.q_array)

    def sort(self, variable, in_place=True) -> np.ndarray:
        """Sort ParticleArray in place according to the chosen key.

        :param variable: One of "x", "px", "y", "py", "tau", "p" or one of the other
        macroparticle properties (e.g. pz, momenta, etc.)  with which to sort the
        macroparticle array by.

        """

        try:
            member = getattr(self, variable)
        except AttributeError:
            pass

        try:
            indices = member().argsort()
        except TypeError:
            try:
                indices = member.argsort()
            except TypeError:
                raise ValueError(f"Unknown variable name for ParticleArray: {variable}")

        if in_place:
            self.rparticles = self.rparticles[..., indices]
            self.q_array = self.q_array[indices]

        return indices

    def thin_out(self, nth=10, n0=0):
        """
        Method to thin out the particle array in n-th times. Means every n-th particle will be saved in new Particle array

        :param nth: 10, every n-th particle will be taken to new Particle array
        :param n0: start from n0 particle
        :return: New ParticleArray
        """
        nth = int(nth)
        if nth <= 1:
            print("Nothing to do. nth number must be bigger 1")
            return self
        if nth > self.n:
            raise ValueError("nth number is bigger of particles number")
        if n0 > self.n:
            raise ValueError("n0 number is bigger of particles number")
        n = int((self.n - n0) / nth)
        if n < 1:
            raise ValueError("Number of particles in new ParticleArray is less then 1")

        n_end = n0 + nth * n
        p = ParticleArray(n)
        p.rparticles[:, :] = self.rparticles[:, n0:n_end:nth]
        p.q_array[:] = self.q_array[n0:n_end:nth] * nth
        p.s = self.s
        p.E = self.E
        return p

    def rm_particle(self, index):
        """
        Method removes a "bad" particle with particular "index".
        :param index:
        :return:
        """
        _logger.warning("rm_particle is deprecated and will be removed in the future. Use delete_particles instead.")
        self.rparticles = np.delete(self.rparticles, index, 1)
        self.q_array = np.delete(self.q_array, index, 0)

    def rescale2energy(self, energy):
        """
        Method to rescale beam coordinates with new energy

        :param energy: new energy
        :return:
        """
        Ei = self.E
        betai = np.sqrt(1 - (glb.m_e_GeV / Ei) ** 2)
        Ef = energy
        betaf = np.sqrt(1 - (glb.m_e_GeV / Ef) ** 2)

        self.E = Ef
        Rnn = Ei * betai / (Ef * betaf)
        self.px()[:] = Rnn * self.px()[:]
        self.py()[:] = Rnn * self.py()[:]
        self.p()[:] = Rnn * self.p()[:]

    def __str__(self):
        val = "ParticleArray: \n"
        val += "Ref. energy : " + str(np.round(self.E, 4)) + " GeV \n"
        val += "Ave. energy : " + str(np.around(self.E * (1 + np.mean(self.p())), 4)) + " GeV \n"
        val += "std(x)      : " + str(np.round(np.std(self.x()) * 1e3, 3)) + " mm\n"
        val += "std(px)     : " + str(np.round(np.std(self.px()) * 1e3, 3)) + " mrad\n"
        val += "std(y)      : " + str(np.round(np.std(self.y()) * 1e3, 3)) + " mm\n"
        val += "std(py)     : " + str(np.round(np.std(self.py()) * 1e3, 3)) + " mrad\n"
        val += "std(p)      : " + str(np.round(np.std(self.p()), 4)) + "\n"
        val += "std(tau)    : " + str(np.round(np.std(self.tau()) * 1e3, 3)) + " mm\n"
        val += "Charge      : " + str(np.around(np.sum(self.q_array) * 1e9, 4)) + " nC \n"
        val += "s pos       : " + str(self.s) + " m \n"
        val += "n particles : " + str(self.n) + "\n"
        return val

    def __len__(self):
        return self.size()

    def delete_particles(self, inds, record=True):
        """
        Deletes particles from the particle array via index.
        :param inds: Indices that will be removed from the particle array.
        :param record: If record is true the deleted particles will be saved in self.lost_particle_recorder.
        :return:
        """
        if record:
            self.lost_particle_recorder.add(inds, self.s)
        self.rparticles = np.delete(self.rparticles, inds, 1)
        self.q_array = np.delete(self.q_array, inds, 0)

    def copy(self) -> TypeParticleArray:
        """Return a copy of this ParticleArray instance."""
        return deepcopy(self)

    def get_twiss(self, tws_i=None, bounds=None, slice=None, auto_disp=False):
        """
        Calculate Twiss parameters from the ParticleArray.

        This method computes the statistical beam parameters (Twiss) from the particle distribution,
        optionally using dispersion correction, bounding filters, and a reference slice definition.

        Parameters:
        - tws_i : Twiss, optional
            Design Twiss parameters used for dispersion correction. Defaults to None (no correction).
        - bounds : list, optional
            Specifies the region of interest as [left_bound, right_bound], in units of std(tau).
            Only particles within these longitudinal bounds are considered.
        - slice : str or None, optional
            Defines how to choose the reference slice when `bounds` is used:
            - None (default): Uses the central slice at z0 = mean(tau).
            - "Imax": Uses the slice where the current is maximal.
         - auto_disp : bool, optional
            If True and tws_i is None, estimate and subtract linear dispersion from the statistics of the particle array.
            Default is False.

        Returns:
        - Twiss
            The Twiss parameters computed from the filtered particle array.
        """
        tws = analysis.get_envelope(self, tws_i=tws_i, bounds=bounds, slice=slice, auto_disp=auto_disp)
        return tws

    def get_twiss_from_slice(self, slice="Imax", nparts_in_slice=5000, smooth_param=0.05, filter_base=2, filter_iter=2):
        """
        Function calculates twiss parameters in a beam slice

        :param parray: ParticleArray
        :param slice: "Imax" or "Emax" or center of bunch
        :param nparts_in_slice: 5000, number of particles in the slice (in moving window)
        :param smooth_param: 0.01, smoothing parameters to calculate the beam current: smooth_param = m_std * np.std(p_array.tau())
        :param filter_base: 2, filter parameter in the func: simple_filter
        :param filter_iter: 2, filter parameter in the func: simple_filter
        :return: Twiss
        """
        tws = core.Twiss()
        slice_params = analysis.global_slice_analysis(self, nparts_in_slice=nparts_in_slice, smooth_param=smooth_param,
                                             filter_base=filter_base, filter_iter=filter_iter)

        if slice == "Imax":
            ind0 = np.argmax(slice_params.I)
        elif slice == "Emax":
            ind0 = np.argmax(slice_params.me)
        else:
            ind0 = np.argsort(np.abs(slice_params.s))[0]

        tws.beta_x = slice_params.beta_x[ind0]
        tws.alpha_x = slice_params.alpha_x[ind0]
        tws.beta_y = slice_params.beta_y[ind0]
        tws.alpha_y = slice_params.alpha_y[ind0]
        tws.emit_x = slice_params.ex[ind0]
        tws.emit_y = slice_params.ey[ind0]
        tws.pp = (slice_params.se[ind0] * 1e-9 / self.E) ** 2
        tws.E = self.E
        tws.s = self.s

        return tws

    def I(self, num_bins=None):
        """
        simple function to calculate current profile form the beam distribution.
        :return: np.array(Nx2), where s = B[:, 0], I = B[:, 1]
        """
        if num_bins is None:
            sigma = np.std(self.tau()) / 10.
            q0 = np.sum(self.q_array)
            relgamma = self.E / glb.m_e_GeV
            relbeta = np.sqrt(1 - relgamma ** -2) if relgamma != 0 else 1.
            v = relbeta * glb.speed_of_light
            B = analysis.s_to_cur(self.tau(), sigma, q0, v)
        else:
            s, I = analysis.get_current(self, num_bins=num_bins)
            B = np.column_stack((s, I))
        return B

    def quietify_tau(self, inv_cdf=None, quiet_method="stratified", quiet_alpha=1.0,
                     quiet_jitter_fraction=0.2, noise_Ne=None,
                     kmin=2 * np.pi / 1e-5, kmax=2 * np.pi / 1e-6, rng=None):
        """
        Quietify the longitudinal coordinate (tau).

        Parameters
        ----------
        inv_cdf : callable or None
            Inverse CDF function for tau distribution. If None, an empirical
            inverse CDF is built from the current tau histogram.
        quiet_method : str
            One of {"stratified","sobol","jittered"}.
        quiet_alpha : float
            Blend factor (0=keep random, 1=fully quiet).
        quiet_jitter_fraction : float
            Jitter fraction if method="jittered".
        noise_Ne : int or None
            If set, inject spectral noise calibrated to 1/sqrt(N_e) in [kmin,kmax].
        kmin, kmax : float
            Band for spectral noise injection (in 1/m).
        rng : np.random.Generator or None
            RNG for reproducibility.

        """
        if rng is None:
            rng = np.random.default_rng()
        tau = self.tau()
        if inv_cdf is None:
            inv_cdf = noise.make_inverse_cdf_from_samples(tau, bins=500)

        tau = noise.quietify_1d(
            tau, inverse_cdf=inv_cdf,
            method=quiet_method, alpha=quiet_alpha,
            jitter_fraction=quiet_jitter_fraction, rng=rng
        )

        if noise_Ne is not None:
            tau = noise.inject_spectral_noise_tau(
                tau, N_e=noise_Ne,
                kmin=kmin, kmax=kmax, M=8192, rng=rng
            )
        self.rparticles[4][:] = tau

    def quietify_transverse(self, coords=("x", "y", "px", "py", "dp"),
                            quiet_method="stratified", quiet_alpha=1.0,
                            quiet_jitter_fraction=0.2, rng=None):
        """
        Quietify transverse (and momentum) coordinates assuming Gaussian distributions.

        Parameters
        ----------
        coords : tuple of str
            Subset of {"x","y","px","py","dp"} to quietify. Default = all.
        quiet_method, quiet_alpha, quiet_jitter_fraction, rng : see quietify_tau.
        """
        if rng is None:
            rng = np.random.default_rng()

        coord_map = {
            "x": (0, self.x().std()),
            "px": (1, self.px().std()),
            "y": (2, self.y().std()),
            "py": (3, self.py().std()),
            "dp": (5, self.p().std()),
        }

        for name, (idx, sigma) in coord_map.items():
            if name in coords:
                arr = self.rparticles[idx][:]
                inv_cdf_gauss = lambda u, s=sigma: s * noise.invnorm(u)
                self.rparticles[idx][:] = noise.quietify_1d(
                    arr, inverse_cdf=inv_cdf_gauss,
                    method=quiet_method, alpha=quiet_alpha,
                    jitter_fraction=quiet_jitter_fraction, rng=rng
                )

class Particle:
    """
    particle
    to be used for tracking
    """

    def __init__(self, x=0.0, y=0.0, px=0.0, py=0.0, s=0.0, p=0.0, tau=0.0, E=0.0, q=0.0):
        self.x = x
        self.y = y
        self.px = px  # horizontal (generalized) momentum
        self.py = py  # vertical (generalized) momentum
        self.p = p  # longitudinal momentum
        self.s = s
        self.tau = tau  # time-like coordinate wrt reference particle in the bunch (e.g phase)
        self.E = E  # energy
        self.q = q  # charge in C

    def __str__(self):
        val = ""
        val = val + "x = " + str(self.x) + "\n"
        val = val + "px = " + str(self.px) + "\n"
        val = val + "y = " + str(self.y) + "\n"
        val = val + "py = " + str(self.py) + "\n"
        val = val + "tau = " + str(self.tau) + "\n"
        val = val + "p = " + str(self.p) + "\n"
        val = val + "E = " + str(self.E) + "\n"
        val = val + "s = " + str(self.s)
        return val

    def multiply_with_tm(self, tm: 'TransferMap', length):
        tm.apply(self)
        return deepcopy(self)


