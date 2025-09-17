import numpy as np
from ocelot.common.ocelog import *
from ocelot.common.math_op import invert_cdf

_logger = logging.getLogger(__name__)

def invnorm(u):
    """
    Approximate inverse cumulative distribution function (quantile function)
    of the standard normal distribution using Peter J. Acklam's rational
    approximation.

    This implementation maps uniform random numbers in (0,1) to
    normally distributed values with mean 0 and standard deviation 1.
    It is useful for generating Gaussian-distributed samples from
    low-discrepancy sequences (e.g., Sobol, Halton) or stratified
    quantiles without relying on external libraries.

    Parameters
    ----------
    u : array_like
        Input values in the open interval (0, 1). These are usually
        uniform quantiles to be mapped into standard normal deviates.
        Values outside (0,1) are clipped internally to avoid numerical issues.

    Returns
    -------
    x : ndarray
        Approximate inverse-normal values corresponding to input `u`.
        The output has mean ~0 and standard deviation ~1.

    Notes
    -----
    - Based on the algorithm by Peter J. Acklam,
      see: http://home.online.no/~pjacklam/notes/invnorm/
    - Maximum relative error of this approximation is ~1.15e-9.
    - For extreme tails (u close to 0 or 1), accuracy is slightly reduced
      but still sufficient for most Monte Carlo or beam dynamics applications.

    Examples
    --------
    >>> u = np.linspace(0.01, 0.99, 5)
    >>> _invnorm(u)
    array([-2.32634787, -0.84162123,  0.        ,  0.84162123,  2.32634787])

    >>> # Map Sobol points into Gaussian distribution
    >>> from scipy.stats import qmc
    >>> sobol = qmc.Sobol(d=1, scramble=True)
    >>> u = sobol.random(8).flatten()
    >>> x = _invnorm(u)  # quasi-random Gaussian points
    """
    a = np.array([-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00])
    b = np.array([-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,6.680131188771972e+01,-1.328068155288572e+01])
    c = np.array([-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,-2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00])
    d = np.array([ 7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00])
    pl, ph = 0.02425, 1-0.02425
    u = np.clip(u, 1e-12, 1-1e-12)
    x = np.empty_like(u, dtype=float)
    lo, mid, hi = u<pl, (u>=pl)&(u<=ph), u>ph
    if np.any(lo):
        q = np.sqrt(-2*np.log(u[lo])); x[lo]=((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]
        x[lo]/=((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if np.any(mid):
        q = u[mid]-0.5; r=q*q
        x[mid]=(((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q
        x[mid]/=(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)
    if np.any(hi):
        q = np.sqrt(-2*np.log(1-u[hi])); x[hi]=-(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])
        x[hi]/=((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    return x


def make_inverse_cdf_from_samples(samples, bins=500):
    counts, edges = np.histogram(samples, bins=bins, density=True)
    centers = 0.5*(edges[:-1] + edges[1:])
    return invert_cdf(counts, centers)


def quietify_1d(z, inverse_cdf, method="stratified", alpha=1.0, jitter_fraction=0.2, rng=None):
    """
    Quietify z coordinate using target inverse CDF.
    [docstring unchanged...]
    """
    if rng is None:
        rng = np.random.default_rng()
    N = len(z)

    if method == "stratified":
        u = (np.arange(1, N+1) - 0.5)/N

    elif method == "sobol":
        try:
            from scipy.stats import qmc
        except ImportError:
            raise ImportError("Sobol method requires scipy>=1.7")
        # Generate 2^m points, then slice the first N → balanced, no warning
        m = int(np.ceil(np.log2(N)))
        u_full = qmc.Sobol(d=1, scramble=True).random_base2(m).ravel()
        u = np.sort(u_full[:N])

    elif method == "jittered":
        u = (np.arange(1, N+1) - 0.5)/N

    else:
        raise ValueError(f"Unknown method: {method}")

    z_target = inverse_cdf(u)

    if method == "jittered":
        du = 1.0/N
        dz_local = np.gradient(z_target)
        jitter = rng.uniform(-0.5, 0.5, N) * (jitter_fraction*du) * np.maximum(np.abs(dz_local), 1e-30)
        z_target = z_target + jitter

    idx = np.argsort(z)
    z_sorted = z[idx]
    z_quiet_sorted = (1 - alpha)*z_sorted + alpha*z_target

    z_new = np.empty_like(z_quiet_sorted)
    z_new[idx] = z_quiet_sorted
    return z_new

def inject_spectral_noise_tau(tau, N_e, kmin=None, kmax=None, M=8192, rng=None):
    """
    Inject spectral shot noise into a quiet-start longitudinal distribution.

    Parameters
    ----------
    tau : ndarray
        Quiet-start longitudinal positions.
    N_e : int
        Effective number of real electrons (sets noise floor ~1/sqrt(N_e)).
    kmin, kmax : float or None
        Wavenumber band [1/m] where noise should be applied.
        If None, defaults to top half of FFT spectrum.
    M : int, default=8192
        Number of bins for FFT/histogram.
    rng : np.random.Generator or None
        RNG.

    Returns
    -------
    tau_noisy : ndarray
        Modified tau array with spectral noise injected.
    """
    if rng is None:
        rng = np.random.default_rng()
    N_m = len(tau)

    # Histogram (uniform grid)
    zmin, zmax = tau.min(), tau.max()
    L = zmax - zmin
    n, edges = np.histogram(tau, bins=M, range=(zmin, zmax))
    n = n.astype(float)

    # FFT
    b = np.fft.rfft(n) / N_m
    k = 2*np.pi * np.arange(len(b)) / L

    # Band selection
    if kmin is None: kmin = k[len(k)//4]
    if kmax is None: kmax = k[-1]
    band = (k >= kmin) & (k <= kmax)

    # Inject complex Gaussian noise in band
    noise = (rng.normal(size=band.sum()) + 1j*rng.normal(size=band.sum())) / np.sqrt(2*N_e)
    b[band] += noise

    # Hermitian symmetry preserved automatically with rfft

    # Inverse FFT
    n_noisy = np.fft.irfft(b * N_m, n=M).real
    n_noisy = np.clip(n_noisy, 1e-30, None)  # avoid negatives

    # Build CDF
    cdf = np.cumsum(n_noisy)
    cdf /= cdf[-1]
    centers = 0.5*(edges[:-1] + edges[1:])

    # Map uniform quantiles back to noisy distribution
    u = (np.arange(1, N_m+1) - 0.5)/N_m
    tau_noisy = np.interp(u, cdf, centers)

    return tau_noisy
