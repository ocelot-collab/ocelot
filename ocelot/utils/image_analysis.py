import numpy as np
from scipy import ndimage
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter, uniform_filter
from scipy.optimize import curve_fit


def simple_filter_and_mask(image, sigma=5, threshold=0.0):
    """
    The function returns the same image but with the background cut (masked).
    First, the function filters the original image with a Gaussian filter (sigma=5 by default).
    Second, a mask is generated corresponding to areas where the intensity of the filtered image is less than the 'threshold'.
    Lastly, the mask is applied to the original image by zeroing out the masked area.

    :param image: 2D numpy array - particle density distribution
    :param sigma: gaussian filter 5 pixes by default
    :param threshold: below threshold the area will be set to zero
    :return: masked image with the same sizes
    """

    Hg = ndimage.gaussian_filter(image, sigma=sigma, truncate=2)
    inds_mask_neg = ((Hg - np.max(Hg) * threshold) < 0).nonzero()
    for i in range(10):
        img_neg = np.copy(image)
        img_neg[:] = 1
        img_neg[inds_mask_neg] = 0.0
        proj = np.sum(img_neg, axis=0)
        if np.all(proj[0:20] == 0) and np.all(proj[-20:] == 0):
            break
        threshold += 0.01
        print(f"threshold = {threshold}")
        inds_mask_neg = ((Hg - np.max(Hg) * threshold) < 0).nonzero()

    image[inds_mask_neg] = 0.0

    return image


def simple_filter_and_mask_v2(image, sigma=5, threshold=0.0, max_iterations=10):
    """
    The function returns the same image but with the background cut (masked).
    1) The function filters the original image with a Gaussian filter (sigma=5 by default).
    2) A mask is generated for areas where intensity < threshold * max(Hg).
    3) threshold is iteratively adjusted so that the first and last 20 columns become fully masked.
    4) The mask is then applied in-place (zeroing out masked areas) to the original image.

    :param image: 2D numpy array - particle density distribution
    :param sigma: Gaussian filter width (5 pixels by default)
    :param threshold: Initial threshold fraction of the max intensity in the filtered image
    :param max_iterations: Maximum number of increments to find a threshold that masks edges
    :return: The modified (in-place) masked image
    """

    # 1) Filter the image
    Hg = ndimage.gaussian_filter(image, sigma=sigma, truncate=2)
    max_val = np.max(Hg)

    # 2) Iteratively adjust threshold until the edges are fully masked or we hit max_iterations
    for _ in range(max_iterations):
        mask = (Hg < threshold * max_val)
        # Check if the first 20 and last 20 columns are fully masked
        if mask[:, :20].all() and mask[:, -20:].all():
            break
        threshold += 0.01
        # Optionally: print the updated threshold if you want to track it
        # print(f"threshold = {threshold}")

    # 3) Create the final mask after adjustments
    mask = (Hg < threshold * max_val)

    # 4) Zero out the masked area in-place
    image[mask] = 0.0

    return image


def get_image_proj(image):
    """
    function gets projection on horizontal asis
    :param image: 2D numpy array
    :return: 1D numpy array - projection on horizontal axis
    """
    proj = np.sum(image, axis=0)
    return proj


def find_peaks_regions(y):
    """
    Finds peaks in the 1D array y and returns the indices of the regions around the peaks.
    """
    # Find peaks
    peaks, properties = find_peaks(y, prominence=0.05 * np.max(y))

    # You can adjust 'prominence' to suit your data
    # Alternatively, use 'height', 'distance', 'width', etc.

    regions = []
    for peak in peaks:
        # Find the left and right bases of the peak
        left_base = properties['left_bases'][np.where(peaks == peak)[0][0]]
        right_base = properties['right_bases'][np.where(peaks == peak)[0][0]]
        regions.append((left_base, right_base))
    #print("find_peaks_regions ", regions, len(regions))
    if len(regions) == 1:
        return regions[0]
    elif len(regions) == 0:
        return regions
    dists = [r[1] - r[0] for r in regions]
    ind = np.argmax(dists)
    return regions[ind]


def roi_1d(y, **kwargs):
    threshold = kwargs.get("threshold", 0)
    lmargin = kwargs.get("lmargin", 0)
    rmargin = kwargs.get("rmargin", 0)
    #regions = roi_1d_simple(y, threshold=threshold, lmargin=lmargin, rmargin=rmargin)
    regions = find_peaks_regions(y)
    return regions


def crop_1d(x, y, threshold=0.01, lmargin=0, rmargin=0):
    """
    function crops current distribution below thresholds + margins
    :param x:
    :param y:
    :param threshold: find indices when yi == max(y)*threshold
    :param lmargin: add number of pixes (indices) on the left side
    :param rmargin: add number of pixes (indices) on the right side
    :return: x[indx1:indx2], y[indx1:indx2]
    """
    indx1, indx2 = roi_1d(y, threshold=threshold, lmargin=lmargin, rmargin=rmargin)
    x_new = x[indx1:indx2]
    y_new = y[indx1:indx2]
    return x_new, y_new


def roi_2d(image, threshold=0.01, hor_margin=[5, 5], ver_margin=[5, 5]):
    """
    function finds indices of image below threshold.
    :param image:
    :param threshold:
    :param hor_margin: [left margin, right margin]
    :param ver_margin: [bottoь margin, top margin]
    :return:
    """
    h_proj = get_image_proj(image)
    h_proj = gaussian_filter(h_proj, sigma=2)
    #h_proj = uniform_filter(h_proj, size=2)
    ix1, ix2 = roi_1d(y=h_proj, threshold=threshold, lmargin=hor_margin[0], rmargin=hor_margin[1])

    v_proj = get_image_proj(image.T)
    v_proj = gaussian_filter(v_proj, sigma=2)
    #v_proj = uniform_filter(v_proj, size=2)
    iy1, iy2 = roi_1d(y=v_proj, threshold=threshold, lmargin=ver_margin[0], rmargin=ver_margin[1])

    return ix1, ix2, iy1, iy2


def crop_2d(image, threshold=0.01, hor_margin=[5, 5], ver_margin=[5, 5]):
    """
    function crop image

    :param image:
    :param threshold:
    :param hor_margin:
    :param ver_margin:
    :return:
    """
    ix1, ix2, iy1, iy2 = roi_2d(image, threshold=threshold, hor_margin=hor_margin, ver_margin=ver_margin)
    #print("crop = ", ix1, ix2, iy1, iy2)
    return image[iy1:iy2, ix1:ix2]


def gauss_fit(x, y, sigma_estm=None, mu_estm=0.0):
    # TODO: create unit tests
    """
    Fit a Gaussian function to the given data.

    Parameters
    ----------
    x : array_like
        Independent variable data (e.g., positions).
    y : array_like
        Dependent variable data (e.g., intensities at positions x).
    sigma_estm : float, optional
        Initial estimate of the standard deviation (sigma) of the Gaussian.
        If None, it is estimated as one-fourth of the x-range.
    mu_estm : float, optional
        Initial estimate of the mean (mu) of the Gaussian. Default is 0.0.

    Returns
    -------
    A : float
        Amplitude of the fitted Gaussian function.
    mu : float
        Mean (center position) of the fitted Gaussian function.
    sigma : float
        Standard deviation (spread or width) of the fitted Gaussian function.

    Notes
    -----
    The Gaussian function is defined as:
        f(x) = A * exp(- (x - mu)^2 / (2 * sigma^2))
    where:
        - A is the amplitude,
        - mu is the mean,
        - sigma is the standard deviation.
    """
    if sigma_estm is None:
        sigma_estm = (np.max(x) - np.min(x)) / 4.0

    def gauss(x, A, mu, sigma):
        return A * np.exp(- (x - mu) ** 2 / (2.0 * sigma ** 2))

    # Initial guess for the fitting coefficients: [A, mu, sigma]
    p0 = [np.max(y), mu_estm, sigma_estm]

    coeff, _ = curve_fit(gauss, x, y, p0=p0)
    A, mu, sigma = coeff
    return A, mu, sigma


def get_slice_parameters(image, gauss_fit_enabled=True):
    # TODO: create unit tests
    """
    Calculate slice parameters of an image spot, optionally fitting Gaussian profiles.

    Parameters
    ----------
    image : 2D array_like
        The input image containing the spot (intensity values).
    gauss_fit_enabled : bool, optional
        If True, perform Gaussian fitting on each column (slice) of the image.
        Default is True.

    Returns
    -------
    x : ndarray
        Array of x-axis coordinates (column indices).
    y : ndarray
        Array of y-axis coordinates (row indices).
    proj_x : ndarray
        Projection of the image along the y-axis (sum over rows for each column).
    proj_y : ndarray
        Projection of the image along the x-axis (sum over columns for each row).
    slice_energy : ndarray
        Weighted mean position (energy) along the y-axis for each column.
    slice_energy_spread : ndarray
        Standard deviation (spread) of the energy along the y-axis for each column.
    x_sigma : float
        Standard deviation of the spot along the x-axis.

    Notes
    -----
    - The function computes the centroid and spread of the intensity distribution
      for each column (x-axis) of the image.
    - If `gauss_fit_enabled` is True, it refines the centroid and spread by fitting
      a Gaussian function to the intensity profile of each column.
    """
    ny, nx = image.shape
    x = np.arange(nx)
    y = np.arange(ny)

    # Projections along x and y axes
    proj_x = np.sum(image, axis=0)  # Sum over rows (y-axis) for each column
    proj_y = np.sum(image, axis=1)  # Sum over columns (x-axis) for each row

    total_intensity = np.sum(proj_y)

    # Mean position along x-axis
    x_mean = np.sum(x * proj_x) / total_intensity

    # Standard deviation along x-axis
    x_variance = np.sum((x - x_mean) ** 2 * proj_x) / total_intensity
    x_sigma = np.sqrt(x_variance)

    # Weighted mean position (slice energy) along y-axis for each column
    slice_energy = np.sum(image * y[:, np.newaxis], axis=0) / proj_x

    # Compute the squared differences from the mean for each column
    y_diff_squared = (y[:, np.newaxis] - slice_energy) ** 2

    # Weighted sum of squared differences (variance) for each column
    y_variance = np.sum(image * y_diff_squared, axis=0) / proj_x

    # Standard deviation (spread) along y-axis for each column
    slice_energy_spread = np.sqrt(y_variance)

    if gauss_fit_enabled:
        # Refine slice energy and spread by Gaussian fitting
        for i in range(nx):
            try:
                A, mu, sigma = gauss_fit(
                    y,
                    image[:, i],
                    sigma_estm=slice_energy_spread[i],
                    mu_estm=slice_energy[i]
                )
                slice_energy[i] = mu
                slice_energy_spread[i] = sigma
            except RuntimeError as e:
                print(f"Gaussian fit failed for slice {i}: {e}")
            except Exception as e:
                print(f"An error occurred for slice {i}: {e}")

    return x, y, proj_x, proj_y, slice_energy, slice_energy_spread, abs(x_sigma)
