import numpy as np
from scipy import ndimage
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter, median_filter
from scipy.optimize import curve_fit


def simple_filter_and_mask(image, sigma=5, threshold=0.0):
    """
    Applies a Gaussian filter to the image and masks (sets to zero) pixels below a dynamic threshold.

    Parameters:
        image (np.ndarray): Input image.
        sigma (float): Standard deviation for Gaussian kernel.
        threshold (float): Initial relative threshold for masking.

    Returns:
        np.ndarray: Image with masked low-intensity background.
    """
    Hg = ndimage.gaussian_filter(image, sigma=sigma, truncate=2)
    inds_mask_neg = ((Hg - np.max(Hg) * threshold) < 0).nonzero()
    for i in range(10):
        img_neg = np.ones_like(image)
        img_neg[inds_mask_neg] = 0.0
        proj = np.sum(img_neg, axis=0)
        if np.all(proj[:20] == 0) and np.all(proj[-20:] == 0):
            break
        threshold += 0.01
        inds_mask_neg = ((Hg - np.max(Hg) * threshold) < 0).nonzero()
    image[inds_mask_neg] = 0.0
    return image


def simple_filter_and_mask_v2(image, sigma=5, threshold=0.0, max_iterations=10):
    """
    Enhanced version of the filtering and masking function with iterative threshold adjustment.

    Parameters:
        image (np.ndarray): Input image.
        sigma (float): Standard deviation for Gaussian filter.
        threshold (float): Initial relative threshold.
        max_iterations (int): Number of iterations to adjust the threshold.

    Returns:
        np.ndarray: Image with background masked out.
    """
    Hg = ndimage.gaussian_filter(image, sigma=sigma, truncate=2)
    max_val = np.max(Hg)
    for _ in range(max_iterations):
        mask = (Hg < threshold * max_val)
        if mask[:, :20].all() and mask[:, -20:].all():
            break
        threshold += 0.01
    image[Hg < threshold * max_val] = 0.0
    return image


def get_image_proj(image):
    """
    Computes the horizontal projection (sum along vertical axis) of the image.

    Parameters:
        image (np.ndarray): Input image.

    Returns:
        np.ndarray: 1D array of summed pixel values along columns.
    """
    return np.sum(image, axis=0)


def find_widest_peak_region(y, threshold=0.05):
    """
    Finds prominent peaks and returns the widest base region among them.

    The function detects all prominent peaks in the 1D input signal `y` using scipy's `find_peaks` (with a prominence threshold of 5% of the maximum signal).
    For each detected peak, the function determines the left and right bases (the bounds where the peak starts and ends).
    If only one peak is found, its base region is returned as a tuple.
    If multiple peaks are found, the region corresponding to the *widest* base (largest distance between left and right) is returned.
    If no peaks are found, an empty list is returned.

    Parameters:
        y (np.ndarray): 1D input signal.

    Returns:
        tuple or list: (left_index, right_index) of the widest region if peaks are found,
                       or an empty list if none are detected.

    Note:
        The region is chosen by the widest base, not the highest peak.
    """
    peaks, properties = find_peaks(y, prominence=threshold * np.max(y))
    regions = []
    for peak in peaks:
        left = properties['left_bases'][np.where(peaks == peak)[0][0]]
        right = properties['right_bases'][np.where(peaks == peak)[0][0]]
        regions.append((left, right))
    if len(regions) == 1:
        return regions[0]
    elif len(regions) == 0:
        return 0, len(y)
    return regions[np.argmax([r[1] - r[0] for r in regions])]


def roi_1d(y, threshold=0.05, lmargin=0, rmargin=0, **kwargs):
    """
    Returns region of interest (ROI) from a 1D signal using peak detection.

    Parameters:
        y (np.ndarray): 1D signal.
        **kwargs: Additional arguments passed to `find_peaks_regions`.

    Returns:
        tuple: Indices defining the ROI.
    """
    indx1, indx2 = find_widest_peak_region(y, threshold=threshold)
    indx1 = max(0, indx1 - lmargin)
    indx2 = min(len(y), indx2 + 1 + rmargin)
    return indx1, indx2


def crop_1d(x, y, threshold=0.01, lmargin=0, rmargin=0):
    """
    Crops 1D data based on signal region of interest.

    Parameters:
        x (np.ndarray): X-axis values.
        y (np.ndarray): Y-axis values (signal).
        threshold (float): Threshold for ROI detection.
        lmargin (int): Left margin to extend ROI.
        rmargin (int): Right margin to extend ROI.

    Returns:
        tuple: Cropped (x, y) arrays.
    """
    indx1, indx2 = roi_1d(y, threshold=threshold, lmargin=lmargin, rmargin=rmargin)
    return x[indx1:indx2], y[indx1:indx2]


def roi_2d(image, threshold=0.01, hor_margin=None, ver_margin=None, gauss_filter_sigma=2):
    """
    Computes 2D region of interest using horizontal and vertical projections.

    Parameters:
        image (np.ndarray): 2D input image.
        threshold (float): Threshold for ROI detection.
        hor_margin (list or tuple): [left, right] horizontal margins. Default [5, 5].
        ver_margin (list or tuple): [top, bottom] vertical margins. Default [5, 5].
        gauss_filter_sigma (float): Standard deviation for Gaussian smoothing of projections.

    Returns:
        tuple: Indices (ix1, ix2, iy1, iy2) for cropping the image as image[iy1:iy2, ix1:ix2].
    """
    if hor_margin is None:
        hor_margin = [5, 5]
    if ver_margin is None:
        ver_margin = [5, 5]

    h_proj = get_image_proj(image)
    if gauss_filter_sigma:
        h_proj = gaussian_filter(h_proj, sigma=gauss_filter_sigma)
    ix1, ix2 = roi_1d(h_proj, threshold=threshold, lmargin=hor_margin[0], rmargin=hor_margin[1])
    v_proj = get_image_proj(image.T)
    if gauss_filter_sigma:
        v_proj = gaussian_filter(v_proj, sigma=gauss_filter_sigma)
    iy1, iy2 = roi_1d(v_proj, threshold=threshold, lmargin=ver_margin[0], rmargin=ver_margin[1])

    return ix1, ix2, iy1, iy2


def crop_2d(image, threshold=0.01, hor_margin=[5, 5], ver_margin=[5, 5], gauss_filter_sigma=2):
    """
    Crops 2D image to region of interest based on projections.

    Parameters:
        image (np.ndarray): Input image.
        threshold (float): Threshold for ROI detection.
        hor_margin (list): Horizontal margins for cropping.
        ver_margin (list): Vertical margins for cropping.

    Returns:
        np.ndarray: Cropped image.
    """
    ix1, ix2, iy1, iy2 = roi_2d(image, threshold=threshold, hor_margin=hor_margin, ver_margin=ver_margin, gauss_filter_sigma=gauss_filter_sigma)
    return image[iy1:iy2, ix1:ix2]


def gauss_fit(x, y, sigma_estm=None, mu_estm=0.0):
    """
    Fits a Gaussian function to 1D data.

    Parameters:
        x (np.ndarray): X values.
        y (np.ndarray): Y values.
        sigma_estm (float): Initial estimate for standard deviation.
        mu_estm (float): Initial estimate for mean.

    Returns:
        list: Fitted parameters [A, mu, sigma].
    """
    if sigma_estm is None:
        sigma_estm = (np.max(x) - np.min(x)) / 4.0

    def gauss(x, A, mu, sigma):
        return A * np.exp(- (x - mu) ** 2 / (2.0 * sigma ** 2))

    p0 = [np.max(y), mu_estm, sigma_estm]
    non_zero_count = np.count_nonzero(y)
    if non_zero_count < 3:
        return p0
    try:
        coeff, _ = curve_fit(gauss, x, y, p0=p0, maxfev=500)
        return coeff
    except RuntimeError as e:
        print(f"Fit failed: {e}")
        return p0  # fallback if fitting fails


def remove_column_spikes_vectorized(image, threshold=5.0):
    """
    Removes isolated vertical spikes ("hot pixels") from each column of a 2D image.

    A spike is defined as a pixel that deviates from the average of its immediate
    vertical neighbors (above and below) by more than a specified threshold.

    Parameters
    ----------
    image : np.ndarray
        2D input image (height × width), typically with columns representing a spatial or temporal axis.
    threshold : float, optional
        Multiplier for detecting outliers. A pixel is replaced if its deviation from the
        vertical average exceeds `threshold * |local average|`. Default is 5.0.

    Returns
    -------
    cleaned : np.ndarray
        The filtered image with vertical spike noise removed.

    Notes
    -----
    - The first and last rows are not modified, as they lack two neighbors.
    - This function is vectorized and much faster than a nested loop implementation.

    Examples
    --------
    >>> filtered_image = remove_column_spikes_vectorized(raw_image, threshold=4.0)
    """
    image = np.asarray(image)
    cleaned = image.copy()

    # Shifted views along vertical (row) axis
    top = image[:-2, :]
    center = image[1:-1, :]
    bottom = image[2:, :]

    # Local average of top and bottom neighbors
    local_avg = 0.5 * (top + bottom)

    # Compute deviation
    deviation = np.abs(center - local_avg)
    mask = deviation > threshold * np.abs(local_avg)

    # Apply correction
    cleaned[1:-1, :][mask] = local_avg[mask]

    return cleaned


def get_slice_parameters(image, gauss_fit_enabled=True):
    """
    Calculates image slice parameters like mean, sigma, and energy distribution.

    Parameters
    ----------
    image : np.ndarray
        2D image (e.g., from TDS or spectrometer).
    gauss_fit_enabled : bool, optional
        If True, fit Gaussian to each column for better precision.

    Returns
    -------
    x, y : np.ndarray
        Axis vectors.
    proj_x, proj_y : np.ndarray
        Projections along each axis.
    slice_energy : np.ndarray
        Mean position of slices along y.
    slice_energy_spread : np.ndarray
        RMS spread of each slice along y.
    x_sigma : float
        Global RMS size along x.
    """
    image = remove_column_spikes_vectorized(image, threshold=5.0)
    ny, nx = image.shape
    x = np.arange(nx)
    y = np.arange(ny)
    proj_x = np.sum(image, axis=0)
    safe_proj_x = np.where(proj_x == 0, 1e-12, proj_x)
    proj_y = np.sum(image, axis=1)
    total_intensity = np.sum(proj_y)
    x_mean = np.sum(x * proj_x) / total_intensity
    x_variance = np.sum((x - x_mean) ** 2 * proj_x) / total_intensity
    x_sigma = np.sqrt(x_variance)

    slice_energy = np.sum(image * y[:, np.newaxis], axis=0) / safe_proj_x
    y_diff_squared = (y[:, np.newaxis] - slice_energy) ** 2
    y_variance = np.sum(image * y_diff_squared, axis=0) / safe_proj_x
    slice_energy_spread = np.sqrt(y_variance)

    if gauss_fit_enabled:
        for i in range(nx):
            A, mu, sigma = gauss_fit(y, image[:, i], sigma_estm=slice_energy_spread[i], mu_estm=slice_energy[i])
            slice_energy[i] = mu
            slice_energy_spread[i] = sigma
    return x, y, proj_x, proj_y, slice_energy, slice_energy_spread, abs(x_sigma)


def center_images(images):
    """
    Centers a list of images based on their center of mass.

    Parameters:
        images (list of np.ndarray): List of 2D images.

    Returns:
        list of np.ndarray: Centered images.
    """
    if len(images) <= 1:
        return images

    jj, ii = [], []
    nyy, nxx = [], []
    for img in images:
        j, i = ndimage.center_of_mass(img)
        ny, nx = img.shape
        jj.append(int(j))
        ii.append(int(i))
        nyy.append(ny)
        nxx.append(nx)

    j0 = np.array(jj) - np.min(jj)
    j1 = np.array(jj) + np.min(np.array(nyy) - np.array(jj))
    i0 = np.array(ii) - np.min(ii)
    i1 = np.array(ii) + np.min(np.array(nxx) - np.array(ii))

    return [images[i][j0[i]:j1[i], i0[i]:i1[i]] for i in range(len(images))]


def averaging_images(images):
    """
    Averages a list of images.

    Parameters:
        images (list of np.ndarray): List of 2D images.

    Returns:
        np.ndarray: Averaged image.
    """
    return np.mean(images, axis=0) if len(images) > 1 else images[0]


def shearing(image, shear):
    """
    Applies horizontal shearing transformation to an image.

    Parameters:
        image (np.ndarray): Input image.
        shear (float): Shear angle in radians.

    Returns:
        np.ndarray: Sheared image.
    """
    X, Y = np.meshgrid(np.arange(image.shape[1]), np.arange(image.shape[0]))
    coords = np.vstack((X.flatten(), Y.flatten()))
    transform = np.array([[1, np.tan(shear)], [0, 1]])
    new_coords = np.dot(transform, coords)
    return ndimage.map_coordinates(image, coordinates=[new_coords[1], new_coords[0]], order=1).reshape(image.shape)