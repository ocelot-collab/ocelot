"""
Collection of function for dechirper LPS reconstrucion

S.Tomin, 07.2024, DESY
"""

import numpy as np
from scipy import ndimage, interpolate
from scipy.integrate import solve_ivp
from ocelot.common.math_op import invert_cdf
from ocelot.common.globals import *
from ocelot.cpbd.beam import s_to_cur, generate_parray
from ocelot.cpbd.wake3D import Wake, WakeTableDechirperOffAxis
from scipy.optimize import fmin, curve_fit
from scipy.ndimage import gaussian_filter, uniform_filter
import cv2
from scipy.ndimage import label
from scipy.signal import find_peaks


def image2distrib(image, n_particles=200000):
    """
    The function returns the X and Y sample distribution equivalent to the 2D image density distribution.

    _______
    Example:
    x_distr, y_distr = image2distrib(img, n_particles=2000000)
    show_density(x_distr, y_distr)
    -------

    :param image: 2D numpy array - particle density distribution
    :param n_particles: number of particles in distribution
    :return: x and y distribution but number of particles maybe less than n_particles
    """
    nx, ny = np.shape(image)
    npart_in_slice = (np.sum(image, axis=0) * n_particles / np.sum(image)).astype(int)

    y_distr = np.array([])
    x_distr = np.array([])
    for i in range(ny):
        inv_cdf = invert_cdf(y=image[:, i], x=np.arange(nx))
        if np.isnan(inv_cdf(0.5)):
            continue
        y_distr = np.append(y_distr, (inv_cdf(np.random.rand(npart_in_slice[i]))))
        x_distr = np.append(x_distr, np.random.rand(npart_in_slice[i]) - 0.5 + i)
    return x_distr, y_distr


def image2distrib_v2(image, n_particles=200000):
    """
    The function returns the X and Y sample distribution equivalent to the 2D image density distribution.

    _______
    Example:
    x_distr, y_distr = image2distrib(img, n_particles=2000000)
    show_density(x_distr, y_distr)
    -------

    :param image: 2D numpy array - particle density distribution
    :param n_particles: number of particles in distribution
    :return: x and y distribution but number of particles maybe less than n_particles
    """
    nx, ny = np.shape(image)
    proj = get_image_proj(image)
    inv_cdf = invert_cdf(y=proj, x=np.arange(len(proj)))
    x_distr = inv_cdf(np.random.rand(n_particles))

    y_distr = np.zeros(len(x_distr))
    for i in range(ny):
        inv_cdf = invert_cdf(y=image[:, i], x=np.arange(nx))
        idx = np.where((i <= x_distr) & (x_distr < i + 1))[0]
        n_slice = len(idx)
        if np.isnan(inv_cdf(0.5)):
            continue
        y_distr[idx] = inv_cdf(np.random.rand(n_slice))
    return x_distr, y_distr


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


def current_processing(x, y, nbins=250, threshold=0.01, normilize=True):
    """
    The function processes the beam current (or density distribution or screen projection).
    First, it crops the current by setting an appropriate threshold.
    Second, a new X-axis is created corresponding to the newly cropped area and length equal to nbins.
    The current in new points will be obtained by the linear interpolation
    Third, the current is normalized to achieve an integral equal to 1.

    :param x: horizontal axis
    :param y: current or density or screen projection
    :param nbins: length of new current array
    :param threshold: below threshold the area will be cut
    :return: x_new, y_new - with length equal  nbins
    """

    if threshold > 0:
        zero_crossings = np.where(np.diff(np.sign(y - max(y) * threshold)))[0]
        indx1 = zero_crossings[0]
        indx2 = zero_crossings[-1]
        x = x[indx1:indx2]
        y = y[indx1:indx2]

    x_proj = np.linspace(x[0], x[-1], num=nbins)
    y_proj = np.interp(x_proj, x, y)

    y_proj[0] = 0
    y_proj[-1] = 0
    if normilize:
        y_proj = y_proj / np.trapz(y_proj, x_proj)

    x_proj -= x_proj[0]
    return x_proj, y_proj


def direct_transform(x_crisp, y_crisp, x_proj_img, y_proj_img):
    """
    Function finds a transformation to translate one density distribution (CRISP current) to another (screen projection).

    :param x_crisp: CRISP current coordinates
    :param y_crisp: CRISP current
    :param x_proj_img: screen projection coordinates
    :param y_proj_img: screen projection
    :return: x_transform, y_transform - function to transform one distribution to another
    """

    dx_crisp = x_crisp[1] - x_crisp[0]
    dx_proj = x_proj_img[1] - x_proj_img[0]

    crisp_int = np.cumsum((y_crisp[1:] + y_crisp[:-1]) / 2 * dx_crisp)
    img_proj_int = np.cumsum((y_proj_img[1:] + y_proj_img[:-1]) / 2 * dx_proj)

    # Transformation
    # we will go from the end of the CRISP current and will move from the right to the left.
    x_transf = [0]
    y_ytansf = [0]
    n = len(x_crisp)
    for i, Ii in enumerate(crisp_int):
        x_transf.append(x_crisp[i])

        indx = np.argmin(np.abs(img_proj_int - Ii))
        dx_proj = x_proj_img[indx]
        y_ytansf.append(dx_proj)
    return x_transf, y_ytansf


def inverse_transform(x_transf, y_ytansf):
    """
    function invert transformation and return inverted transform function

    :param x_transf:
    :param y_ytansf:
    :return: f_inverse function
    """
    f_inverse_transform = interpolate.interp1d(y_ytansf, x_transf)
    return f_inverse_transform


def get_image_proj(image):
    """
    function gets projection on horizontal asis
    :param image: 2D numpy array
    :return: 1D numpy array - projection on horizontal axis
    """
    proj = np.sum(image, axis=0)
    return proj


def roi_1d_simple(y, threshold=0.01, lmargin=0, rmargin=0):
    """
    function returns indices of y array when y values are higher thresholds.

    :param y: current like 1D array
    :param threshold: find indices when yi == max(y)*threshold
    :param lmargin: add number of pixes (indices) on the left side
    :param rmargin: add number of pixes (indices) on the right side
    :return: indx1, indx2
    """
    zero_crossings = np.where(np.diff(np.sign(y - max(y) * threshold)))[0]
    print("zero_crossings = ", zero_crossings)
    if len(zero_crossings) == 1:
        zero_crossings = np.insert(zero_crossings, 0, 0)
    indx1 = zero_crossings[0]
    indx2 = zero_crossings[1]
    indx1 = max([0, indx1 - lmargin])
    indx2 = min([len(y), indx2 + rmargin])
    return indx1, indx2
    #print("roi_1d = ", y - max(y) * threshold)
    #print("roi_1d = diff", np.diff(np.sign(y - max(y) * threshold)))
    #print("roi_1d = where", np.where(np.diff(np.sign(y - max(y) * threshold))))


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
    print("find_peaks_regions ", regions, len(regions))
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
    print("crop = ", ix1, ix2, iy1, iy2)
    return image[iy1:iy2, ix1:ix2]


def dchirper_recon(img, t_crisp, y_crisp, img_mask_thresh=0.03, img_sigma=5, img_crop_thresh=0.01, crisp_thresh=0.01,
                   n_particles=200000):
    """
    main function transforms streaked beam image from nonlinear passive streaker to linear streak using as reference
    CRISP current reconstruction

    :param img: image
    :param t_crisp: time CRISP in [fm]
    :param y_crisp: current CRISP in [A]
    :param img_mask_thresh: to mask image background
    :param img_sigma: gaussian filter 5 pixes by default
    :param img_crop_thresh: to crop image
    :param crisp_thresh: crop CRISP current
    :param n_particles: number of particle for distribution
    :return: (x_distrib_transformed, y_distrib), (x_crisp, y_crisp), image_processed
    """
    img = simple_filter_and_mask(img, sigma=img_sigma, threshold=img_mask_thresh)
    img = crop_2d(img, threshold=img_crop_thresh, hor_margin=[0, 0], ver_margin=[0, 0])
    proj = get_image_proj(img)

    x_proj_img, y_proj_img = current_processing(x=np.arange(len(proj)), y=proj, nbins=2000, threshold=0)

    x_crisp, y_crisp = current_processing(x=t_crisp * 1e-15 * speed_of_light * 1e6, y=y_crisp, nbins=250,
                                          threshold=crisp_thresh)

    x_transf, y_ytansf = direct_transform(x_crisp, y_crisp, x_proj_img, y_proj_img)

    inv_transf_x = np.linspace(y_ytansf[0], y_ytansf[-1], num=400)
    f_inverse_transform = inverse_transform(x_transf, y_ytansf)

    # generate distribution out of image
    x_distr, y_distr = image2distrib(img, n_particles=n_particles)

    indeces = np.where(np.logical_and(x_distr >= np.min(inv_transf_x), x_distr <= np.max(inv_transf_x)))

    x_img_dist_flt = x_distr[indeces]
    y_img_dist_flt = y_distr[indeces]

    x_distr_transformed = f_inverse_transform(x_img_dist_flt)

    crisp_proc = (x_crisp, y_crisp)
    transformed_distr = (x_distr_transformed, y_img_dist_flt)
    return transformed_distr, crisp_proc, img, x_proj_img, f_inverse_transform


def find_solution_with_rk(x_crisp, y_crisp, x_proj_img, y_proj_img, der0, num=500):
    """

    :param x_crisp: time CRISP in [um]
    :param y_crisp: normalised current CRISP, I(x)/Integral(I(x))
    :param x_proj_img: x projection in [px]
    :param y_proj_img: normalised image projection I(x)/Integral(I(x))
    :param der0: initial derivative
    :param num: number of points
    :return: trasformation x and y(x)
    """
    #  fill_value="extrapolate"
    #f = interpolate.interp1d(x_crisp, y_crisp, fill_value="extrapolate")
    #g = interpolate.interp1d(x_proj_img, y_proj_img, fill_value="extrapolate")

    f = interpolate.interp1d(x_crisp, y_crisp, bounds_error=False, fill_value=(0, 0))
    g = interpolate.interp1d(x_proj_img, y_proj_img, bounds_error=False, fill_value=(0, 0))

    #g = interpolate.CubicSpline(x_proj_img, y_proj_img, bc_type='natural')
    #f = interpolate.CubicSpline(x_crisp, y_crisp, bc_type='natural')
    def func(y, x, f, g):
        """
        f(x) - crisp current
        g(y) - image projection

        dy/dt = func(y, t, ...)
        """
        dxdy = g(y) / f(x) if f(x) != 0 else der0
        return dxdy

    y0 = [0]
    x = np.linspace(x_proj_img[0], x_proj_img[-1], num=num)

    sol = solve_ivp(func, [x[0], x[-1]], y0, t_eval=x, args=(f, g), method="RK23")
    return sol.t, sol.y[0]

def get_largest_contour_roi(contours):
    largest_contour = max(contours, key=cv2.contourArea)
    x, y, w, h = cv2.boundingRect(largest_contour)
    return x, y, w, h

def improved_roi_detection(image):
    img_blur = cv2.GaussianBlur(image, (5, 5), 0)
    _, img_thresh = cv2.threshold(
        img_blur.astype(np.uint8), 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
    kernel = np.ones((5, 5), np.uint8)
    img_close = cv2.morphologyEx(img_thresh, cv2.MORPH_CLOSE, kernel)
    contours, _ = cv2.findContours(
        img_close, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if contours:
        x, y, w, h = get_largest_contour_roi(contours)
        roi = image[y:y+h, x:x+w]
        return roi
    else:
        return None


def detect_roi_with_opencv(image):
    # Normalize and convert to uint8
    image_normalized = cv2.normalize(image, None, 0, 255, cv2.NORM_MINMAX)
    image_uint8 = image_normalized.astype(np.uint8)

    # Apply Gaussian Blur
    img_blur = cv2.GaussianBlur(image_uint8, (5, 5), 0)

    # Apply Otsu's Thresholding
    _, img_thresh = cv2.threshold(
        img_blur, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)

    # Morphological Operations
    kernel = np.ones((3, 3), np.uint8)
    img_morph = cv2.morphologyEx(img_thresh, cv2.MORPH_CLOSE, kernel)
    img_morph = cv2.morphologyEx(img_morph, cv2.MORPH_OPEN, kernel)

    # Find Contours
    contours, _ = cv2.findContours(
        img_morph, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    if contours:
        # Assume the largest contour is the ROI
        largest_contour = max(contours, key=cv2.contourArea)
        # Get bounding rectangle
        x, y, w, h = cv2.boundingRect(largest_contour)
        # Extract ROI
        roi = image[y:y+h, x:x+w]
        return roi, (x, y, w, h)
    else:
        return None, None


def dchirper_recon_RK(img, t_crisp, y_crisp, img_mask_thresh=0.03, img_sigma=3, img_crop_thresh=0.01, crisp_thresh=0.01,
                      n_particles=200000, is_image_processed=False):
    """
    main function transforms streaked beam image from nonlinear passive streaker to linear streak using beam current
     from CRISP reconstruction as the reference



    :param img: image
    :param t_crisp: time CRISP in [fm]
    :param y_crisp: current CRISP in [A]
    :param img_mask_thresh: to mask image background
    :param img_sigma: gaussian filter 5 pixes by default
    :param img_crop_thresh: to crop image
    :param crisp_thresh: crop CRISP current
    :param n_particles: number of particle for distribution
    :param is_image_processed: if True image was processed before, if False it will be processed here
    :return: (x_distrib_transformed, y_distrib), (x_crisp, y_crisp), image_processed, (x_transf_curve, y_transf_curve)
    """
    if not is_image_processed:
        img = simple_filter_and_mask(img, sigma=img_sigma, threshold=img_mask_thresh)

        img = crop_2d(img, threshold=img_crop_thresh, hor_margin=[5, 5], ver_margin=[5, 5])

    proj = get_image_proj(img)

    x_proj_img, y_proj_img = current_processing(x=np.arange(len(proj)), y=proj, nbins=750, threshold=0)

    x_crisp, y_crisp = current_processing(x=t_crisp * 1e-15 * speed_of_light * 1e6, y=y_crisp, nbins=250,
                                          threshold=crisp_thresh)

    # find roughly inital derivative
    intg_crisp = (y_crisp[0] + y_crisp[1]) / 2 * (x_crisp[1] - x_crisp[0])
    intg_img_proj = (y_proj_img[0:20] + y_proj_img[1:21]) / 2 * (x_proj_img[1] - x_proj_img[0])
    indx_min = np.argmin(np.abs(intg_img_proj - intg_crisp))
    der = y_proj_img[indx_min + 1] / y_crisp[1]

    x, y = find_solution_with_rk(x_crisp, y_crisp, x_proj_img, y_proj_img, der0=der, num=500)
    f_transform = interpolate.interp1d(x, y, fill_value="extrapolate")
    # generate distribution out of image
    x_distr, y_distr = image2distrib(img, n_particles=n_particles)

    x_distr_transformed = f_transform(x_distr)

    crisp_proc = np.array([x_crisp, y_crisp])
    transformed_distr = (x_distr_transformed, y_distr)
    return transformed_distr, crisp_proc, img, x_proj_img, f_transform


def get_hist(x, nbins=100):
    hist, bin_edges = np.histogram(x, bins=nbins, weights=None)
    bin_edges = (bin_edges[:-1] + bin_edges[1:]) / 2.
    return bin_edges, hist


def direct_transform_backwards(x_crisp, y_crisp, x_proj_img, y_proj_img):
    # for simplicity we invert both currents and caclulate integrals of both currents
    y_crisp_inv = y_crisp[::-1]
    y_proj_inv = y_proj_img[::-1]

    dx_crisp = x_crisp[1] - x_crisp[0]
    dx_proj = x_proj_img[1] - x_proj_img[0]

    crisp_int = np.cumsum((y_crisp_inv[1:] + y_crisp_inv[:-1]) / 2 * dx_crisp)
    img_proj_int = np.cumsum((y_proj_inv[1:] + y_proj_inv[:-1]) / 2 * dx_proj)

    # Transformation
    # we will go from the end of the CRISP current and will move from the right to the left.
    x_transf = [x_crisp[-1]]
    y_ytansf = [0]
    n = len(x_crisp)
    for i, Ii in enumerate(crisp_int):
        x_transf.append(x_crisp[n - 1 - i])

        indx = np.argmin(np.abs(img_proj_int - Ii))
        dx_proj = x_proj_img[-1] - x_proj_img[len(x_proj_img) - 1 - indx]
        y_ytansf.append(dx_proj)
    return x_transf, y_ytansf


def get_img_from_distr(x, y, nbins_x=800, nbins_y=500):
    """
    function generate image from the particle distribution and also return px sizes.
    Can be useful to simulate reconstruction process.

    :param x: particle coordinates in horizontal direction
    :param y: particle coordinates in vertical direction
    :param nbins_x: number of bins in horizontal direction
    :param nbins_y: number of bins in vertical direction
    :return: img (2d numpy array), x px size, y px size
    """
    H, xedges, yedges = np.histogram2d(x=x, y=y, bins=(nbins_x, nbins_y))
    img = H.T
    y_px_size = (yedges[-1] - yedges[0]) / nbins_y
    x_px_size = (xedges[-1] - xedges[0]) / nbins_x
    return img, x_px_size, y_px_size


def track_1D(distance, R34, s_coord, y_scr_arr, charge=250e-12, energy=14, beta_y=5, emit=1e-6):
    """
    function transforms s coordinates obtained to y coordinates using R34 and dipole component of the corrugated plate wake.
    it is also does convolution of the projection curve with gaussian curve with sigma sqrt(beta8emit/gamma)

    :param distance: distance between passive streaker and corrugated plate
    :param R34: between passive streaker and screen
    :param s_coord: numpy array with s coordinates
    :param y_scr_arr: screen mesh, numpy array where projection calculated
    :param charge: beam charge in C
    :param energy: beam energy in GeV
    :param beta_y: vertical betatron function in m
    :param emit: normilized vertical emittance
    :return: screen projection calculated in y_scr_arr
    """

    wk_tv_kick = WakeTableDechirperOffAxis(b=distance,  # distance from the plate in [m]
                                           a=0.01,  # half gap between plates in [m]
                                           width=0.02,  # width of the corrugated structure in [m]
                                           t=0.25 * 1e-3,  # longitudinal gap in [m]
                                           p=0.5 * 1e-3,  # period of corrugation in [m]
                                           length=5.,  # length of the corrugated structure in [m]
                                           sigma=60e-6,  # characteristic (rms) longitudinal beam size in [m]
                                           orient="horz")  # "horz" or "vert" plate orientation

    ws = Wake()
    ws.w_sampling = 500
    ws.wake_table = wk_tv_kick
    ws.prepare(lat=None)
    B = s_to_cur(s_coord, sigma=0.1 * np.std(s_coord), q0=charge, v=speed_of_light)

    x, Wd = ws.get_dipole_wake(current_profile=B)
    f_transform = interpolate.interp1d(x, R34 * Wd / (energy * 1e9), fill_value="extrapolate")

    y_coord = f_transform(s_coord)

    bin_edges, hist = get_hist(y_coord, nbins=500)

    if emit is not None or beta_y is not None:
        # gaussian filter
        dx = bin_edges[1] - bin_edges[0]
        sigma_y = np.sqrt(emit / energy * m_e_GeV * beta_y)

        gx = np.arange(-3 * sigma_y, 3 * sigma_y, dx)
        gaussian = np.exp(-(gx / sigma_y) ** 2 / 2)

        hist = np.convolve(hist, gaussian, mode="same")
    # normalization
    g_track = interpolate.interp1d(bin_edges, hist, fill_value=(0, 0), bounds_error=False)

    g_tr = g_track(y_scr_arr) / np.trapz(g_track(y_scr_arr), y_scr_arr)

    return g_tr


def track_6D(distance, R, s_coord, p_coord, tws_ws, y_scr_arr, charge):
    """
    function transforms s coordinates obtained to y coordinates using R34 and dipole component of the corrugated plate wake.
    it is also does convolution of the projection curve with gaussian curve with sigma sqrt(beta8emit/gamma)

    :param distance: distance between passive streaker and corrugated plate
    :param R: matrix between passive streaker and screen
    :param s_coord: numpy array with s coordinates in [m]
    :param p_coord: numpy array with s coordinates in [MeV]
    :param tws_ws: Twiss parameters at dechirper position. Twiss object should contain emittances
    :param y_scr_arr: screen mesh, numpy array where projection calculated
    :param charge: beam charge in C
    :return: screen projection calculated in y_scr_arr
    """

    wk_tv_kick = WakeTableDechirperOffAxis(b=distance,  # distance from the plate in [m]
                                           a=0.01,  # half gap between plates in [m]
                                           width=0.02,  # width of the corrugated structure in [m]
                                           t=0.25 * 1e-3,  # longitudinal gap in [m]
                                           p=0.5 * 1e-3,  # period of corrugation in [m]
                                           length=5.,  # length of the corrugated structure in [m]
                                           sigma=60e-6,  # characteristic (rms) longitudinal beam size in [m]
                                           orient="horz")  # "horz" or "vert" plate orientation

    ws = Wake()
    ws.w_sampling = 500
    ws.wake_table = wk_tv_kick
    ws.prepare(lat=None)
    B = s_to_cur(s_coord, sigma=0.1 * np.std(s_coord), q0=charge, v=speed_of_light)

    x, Wd = ws.get_dipole_wake(current_profile=B)

    p_array = generate_parray(tws=tws_ws, charge=charge, nparticles=len(s_coord), chirp=0.0, energy=tws_ws.E,
                              sigma_tau=np.std(s_coord))

    p_array.tau()[:] = s_coord[:]
    p_array.p()[:] = p_coord[:] / (tws_ws.E * 1000)

    f_transform = interpolate.interp1d(x, Wd / (tws_ws.E * 1e9), fill_value="extrapolate")
    p_array.py()[:] += f_transform(s_coord)[:]
    a = np.dot(R, p_array.rparticles)
    p_array.rparticles[:] = a[:]

    bin_edges, hist = get_hist(-p_array.y(), nbins=250)

    # normalization
    g_track = interpolate.interp1d(bin_edges - bin_edges[0], hist, fill_value=(0, 0), bounds_error=False)
    g_tr = g_track(y_scr_arr) / np.trapz(g_track(y_scr_arr), y_scr_arr)

    return g_tr


def find_distance_to_corrugated_plate_6D(s_coord, p_coord, img_proc,
                                         x_px_size, tws_ws_start, R, charge=250e-12, align_peaks=False):
    """

    :param s_coord: s coordinates of the beam before corrugated plate [m]
    :param p_coord: energy coordinates of the beam before corrugated plate [MeV]
    :param img_proc:
    :param x_px_size:
    :param tws_ws_start:
    :param R:
    :param charge:
    :param align_peaks:
    :return:
    """
    img_proj = get_image_proj(img_proc)
    y_proj = np.arange(len(img_proj)) * x_px_size

    y_scr_array = np.linspace(0, y_proj[-1], num=500)

    def err_func(distance, img_proc, x_px_size, y_scr_arr):
        if distance < 50e-6 or distance > 0.002:
            return 100000
        proj = get_image_proj(img_proc)
        g_img = interpolate.interp1d(np.arange(len(proj)) * x_px_size, proj, fill_value=(0, 0), bounds_error=False)
        g_img_n = g_img(y_scr_arr) / np.trapz(g_img(y_scr_arr), y_scr_arr)

        g_track = track_6D(distance, R=R, s_coord=s_coord, p_coord=p_coord, tws_ws=tws_ws_start,
                           y_scr_arr=y_scr_array, charge=charge)
        if align_peaks:
            # Find peaks
            peak_im = np.argmax(g_img_n)
            peak_tr = np.argmax(g_track)

            # Calculate shift
            shift = peak_im - peak_tr
            # interpolate for sift
            g_tr_shift_func = interpolate.interp1d(np.arange(len(g_track)) + shift, g_track, fill_value=(0, 0),
                                                   bounds_error=False)

            # get new array
            i = np.arange(len(g_img_n))
            g_track = g_tr_shift_func(i)

        res = np.sqrt(np.sum((g_track - g_img_n) ** 2))  #+ 10000*np.abs(max(g_track) - max(g_img_n))
        return res

    res = fmin(err_func, x0=800e-6, args=(img_proc, x_px_size, y_scr_array), xtol=0.1, ftol=0.1, maxiter=50)

    return res[0]


def find_distance_to_corrugated_plate(s_coord, img_proc, x_px_size, R34=35.33, charge=250e-12, energy=14, beta_y=5,
                                      emit=1e-6, align_peaks=True):
    """
    function finds distance between corrugated plate and the beam using simplex optimization function.
    Objective function is difference between screen projections of tracked beam image on the screen and measured screen image.

    :param s_coord: s coordinates of the beam before corrugated plate
    :param img_proc: screen image processed (cropped)
    :param x_px_size: pixel size of the screen in time direction
    :param R34: between passive streaker and screen
    :param charge: beam charge in C
    :param energy: beam energy in GeV
    :param beta_y: vertical betatron function in m
    :param emit: normilized vertical emittance
    :param align_peaks: if True that image projection and projection finded from tracking will be aligned wrt peaks
    :return: distance to the plate.
    """
    img_proj = get_image_proj(img_proc)
    y_proj = np.arange(len(img_proj)) * x_px_size

    y_scr_array = np.linspace(0, y_proj[-1], num=500)

    def err_func(distance, img_proc, x_px_size, y_scr_arr):
        if distance < 50e-6 or distance > 0.002:
            return 100000
        proj = get_image_proj(img_proc)
        g_img = interpolate.interp1d(np.arange(len(proj)) * x_px_size, proj, fill_value=(0, 0), bounds_error=False)
        g_img_n = g_img(y_scr_arr) / np.trapz(g_img(y_scr_arr), y_scr_arr)
        g_track = track_1D(distance, R34=R34, s_coord=s_coord * 1e-6, y_scr_arr=y_scr_array, charge=charge,
                           energy=energy,
                           beta_y=beta_y, emit=emit)
        print(distance, np.sqrt(np.sum((g_track - g_img_n) ** 2)), 1000 * np.abs(max(g_track) - max(g_img_n)))
        if align_peaks:
            # Find peaks
            peak_im = np.argmax(g_img_n)
            peak_tr = np.argmax(g_track)

            # Calculate shift
            shift = peak_im - peak_tr
            # interpolate for sift
            g_tr_shift_func = interpolate.interp1d(np.arange(len(g_track)) + shift, g_track, fill_value=(0, 0),
                                                   bounds_error=False)

            # get new array
            i = np.arange(len(g_img_n))
            g_track = g_tr_shift_func(i)

        res = np.sqrt(np.sum((g_track - g_img_n) ** 2)) + 1000 * np.abs(max(g_track) - max(g_img_n))
        return res

    res = fmin(err_func, x0=500e-6, args=(img_proc, x_px_size, y_scr_array), xtol=0.01, ftol=0.01, maxiter=100)

    print("res = ", res)
    return res[0]


def calculate_long_wake(s_coord, dist):
    """
    A simple function to calculate longitudinal wake using dechirper wake table and current profile

    :param s_coord: longitudinal distribution  [m]
    :param dist: distance between beam and plate
    :return: current, wake
    """
    wk_tv_kick = WakeTableDechirperOffAxis(b=dist,  # distance from the plate in [m]
                                           a=0.01,  # half gap between plates in [m]
                                           width=0.02,  # width of the corrugated structure in [m]
                                           t=0.25 * 1e-3,  # longitudinal gap in [m]
                                           p=0.5 * 1e-3,  # period of corrugation in [m]
                                           length=5.,  # length of the corrugated structure in [m]
                                           sigma=60e-6,  # characteristic (rms) longitudinal beam size in [m]
                                           orient="horz")  # "horz" or "vert" plate orientation

    ws = Wake()
    # w_sampling - defines the number of the equidistant sampling points for the one-dimensional
    # wake coefficients in the Taylor expansion of the 3D wake function.
    ws.w_sampling = 500
    ws.wake_table = wk_tv_kick

    ws.prepare(lat=None)

    B = s_to_cur(s_coord, sigma=0.1 * np.std(s_coord), q0=250e-12, v=speed_of_light)

    x_w, Wl = ws.get_long_wake(current_profile=B)

    return (B[:, 0], B[:, 1]), (x_w, Wl)


def subtract_long_wake(s_coord, p_coord, dist):
    """
    Function to subtract energy chirp induced by dechirper from particle distribution

    :param s_coord: longitudinal distribution in [m]
    :param p_coord: energy distribution in [MeV]
    :param dist: distance between beam and corrugated plate in [m]
    :return: distribution after energy chirp subtraction
    """

    (x, I), (x_w, Wl) = calculate_long_wake(s_coord, dist)

    energy_kick = interpolate.interp1d(x_w, Wl * 1e-6, fill_value="extrapolate")  # wake originally is in [eV]
    y = p_coord - energy_kick(s_coord)

    return s_coord, y


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
