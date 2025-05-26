import numpy as np
import pytest

from ocelot.utils.image_analysis import (
    simple_filter_and_mask_v2,
    get_image_proj,
    find_widest_peak_region,
    roi_1d,
    crop_1d,
    roi_2d,
    crop_2d,
    gauss_fit,
    remove_column_spikes_vectorized,
    get_slice_parameters,
    center_images,
    averaging_images,
)

def test_simple_filter_and_mask_v2():
    img = np.zeros((50, 50))
    img[20:30, 20:30] = 10
    img_copy = img.copy()
    masked = simple_filter_and_mask_v2(img_copy, sigma=2, threshold=0.5, max_iterations=3)
    # Most of the background should be zero
    assert np.count_nonzero(masked) <= 100
    # The central square should not be zero
    assert np.max(masked) > 0

def test_get_image_proj():
    img = np.ones((4, 4))
    proj = get_image_proj(img)
    np.testing.assert_array_equal(proj, [4, 4, 4, 4])

def test_find_widest_peak_region():
    y = np.zeros(50)
    y[10] = 3
    y[20] = 1
    y[21] = 0.5
    y[30] = 2
    region = find_widest_peak_region(y)
    assert isinstance(region, tuple)
    assert region[0] <= 19 and region[1] >= 22

def test_roi_1d():
    y = np.zeros(50)
    y[20:25] = 1
    roi = roi_1d(y)
    assert isinstance(roi, tuple)
    assert roi[1] > roi[0]

def test_crop_1d():
    x = np.arange(10)
    y = np.zeros(10)
    y[3:7] = 1

    # Default margins: should include both zeros at edges if implemented as discussed
    x_crop, y_crop = crop_1d(x, y)
    assert x_crop.shape == y_crop.shape
    # Should include left and right zeros if roi returns base [2, 7]
    assert np.all(y_crop[1:-1] == 1)
    assert y_crop[0] == 0 and y_crop[-1] == 0
    assert len(x_crop) == 6  # indices 2 to 7 (inclusive): 2,3,4,5,6,7

    # With left margin
    x_crop2, y_crop2 = crop_1d(x, y, lmargin=1)
    assert x_crop2[0] == 1  # extended to the left by 1
    assert y_crop2[0] == 0

    # With right margin
    x_crop3, y_crop3 = crop_1d(x, y, rmargin=1)
    assert x_crop3[-1] == 8  # extended to the right by 1
    assert y_crop3[-1] == 0

    # With both margins, test cropping does not go out of bounds
    x_crop4, y_crop4 = crop_1d(x, y, lmargin=5, rmargin=5)
    assert x_crop4[0] == 0 and x_crop4[-1] == 9  # capped at array bounds

    # All output slices are same length and correct shape
    assert all(a.shape == b.shape for a, b in [
        (x_crop, y_crop),
        (x_crop2, y_crop2),
        (x_crop3, y_crop3),
        (x_crop4, y_crop4),
    ])


def test_roi_2d():
    # Case 1: Central block in image
    img = np.zeros((20, 20))
    img[5:15, 7:13] = 1  # central block (10 rows, 6 columns)
    ix1, ix2, iy1, iy2 = roi_2d(img, threshold=0.1, hor_margin=[0, 0], ver_margin=[0, 0], gauss_filter_sigma=0)
    # Block is from x=7 to 12 (inclusive), y=5 to 14 (inclusive), so indices should be [7,13) and [5,15)
    assert (ix1, ix2) == (6, 14)
    assert (iy1, iy2) == (4, 16)
    cropped = img[iy1:iy2, ix1:ix2]
    assert cropped.shape == (12, 8)
    assert np.all(cropped[1:-1, 1:-1] == 1)


    # Case 2: Left/top shifted block
    img = np.zeros((20, 20))
    img[5:15, 7:13] = 1  # central block (10 rows, 6 columns)
    ix1, ix2, iy1, iy2 = roi_2d(img, threshold=0.1, hor_margin=[2, 3], ver_margin=[3, 4], gauss_filter_sigma=0)
    # Block is from x=7 to 12 (inclusive), y=5 to 14 (inclusive), so indices should be [7,13) and [5,15)
    assert (ix1, ix2) == (4, 17)
    assert (iy1, iy2) == (1, 20)
    cropped = img[iy1:iy2, ix1:ix2]
    assert cropped.shape == (19, 13)

    # Case 3: Add margins
    img3 = np.zeros((30, 30))
    img3[10:20, 15:22] = 3
    ix1, ix2, iy1, iy2 = roi_2d(img3, threshold=0.1, hor_margin=[2, 3], ver_margin=[1, 4], gauss_filter_sigma=0)
    assert ix1 == 12
    assert ix2 == 26
    assert iy1 == 8
    assert iy2 == 25
    cropped3 = img3[iy1:iy2, ix1:ix2]
    # Should be larger due to margins
    assert cropped3.shape == (17, 14)
    # Central region should be all 3
    assert np.all(cropped3[2:-5, 3:-4] == 3)

    # Case 4: All zeros input should cover the whole image
    img4 = np.zeros((15, 18))
    ix1, ix2, iy1, iy2 = roi_2d(img4, threshold=0.1, hor_margin=[0,0], ver_margin=[0,0])
    assert ix1 == 0 and ix2 == 18
    assert iy1 == 0 and iy2 == 15


def test_crop_2d():
    img = np.zeros((20, 20))
    img[5:15, 6:16] = 10
    cropped = crop_2d(img, hor_margin=[0,0], ver_margin=[0,0], gauss_filter_sigma=0)
    # The result should be square and match the nonzero region
    assert cropped.shape[0] == cropped.shape[1]
    assert np.all(cropped[1:-1, 1:-1] == 10)

def test_gauss_fit():
    x = np.linspace(-5, 5, 100)
    y = np.exp(-(x-1)**2 / (2*2.0**2))
    params = gauss_fit(x, y, sigma_estm=2.0, mu_estm=1.0)
    A, mu, sigma = params
    assert np.isclose(mu, 1.0, atol=0.1)
    assert np.isclose(sigma, 2.0, atol=0.1)

def test_remove_column_spikes_vectorized():
    img = np.ones((10, 10))
    img[5, 3] = 100
    cleaned = remove_column_spikes_vectorized(img, threshold=4)
    assert abs(cleaned[5, 3] - 1) < 1e-6  # The spike should be removed

def test_get_slice_parameters():
    img = np.zeros((20, 10))
    img[10, :] = 5
    img[11, :] = 4
    result = get_slice_parameters(img)
    x, y, proj_x, proj_y, slice_energy, slice_energy_spread, x_sigma = result
    assert len(x) == 10 and len(y) == 20
    assert proj_x.shape == (10,)
    assert proj_y.shape == (20,)
    assert slice_energy.shape == (10,)
    assert x_sigma > 0

def test_center_images():
    img1 = np.zeros((10, 10))
    img1[4, 4] = 1
    img2 = np.zeros((10, 10))
    img2[5, 5] = 1
    centered = center_images([img1, img2])
    assert len(centered) == 2
    # Check if the mass is near the center in both images
    for cimg in centered:
        cy, cx = np.unravel_index(np.argmax(cimg), cimg.shape)
        assert abs(cy - cimg.shape[0]//2) <= 1
        assert abs(cx - cimg.shape[1]//2) <= 1

def test_averaging_images():
    img1 = np.ones((5, 5))
    img2 = np.zeros((5, 5))
    avg = averaging_images([img1, img2])
    assert np.allclose(avg, 0.5)

    avg_single = averaging_images([img1])
    assert np.allclose(avg_single, img1)