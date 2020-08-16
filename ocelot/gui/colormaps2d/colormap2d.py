"""
This script contains the colormap2d functions.

M Meschede 2016
https://github.com/MMesch/cmap_builder
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import map_coordinates

def data2d_to_rgb(data, cmap2d='brightwheel', huenorm=None, huevmin=None,
                  huevmax=None, fill_value=(1., 1., 1.),
                  lightnorm=None, lightvmin=None, lightvmax=None):
    """
    Map 2 parameter 2D data array to rgb values.

    :param data: numpy array with shape (2, nwidth, nheight). The first index
                 corresponds to the hue and the second to the lightness of the
                 colors.
    :param cmap2d: either:
                   numpy array with shape (nwidth, nheight, 4) that contains
                   the 4 rgba values in hue (width) and lightness (height).
                   Can be obtained by a call to get_cmap2d(name).
                   or:
                   name where name is one of the following strings:
                   'brightwheel', 'darkwheel', 'hardwheel', 'newwheel',
                   'smoothwheel', 'wheel'
    :param huenorm: a plt.Normalize() instance that normalizes the hue values.
    :param huevmin: the minimum of the huevalues. Only used if huenorm=None.
    :param huevmax: the maximum of the huevalues. Only used if huenorm=None.
    :param fill_value: e.g. (0., 0., 0.) for values that are not finite as nan
    :param lightnorm: a plt.Normalize() instance that normalizes the lightness
                      values.
    :param lightvmin: the minimum of the lightness values.
                      Only used if lightnorm=None.
    :param lightvmax: the maximum of the lightness values.
                      Only used if lightnorm=None.
    :returns: numpy array with shape (3, nwidth, nheight) that can be plotted
              with imshow
    """
    # make copy of data that ultimately stores the data value in colormap index
    # coordinates (hence the name 'idata')
    data_dim, nrows, ncols = data.shape
    idata = np.copy(data)
    mask = np.all(np.isfinite(idata), axis=0)
    idata[:, ~mask] = 0.

    # normalize data if required
    if huenorm is None:
        if huevmin is None:
            huevmin = np.nanmin(data[0])
        if huevmax is None:
            huevmax = np.nanmax(data[0])
        huenorm = plt.Normalize(huevmin, huevmax, clip=True)

    if lightnorm is None:
        if lightvmin is None:
            lightvmin = np.nanmin(data[1])
        if lightvmax is None:
            lightvmax = np.nanmax(data[1])
        lightnorm = plt.Normalize(lightvmin, lightvmax, clip=True)

    idata[0] = huenorm(idata[0])
    idata[1] = lightnorm(idata[1])

    # get colormap if it is not given by numpy array
    if not isinstance(cmap2d, np.ndarray):
        cmap2d = get_cmap2d(cmap2d)

    # multiply data that is normalized to [0-1] with the colormap shape to
    # map it to the corresponding colormap index. The values are then between:
    # [0 - nhue_cmap], [0 - nlight_cmap]
    idata[0] *= cmap2d.shape[0]
    idata[1] *= cmap2d.shape[1]

    # now search for the index in the colormap for each color independently
    r = map_coordinates(cmap2d[:, :, 0], idata, order=1, mode='nearest')
    g = map_coordinates(cmap2d[:, :, 1], idata, order=1, mode='nearest')
    b = map_coordinates(cmap2d[:, :, 2], idata, order=1, mode='nearest')

    # assemble the values in a [3, nrows, ncols] array and transpose it
    # to [ncols, nrows, 3] which can be used by matplotlib imshow
    rgb = np.array([r, g, b])
    rgb[:, ~mask] = np.array(fill_value)[:, None]
    rgb = rgb.reshape(3, nrows, ncols).transpose(1, 2, 0)
    return rgb


def imshow2d(data, ax=None, cmap2d='brightwheel', huenorm=None, huevmin=None,
             huevmax=None, lightnorm=None, lightvmin=None, lightvmax=None,
             **kwargs):
    """
    Plot 2 parameter 2D data array to current axis.

    :param data: numpy array with shape (2, nwidth, nheight). The first index
                 corresponds to the hue and the second to the lightness of the
                 colors.
    :param ax: a matplotlib axis instance.
    :param cmap: either:
                 numpy array with shape (nwidth, nheight, 4) that contains
                 the 4 rgba values in hue (width) and lightness (height).
                 Can be obtained by a call to get_cmap2d(name).
                 or:
                 name where name is one of the following strings:
                 'brightwheel', 'darkwheel', 'hardwheel', 'newwheel',
                 'smoothwheel', 'wheel'
    :param huenorm: a plt.Normalize() instance that normalizes the hue values.
    :param huevmin: the minimum of the huevalues. Only used if huenorm=None.
    :param huevmax: the maximum of the huevalues. Only used if huenorm=None.
    :param lightnorm: a plt.Normalize() instance that normalizes the lightness
                      values.
    :param lightvmin: the minimum of the lightness values.
                      Only used if lightnorm=None.
    :param lightvmax: the maximum of the lightness values.
                      Only used if lightnorm=None.
    :param **kwargs: remaining kwargs are passed to plt.imshow()
    """
    if ax is None:
        ax = plt.gca()
    rgb_data = data2d_to_rgb(data, cmap2d=cmap2d,
                             huenorm=huenorm, huevmin=huevmin,
                             huevmax=huevmax, lightnorm=lightnorm,
                             lightvmin=lightvmin, lightvmax=lightvmax)
    im = ax.imshow(rgb_data, **kwargs)
    return im


def get_cmap2d(name):
    """
    Return 2d colormap as [nwidth, nheight, 4] numpy array.

    :param name: can be one of the following strings:
                 'brightwheel', 'darkwheel', 'hardwheel', 'newwheel',
                 'smoothwheel', 'wheel'
    """
    package_dir = os.path.dirname(__file__)
    cmap_file = name + '.npy'
    cmap_path = os.path.join(package_dir, cmap_file)
    cmap2d = np.load(cmap_path)
    return cmap2d