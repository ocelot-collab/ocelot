from ocelot.optics.wave import *
from ocelot.rad.optics_elements import *
from ocelot.rad.optics_line import *
from ocelot.common.ocelog import *
import numpy as np
import copy
_logger = logging.getLogger(__name__)


def get_transfer_function(element):
    element.mask = Mask()

    if element.__class__ is ApertureRect:
        mask = RectMask()
        mask.lx = element.lx
        mask.ly = element.ly
        element.mask = mask
    elif element.__class__ is FreeSpace:
        element.mask = DriftMask(element.l, element.mx, element.my)
    elif element.__class__ is ThinLens:
        element.mask = LensMask(element.fx, element.fy)


class Mask:
    def __init__(self):
        pass

    def apply(self, dfl):
        return dfl

    def get_mask(self, dfl):
        pass

    def __mul__(self, other):
        m = copy.deepcopy(self)
        if other.__class__ in [self] and self.mask is not None and other.mask is not None:
            m.mask = self.mask * other.mask
            return m


class RectMask(Mask):
    def __init__(self):
        Mask.__init__(self)
        self.lx = np.inf
        self.ly = np.inf
        self.mask = None

    def apply(self, dfl):
        if self.mask is None:
            self.get_mask(dfl)
        mask_idx = np.where(self.mask == 0)

        # dfl_out = deepcopy(dfl)
        dfl_energy_orig = dfl.E()
        dfl.fld[:, mask_idx[0], mask_idx[1]] = 0

        if dfl_energy_orig == 0:
            _logger.warn(ind_str + 'dfl_energy_orig = 0')
        elif dfl.E() == 0:
            _logger.warn(ind_str + 'done, %.2f%% energy lost' % (100))
        else:
            _logger.info(ind_str + 'done, %.2f%% energy lost' % ((dfl_energy_orig - dfl.E()) / dfl_energy_orig * 100))
        # tmp_fld = dfl.fld[:,idx_x1:idx_x2,idx_y1:idx_y2]
        return dfl

    def get_mask(self, dfl):
        """
        model rectangular aperture to the radaition in either domain
        """
        _logger.info('applying square aperture to dfl')

        if np.size(self.lx) == 1:
            self.lx = [-self.lx / 2, self.lx / 2]
        if np.size(self.ly) == 1:
            self.ly = [-self.ly / 2, self.ly / 2]
        _logger.debug(ind_str + 'ap_x = {}'.format(self.lx))
        _logger.debug(ind_str + 'ap_y = {}'.format(self.ly ))

        idx_x = np.where((dfl.scale_x() >= self.lx[0]) & (dfl.scale_x() <= self.lx[1]))[0]
        idx_x1 = idx_x[0]
        idx_x2 = idx_x[-1]

        idx_y = np.where((dfl.scale_y() >= self.ly [0]) & (dfl.scale_y() <= self.ly [1]))[0]
        idx_y1 = idx_y[0]
        idx_y2 = idx_y[-1]

        _logger.debug(ind_str + 'idx_x = {}-{}'.format(idx_x1, idx_x2))
        _logger.debug(ind_str + 'idx_y = {}-{}'.format(idx_y1, idx_y2))

        self.mask = np.zeros_like(dfl.fld[0, :, :])
        self.mask[idx_y1:idx_y2, idx_x1:idx_x2] = 1
        return self.mask

    def __mul__(self, other):
        print("mul")
        m = RectMask()
        if other.__class__ in [RectMask, ] and self.mask is not None and other.mask is not None:
            m.mask = self.mask * other.mask
            return m


class DriftMask(Mask):

    def __init__(self, l, mx, my):
        Mask.__init__(self)
        self.l = l
        self.mx = mx
        self.my = my

    def apply(self, dfl):
        dfl.prop_m(z=self.l, m=(self.mx, self.my))
        return dfl


class LensMask(Mask):

    def __init__(self, fx, fy):
        Mask.__init__(self)
        self.fx = fx
        self.fy = fy

    def apply(self, dfl):
        print("LENS APPLY fx_{}".format(self.fx))
        dfl.curve_wavefront(r=self.fx, plane='x', domain_z=None)
        return dfl
