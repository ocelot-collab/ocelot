"""
module for misalignment error generation
Sergey Tomin
"""
import copy
from numbers import Integral

from numpy.random import normal
import scipy.stats as stats
import numpy as np
from ocelot.cpbd.elements import Sextupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice


class Errors:
    def __init__(self):
        self.sigma_x = 100e-6
        self.sigma_y = 100e-6


def tgauss(mu=0, sigma=1, trunc=2, num=1):
    if sigma == 0:
        return mu if num == 1 else np.full(num, mu)
    lower = -sigma*trunc
    upper = sigma*trunc
    X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
    if num == 1:
        x = X.rvs(1)[0]
    else:
        x = X.rvs(num)
    return x


def tgauss_old(mu=0, sigma=1, trunc=2):
    if sigma == 0:
        return mu
    err = normal(mu, sigma, size=50)
    err_ok = err[abs(err) < trunc*sigma]
    return err_ok[0]


def create_copy(lattice, nsuperperiods):

    lat_copy_seg = []
    for i in range(nsuperperiods):

        for n, elem in enumerate(lattice.sequence):
            lat_copy_seg.append(copy.deepcopy(elem))
            if lat_copy_seg[-1].id is None:
                lat_copy_seg[-1].id = lat_copy_seg[-1].type[:2] + "_sp" + str(i) + "_" + str(n)
            else:
                if elem.__class__ != Sextupole:
                    lat_copy_seg[-1].id = lat_copy_seg[-1].id + "_sp" + str(i) + "_" + str(n)
    return MagneticLattice(lat_copy_seg, method=lattice.method)


def _sample_alignment_errors(lattice, er_list, shared_groups):
    sigmas = {}
    for family, settings in er_list.items():
        offset, roll = settings["offset"], settings["dtilt"]
        if not (np.isfinite(offset) and offset >= 0 and np.isfinite(roll) and roll >= 0):
            raise ValueError("Alignment RMS values must be finite and nonnegative")
        sigmas[family] = (offset, roll)

    sequence = lattice.sequence
    leaders = {}
    for group in shared_groups:
        leader = None
        for i in group:
            if not isinstance(i, Integral) or not 0 <= i < len(sequence):
                raise ValueError(f"Invalid slice index {i!r}")
            if i in leaders:
                raise ValueError(f"Slice index {i} belongs to more than one group or occurs twice")
            if type(sequence[i]) not in sigmas:
                raise ValueError(f"Slice index {i} is outside the configured magnet families")
            if leader is None:
                leader = i
            elif sigmas[type(sequence[i])] != sigmas[type(sequence[leader])]:
                raise ValueError("Slices of one magnet must have the same alignment RMS values")
            leaders[i] = leader

    records = [None] * len(sequence)
    for i, elem in enumerate(sequence):
        if type(elem) in sigmas:
            leader = leaders.get(i, i)
            if records[leader] is None:
                offset, roll = sigmas[type(elem)]
                records[leader] = dict(dx=tgauss(sigma=offset), dy=tgauss(sigma=offset),
                                       dtilt=tgauss(sigma=roll))
            if i != leader:
                records[i] = records[leader].copy()
    return records


def apply_alignment_errors(lattice, records):
    """Apply absolute alignment values by sequence index and return the same lattice.

    Supply one entry per occurrence in the original sequence order. Each entry
    is either None (keep its alignment) or a dict with dx, dy [m] and dtilt
    [rad]. If a reused element needs different alignments, the affected
    occurrences receive independent copies, including their magnet strengths.
    """
    if len(records) != len(lattice.sequence):
        raise ValueError("Expected one alignment record per lattice occurrence")

    occurrences = {}
    for i, (elem, record) in enumerate(zip(lattice.sequence, records)):
        values = None
        if record is not None:
            if not hasattr(elem.element, "dtilt"):
                raise ValueError(f"{type(elem).__name__} at index {i} does not support mechanical roll")
            values = tuple(float(record[key]) for key in ("dx", "dy", "dtilt"))
            if not np.isfinite(values).all():
                raise ValueError(f"Nonfinite alignment at index {i}")
        occurrences.setdefault(elem, {}).setdefault(values, []).append(i)

    assignments = []
    for elem, groups in occurrences.items():
        # Keep the original for untouched occurrences, or for the first draw.
        original = None if None in groups else next(iter(groups))
        for values, indices in groups.items():
            if values is None:
                continue
            target = elem if values == original else copy.deepcopy(elem)
            assignments.append((target, values, indices))

    for target, values, indices in assignments:
        target.dx, target.dy, target.dtilt = values
        for i in indices:
            lattice.sequence[i] = target
    # Alignment setters invalidate the affected maps, including first-order optics.
    lattice._refresh_element_lookup()
    return lattice


def errors_seed(lattice, er_list, *, shared_groups=()):
    """Seed magnet alignments in place and return (lattice, records).

    er_list maps exact element classes to {"offset": sigma_xy, "dtilt": sigma_roll}.
    Each transverse offset [m] and the mechanical roll [rad] are independent
    zero-mean Gaussian draws truncated at two sigma. Zero sigma sets that
    component to zero. Supported magnetic families include quadrupoles,
    sextupoles, octupoles, bends, XYQuadrupole and thin Multipole. Correctors
    have a separate horizontal/vertical kick model and cannot be seeded here.

    Occurrences are independent by default. Pass groups of zero-based sequence
    indices, e.g. shared_groups=[(0, 2)], to give those occurrences the same
    three values. Each bend occurrence still rolls about its own entrance.
    To slice one rigidly rolled bend, use get_section_tms on a single bend.
    records[i] is the alignment dict for sequence[i], or None for an
    untargeted occurrence. Replay it with apply_alignment_errors.
    """
    records = _sample_alignment_errors(lattice, er_list, shared_groups)
    return apply_alignment_errors(lattice, records), records


if __name__ == "__main__":
    for i in range(1):
        rv = tgauss()
        if abs(rv) > 2.99:
            print(rv)
