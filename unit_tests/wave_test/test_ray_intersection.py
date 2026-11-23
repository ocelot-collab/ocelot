"""Plane intersections retain the surface frame and the nearest accepted hit."""
import numpy as np

from ocelot.optics.elements import Detector
from ocelot.optics.ray import Ray, find_intersections, trace


def test_tilted_detector_and_closer_aperture_miss():
    normal = np.array([0., .6, .8])
    along = np.array([0., .8, -.6])
    centre = np.array([0., 0., 1.])
    hit = centre + .2 * along
    ray = Ray(r0=hit - normal, k=normal)
    detector = Detector(r=centre, no=normal, size=[1., .3, 0], nx=5, ny=12)
    # This closer plane intersects the ray outside its active aperture.
    miss = Detector(r=hit - .5 * normal + [2., 0., 0.], no=normal, size=[.1, .1, 0])
    distance, obj, local, returned_normal = find_intersections(ray, lambda: [detector, miss])
    assert obj is detector
    np.testing.assert_allclose(distance, 1.)
    np.testing.assert_allclose(local, [0., .2, 0.], atol=1e-15)
    np.testing.assert_allclose(returned_normal, normal)
    trace(ray, lambda: [detector, miss])
    assert detector.matrix.sum() == detector.matrix[2, 10] == 1
    detector.clear()
    assert detector.matrix.sum() == 0
