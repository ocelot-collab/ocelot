"""Magnetic field-map adapters preserve Cartesian components and coordinates."""

import numpy as np
import pytest

from ocelot.cpbd.field_map import FieldMap, field_map_to_field_func


@pytest.mark.parametrize('units,scale', [('mm', 1000.), ('m', 1.)])
def test_rectangular_map_interpolates_all_three_components(units, scale):
    field_map = FieldMap(None, xyz_units=units)
    field_map.x_arr = np.array([-.02, 0., .02]) * scale
    field_map.y_arr = np.array([-.01, .01]) * scale
    field_map.z_arr = np.array([.1, .2, .3]) * scale
    z, y, x = np.meshgrid(field_map.z_arr / scale - .1,
                          field_map.y_arr / scale, field_map.x_arr / scale, indexing='ij')
    field_map.Bx_arr = (x + 2*y + 3*z).ravel()
    field_map.By_arr = (4*x - y + z).ravel()
    field_map.Bz_arr = (x - y - z).ravel()
    field = field_map_to_field_func(field_map)
    x, y, z = .005, .002, .15
    np.testing.assert_allclose(field(x, y, z), [x+2*y+3*z, 4*x-y+z, x-y-z], atol=1e-15)
    values = field(np.array([x, x]), y, np.array([z, z]))
    np.testing.assert_allclose(values, np.repeat(np.array(field(x, y, z))[:, None], 2, axis=1))


def test_repeated_short_line_map_preserves_components_and_period():
    field_map = FieldMap(None, xyz_units='m')
    field_map.z_arr = np.array([.2, .3])
    field_map.Bx_arr = np.array([1., 2.])
    field_map.By_arr = np.array([3., 4.])
    field_map.Bz_arr = np.array([5., 6.])
    field_map.field_file_rep = 2
    field = field_map_to_field_func(field_map)
    np.testing.assert_allclose(field(0., 0., .05), [1.5, 3.5, 5.5])
    np.testing.assert_allclose(field(0., 0., .15), [1.5, 3.5, 5.5])
    np.testing.assert_allclose(field(0., 0., 2*(.3-.2)), [2., 4., 6.])


def test_six_column_map_reads_and_interpolates_unordered_grid(tmp_path):
    path = tmp_path / 'grid.txt'
    z, y, x = np.meshgrid([10., 20., 30.], [-2., 2.], [-1., 0., 1.], indexing='ij')
    rows = np.column_stack([a.ravel() for a in (x, y, z, x+y+z, x-y+z, x+y-z)])
    np.random.default_rng(31).shuffle(rows)
    np.savetxt(path, rows)
    field = field_map_to_field_func(FieldMap(str(path)))
    # File coordinates are in mm; callback Z is relative to map entrance.
    np.testing.assert_allclose(field(.0005, .001, .015), [26.5, 24.5, -23.5], atol=1e-12)
