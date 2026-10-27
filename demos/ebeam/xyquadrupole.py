"""Linear optics and survey for an offset quadrupole with angled entry.

Run from the repository root: python demos/ebeam/xyquadrupole.py
The weak-curvature model requires small reference offsets and angles.
"""
import numpy as np
import ocelot as ocl


# Offsets are in metres in the unrolled magnet frame. Entrance directions
# are dX/ds, dY/ds (unit-tangent components, approximately angles in radians).
# l is reference arc length, not the straight mechanical magnet length.
quad = ocl.XYQuadrupole(l=1.0, k1=1.0, x_offs=0.02, y_offs=0.03,
                       x_angle=0.01, y_angle=-0.005, tilt=0.2, eid="XYQ")
lat = ocl.MagneticLattice([quad, ocl.Drift(l=0.5, eid="D")])
energy = 0.13  # GeV
matrix = quad.R(energy)[0]
np.set_printoptions(precision=6, suppress=True)
print("Exit dispersion [Dx, Dxp, Dy, Dyp]:", matrix[:4, 5])
print("R56 in Ocelot's tau convention [m]:", matrix[4, 5])

# Two slices use the same reference and reproduce the full-element map.
first = quad.get_section_tms(0.5)[0].get_params(energy).get_rotated_R()
second = quad.get_section_tms(0.5, start_l=0.5)[0].get_params(energy).get_rotated_R()
print("Full map versus two slices, max difference:", np.max(np.abs(second @ first - matrix)))

tws0 = ocl.Twiss(beta_x=10., beta_y=12., E=energy)
optics = ocl.twiss(lat, tws0, nPoints=31)
print("Lattice exit beta_x, beta_y [m]:", optics[-1].beta_x, optics[-1].beta_y)

# Survey coordinates start in the incoming reference frame. The entrance
# angles describe its direction relative to the magnet axis; do not add
# them again to the lattice's initial survey direction.
midpoints, endpoints = lat.survey()
for element, endpoint in zip(lat.sequence, endpoints[1:]):
    print(element.id, "reference exit [m]:", endpoint["r_end"])
    print("  Exit unit tangent:", endpoint["W"][:, 2])
