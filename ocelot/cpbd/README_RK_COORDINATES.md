# Runge–Kutta tracking and coordinates

Implementation status, 2026-10-23. This replaces the earlier analysis-only note.
The changes keep the existing transfer-map classes and the Cartesian radiation
trajectory format. No fringe-field model or vector-potential conversion is added.

## Coordinate conventions

| Interface | Transverse quantities | Observation surface |
| --- | --- | --- |
| `rk_track_in_field` | Cartesian positions and **slopes** `dX/dZ`, `dY/dZ` | Fixed Cartesian Z |
| `RungeKuttaGlobalTM`, legacy `RungeKuttaTM` | Cartesian positions and **mechanical momenta / p0** | Fixed Cartesian Z |
| `RungeKuttaOcelotTM` | Positions and **mechanical momenta / p0** in transported Bishop axes | Plane normal to the bent reference at arc length s |

The wrappers accept `ParticleArray` momentum coordinates. They no longer copy
momenta into slope slots or return slopes as momenta. Scripts that deliberately
used the old erroneous slope interpretation of these slots must convert them.
The low-level slope interface and nine-row SR output remain available unchanged.

`delta=(E-E0)/(p0*c)`, where E0 is total reference energy. The longitudinal
coordinate is `tau=ct-s/beta0` for the relative map and `tau=ct-Z/beta0` for the
fixed-frame map, with the incoming offset preserved. Global timing is relative
to the straight field-frame Z axis, not to an internally subtracted bent orbit.

Inside fields these mechanical momenta are not generally canonical momenta.
The comparison code's transverse vector-potential terms are intentionally omitted.

## Implemented corrections

### Particle equations, rigidity and timing

[`high_order.py`](high_order.py) still uses RK4. Internally it now advances
normalized mechanical momenta, with explicit scalar equations:

```
PZ = sqrt(1 + 2*delta/beta0 + delta**2 - PX**2 - PY**2)
dX/dZ  = PX/PZ
dPX/dZ = (PY/PZ)*bZ - bY
dY/dZ  = PY/PZ
dPY/dZ = bX - (PX/PZ)*bZ
dpath/dZ = sqrt(1 + 2*delta/beta0 + delta**2)/PZ
```

Here `b = B / (B*rho)_reference`, retaining Ocelot's existing field-sign
convention. This uses magnetic rigidity rather than the previous 1/E
approximation. For a static magnetic field each particle has constant speed,
so elapsed `ct=path/beta_particle`. Path is integrated at the same RK4 stages.

`rk_track_in_field(..., return_path=True)` optionally returns `(trajectory, paths)`;
its default return is still the historical `(9*N, n_particles)` trajectory.
Global tracking through `rk_field` (`RungeKuttaGlobalTM` and legacy `RungeKuttaTM`)
now stores only the current particle state and accumulated path. Its temporary
array memory scales with particle count, independently of the number of steps.
It uses the same RK4 step directly in mechanical momenta, avoiding both full
trajectory storage and unnecessary conversions through slopes. SR/CSR callers
of `rk_track_in_field` still receive every trajectory and field sample.

For 10,000 particles and 3,000 points, the former global implementation allocated
2.4 GB for trajectory and path arrays alone. A `tracemalloc` measurement of the
new global path used about 3.4 MB of temporary allocations in an undulator example
(excluding the input array and prior field compilation). Endpoint comparisons,
step-count memory checks, input validation and empty-beam tests cover this change.
Roundoff differences from removing slope conversions do not require new snapshots.
Focused verification of this memory change: 93 tests passed and 5 were skipped
across the coordinate tests, RK tracking regressions and complete SR suite.

The wrappers use the RK path integral for timing. The spurious addition of
`s_start` to tau is removed. Invalid energy, grids and non-forward momenta raise.

### Normal-plane crossings and slices

`rk_field_to_ocelot` computes the reference from the same element entrance for
every call. A slice uses the entrance and exit positions and Bishop frames at
its absolute arc coordinates. Incoming local positions/momenta are transformed
back to the original Cartesian field frame before integration.

Particles start at their own Cartesian Z on the entrance normal plane. At the
exit, a crossing iteration advances position, momentum and path through the
field together. It does not approximate that last motion by a straight tangent.

The maximum step for relative tracking and its reference is based on the full
element length, `l/(npoints-1)`, including for slices. Reference data are recomputed,
not cached on shared element objects. This costs more than a cache but avoids
stale references and occurrence-dependent mutable state. The low-level/global
interfaces retain N samples over their requested Cartesian interval.

### XYQuadrupole CSR trajectory

[`CSR.prepare`](csr.py) no longer substitutes a uniform dipole for XYQuadrupole:

- Without RK it samples the same varying-curvature analytical reference as the
  XYQuadrupole survey; it no longer falls through to a straight line.
- With RK it follows the actual ideal quadrupole field, including entrance
  offsets, entrance directions (`x_angle`, `y_angle`) and roll, using the
  specified CSR energy. Positions and all three field components are
  transformed with the minimal-rotation entrance frame from Appendix B of
  the updated paper. The analytical branch includes these angles in its
  varying curvatures as well.

The segment is joined to the preceding reference position and transported axes.
CSR's existing `(s,X,Y,Z,tx,ty,tz)` format and wake calculation are retained.
The exact RK reference and analytical weak-curvature reference are distinct
models. XYQuadrupole itself still supports analytical `TransferMap` only.

### Radiation trajectories and returned particles

[`radiation_py.py`](../rad/radiation_py.py) still consumes Cartesian trajectories.
Undulator entrance/exit momenta now convert to/from slopes explicitly. Its
`ParticleArray` endpoint retains incoming timing, adds flight time consistently,
and advances s. Coherent radiation preserves the original bunch arrival offsets.
`BeamTraject` can retain that endpoint and restore it without reconstructing a
reference from the mean bunch path or accumulating the path twice.

Standalone `BeamTraject(...).p_array_end(entrance_array)` without a saved endpoint
also has corrected momentum and timing conventions. That legacy reconstruction
uses sampled-path quadrature and assumes a fixed frame and constant reference
energy; it is less accurate than the retained RK endpoint.

Ordinary, non-undulator trajectory sampling now includes the true entrance and
exit, uses N-1 advances for N samples, converts momenta to slopes, and no longer
uses `tau+Z` as Cartesian Z. When the existing aggregate undulator energy-loss
option is used, endpoint energy and momentum normalization are updated together.
Radiation gamma uses each segment's energy deviation, not the final deviation
for every preceding segment. Incoherent SR accumulation converts each segment's
field amplitude to a common gamma normalization when energy changes.

## Validation and examples

- `unit_tests/cpbd/elements/test_runge_kutta_coordinates.py`: exact finite-beta
  drift, off-energy circular motion, nonzero start, initial arrival offsets,
  full versus uneven slices, repeated calls, and an independent DOP853 Lorentz
  calculation with a plane event and longitudinal magnetic field. A planar
  undulator check also compares flight time with quadrature of its exact orbit
  for negative, zero and positive energy deviations.
- `unit_tests/cpbd/elements/test_xyquadrupole_csr.py`: analytical CSR versus
  survey, preceding tilted bend, both k1 signs, roll, and independent Lorentz
  reference checks for the RK trajectory.
- `unit_tests/sr_test/test_radiation_helpers.py`: SR endpoint versus fixed-frame
  tracking, standalone conversion, sampling, timing and energy-loss conventions.
- Existing RK, XYQuadrupole, architecture and radiation regressions also run.

Initial targeted verification: 376 CPBD tests passed; 52 SR tests passed with
4 skipped; 9 laser-heater tests passed with 1 skipped; the existing CSR arc-versus-RK
tracking check passed. There are now 52 new parametrized regression cases across the
new coordinate/CSR tests and the added radiation helper tests. All eight tutorial
code cells executed successfully, the plots were inspected, and both updated
public Markdown pages compiled as MDX.

A separate read-only comparison with the supplied `inspect_methods` code gives
maximum six-coordinate differences about `2.2e-15` at beta0=0.8 and 0.999999,
using 501 points. Its canonical outputs were converted to mechanical momenta
using its own known potential, and sigma was changed to -tau. Uneven slices
agreed with the full map to about `1.5e-14` in that example.

The SR and laser-heater CSR snapshots change with the corrected rigidity, integration variables,
sampling, momentum and arrival conventions. They are regenerated after checking
the original implementation against the old references and checking the new
undulator trajectory against independent DOP853 integration. Doubling the grid
reduced its position and slope errors by about 16, as expected for RK4. Spectrum
convergence was checked separately; saved snapshots are regressions at their
specified accuracy, not claims of absolute spectral accuracy. At 130 MeV the
rigidity correction shifts the laser-heater CSR reference by about 0.23 micrometres
and its tangent by about 0.77 microradians; its dependent CSR beam snapshots are
updated as well. The old implementation was checked against those old references.

The full-suite follow-up identified six more outdated snapshots: the two
`csr_ex` undulator cases, three `rk_track` chirp cases and `section_track`'s CSR
case. All six passed with the old RK formulas restored in memory. Their
references are now updated without relaxing tolerances. For the zero-chirp
600 MeV, K=44.81 undulator, global tau increases by 1.31185625 mm relative to
the old result: this is the extra flight time along the oscillating trajectory,
previously subtracted as a reference path. Independent exact-orbit quadrature
checks the new value. The CSR cases change through their RK reference orbits;
their largest particle energy-deviation differences are about 3.8e-7, 1.6e-7
and 4.3e-8, respectively. The chirp test also had a pre-existing branch mistake:
its positive-chirp case reset the chirp to zero. It now preserves positive chirp,
so all three cases exercise distinct inputs.

After these updates, `python -m pytest unit_tests -q` completed with **852 passed,
59 skipped and one warning** (the beam-matching test reaches its optimizer's
function-evaluation limit). All six previously failing cases pass; the three
additional cases check the planar undulator flight time described above.

The executed notebook [`18_runge_kutta_tracking.ipynb`](../../demos/docs/18_runge_kutta_tracking.ipynb)
and its generated public tutorial cover coordinate choice, slicing, offsets,
finite-beta drift, the correct fixed-Z circular solution, a numerical matrix
check, and SR trajectories/spectra. The matrix check uses existing tracking;
there is no new automatic optics API.

## Work deliberately remaining

1. **Canonical conversion inside fields.** Supply compatible element-specific
   vector potentials/gauges before calling the mechanical output canonical.
2. **Fringes and edge maps.** A supplied physical field or an explicit element
   edge model is required. Coordinate conversion cannot create missing forces.
3. **Arbitrary-field optics.** `first_order_tms` and `twiss` retain the element's
   nominal analytical model. The tutorial shows finite differencing of the
   actual RK map with a fixed reference. Automatic numerical optics, coupled
   optics and canonical interior-plane validation are separate work.
4. **Nominal survey versus actual orbit.** The installed layout must remain
   distinct from the energy-dependent reference through a supplied field.
   No automatic field-derived replacement of `MagneticLattice.survey` is made.
5. **General laboratory assembly for SR/CSR.** SR's mixed-element sampling still
   assumes a straight layout and samples nominal maps in non-undulator sections.
   General placement of arbitrarily oriented fields,
   conversion of local-map samples into a common laboratory frame, and an SR
   endpoint relative to a curved reference remain to be implemented. The CSR
   correction here is specifically its XYQuadrupole branch and frame continuity;
   other existing CSR field-placement approximations are not a general survey.
6. **Integrator scope.** Static B fields and forward laboratory-Z motion only;
   electric/time-dependent fields and trajectories turning back in Z need
   additional physics or another independent variable. RK4 is not an exactly
   symplectic integrator. Accurate particle timing does not make every SR
   radiation formula exact at finite beta.

Research reference:
`/Users/tomins/Nextcloud/publications/2025-10-15_XYQuadrupole/theory_analysis/inspect_methods`,
especially `lorentz.py` and `frames.py`.
