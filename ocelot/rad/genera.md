# Radiation module notes

This document describes OCELOT's native Python synchrotron-radiation (SR)
implementation in:

- `ocelot/rad/radiation_py.py`
- `ocelot/rad/screen.py`
- `ocelot/rad/spline_py.py`

It records current behavior and unresolved problems. Implementation history and
detailed optimization logs belong in version control, not here.

## Tutorials

The following notebooks are usage examples. They are not part of the SR module
implementation:

- `demos/ipython_tutorials/pfs_1_synchrotron_radiation.ipynb` introduces
  spontaneous-radiation spectra, spatial distributions, magnetic-field maps,
  and custom magnetic fields.
- `demos/ipython_tutorials/9_thz_source.ipynb` combines accelerator tracking
  with coherent THz-radiation calculation from a macroparticle bunch.

## Calculation flow

1. `track4rad_beam` tracks particles through the lattice and creates trajectory
   segments.
2. `traj2motion` converts a segment to millimetres, inserts three-point
   Gauss-Legendre nodes, interpolates coordinates and fields, and integrates
   the squared transverse slopes.
3. `radiation_py` integrates the complex horizontal and vertical fields on a
   `Screen`.
4. `calculate_radiation` or `coherent_radiation` converts the fields to photon
   distributions.

A raw trajectory has shape `(9 * n_points, n_particles)`. Each nine-row block
contains:

```text
x, x', y, y', z, relative momentum, Bx, By, Bz
```

## Public calculations

### `calculate_radiation`

Calculates conventional incoherent synchrotron radiation from a `Beam`.
`Beam.I` is interpreted as current in amperes. The photon distributions are
photon flux:

```text
photons / second / mm^2 / (10^-3 relative bandwidth)
```

The supplied `Screen` is mutated and returned.

### `coherent_radiation`

Calculates coherent radiation from a finite macroparticle bunch. Each
macroparticle field is weighted by its represented electron count:

```python
n_e = p_array.q_array[i] / q_e
```

Complex fields are summed before calculating intensity. No repetition rate or
beam current is used, so the photon distributions are:

```text
photons / bunch / mm^2 / (10^-3 relative bandwidth)
```

The shared plotting code may display `ph/sec`. For coherent radiation this is
numerically equivalent to assuming one bunch per second. For repetition rate
`f_rep` in hertz:

```python
photons_per_second = photons_per_bunch * f_rep
```

Both the supplied `Screen` and `ParticleArray` are mutated. Initial `tau`
coordinates define the coherent phase; final tracked coordinates are written
back to the particle array.

### Integrating a spectrum

`10^-3 BW` means:

```text
delta_E / E = abs(delta_lambda) / lambda = 10^-3
```

It is not the spacing between adjacent energy samples. For a one-point
spectrum, photons in the selected energy range are:

```python
E = screen.Eph
S = screen.Total
N = np.trapezoid(S / (1e-3 * E), E)
```

For coherent radiation, `N` remains photons per bunch per square millimetre.
For incoherent radiation, it remains photons per second per square millimetre.
Obtaining the total photon number also requires integration over screen area
or solid angle.

## Units

- Beamline and particle coordinates: metres.
- Internal `Motion` and screen integration coordinates: millimetres.
- Photon energy: eV.
- Electron/reference energy: GeV.
- Magnetic field: tesla.
- Stored screen fields: flattened logical order `(energy, y, x)`.

## Trajectory sampling

For an undulator, an explicit `npoints` value is the exact number of
Runge-Kutta trajectory points, including endpoints. It overrides the automatic
estimate and `accuracy`. The value must be an integer of at least four.

The radiation module reads `npoints` first from the active Runge-Kutta
transformation and then from the undulator constructor configuration. This is
required because radiation tracking calls `rk_track_in_field` directly rather
than using the normal lattice transformation path.

Non-undulator sections before, between, and after undulators are included as
separate trajectory segments.

## Current implementation status

- Coherent multi-segment fields are accumulated once per segment.
- `Undulator` and `ParticleArray` subclasses are supported.
- Explicit undulator `npoints` is honored.
- Trailing non-undulator sections are tracked.
- Repeated array appends and avoidable full-screen deep copies were removed.
- One field workspace is reused across coherent macroparticles.
- Screen field arrays remain views of `memory_screen` after phase rebuilding.
- Spline integration and long one-point spectra use specialized Numba paths
  when Numba is available.
- A non-trackable lattice raises a descriptive `ValueError`.

## Open problems

### General curved beamline assembly

Non-undulator sampling now stores the true entrance and exit, using N-1 advances
for N samples. It converts momenta to slopes and no longer adds tau to Cartesian
Z. Undulator tracking retains flight time and the actual endpoint ParticleArray;
SR no longer infers its reference from the bunch's mean path.

General placement of curved local-map samples in a common laboratory frame
remains to be implemented. The current mixed-element SR path assumes a straight
layout. See [RK coordinate status](../cpbd/README_RK_COORDINATES.md).

### Nonuniform input to `x2xgaus`

`x2xgaus` uses the first interval length for every interval. Current generated
trajectories are intended to be uniform, but the helper neither validates that
assumption nor supports a nonuniform grid.

### Random-number control

Quantum diffusion uses NumPy's global random state. Accepting an optional
`numpy.random.Generator` would allow reproducible calculations without global
seeding.

### Input validation

`accuracy` should be validated as finite and positive. Very short or malformed
trajectory arrays should also produce descriptive exceptions before cubic
interpolation.

### Zero-current beam behavior

`calculate_radiation` currently changes `ebeam.I` to `0.1 A` when it is zero.
This input mutation is surprising. A future change should either require a
positive current or use a documented local fallback without modifying the
beam.

### Public API boundary

`ocelot.rad` imports `radiation_py` with `*`, while the module has no explicit
`__all__`. Internal helpers are therefore accidentally public. Define the
supported API before renaming or removing legacy helpers.

### Remaining large-scale performance work

The expensive coherent path still performs spline interpolation separately
for every particle. Batching trajectory preprocessing across particles is the
main remaining structural optimization, but it requires benchmarks and
numerical comparison of fields and photon distributions.

## Verification

Run the synchrotron-radiation tests with:

```bash
PYTHONPATH=. pytest unit_tests/sr_test -q
```

Performance changes should be measured after Numba warm-up and checked for
spectrum, spatial-screen, and full-3D cases.

## References

- S. Tomin and G. Geloni, "Synchrotron Radiation Module in OCELOT Toolkit,"
  Proceedings of IPAC 2019, WEPTS017,
  [doi:10.18429/JACoW-IPAC2019-WEPTS017](https://doi.org/10.18429/JACoW-IPAC2019-WEPTS017).
  The paper describes the Python Runge-Kutta trajectory solver, the
  frequency-domain SR solver, coherent macroparticle calculations, radiation
  energy loss, and quantum diffusion.
- G. Geloni, T. Tanikawa, and S. Tomin, "Dynamical effects on superradiant THz
  emission from an undulator," Journal of Synchrotron Radiation 26 (2019),
  737-749,
  [doi:10.1107/S1600577519002509](https://doi.org/10.1107/S1600577519002509).
