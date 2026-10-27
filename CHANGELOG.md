# Changelog

All notable changes to OCELOT will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning is date based: `YY.MM.patch`.

## [Unreleased]

### Added

- Analytical XYQuadrupole transfer maps and survey support, with CSR reference
  trajectories and physics regression tests. Entrance directions `x_angle`
  and `y_angle` extend the map, arrival time, survey and CSR field geometry
  according to the updated paper; both default to zero for compatibility.
- ExactDriftTM for finite-angle, finite-energy drift tracking.
- Matcher phase-advance targets and explicit lattice-occurrence resolution for
  repeated element instances.
- Release commands in `tasks.py` with explicit PyPI/Anaconda upload flags,
  clean tagged-source checks, and a restored release checklist.
- Demo execution reports, per-example logs, copied input trees and timeouts.

### Changed

- Separate `periodic_twiss()` from propagation with `twiss()`; unstable periodic
  solutions raise `UnstableLatticeError`. Update callers that relied on implicit
  periodic initialization.
- Use lazy public imports to reduce startup costs.
- RK tracking uses mechanical normalized momenta, corrected finite-beta rigidity
  and flight time, and a common entrance reference for curvilinear slices.
  Global tau is relative to the straight Cartesian Z axis.
- Global RK tracking retains only the current state and path, avoiding full
  trajectory storage for large bunches; SR keeps its sampled trajectories.
- Update RK tutorials, coordinate documentation and affected regression references.
- Development/release scripts are consolidated in `tasks.py`. Local preparation
  no longer pushes or uploads automatically, and failed checks stop the command.
- Raised the minimum supported Python version from 3.10 to 3.11 in the Python
  package and conda recipe, matching the pandas 3 and scikit-learn 1.9 minimums.
- Added scikit-learn and TFS-Pandas dependencies and synchronized their conda
  requirements with the Python package.
- Added `mpi` and `moga` installation extras. Kept the `openpmd` extra for
  openPMD I/O and expanded `pmd` to include both openPMD and openpmd-beamphysics.
- Replaced the obsolete `pmd-beamphysics` distribution with the maintained
  `openpmd-beamphysics` package and its `beamphysics` import, fixing PMD file
  writing with NumPy 2.
- Added separate script and notebook exclusions to the demo runner.

### Fixed

- Cavity and TWCavity longitudinal maps now consistently use the time-like
  coordinate `tau = c * (t - t_ref)`: corrected RF phase and energy-kick beta
  factors, matched the tracking Jacobian to the linear map, and stabilized
  first- and second-order flight-time coefficients at either zero crossing
  and at vanishing voltage. Updated the R-matrix notebook derivation.
- Keep Genesis beam-rematching Twiss parameters real when NumPy returns complex
  eigenpair arrays with zero imaginary parts, preventing CPBD phase errors on
  subsequent propagation.
- SR momentum/slope conversion, endpoint timing and restoration, trajectory
  sampling, and normalization when radiation energy loss changes beam energy.
- MAD second-order map integrals and singular limits; Twiss phase advances
  between pi and 2*pi; cumulative particle propagation in `trace_z`.
- Matcher matrix energies after acceleration and success reporting when targets
  remain unmet; preserve metadata and loss tracking when extracting particles.
- Release commits with staged differences and Conda source-tag version updates.
- Kept extracted and loaded response matrices writable with pandas 3 so callers
  can continue editing and injecting response-matrix values.
- Corrected the longitudinal-coordinate sign when importing PMD particle groups,
  preserving bunch order through an export/import round trip.
- Added an absolute tolerance of `1e-15` to the CSR particle and space-charge
  Twiss reference comparisons so roundoff near zero does not fail relative
  checks. Reference data and simulation algorithms are unchanged.
- Documented a Numba `workqueue` workaround for conflicting OpenMP runtimes in
  mixed conda/pip environments on macOS.

### Removed

- Removed the unused `docs` extra and Sphinx dependencies. Public documentation
  is built with Docusaurus in the website repository.

## [26.06.1] - 2026-06-08

### Added

- Added `scripts/publish_release.py` for explicit PyPI and Anaconda artifact
  publishing after release preparation.
- Added `scripts/RELEASE.md` with the recommended release checklist and command
  order.

### Changed

- Made global lattice transfer-map requests permissive and silent when an
  element family does not support the requested active transformation. The
  element now falls back to its `default_tm`, while explicit unsupported
  family-specific requests still raise an error.
- Updated element architecture documentation to describe explicit versus global
  transfer-map selection.
- Ignored generated `*.egg-info/` package metadata directories.
- Cleaned up selected electron-beam demo scripts.

### Fixed

- Fixed the conda release build script metadata for the `ocelot-collab`
  package.
- Removed repeated global transfer-map fallback warnings for cavities during
  common optics and tracking calls such as `twiss`.

## [26.06.0] - 2026-06-07

### Added

- Added Runge-Kutta transfer maps for arbitrary magnetic fields, including the
  OCELOT and global-frame RK variants, constructor parameters, field tests, and
  tutorial material.
- Added element offset support in constructors and tests covering offset-aware
  RK tracking.
- Added the object-oriented matcher module, matcher tests, SVD weights, and
  ParticleArray Twiss matching support.
- Added BBA utilities, response-matrix caching, per-measurement launch modes,
  tutorial material, and unit tests.
- Added LongWake and LinLongWake tracking support with tests and reference data.
- Added lattice layout and survey tooling, including 2D/3D plotting helpers,
  survey demos, and survey tests.
- Added architecture contract tests for element transfer-map selection,
  slicing, wrapper atoms, matrix maps, energy gain, and undulator behavior.
- Added legacy matching tutorial notebooks and updated THz source, R-matrix,
  compression, and high-resolution optics tutorials.
- Added release-preparation automation in `scripts/prepare_release.py`.

### Changed

- Refactored `ocelot.cpbd.beam` into a package with dedicated modules for core
  beam objects, particles, generation, analysis, utilities, and noise.
- Clarified and documented CPBD element and transfer-map contracts, parameter
  containers, wrapper APIs, slice behavior, and supported methods.
- Updated lattice serializers and adaptors to use the public wrapper API.
- Cleaned up default element transfer-map selection rules and matrix slice atom
  naming.
- Made `numexpr`, `pyfftw`, and `numba` required package dependencies in
  `setup.py`; the public installation docs previously listed them as optional
  speed-up dependencies.
- Declared Python 3.10+ package support.

### Fixed

- Fixed sliced Runge-Kutta tracking behavior and revised the RK tracking
  tutorial.
- Fixed RBend edge angles after angle updates and deduplicated first-order bend
  edge matrix calculations.
- Fixed Twiss emittance handling when energy is not provided.
- Fixed adaptor attribute access and TFS import behavior.
- Fixed MPI imports so MPI is imported only when used.
- Fixed cavity second-order terms for the phase-90-degree removable singularity
  and negative-gain fallback behavior.
- Fixed matching edge cases, including Drift length variables and delta
  constraints.
- Fixed floating-point handling reported in issue #293.

### Removed

- Removed the monolithic `ocelot/cpbd/beam.py` implementation in favor of the
  new `ocelot.cpbd.beam` package structure.

[Unreleased]: https://github.com/ocelot-collab/ocelot/compare/v26.06.1...HEAD
[26.06.1]: https://github.com/ocelot-collab/ocelot/compare/v26.06.0...v26.06.1
[26.06.0]: https://github.com/ocelot-collab/ocelot/compare/v25.07.0...v26.06.0
