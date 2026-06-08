# Changelog

All notable changes to OCELOT will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning is date based: `YY.MM.patch`.

## [Unreleased]

- None yet.

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
