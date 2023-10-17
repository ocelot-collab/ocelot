# Changelog

All notable changes to OCELOT will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
~and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).~ (Not yet, maybe it won't.  Versioning is based on year.month.day.

## [Unreleased]

### Changed

- The `CSR` class now gives an explicit `logging` `error` when failing in `arcline` and
  a more useful explanation.  This error often occurs because the user has mixed CSR for
  vertical bends with horizontal bends.
- `Navigator.add_physics_process` should now be a lot faster by doing fewer unnecessary deepcopies.

### Added

- New properties on `ParticleArrray`, `beta` and `gamma`, providing the relativistic
  properties of each particle.
- New method `sort` on `ParticleArray` for possibly sorting in place with respect to one
  of the properties or functions (e.g. `beta`, `x`, `p0c`, etc.) returning the indices
  that sort the `ParticleArray`.
- Implementation of `__len__` for `ParticleArray` which simply returns the number of
  particles it contains.  Identical in function to `size`.
- `get_envelope` now also calculates the normalised emittances `emit_xn` and `emit_yn`,
  not just the geometric ones, as was previously the case.
- Added `overwrite_progress` kwarg to `cpbd.track.track` function which allows
  tracking progress to optionally be written on new lines rather than
  overwritten repeatedly using carriage returns. This is particularly useful
  when logging the output of a running ocelot simulation to avoid extremely long
  lines and unreadable log files.
- Inactive processes in `Navigator` instances are now stored in the
  `inactive_proccses` attribute. This way processes attached to a `Navigator`
  remain accessible in that for the lifetime of the Navigator.
- New `ParameterScanner` class in `cpbd.track` for scanning arbitrary parameters
  in parallel (either using multiprocessing mpi4py). For example one might scan
  different compression schemes. Results are compiled into a single hdf5 file.
- Added new method `Navigator.add_physics_processes` for adding multiple physics processes at the same time.  This will be a lot fastwe when lots of physics processes are to be added.

- `Navigator.jump_to` method allowing a `Navigator` instance to jump
  to arbitrary points in z along the magnetic lattice.  Useful as it
  does not require modifying the underling beamling to achieve
  equivalent behaviour.
- new `__repr__` methods for some common physics processes: `CSR`, `SmoothBeam`, `SpaceCharge` and `WakeTable`.
  

### Fixed

- None yet.

### Removed

- None yet.

### Deprecated

- None yet.

[unreleased]: https://github.com/ocelot-collab/ocelot/compare/dev_2021..HEAD
