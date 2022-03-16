# Changelog

All notable changes to OCELOT will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
~and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).~ (Not yet, maybe it won't.  Versioning is based on year.month.day.

## [Unreleased]

### Changed

- The `CSR` class now gives an explicit `logging` `error` when failing in `arcline` and
  a more useful explanation.  This error often occurs because the user has mixed CSR for
  vertical bends with horizontal bends.

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

### Fixed

- Fix rare crash in CSR process.

### Removed

- None yet.

### Deprecated

- None yet.

[unreleased]: https://github.com/ocelot-collab/ocelot/compare/dev_2021..HEAD
