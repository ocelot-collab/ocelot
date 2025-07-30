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
- insert_markers_by_predicate now has optional arbitrary name suffixes for the markers
- Twiss.from_series allows for Twiss instances to be made from pd.Series (or a dictionary for example).
- ParticleArray slicing: `parray[10:100]` returns a ParticleArray instance with slices rparticles and q_array.
- `remove_coupler_kick` convenience method for Cavity class
- Integrated quad strengths for `k1` and `k2` convenience getters/setters (`k1l` and `k2l).
- Added new method `Navigator.add_physics_processes` for adding multiple physics processes at the same time.  This will be a lot fastwe when lots of physics processes are to be added.
>>>>>>> 6fdd98a24d013057498abd8e9ccf734939ef30e7
- `Navigator.jump_to` method allowing a `Navigator` instance to jump
  to arbitrary points in z along the magnetic lattice.  Useful as it
  does not require modifying the underling beamling to achieve
  equivalent behaviour.
- new `__repr__` methods for some common physics processes: `CSR`, `SmoothBeam`, `SpaceCharge` and `WakeTable`.
- `SmoothBeam` mslice attribute can now be set in the `__init__`.
- `extract_slice` method on SliceParameters instances
- modified the 'MagneticLattice.transfer_maps()' method. When using the 'output_at_each_step=True' flag,
  the method returns (Bs, Rs, Ts, S), where 'S' is a list of coordinates after each transfer map.
- New `WakeTable3` and `Wake3` class for third order Taylor expansion of wakefield tracking.
- New wake table for the wakefields of parallel plate structure based on analytical results, `WakeTableParallelPlate_origin`, `WakeTableParallelPlate`, 
  `WakeTableParallelPlate3_origin`, `WakeTableParallelPlate3`. The waketables with `origin` use point charge wake at the vicinity of the drive particle, which is the upper limit for wake. 
  The other two waketables contain exponential decay term. The waketables with `3` provide third order Taylor expansion and should be used with `Wake3`.

- added ParticleArray.get_twiss() method which calls get_envelop()
- Twiss recalculates gamma_x/y and emit_x/y automatically
- added link between quads to match function, e.g. vars = [{QF: 1.0, QD: -1.0}],
- Added dispersion auto correction in get_envelope function, it calculates dispersion from statistics and correct it if needed
- added image analysis functions to utils 
- added LSC for undulators 
- added in LSC feature to control size of the slice for transverse beam sizes calculation
- added new feature to show_density(). now it can plot particles with scatter plot and color of the individual particle will correspond to density - slower but maybe can be visual pleasing for presentations
  
>>>>>>> 6fdd98a24d013057498abd8e9ccf734939ef30e7

### Fixed

- fixed bugs in matching function if Drift length is in list of variables. 
- fixed bug with match function delta constrain
- fixed geometrical angle approximation in get_envelope function 
- fixed bug with navigator reset_position() method. now it also call physics process .prepare() method

### Removed

- None yet.

### Deprecated

- None yet.

[unreleased]: https://github.com/ocelot-collab/ocelot/compare/dev_2021..HEAD
