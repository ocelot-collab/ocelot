# Agent Guide for Ocelot

This file is the starting point for coding agents working in this repository.
It is intentionally compact: use it to find the right modules, examples, and
tests before making changes.

## Project Context

Ocelot is a Python toolkit for accelerator and photon simulations, especially
FELs, storage rings, and transport lines. Public user documentation lives at:

- Website: https://www.ocelot-collab.com
- Documentation: https://www.ocelot-collab.com/docs/docu/intro/
- Tutorials: https://www.ocelot-collab.com/docs/tutorial/intro/
- Local documentation checkout:
  `/Users/tomins/Nextcloud/DESY/repository/ocelot-collab.github.io`

Use that local checkout whenever a source change requires public documentation
or generated tutorial Markdown updates; the user does not need to provide its
path again.

The online documentation is useful for workflow orientation, but source code,
docstrings, demos, and tests are the most reliable implementation references.
If docs and code disagree, follow the code and existing tests.

## Repository Map

- `ocelot/`: importable package.
- `ocelot/cpbd/`: charged-particle beam dynamics: lattices, elements,
  tracking, optics, matching, collective effects, wakefields, and physics
  processes.
- `ocelot/cpbd/elements/`: public accelerator element wrappers and atom
  classes. This is the main entry point for `Drift`, `Quadrupole`, `Bend`,
  `Cavity`, `Undulator`, monitors, correctors, and apertures.
- `ocelot/cpbd/transformations/`: transfer-map and tracking-method classes.
- `ocelot/cpbd/tm_params/`: typed parameter containers passed from elements to
  transformations.
- `ocelot/cpbd/beam/`: `Twiss`, `Beam`, `Particle`, `ParticleArray`, beam
  generation, and beam analysis helpers.
- `ocelot/rad/`: synchrotron/FEL radiation calculations.
- `ocelot/optics/`: photon optics and wavefront utilities.
- `ocelot/adaptors/`: import/export adapters for external tools and file
  formats.
- `ocelot/gui/`: plotting and GUI helpers. Keep GUI imports out of core
  simulation paths.
- `demos/`: runnable examples. `demos/ipython_tutorials/` mirrors the public
  tutorial workflows; `demos/ebeam/` and `demos/sr/` are good source-level
  examples.
- `unit_tests/`: regression and architecture tests. New behavior should usually
  get a focused test here.

## Accelerator Workflow Pointers

For a simple lattice or optics task, start with:

- Elements: `ocelot.cpbd.elements`
- Lattice container: `ocelot.cpbd.magnetic_lattice.MagneticLattice`
- Linear optics: `ocelot.cpbd.optics.twiss` and
  `ocelot.cpbd.optics.periodic_twiss`
- Tracking: `ocelot.cpbd.track.track`
- Beam objects: `ocelot.cpbd.beam.Twiss`, `ParticleArray`, `generate_parray`
- Navigation and physics process scheduling:
  `ocelot.cpbd.navi.Navigator` and `ocelot.cpbd.physics_proc`

Minimal example shape:

```python
from ocelot.cpbd.elements import Drift, Quadrupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.optics import periodic_twiss

d = Drift(l=0.5)
qf = Quadrupole(l=0.2, k1=0.3)
qdh = Quadrupole(l=0.1, k1=-0.3)
cell = (qdh, d, d, qf, d, d, qdh)
lat = MagneticLattice(cell)
tws = periodic_twiss(lat)
```

For workflow examples, prefer these before inventing new patterns:

- Linear optics and lattice design: `demos/ipython_tutorials/1_introduction.ipynb`,
  `demos/ipython_tutorials/7_lattice_design.ipynb`, `demos/ebeam/dba.py`
- Tracking and Runge-Kutta examples: `demos/ipython_tutorials/2_tracking.ipynb`,
  `demos/docs/18_runge_kutta_tracking.ipynb`, `demos/ebeam/rk_vs_matrix.py`
- Space charge, wake, CSR, and laser heater workflows:
  `demos/ipython_tutorials/3_space_charge.ipynb`,
  `demos/ipython_tutorials/4_wake.ipynb`,
  `demos/ipython_tutorials/5_CSR.ipynb`,
  `demos/ipython_tutorials/8_laser_heater.ipynb`
- Synchrotron radiation and photon field simulations:
  `demos/ipython_tutorials/pfs_1_synchrotron_radiation.ipynb`,
  `demos/sr/`, `demos/optics/`

## Architecture Notes

The CPBD element implementation uses a wrapper/atom/parameter/transformation
structure:

1. Public wrapper: user-facing element class, usually in
   `ocelot/cpbd/elements/*.py`.
2. Atom: physics state and `create_*_params(...)` hooks, often named
   `*_atom.py`.
3. TMParams: data objects in `ocelot/cpbd/tm_params/`.
4. Transformation: tracking algorithm in `ocelot/cpbd/transformations/`.

When adding or changing an element, check the architecture-contract tests in
`unit_tests/cpbd/architecture_contract/`. Preserve both the active tracking
method path and the first-order optics path unless the existing contract says
otherwise.

## Element Identity and Lattice Occurrences

The same element object may intentionally appear more than once in
`MagneticLattice.sequence`, for example to share a magnet strength or reuse
equal-length drifts. Repeated instances are valid for ordered propagation, but
an element object alone does not identify one lattice position.

For location-sensitive code:

- Compare element occurrences by identity (`is`), not by `id` strings or
  physics parameters. Element IDs are not required to be unique.
- Use `MagneticLattice.find_element_indices(element)` to obtain all occurrence
  indices.
- Use `MagneticLattice.resolve_element_index(element)` when the caller requires
  a unique occurrence. It raises for absent or repeated instances. Pass the
  zero-based `occurrence` argument only in APIs that explicitly support it.
- Do not use `sequence.index(element)` or map an element object directly to one
  position-dependent value without first proving uniqueness.
- Keep repeated instances legal at lattice construction. Validate ambiguity at
  the API boundary that needs a unique position.

`twiss(..., attach2elem=True)` and `periodic_twiss(..., attach2elem=True)`
therefore require every attached element to be unique. To attach only selected
unique elements when other lattice objects are reused, pass an iterable such as
`attach2elem=[q1, b1]`. Matcher state uses
`state.twiss_at(element, occurrence=n)` for deliberate repeated occurrences.

## Twiss API

Keep propagation and periodic-solution calculation explicit:

- `twiss(lattice, tws0, ...)` propagates supplied initial Twiss parameters and
  requires positive `beta_x` and `beta_y`.
- `periodic_twiss(lattice, tws0=None, ...)` calculates the initial periodic
  solution and propagates it through the lattice. The optional seed supplies
  values such as energy and emittance.
- An unavailable periodic solution raises `UnstableLatticeError`; it does not
  log a warning or return `None`. Optimization code that deliberately explores
  unstable lattices should catch this exception and convert it to a penalty.
- `MagneticLattice.periodic_twiss(tws=None)` remains the lower-level API when
  only the periodic initial `Twiss` object is required.

## Import Guidance

`from ocelot import *` is a legacy tutorial style. Keep the public facade
working for existing user scripts, but do not add it to new source code, tests,
demos, generated lattice files, or documentation examples.

For user-facing scripts and tutorials, prefer the lazy root facade through a
namespace alias:

```python
import ocelot as ocl

d = ocl.Drift(l=0.5)
qf = ocl.Quadrupole(l=0.2, k1=0.3)
qdh = ocl.Quadrupole(l=0.1, k1=-0.3)
lat = ocl.MagneticLattice((qdh, d, d, qf, d, d, qdh))
tws = ocl.periodic_twiss(lat)
```

This gives users one namespace to remember without polluting the global
namespace or forcing every root export to load at startup.

For library code, tests, and agent-authored changes, prefer explicit submodule
imports so ownership and dependencies are clear:

```python
from ocelot.cpbd.elements import Drift, Quadrupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.optics import UnstableLatticeError, periodic_twiss, twiss
from ocelot.cpbd.beam import Twiss, ParticleArray, generate_parray
```

The root `ocelot` package and these CPBD package facades are intentionally
lazy-loaded:

- `ocelot/__init__.py`
- `ocelot/cpbd/beam/__init__.py`
- `ocelot/cpbd/elements/__init__.py`
- `ocelot/cpbd/tm_params/__init__.py`
- `ocelot/cpbd/transformations/__init__.py`

When adding a public facade name there, add it to `__all__` and the local
lazy-export map instead of importing the implementation at module import time.
This keeps `import ocelot`, narrow CPBD imports, and short-lived agent scripts
fast while preserving tutorial-facing names.

Avoid adding heavyweight imports to:

- `ocelot/__init__.py`
- `ocelot/cpbd/__init__.py`
- package `__init__.py` files that are imported by core workflows

Importing plotting, GUI, HDF5, pandas-heavy analysis, or optional acceleration
libraries at package import time makes every simulation startup slower. Prefer
function-local imports when the dependency is only needed by a specific feature.
This also applies to SciPy subpackages and `numba`: defer them until the
calculation that needs them.

## Testing

Useful focused commands:

```bash
python -m pytest unit_tests/cpbd -q
python -m pytest unit_tests/cpbd/architecture_contract -q
python -m pytest unit_tests/ebeam_test/dba/dba_test.py -q
python -m pytest unit_tests/sr_test -q
```

For import-related changes, measure fresh interpreter startup, not repeated
imports in one process:

```bash
python -X importtime -c "import ocelot"
python -c "import subprocess, sys, time; t=time.perf_counter(); subprocess.run([sys.executable, '-c', 'import ocelot']); print(time.perf_counter() - t)"
```

## Change Discipline

- Keep public APIs and tutorial-facing names stable unless the task explicitly
  asks for a breaking change.
- Prefer small changes with targeted tests over broad refactors.
- Check existing demos/tests for the workflow before introducing a new helper.
- Do not mix formatting-only churn with behavior changes.
