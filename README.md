Accelerator, radiation and x-ray optics simulation framework

# Ocelot

Ocelot is an open-source multiphysics simulation toolkit for accelerator physics.
It is designed to model Free Electron Lasers (FELs), storage rings, and transport lines using a modular and scriptable Python interface.

- **Website**: [https://www.ocelot-collab.com](https://www.ocelot-collab.com)
- **Documentation**: [https://www.ocelot-collab.com/docs/docu/intro](https://www.ocelot-collab.com/docs/docu/intro)
- **Tutorials**: [https://www.ocelot-collab.com/docs/tutorial/intro](https://www.ocelot-collab.com/docs/tutorial/intro)
- **Source code**: [https://github.com/ocelot-collab/ocelot](https://github.com/ocelot-collab/ocelot)
- **Bug reports**: [https://github.com/ocelot-collab/ocelot/issues](https://github.com/ocelot-collab/ocelot/issues)
- **License**: [GPL-3.0 license](https://github.com/ocelot-collab/ocelot/blob/master/LICENSE)
- **Agent guide**: [AGENTS.md](./AGENTS.md)

---

## Features

Ocelot provides:

- A modular framework for beam dynamics simulations (tracking, optics, matching)
- Physics processes including:
  - Space charge
  - Coherent synchrotron radiation (CSR)
  - Wakefields
  - and many more
- A framework for FEL-related studies and synchrotron calculation
- Jupyter-based interactive tutorials for education and development


---

## Getting Started

For requirements and installation instructions, see the official guide:
👉 [**Installation & Setup**](https://www.ocelot-collab.com/docs/docu/intro)

To explore tutorials, visit:
👉 [**Tutorial Overview**](https://www.ocelot-collab.com/docs/tutorial/intro)
👉 [**Student-Friendly Introduction**](https://www.ocelot-collab.com/docs/tutorial/tutorial-beam-dynamics/for_students)

## Import and Startup Notes

Legacy tutorial code using `from ocelot import *` remains supported, but it is
not recommended for new scripts, tests, demos, or generated lattice files.

For user scripts and tutorials, prefer:

```python
import ocelot as ocl

d = ocl.Drift(l=0.5)
qf = ocl.Quadrupole(l=0.2, k1=0.3)
qdh = ocl.Quadrupole(l=0.1, k1=-0.3)
lat = ocl.MagneticLattice((qdh, d, d, qf, d, d, qdh))
tws = ocl.periodic_twiss(lat)
```

This keeps one memorable namespace while avoiding star-import side effects. For
library code, tests, and agent-authored changes, prefer explicit imports from
the submodule that owns the API:

```python
from ocelot.cpbd.elements import Drift, Quadrupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.optics import periodic_twiss, twiss
from ocelot.cpbd.beam import Twiss, ParticleArray, generate_parray
```

The root facade is curated and lazy-loaded, so short-lived scripts can start
quickly without importing plotting, radiation, pandas, SciPy, or optional
acceleration modules until those features are actually used.

`twiss(lat, tws0)` propagates explicitly supplied initial Twiss parameters.
Use `periodic_twiss(lat, tws0=None)` when the lattice should determine the
initial periodic optics. An unstable lattice raises `UnstableLatticeError`
instead of logging a warning and returning `None`:

```python
from ocelot.cpbd.optics import UnstableLatticeError, periodic_twiss

try:
    tws = periodic_twiss(lat)
except UnstableLatticeError as exc:
    print(f"Periodic optics are unavailable: {exc}")
```

## Repeated Elements and Position Lookups

A lattice may reuse the same element object at several positions. This is
useful when those placements intentionally share parameters, but the object by
itself then does not identify one position:

```python
import ocelot as ocl

d = ocl.Drift(l=1.0, eid="D")
q = ocl.Quadrupole(l=0.3, k1=1.0, eid="Q")
lat = ocl.MagneticLattice((d, q, d))

lat.find_element_indices(d)                 # [0, 2]
lat.resolve_element_index(d, occurrence=0) # 0
lat.resolve_element_index(d, occurrence=1) # 2
```

APIs that accept only an element object and require one lattice position raise
`ValueError` when that instance is repeated instead of silently choosing its
first or last occurrence. Element IDs are not used for this check because IDs
need not be unique.

Attaching calculated Twiss parameters directly to elements remains available
as a convenience for small scripts. Attach only selected unique elements when
the lattice contains reused objects:

```python
tws0 = ocl.Twiss(beta_x=10.0, beta_y=12.0)
tws = ocl.twiss(lat, tws0, attach2elem=[q])
print(q.tws.beta_x)  # Twiss at the exit of the unique q occurrence
```

`attach2elem=True` requires every element instance in the lattice to be unique.
A selected repeated instance also raises because one `element.tws` attribute
cannot represent several locations. Matcher results can address deliberate
repetitions explicitly with `state.twiss_at(element, occurrence=n)`.

---

## Core Modules & API Reference

Ocelot's core functionality is organized into key modules:

### 📘 Lattice Elements & Design
- [**Elements**](./ocelot/cpbd/elements/) - Lattice element definitions (dipoles, quadrupoles, cavities, etc.)
  - `OpticElement`, `Element`, `Magnet` - Base classes with inheritance hierarchy
  - `Drift`, `Bend`, `RBend`, `SBend`, `Quadrupole`, `Sextupole`, `Octupole`
  - `Cavity`, `TDCavity`, `Undulator`, `Marker`, `Monitor`, `Aperture`

### 🔬 Beam Physics & Tracking
- [**Beam**](./ocelot/cpbd/beam/) - Beam models and particle arrays
  - `Beam`, `Twiss` - Beam envelopes and Twiss parameters
  - `ParticleArray`, `Particle` - Individual particle tracking
  - `generate_parray()`, `ellipse_from_twiss()` - Beam generation utilities

### 🎯 Tracking & Optics
- [**Tracking**](./ocelot/cpbd/track.py) - Particle and beam tracking algorithms
- [**Optics**](./ocelot/cpbd/) - Optics calculations and transfer maps
  - `MagneticLattice` - Core lattice object for simulations
  - `Navigator` - Lattice navigation utilities

### ⚡ Physics Processes
- [**Space Charge**](./ocelot/cpbd/sc.py) - Space charge effects
- [**CSR**](./ocelot/cpbd/csr.py) - Coherent synchrotron radiation
- [**Wake Fields**](./ocelot/cpbd/wake3D.py) - Longitudinal and transverse wake effects
- [**Physics Processes**](./ocelot/cpbd/physics_proc.py) - Extensible framework for custom effects

### 🔧 Advanced Features
- [**Chromaticity Compensation**](./ocelot/cpbd/chromaticity.py)
- [**Beam Matching**](./ocelot/cpbd/match.py) - Automated optics matching
- [**Beam Parameter Calculations**](./ocelot/cpbd/beam_params.py)

---

## We welcome feedback, contributions, and new ideas



- [How to create a pull request](https://www.ocelot-collab.com/docs/docu/how-to/pull_request)
- [How to add unit tests](https://www.ocelot-collab.com/docs/docu/how-to/unit_test)
- [How to create your own Physics Process class](https://www.ocelot-collab.com/docs/docu/how-to/phys_proc)

**Disclaimer:** The OCELOT code comes with absolutely NO warranty. The authors of the OCELOT do not take any responsibility for any damage to equipments or personnel injury that may result from the use of the code.
