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
