# Autodiff-Friendly Linear Transfer Matrices

## Goal

Add a minimal autodiff-friendly path for linear transfer matrices in Ocelot without changing the normal lattice-building workflow and without rewriting CPBD into JAX.

The intended user workflow stays the same:

```python
from ocelot import Drift, Quadrupole, Bend, MagneticLattice

q = Quadrupole(l=0.3, k1=1.2)
b = Bend(l=1.0, angle=0.1, e1=0.05, e2=0.05)
d = Drift(l=0.8)

cell = [q, d, b, d]
lat = MagneticLattice(cell)
```

The new work is about exposing the linear `R` matrices behind these elements in a way that can later be used with `jax.numpy` and autodiff.

Current scope:

- `Drift`
- `Quadrupole`
- `Bend` / `SBend` / `RBend`
- `Cavity`

## Non-Goal

This project does **not** try to make the full current Ocelot runtime differentiable.

In particular, this work does not try to autodiff through:

- `TransferMap`
- `MagneticLattice`
- `track()`
- `match.py`
- `matcher.py`

Those layers currently rely on mutable objects, caches, and the standard Ocelot execution model.

## Current Direction

The shared linear matrix formulas are being centralized in:

- [r_matrix.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/r_matrix.py)

This file now contains:

- backend-aware linear magnet-body matrix construction
- shared bend-edge matrix construction
- one place for the first-order linear physics formulas

The key design rule is:

> Define the physics formulas once, and let both NumPy and `jax.numpy` call the same low-level implementation.

## What Is Already Implemented

This is no longer only a design sketch. The branch already contains a working first implementation.

Implemented pieces:

- shared backend-aware first-order matrix helpers in [r_matrix.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/r_matrix.py)
- public element-level linear matrix API in [optic_element.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/elements/optic_element.py)
- atom-level support for magnets, bends, and cavities
- pure backend-aware Twiss helpers in [twiss_linear.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/twiss_linear.py)
- reuse of the same Twiss formulas from the legacy wrapper path in [core.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/beam/core.py) and [optics.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/optics.py)
- two tutorial notebooks in `demos/ipython_tutorials`

In practical terms, the branch now supports:

- element-level `R(parameter)` access without mutating the live lattice
- boundary-to-boundary linear Twiss tracing from explicit initial Twiss state
- JAX `grad(...)` through several real linear optics examples
- cavity participation in the pure linear Twiss path

## Implemented Files And Responsibilities

### Low-level linear physics

- [r_matrix.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/r_matrix.py)

This file now contains the shared first-order kernels:

- `linear_magnet_matrix(...)`
- `bend_edge_matrix(...)`
- `cavity_coupler_edge_matrix(...)`
- `standing_wave_cavity_matrix(...)`
- `rot_mtx(..., xp=...)`

These are the formulas intended to be shared by both the NumPy path and the future JAX path.

### Public element wrappers

- [optic_element.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/elements/optic_element.py)

This wrapper layer now exposes:

- `linear_r(...)`
- `linear_r_blocks(...)`
- `linear_r_full(...)`
- `linear_delta_e_blocks(...)`
- `linear_length_blocks(...)`

The override mechanism works by copying the underlying atom and applying temporary keyword overrides, so a user can evaluate `R(parameter)` without mutating the live element in the lattice.

### Atom-level family support

- [element.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/elements/element.py)
- [magnet.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/elements/magnet.py)
- [bend_atom.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/elements/bend_atom.py)
- [cavity_atom.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/elements/cavity_atom.py)

Current family status:

- `Drift`: supported through the shared magnet-body path
- `Quadrupole`: supported through the shared magnet-body path
- `Bend` / `SBend` / `RBend`: supported for body and edge blocks
- `Cavity`: supported for entrance/body/exit blocks plus block-wise energy gain

### Pure linear Twiss layer

- [twiss_linear.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/twiss_linear.py)

This module currently provides:

- `make_twiss_state(...)`
- `twiss_step_from_r(...)`
- `trace_linear_twiss(...)`
- `periodic_twiss_from_r(...)`

This is the pure linear optics path intended for autodiff experiments.

### Legacy wrapper reuse

- [core.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/beam/core.py)
- [optics.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/optics.py)

The branch also reuses the extracted formulas from the legacy wrapper layer:

- `Twiss.track(...)` now goes through the shared one-step linear helper
- `Twiss.map_x_twiss(...)` now uses the same underlying step logic
- `periodic_twiss(...)` now reuses the shared periodic solver

That means the new work is not isolated in a side module; some of it is already feeding the standard Ocelot optics wrapper path.

## Additional Refactoring Done

Two extra cleanup steps were also done because they help the shared-formula approach:

- bend-edge first-order matrices were centralized in one helper instead of staying duplicated in several places
- cavity linear block handling was moved into the same shared-matrix style as magnets and bends

This matters because the whole point of the autodiff work is to avoid keeping separate NumPy and JAX copies of the same physics formulas.

## Implemented Examples And Notebooks

Two tutorial notebooks were added:

- [matcher_basics.ipynb](/Users/tomins/Nextcloud/DESY/repository/ocelot/demos/ipython_tutorials/matcher_basics.ipynb)
- [matcher_autodiff.ipynb](/Users/tomins/Nextcloud/DESY/repository/ocelot/demos/ipython_tutorials/matcher_autodiff.ipynb)

### `matcher_basics.ipynb`

This notebook repeats simple working matcher-style examples from the matcher documentation using small self-contained lattices:

- vary one drift length to hit a Twiss target
- vary cavity voltage to hit end energy
- target one `R`-matrix element between internal markers
- use one linked knob for two quadrupoles
- use a simple custom objective

### `matcher_autodiff.ipynb`

This notebook shows the same style of problems solved through the new autodiff-friendly linear path:

- drift-length matching through differentiable Twiss propagation
- cavity-voltage matching through differentiable end-energy objective
- linked quadrupole matching through one explicit scalar variable
- custom scalar objective
- a larger four-quadrupole end-Twiss matching problem with before/after beta traces

It also includes one deliberate limitation example:

- direct drift `R12` targeting currently gives `NaN` JAX gradients in the exact drift corner

That notebook is important because it shows both what already works and what still needs low-level cleanup.

## Current Practical Status

What works now:

- `jax.grad(...)` through several element-level and Twiss-level examples
- explicit parameter overrides such as `q.linear_r(k1=..., xp=jnp)`
- explicit lattice objectives built from `trace_linear_twiss(...)`
- simple SciPy optimization loops using exact JAX gradients
- cavity energy changes in the pure linear optics path

What is not done yet:

- full autodiff through `Matcher`
- full autodiff through `track()`
- full autodiff through `TransferMap`
- a pure lattice-segment `R(start, stop, x)` model matching the current `lat.transfer_maps(...)` API
- a dedicated `LinearLatticeModel` object

## Current Limitations

The main current limitations are technical, not conceptual.

### 1. Some wrapper code is not JIT-safe yet

Plain `jax.grad(...)` already works in several examples, but full `jax.jit(...)` over the current wrapper path is not reliable yet.

Reason:

- some wrapper/element paths still contain Python conditionals on traced values
- for example, simple checks such as `if self.l == 0` are acceptable for normal NumPy execution, but not for JAX tracing

So the current milestone should be understood as:

- autodiff-friendly: partly yes
- JIT-friendly: not fully yet

### 2. Exact drift corner still has a bad JAX derivative

The direct drift `R12` example currently exposes a real low-level limitation:

- for an exact drift, the current backend-generic magnet-body helper still passes through a `k -> 0` corner
- the matrix value is correct
- but the JAX derivative can become `NaN`

This is not a high-level API problem. It is a low-level differentiability problem inside the shared magnet matrix helper.

### 3. The current pure path is boundary-based

`trace_linear_twiss(...)` currently traces one state per element boundary and internally applies the block sequence for each element.

That is enough for the current milestone, but it is still a simple pure tracing layer, not a full replacement for the richer CPBD runtime.

## Validation Done So Far

The implemented examples were checked in two ways:

- contract/unit tests for the new linear API and pure Twiss helpers
- direct execution of the tutorial notebook code cells

In particular:

- `matcher_basics.ipynb` was executed cell-by-cell
- `matcher_autodiff.ipynb` was executed cell-by-cell

So the README examples are not only aspirational; they were written against working code paths in the current branch.

## Why Not Reuse `elem.R(energy)` Directly?

Ocelot already exposes:

```python
elem.R(energy)
```

but this is not the best autodiff API.

Reasons:

1. `R()` is part of the existing transfer-map wrapper API.
2. `R()` returns the evaluated Ocelot transfer-map blocks for the current mutable object state.
3. `R()` is tied to the current CPBD object model and cached map machinery.
4. `R()` does not accept a backend argument such as `xp=jax.numpy`.
5. `R()` is less explicit about whether it returns body-only, edge-only, or full composed matrix.

So `R()` should remain the legacy/runtime-facing API, while the autodiff-facing path should be an explicit new API.

## Implemented Element-Level API

The public wrapper API is now:

```python
elem.linear_r(...)
elem.linear_r_blocks(...)
elem.linear_r_full(...)
```

Meaning:

- `linear_r(...)`: rotated main/body linear `R`
- `linear_r_blocks(...)`: rotated list of linear blocks, e.g. entrance/body/exit
- `linear_r_full(...)`: one rotated composed 6x6 linear matrix for the whole element

These methods should be thin wrappers over the shared kernels in `r_matrix.py`.

For energy-changing elements such as `Cavity`, the block API also carries the
same block ordering as the existing first-order TM sequence. This matters for
pure Twiss tracing because the reference energy changes only on the `MAIN`
block.

## Important API Detail

`linear_r(...)` should return a matrix, not a function object.

For autodiff, the callable is the method itself:

```python
R = q.linear_r(k1=1.5, energy=0.0, xp=jax.numpy)
```

That is already a differentiable function call from JAX's point of view.

## Default State vs Explicit Overrides

For ordinary Ocelot use, the element keeps its stored parameters:

```python
q = Quadrupole(l=0.3, k1=1.2)
R0 = q.linear_r()
```

For autodiff and optimization, the same method should allow explicit overrides:

```python
R1 = q.linear_r(k1=1.5)
R2 = q.linear_r(l=0.4, k1=1.5)
```

This avoids mutating the element inside the objective function.

### Drift example

```python
import numpy as np
from ocelot import Drift

d = Drift(l=1.0)

R0 = d.linear_r()
R1 = d.linear_r(l=1.2)

assert np.isclose(R0[0, 1], 1.0)
assert np.isclose(R1[0, 1], 1.2)
assert np.allclose(R0, d.R(0.0)[0])
```

### Quadrupole example

```python
import numpy as np
from ocelot import Quadrupole

q = Quadrupole(l=0.3, k1=1.2)

R0 = q.linear_r()
R1 = q.linear_r(k1=1.5)
R2 = q.linear_r(l=0.4, k1=1.5)

q_ref = Quadrupole(l=0.4, k1=1.5)
assert np.allclose(R2, q_ref.R(0.0)[0])
```

### Bend example

```python
import numpy as np
from ocelot import Bend

b = Bend(l=1.0, angle=0.1, e1=0.05, e2=0.05)

Rb = b.linear_r()                     # body
blocks = b.linear_r_blocks()         # entrance/body/exit
Rf = b.linear_r_full()               # composed full-element R

assert np.allclose(Rb, b.R(0.0)[1])
assert len(blocks) == 3

Rf_manual = blocks[0]
for block in blocks[1:]:
    Rf_manual = block @ Rf_manual

assert np.allclose(Rf, Rf_manual)
```

### Cavity example

```python
import numpy as np
from ocelot import Cavity

c = Cavity(l=0.35, v=0.02, phi=20.0, freq=1.3e9, eid="C1")

blocks = c.linear_r_blocks(energy=0.13)
delta_e_blocks = c.linear_delta_e_blocks()

assert len(blocks) == 3
assert np.allclose(blocks[1], c.linear_r(energy=0.13))
assert np.isclose(sum(delta_e_blocks), c.v * np.cos(np.deg2rad(c.phi)))
```

## Minimal JAX Example

```python
import jax
import jax.numpy as jnp
from ocelot import Quadrupole

q = Quadrupole(l=0.3, k1=1.2)

def objective(k1):
    R = q.linear_r(k1=k1, energy=0.0, xp=jnp)
    return R[0, 1]

grad_k1 = jax.grad(objective)(1.2)
```

This is the intended first milestone: element-level `R(parameter)` support using shared formulas.

## Real Lattice Example

The normal lattice-construction workflow does not change:

```python
import numpy as np
from ocelot import Bend, Drift, MagneticLattice, Quadrupole

q = Quadrupole(l=0.3, k1=1.2)
d = Drift(l=0.8)
b = Bend(l=1.0, angle=0.1, e1=0.05, e2=0.05)
lat = MagneticLattice([q, d, b, d])

energy = 0.0
R_lat = np.eye(6)
for elem in lat.sequence:
    R_lat = elem.linear_r_full(energy=energy) @ R_lat
```

The same idea can later be turned into a pure optimization model by keeping
`lat` as structure and passing explicit override values instead of mutating
element attributes inside the objective.

## Shared Twiss Helpers

The same pattern now also applies to linear Twiss transport:

- shared pure formulas live in [twiss_linear.py](/Users/tomins/Nextcloud/DESY/repository/ocelot/ocelot/cpbd/twiss_linear.py)
- the legacy `Twiss.track(...)` and `periodic_twiss(...)` wrappers reuse those formulas
- a new autodiff-oriented path can call the same helpers directly with `xp=jax.numpy`

The key helpers are:

```python
from ocelot.cpbd.twiss_linear import (
    make_twiss_state,
    periodic_twiss_from_r,
    trace_linear_twiss,
    twiss_step_from_r,
)
```

For energy-changing elements, `trace_linear_twiss(...)` applies the linear
blocks in sequence together with the corresponding block-wise `delta_e`. That
is why a cavity can now participate in the pure Twiss path without reusing the
mutable `TransferMap` runtime.

## Twiss At Element Boundaries

The new tracer works at element boundaries and returns one state per element,
plus the starting state.

```python
from ocelot import Bend, Drift, MagneticLattice, Quadrupole
from ocelot.cpbd.twiss_linear import make_twiss_state, trace_linear_twiss

q = Quadrupole(l=0.3, k1=1.2, eid="Q1")
d = Drift(l=0.8, eid="D1")
b = Bend(l=1.0, angle=0.1, e1=0.05, e2=0.05, eid="B1")
lat = MagneticLattice([q, d, b, d])

tws0 = make_twiss_state(
    beta_x=12.0,
    alpha_x=0.5,
    beta_y=15.0,
    alpha_y=-0.2,
    Dx=0.1,
    Dxp=0.0,
    E=0.0,
)

states = trace_linear_twiss(lat, tws0)

beta_x_start = states[0]["beta_x"]
beta_x_after_q = states[1]["beta_x"]
beta_x_after_b = states[3]["beta_x"]
```

If you want the states together with the element names:

```python
labels = ["START"] + [elem.id for elem in lat.sequence]
for label, state in zip(labels, states):
    print(label, state["s"], state["beta_x"], state["Dx"])
```

## One-Step Twiss Update

You can also update Twiss state through a single explicit matrix:

```python
from ocelot import Quadrupole
from ocelot.cpbd.twiss_linear import make_twiss_state, twiss_step_from_r

q = Quadrupole(l=0.3, k1=1.2)
R = q.linear_r_full()

tws0 = make_twiss_state(beta_x=12.0, alpha_x=0.5, beta_y=15.0, alpha_y=-0.2)
tws1 = twiss_step_from_r(R, tws0, length=q.l)
```

## Periodic Twiss From One-Turn Matrix

For a stable linear lattice, periodic Twiss can be obtained directly from the
one-turn `R` matrix:

```python
from ocelot import Drift, MagneticLattice, Quadrupole
from ocelot.cpbd.twiss_linear import periodic_twiss_from_r

qf = Quadrupole(l=0.2, k1=0.3)
qd = Quadrupole(l=0.1, k1=-0.3)
d = Drift(l=0.5)
lat = MagneticLattice([qd, d, d, d, qf, d, d, d, qd])

R_turn = lat.transfer_maps(energy=0.005)[1]
tws_periodic = periodic_twiss_from_r(R_turn, energy=0.005)

beta_x0 = tws_periodic["beta_x"]
alpha_y0 = tws_periodic["alpha_y"]
```

## Differentiable Twiss Example

This is the first practical autodiff-style pattern:

```python
import jax
import jax.numpy as jnp
from ocelot import Drift, MagneticLattice, Quadrupole
from ocelot.cpbd.twiss_linear import make_twiss_state, trace_linear_twiss

q1 = Quadrupole(l=0.3, k1=1.2, eid="Q1")
q2 = Quadrupole(l=0.3, k1=-0.9, eid="Q2")
d = Drift(l=0.8, eid="D1")
lat = MagneticLattice([q1, d, q2, d])

tws0 = make_twiss_state(beta_x=12.0, alpha_x=0.5, beta_y=15.0, alpha_y=-0.2, xp=jnp)

def objective(k1):
    states = trace_linear_twiss(lat, tws0, xp=jnp, overrides={q1: {"k1": k1}})
    return states[-1]["beta_x"]

grad_beta_x = jax.grad(objective)(1.2)
```

With a cavity in the lattice:

```python
import jax
import jax.numpy as jnp
from ocelot import Cavity, Drift, MagneticLattice, Quadrupole
from ocelot.cpbd.twiss_linear import make_twiss_state, trace_linear_twiss

q1 = Quadrupole(l=0.25, k1=0.8, eid="Q1")
c1 = Cavity(l=0.35, v=0.02, phi=20.0, freq=1.3e9, eid="C1")
d1 = Drift(l=0.5, eid="D1")
lat = MagneticLattice([d1, c1, q1, d1])

tws0 = make_twiss_state(beta_x=11.0, alpha_x=0.4, beta_y=13.0, alpha_y=-0.3, E=0.13, xp=jnp)

def objective(v):
    states = trace_linear_twiss(lat, tws0, xp=jnp, overrides={c1: {"v": v}})
    return states[-1]["beta_x"]

grad_beta_x = jax.grad(objective)(0.02)
```

This keeps:

- the lattice structure in normal Ocelot objects
- the linear optics formulas in one shared implementation
- the differentiated values explicit in the function arguments

It does **not** try to autodiff through the full current `twiss()` runtime.

## Lattice-Level Direction

For future optimization, the preferred design is not to autodiff through the current matcher internals directly.

Instead, build a small pure linear model on top of:

- `lat`
- a list of vary specifications
- a parameter vector `x`
- a backend `xp`

Conceptually:

```python
model = LinearLatticeModel(lat, energy=0.13)
model.vary_element(q, quantity="k1", name="Q1")
model.vary_element(d, quantity="l", name="D1")

x0 = model.x0()
R = model.r_matrix(x0, xp=jax.numpy)
```

This matches the mental model of Ocelot matching:

- keep `lat` as the structural input
- describe "what to shake"
- evaluate a linear model from an explicit parameter vector

## Reusing Matcher Concepts

The newer matcher already exposes a good user-facing idea:

- `vary_element(...)`
- `vary_linked_elements(...)`

That interface style should be reused conceptually.

What should **not** be reused directly for autodiff is the current mutation-based evaluation flow.

The autodiff-friendly version should be:

- pure
- explicit in inputs
- based on a parameter vector
- backend-aware via `xp`

## Naming Notes

`r_func()` is possible, but it is misleading because the method returns a matrix, not a function object.

Better names are those that say:

- this is about the `R` matrix
- this is linear optics
- this is distinct from the existing `R()`

Recommended names:

1. `linear_r`
2. `linear_r_full`
3. `linear_r_blocks`

Other acceptable options:

1. `r_matrix_linear`
2. `r_linear`
3. `linear_matrix`
4. `r_map_linear`
5. `r_eval`

Less recommended:

1. `r_func`
2. `get_jax_r`
3. `get_r`

Why less recommended:

- `r_func` suggests returning a callable
- `get_jax_r` hardcodes one backend into the API
- `get_r` is too close to the existing `R()`

## Recommended Naming Choice

Use:

- `linear_r(...)`
- `linear_r_blocks(...)`
- `linear_r_full(...)`

This is short, explicit, and consistent with the existing `R()` while still clearly distinguishing the new API.

## Current Practical Milestone

The current first milestone is:

1. shared backend-aware first-order formulas in `r_matrix.py`
2. element-level explicit linear `R` access
3. future pure lattice-level `LinearLatticeModel`
4. later optional JAX-based optimization layer

This gives a minimal bridge from current Ocelot usage to autodiff-friendly linear optics, without a major rewrite of the core codebase.
