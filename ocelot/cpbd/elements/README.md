# CPBD Element Architecture: Current Contract and Safer Transition Plan

This note is for developers working inside `ocelot.cpbd`.

It has two goals:

1. explain how the current element machinery really works
2. describe a redesign path that improves clarity without breaking user-visible behavior

The public Ocelot docs already explain how to use elements such as `Bend`, `SBend`, `RBend`, `Cavity`, and `Quadrupole`. What is harder to understand from the outside is the internal contract between wrappers, atoms, parameter containers, and transformations. That contract is the important part to preserve before any architecture change.

## Bottom Line

The current architecture is not fundamentally wrong. In fact, one part of it is quite strong:

- element physics is separated from tracking algorithms
- the same element family can support more than one tracking method
- first-order optics can stay available even when active tracking uses a different method

The main problem is not the split itself. The main problem is that the split is too implicit.

Today a new developer has to discover several important rules by reading code:

- `OpticElement` is not just a thin wrapper; it owns cache invalidation, method selection, slicing, and public compatibility behavior
- many subsystems still expect the wrapper object and often also expect `wrapper.element`
- transformation support is declared inconsistently: sometimes by hook presence, sometimes by wrapper overrides, sometimes by warning-and-fallback behavior

That means the safest first step is not a direct redesign. The safest first step is:

1. document the current contract precisely
2. freeze it with unit tests
3. make capabilities explicit
4. only then try a simpler internal architecture

Users should not notice this work unless a future architecture delivers a clear benefit and preserves current physics results within tolerance.

## The Current Runtime Model

One logical beamline element is usually split into four parts:

1. a public wrapper in `elements/<name>.py`
2. a physics implementation object in `elements/<name>_atom.py`
3. a `TMParams` container in `tm_params/`
4. a tracking transformation in `transformations/`

In short:

- wrapper = public facade and framework behavior
- atom = physics state and element-specific formulas
- `TMParams` = typed data contract between element physics and tracking
- transformation = algorithm that applies the map

This is the central call chain:

```text
MagneticLattice.update_transfer_maps()
    -> element.set_tm(...)
    -> OpticElement._create_tms(...)
    -> Transformation.from_element(...)
    -> atom.create_*_params(...)
    -> TMParams
    -> Transformation.map_function(...)
```

## TM Roles

Before discussing element families, keep these four notions separate.

### 1. Linear optics / Twiss path

For current `OpticElement` wrappers, `first_order_tms` is always built with
`TransferMap`.

That path is used for:

- Twiss / linear optics code
- the always-available first-order cache
- helper accessors such as `R()` and `B()` unless the active method is `SecondTM`

So a family can have a real `TransferMap` optics path even when the wrapper
does not allow `TransferMap` to stay active as the tracking method.

### 2. Active tracking TM

This is `element.tms`, selected through `set_tm(...)` or lattice method
selection.

Examples:

- `TransferMap`
- `SecondTM`
- `KickTM`
- `RungeKuttaTM`
- `CavityTM`
- `MultipoleTM`

### 3. Atom hook surface

Atoms decide what is mechanically buildable through their
`create_*_params(...)` methods.

Examples:

- `Element` gives generic first-order and second-order hooks
- `Magnet` adds generic kick hooks
- `BendAtom` adds edge hooks and Runge-Kutta hooks
- `MultipoleAtom` adds a dedicated `MultipoleTM` hook while still keeping a first-order optics hook

This hook surface may be broader than the wrapper contract.

### 4. Meaning of `supported_tms`

Recommended meaning:

- `supported_tms` should list wrapper-selectable active tracking methods
- it should not be used to encode the always-available `TransferMap` optics path
- if a family keeps a `TransferMap` path only for `first_order_tms`, that should be documented separately

One more important point:

- `SecondTM` is not applied as "`TransferMap` plus a correction map"
- it is its own transformation class with `R`, `B`, and `T`
- many atoms build the `R`/`B` part by calling their own first-order helper internally, but runtime still applies only the `SecondTM`

## Public Element Capability Inventory

This inventory is based on current code behavior, not on an idealized future design.

The columns below deliberately separate:

- the always-available linear optics path
- wrapper-declared / wrapper-selectable active tracking methods
- additional hook-based or internal paths that exist today but should not be confused with the public wrapper contract

### Passive and utility families

| Family | Linear optics / Twiss path | Declared / selectable active TMs | Extra internal or legacy-buildable paths | Notes |
| --- | --- | --- | --- | --- |
| `Aperture` | `TransferMap` via generic `Element` first-order fallback | `TransferMap`, `SecondTM` | none | No dedicated aperture-specific transformation; mainly a holder for aperture metadata and offsets. |
| `Marker` | `TransferMap` via generic `Element` first-order fallback | `TransferMap`, `SecondTM` | none | Zero-length reference point. |
| `Monitor` | `TransferMap` via generic `Element` first-order fallback | `TransferMap`, `SecondTM` | none | Diagnostic state such as `x`, `y`, `x_ref`, `y_ref` lives on the atom. |
| `Matrix` | `TransferMap` from stored `R` / `B` data | `TransferMap`, `SecondTM` | none | Sliced sections currently fall back to a drift-like first-order map rather than slicing the stored matrix. |
| `UnknownElement` | `TransferMap` via inherited `Magnet` first-order hook | no audited `supported_tms` declaration yet | generic wrapper still lets `TransferMap`, `SecondTM`, and `KickTM` activate today | Legacy placeholder family rather than a clearly modeled element type. |
| `Pulse` | n/a | n/a | n/a | Separate helper object for time-dependent kicks; not part of the `OpticElement` wrapper/atom/TM stack. |

### Magnet-like families

| Family | Linear optics / Twiss path | Declared / selectable active TMs | Extra internal or legacy-buildable paths | Notes |
| --- | --- | --- | --- | --- |
| `Drift` | `TransferMap` via `DriftAtom` first-order hook | `TransferMap`, `SecondTM`, `KickTM`, `RungeKuttaTM`, `RungeKuttaTrTM` | none | Useful no-edge reference family. |
| `Quadrupole` | `TransferMap` via inherited `Magnet` first-order hook | `TransferMap`, `SecondTM`, `KickTM` | none | Wrapper also exposes convenience properties `k1l` and `k2l`. |
| `Sextupole` | `TransferMap` via inherited `Magnet` first-order hook | `TransferMap`, `SecondTM`, `KickTM` | none | `KickTM` is the important nonlinear tracking path here. |
| `Octupole` | `TransferMap` via inherited `Magnet` first-order hook | `TransferMap`, `KickTM` | inherited `SecondTM` can still be built on the generic wrapper path today, but it is not declared as audited octupole physics | `k3` only enters the dedicated kick path. |
| `Solenoid` | `TransferMap` via solenoid-specific first-order hook | `TransferMap`, `SecondTM` | none | `SecondTM` currently comes from the generic `Element` second-order fallback, not from a solenoid-specific nonlinear model. |
| `Hcor` / `Vcor` | `TransferMap` via `CorAtom` first-order hook | `TransferMap`, `SecondTM` | none | These are not `Magnet` subclasses, so they do not inherit generic kick hooks. |
| `XYQuadrupole` | `TransferMap` via `XYQuadrupoleAtom` first-order hook | `TransferMap` only | the atom also has `SecondTM` and inherited `KickTM` hook surface, but the wrapper intentionally pins tracking to first order | Important wrapper-level exception. |
| `Bend`, `SBend`, `RBend` | `TransferMap` via bend entrance/main/exit first-order hooks | `TransferMap`, `SecondTM`, `RungeKuttaTM`, `RungeKuttaTrTM` | `KickTM` can still be built mechanically through inherited `Magnet` hooks, but `kick.py` explicitly says it does not work for dipoles, so it is intentionally not declared | Edge maps are part of the real family behavior, not metadata. |

### RF, field-integrated, and special families

| Family | Linear optics / Twiss path | Declared / selectable active TMs | Extra internal or legacy-buildable paths | Notes |
| --- | --- | --- | --- | --- |
| `Cavity` | `TransferMap` via cavity-specific first-order entrance/main/exit hooks | `CavityTM` only | first-order `TransferMap` exists for optics but is not kept as the active tracking method | Most important complex reference family for edge handling and `delta_e`. |
| `TWCavity` | `TransferMap` via traveling-wave cavity first-order hooks | `TWCavityTM` only | first-order `TransferMap` exists for optics but is not kept as the active tracking method | The atom currently warns that this family is unfinished. |
| `TDCavity` | `TransferMap` via `TDCavityAtom` first-order hook | `TransferMap`, `SecondTM` | none | The atom's `additional_tms = [SecondTM]` is informational only; `OpticElement` does not consume it. |
| `Undulator` | `TransferMap` via `UndulatorAtom` first-order hook | `TransferMap`, `RungeKuttaTM`, `RungeKuttaTrTM`, `UndulatorTestTM` | none | `MagneticLattice.update_transfer_maps()` also has special-case length handling when a field map is attached. |
| `Multipole` | `TransferMap` via `MultipoleAtom` first-order hook | `MultipoleTM` only | first-order `TransferMap` exists for optics but the wrapper intentionally pins active tracking to `MultipoleTM` | This first-order path is a linearized multipole optics map, not a pure drift. |

## Layer 1: Public Wrapper

Typical examples:

- `elements/bend.py`
- `elements/cavity.py`
- `elements/quadrupole.py`

Most public classes subclass `OpticElement`.

That wrapper does more than just forward constructor arguments. It is responsible for:

- exposing the legacy public API
- forwarding many attributes to `self.element`
- caching instantiated transformations
- invalidating those caches when physics parameters change
- selecting active tracking methods
- building `ENTRANCE`, `MAIN`, and `EXIT` transformation lists
- building partial element slices for `get_section_tms(...)`

This point is easy to underestimate. The wrapper is currently part of the framework contract, not just syntactic sugar.

## Layer 2: Atom

Typical examples:

- `elements/bend_atom.py`
- `elements/cavity_atom.py`
- `elements/quadrupole_atom.py`

Atoms usually subclass `Element`, `Magnet`, or an existing family atom.

The atom owns:

- physics state such as `l`, `angle`, `k1`, `k2`, `v`, `freq`, `phi`
- geometry offsets and tilt
- edge-specific coefficients when needed
- hook methods that compute transformation parameters

Important hook families include:

- `create_first_order_*_params(...)`
- `create_second_order_*_params(...)`
- `create_cavity_tm_*_params(...)`
- `create_kick_*_params(...)`
- `create_runge_kutta_*_params(...)`
- `create_delta_e(...)`

The atom does not usually track particles directly. It computes the data a transformation needs.

## Layer 3: `TMParams`

Typical examples:

- `tm_params/first_order_params.py`
- `tm_params/second_order_params.py`
- `tm_params/cavity_params.py`
- `tm_params/kick_params.py`
- `tm_params/runge_kutta_params.py`

These objects are the explicit contract between the atom and the transformation.

Examples:

- `FirstOrderParams` carries `R`, `B`, `tilt`
- `SecondOrderParams` adds `T`, `dx`, `dy`
- `CavityParams` adds `v`, `freq`, `phi`
- `KickParams` carries strengths and offsets for kick-style tracking

This part of the design is worth keeping. It is one of the reasons the current system can support multiple tracking algorithms per element family.

## Layer 4: Transformation

Typical examples:

- `transformations/transfer_map.py`
- `transformations/second_order.py`
- `transformations/cavity.py`
- `transformations/kick.py`
- `transformations/runge_kutta.py`

The transformation:

- binds to atom hooks in `from_element(...)`
- requests params for a given energy via `get_params(energy)`
- applies the actual map in `map_function(...)`

Examples:

- `TransferMap` applies `R` and `B`
- `SecondTM` applies `R`, `B`, and `T`
- `CavityTM` uses cavity-specific metadata in addition to the linear map
- `KickTM` computes the kick algorithmically instead of consuming a ready-made matrix

## What Is Good About The Current Design

Before discussing redesign, it is important to preserve the advantages that already exist.

### 1. Physics and algorithm are separated

`CavityAtom` computes cavity-specific matrices and RF parameters. `CavityTM` decides how to apply those parameters to particles. That separation is useful and should remain.

### 2. One element family can support multiple tracking methods

A magnet family can expose first-order, second-order, kick-based, or Runge-Kutta tracking without forcing every algorithm into one monolithic class.

### 3. First-order optics stay available

`OpticElement` always keeps first-order maps available even when the active method is different. This is essential for Twiss, partial tracing, and other linear optics routines.

### 4. Edge handling is uniform

The `ENTRANCE -> MAIN -> EXIT` convention is shared across bends, cavities, and other edge-aware elements. That is a valuable invariant.

## Where The Current Design Is Too Implicit

This is where most developer confusion comes from.

### 1. The wrapper is more than a facade

The wrapper owns important framework behavior:

- lazy map construction
- invalidation of cached maps
- tracking-method switching
- slice construction in `get_section_tms(...)`
- public attribute forwarding

So a redesign cannot simply "remove the wrapper" unless that behavior is moved somewhere explicit and every caller still gets the same semantics.

### 2. Unsupported transformation behavior is inconsistent

Current behavior depends on the element:

- generic `OpticElement.set_tm(...)` relies on transformation construction and fallback
- `Cavity` forces `CavityTM`
- `Multipole` forces `MultipoleTM`
- warnings and fallbacks are not standardized

This is one of the first areas that should be cleaned up before deeper redesign.

### 3. Other subsystems still depend on `wrapper.element`

This is a critical compatibility point.

Several parts of CPBD and related adaptors inspect `element.element` directly, for example:

- lattice serialization in `latticeIO.py`
- some adaptor code
- code paths that assume an `OpticElement` instance rather than a raw physics object

That means a direct one-class element design is not only a local change inside `elements/`. It affects serialization, adaptors, and compatibility surfaces.

### 4. Hook-based capability discovery is hard to read

Today capability is often inferred from whether a method exists:

- if a transformation wants `create_first_order_main_params`, the atom must provide it
- if the element has edges, entrance and exit hooks must also exist
- failure often appears as `AttributeError` or an ad hoc wrapper restriction

This works, but it is hidden and fragile.

### 5. The magic forwarding hides where data lives

`OpticElement.__getattr__` and `__setattr__` are useful for backward compatibility, but they also hide the true storage location of physics state and the exact moment cache invalidation happens.

## Cavity As A Reference Example

`Cavity` is a good example because it exercises several important parts of the architecture.

### Current split

- `Cavity` is the public wrapper
- `CavityAtom` owns RF parameters and coupler-kick coefficients
- `CavityParams` is the contract object
- `CavityTM` owns the longitudinal tracking algorithm

This split is reasonable:

- the wrapper preserves the old public API
- the atom computes `R`, `B`, and cavity-specific metadata
- the transformation implements the nonlinear cavity update

### What `cavity_new.py` shows

`CavityNew` and `DirectOpticElement` are useful as experiments, but they should currently be treated as prototypes, not as a migration target.

They show two useful ideas:

- explicit `default_tm` and `supported_tms`
- a more direct public element object

But they also show two current redesign risks:

1. they duplicate a large part of `OpticElement` instead of extracting shared framework behavior
2. they do not solve the broader compatibility surface that still expects `OpticElement` and `wrapper.element`

So the prototype is valuable as a design probe, but not yet as an architectural answer.

## Important Invariants To Preserve

Any future redesign should preserve the following behavior.

### User-facing invariants

- public classes such as `Bend`, `Cavity`, `Drift`, `Quadrupole`, and `SBend` stay stable
- constructor signatures stay compatible unless a migration plan is provided
- `MagneticLattice(method=...)` continues to work
- numerical optics and tracking results remain unchanged within accepted tolerance

### Framework invariants

- first-order maps remain available even when active tracking uses another method
- edge elements still behave as `ENTRANCE -> MAIN -> EXIT`
- sliced tracking preserves current semantics
- `create_delta_e(...)` scales correctly for sliced maps
- caches are invalidated when tracking-relevant physics parameters change
- the serialization and adaptor surfaces that currently expect `wrapper.element` remain supported until deliberately migrated

## Recommended Transition Plan

The safest transition is staged.

### Phase 0: Freeze the current contract

Do this first:

- tighten this documentation
- add architecture-level unit tests
- list the current per-family behavior for supported transformations

Goal:

- remove guesswork about what the existing system guarantees

### Phase 1: Make capabilities explicit

Each element family should declare, in a predictable place:

- `default_tm`
- `supported_tms`
- what its `TransferMap` / `first_order_tms` optics path means
- `has_edge`

This can be introduced without redesigning the whole hierarchy.

Current policy in the main architecture:

- `default_tm` is the family fallback for active tracking
- `supported_tms`, when declared, should list wrapper-selectable active tracking TMs
- the `TransferMap` optics path used by `first_order_tms` is a separate concept and may exist even when `TransferMap` is not an allowed active TM for the wrapper
- declared TMs must actually build from the atom hooks; otherwise that is treated as a bug in the declaration
- generic wrappers still keep the legacy hook-based path for undeclared requests, with a warning
- pinned wrappers such as `Cavity`, `TWCavity`, `Multipole`, and `XYQuadrupole` intentionally normalize undeclared requests back to `default_tm`

Goal:

- make supported behavior visible without relying on missing-method exceptions
- later remove constructor duplication by letting `tm=None` resolve to the class `default_tm` instead of repeating the same TM in both the `__init__` signature and the class metadata

### Phase 2: Standardize unsupported-TM policy

Current rule:

- declared support that cannot be built raises clearly
- undeclared requests on generic wrappers warn and use the existing hook-based fallback behavior
- undeclared requests on pinned wrappers warn and fall back to `default_tm`

This keeps user-facing behavior conservative while still making bad
declarations fail fast in tests and development.

Goal:

- make method selection readable and testable

### Phase 3: Extract the framework contract before merging classes

Do not start by replacing wrapper-plus-atom everywhere.

Instead, identify the exact framework services currently provided by `OpticElement`:

- forwarding
- cache invalidation
- slicing
- TM creation
- compatibility with current public APIs

Then extract that contract into a reusable base or protocol.

Goal:

- avoid copying `OpticElement` logic into a second architecture

### Phase 4: Pilot only on simple elements

If an experimental direct architecture is still desired, try it first on:

- `Marker`
- `Monitor`
- `Drift`
- `Aperture`

Do not start with `Cavity`, `Bend`, or `Undulator`.

Goal:

- validate the new shape on low-risk elements

### Phase 5: Compare architectures on evidence

Before project-wide migration, compare:

- readability
- amount of duplicated code
- test coverage
- compatibility cost
- numerical behavior

Goal:

- make the migration decision from evidence, not taste

## What To Do If You Add A New Element Today

Until the architecture is redesigned, the safest extension path is still the current one:

1. start from the closest existing family
2. keep the wrapper-plus-atom split
3. reuse existing `TMParams` and transformations when the algorithm is unchanged
4. add new `TMParams` or a new transformation only when the tracking algorithm really changes
5. support first-order hooks even if the active method will usually be something else
6. add edge hooks only when `has_edge=True`
7. add tests before trying to optimize the structure

For new work today, that is lower risk than creating new public direct-style elements.

## Unit Tests Needed Before Redesign

The most useful tests are not only physics regression tests. We also need framework-contract tests.

The best approach is to add one new architecture-focused test module and keep a few element-family tests around it.

Current contract test area:

- `unit_tests/cpbd/architecture_contract/test_wrapper_atom_contract.py`
- `unit_tests/cpbd/architecture_contract/test_tm_selection_contract.py`
- `unit_tests/cpbd/architecture_contract/test_declared_capabilities_contract.py`
- `unit_tests/cpbd/architecture_contract/test_edge_slicing_contract.py`
- `unit_tests/cpbd/architecture_contract/test_energy_gain_contract.py`
- `unit_tests/cpbd/architecture_contract/test_matrix_contract.py`
- `unit_tests/cpbd/architecture_contract/test_undulator_contract.py`
- `unit_tests/cpbd/architecture_contract/test_lattice_compatibility_contract.py`

### A. Core wrapper and cache contract

These tests protect the current public behavior and should be added first.

1. Wrapper reads and writes are forwarded to the atom.
   Example: changing `quad.k1` updates `quad.element.k1`.

2. Changing a tracking-relevant public attribute invalidates both cached active maps and cached first-order maps.
   This is essential for safe refactoring of `__setattr__`.

3. `first_order_tms` always exists even when the active method is `SecondTM`, `KickTM`, or `CavityTM`.

4. `R(...)`, `B(...)`, and `T(...)` use the intended TM list for the current method.

5. `set_tm(...)` behavior is explicit and tested for both supported and unsupported methods.
   This should include current special cases such as `Cavity` and `Multipole`.

### B. Slice and edge contract

These tests are high value because slicing is easy to break during redesign.

1. No-edge element:
   `get_section_tms(...)` returns only one `MAIN` map.

2. Edge element full slice:
   `get_section_tms(start_l=0, delta_l=l)` returns `ENTRANCE`, `MAIN`, `EXIT`.

3. Edge element start slice:
   `get_section_tms(start_l=0, delta_l=l/2)` returns `ENTRANCE`, `MAIN`.

4. Edge element end slice:
   `get_section_tms(start_l=l/2, delta_l=l/2)` returns `MAIN`, `EXIT`.

5. Edge maps are copied, not rescaled, while the main map is rebuilt for the requested slice.

6. `first_order_only=True` is covered explicitly.
   This should be tested carefully because this logic is subtle and easy to make inconsistent.

### C. Energy-gain contract

These tests protect redesigns around cavities and other active elements.

1. Only the `MAIN` map contributes `delta_e`.

2. `create_delta_e(total_length, delta_length)` scales correctly for slices.

3. Full-element `delta_e` equals the sum of slice `delta_e` values within tolerance.

### D. Representative family tests

We do not need one architecture test per element immediately. We do need one representative test per important family.

Recommended coverage:

- `Drift`: simplest no-edge element
- `Quadrupole`: standard magnetic element with wrapper convenience properties
- `Bend` or `SBend`: edge-aware magnetic element
- `Cavity`: edge-aware RF element with custom transformation and energy gain
- `Multipole`: element that restricts the active transformation strongly
- `Undulator` or another nontrivial custom-tracking family when redesign work reaches that area

### E. Cavity-specific contract tests

For future architecture work, `Cavity` deserves dedicated tests beyond `remove_coupler_kick()`.

Recommended additions:

1. Constructor always normalizes unsupported TM choices to `CavityTM`.

2. `len(cavity.tms) == 3` and TM types are `ENTRANCE`, `MAIN`, `EXIT`.

3. Entrance and exit maps do not change reference energy; main map does.

4. `remove_coupler_kick()` changes the entrance and exit maps as expected after cache invalidation.

5. A sliced cavity preserves the same edge semantics as the full cavity.

6. Periodic Twiss and lattice-level behavior remain unchanged for the existing cavity regression cases.

### F. Compatibility tests outside `elements/`

Before migrating away from the wrapper-plus-atom model, add at least a few tests that protect the current external expectations.

Recommended coverage:

1. `MagneticLattice(method=...)` still applies family-specific TM selection correctly.

2. `latticeIO` still serializes representative elements correctly.

3. Any code path that depends on `wrapper.element` continues to work until explicitly migrated.

## If `CavityNew` Stays In The Tree

If the prototype remains in the repository, it should also be tested as a prototype rather than left as an undocumented alternative.

Useful prototype tests:

1. `CavityNew.R(...)`, `B(...)`, and `create_delta_e(...)` match `Cavity` for the same parameters.

2. `CavityNew.get_section_tms(...)` matches the current slice semantics expected from `Cavity`.

3. Any known differences from `Cavity` are documented intentionally rather than discovered accidentally later.

If those tests are not desired, the prototype should stay clearly labeled as experimental and non-public.

## Problems In The Previous Version Of This Note

The earlier README was directionally useful, but it had several weaknesses.

### 1. It repeated itself

There were multiple sections saying similar things about redesign, explicit capabilities, and templates. That made the document longer without making the contract clearer.

### 2. It described the split more cleanly than the code actually behaves

The old note described wrappers mainly as facades, while current code shows that `OpticElement` still owns important framework behavior.

### 3. It moved toward a one-class redesign too quickly

The earlier version mentioned the direct-style prototype as a reasonable next phase before fully freezing the current contract. That is premature given the existing compatibility surface.

### 4. It did not emphasize external compatibility enough

The previous note did not clearly say that other modules still depend on `OpticElement` and `wrapper.element`.

### 5. It mixed extension guidance with redesign advocacy

Those are related topics, but they should not be confused. The current extension path and the future redesign path are not the same thing.

## Recommendation

The right direction is evolutionary, not revolutionary.

The project should move toward:

- explicit capabilities
- consistent TM selection behavior
- stronger architecture tests
- less hidden wrapper magic
- eventually, possibly, a simpler internal model

But the immediate next step should be test-first hardening of the current contract, not replacing the architecture by intuition.
