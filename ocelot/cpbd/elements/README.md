# CPBD Element Architecture: Current Contract and Improvement Plan

This note is for developers working inside `ocelot.cpbd`.

It has two goals:

1. explain how the current element machinery really works
2. describe an improvement path that makes the current architecture clearer without breaking user-visible behavior

The public Ocelot docs already explain how to use elements such as `Bend`, `SBend`, `RBend`, `Cavity`, and `Quadrupole`. What is harder to understand from the outside is the internal contract between wrappers, atoms, parameter containers, and transformations. That contract is the important part to preserve before any architecture change.

In this note:

- `wrapper` means the public element object that users and lattices usually hold, typically an `OpticElement` subclass such as `Quadrupole` or `Cavity`
- `atom` means the internal physics object stored on the wrapper, usually as `wrapper.element`
- `parameter container` or `TMParams` means the structured data object built by the atom and passed to a transformation, for example `FirstOrderParams`, `SecondOrderParams`, or `CavityParams`
- `transformation` means the tracking algorithm object built from the element hooks and `TMParams`, for example `TransferMap`, `SecondTM`, `KickTM`, or `CavityTM`

## Bottom Line

The current architecture already has several strong sides that are worth
keeping:

- element physics lives on atoms, while tracking algorithms live in transformations
- one element family can support more than one active tracking method
- first-order optics can stay available even when active tracking uses another TM
- the public wrapper preserves a stable user-facing API across lattice and tracking code

The main weakness is not the overall split. The main weakness is that too much
of the contract is implicit.

In the current code, a developer still has to discover important rules by
reading implementation details:

- `OpticElement` is part of the real framework behavior, not just a thin facade
- wrapper forwarding hides where physics state actually lives
- TM support and edge behavior have historically been inferred from hook presence, wrapper restrictions, and fallback rules
- some family behavior is clear only after reading both the wrapper and the atom

So the current goal is not to redesign the architecture first. The current goal
is to make the existing architecture clearer, safer, and easier to work with.

That improvement work already started in three ways:

1. the current contract is being documented explicitly in this note
2. architecture-level unit tests now freeze the important behavior
3. wrapper TM declarations and edge behavior are being made more explicit

The next steps should continue in the same direction:

1. finish making per-family capabilities visible
2. keep tightening contract tests around edge and TM behavior
3. simplify local parts of the current design only where that improves clarity without changing physics

Users should not notice these changes except through clearer developer-facing
behavior, safer maintenance, and preserved physics results within tolerance.

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

Every public element wrapper also has a family default active tracking method
(`default_tm`).

- for many families this default is `TransferMap`, meaning the first-order linear map
- some families use a different default, for example `CavityTM`, `TWCavityTM`, or `MultipoleTM`

If the user does not request anything else, the wrapper starts with that family
default.

There are three common ways to request a different active TM:

1. at construction time:

```python
from ocelot.cpbd.elements import Quadrupole, Octupole
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.kick import KickTM

quad = Quadrupole(l=0.4, k1=1.2, tm=SecondTM)
octu = Octupole(l=0.2, k3=3.0, tm=KickTM)
```

2. on an existing element:

```python
from ocelot.cpbd.transformations.second_order import SecondTM

quad.set_tm(SecondTM)
```

3. for a whole lattice or selected element families through `MagneticLattice`:

```python
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.kick import KickTM
from ocelot.cpbd.elements import Octupole

lat = MagneticLattice(cell, method={"global": SecondTM, Octupole: KickTM, "nkick": 5})
```

If `method` is changed after lattice construction, call:

```python
lat.update_transfer_maps()
```

Important:

- not every family allows every TM to stay active
- families such as `Cavity`, `TWCavity`, `Multipole`, and `XYQuadrupole` expose only one active TM by declaring `supported_tms = {default_tm}`
- those single-method families still keep their `TransferMap` optics path in `first_order_tms`, but explicit requests for another active TM raise
- only broad global lattice requests such as `method={"global": SecondTM}` are allowed to warn and fall back to `default_tm`

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

## Edge Contract

`has_edge` is not cosmetic metadata. It changes how CPBD builds both the
always-available optics path and the active tracking path.

- `has_edge = False` means the element builds only one `MAIN` map
- `has_edge = True` means the element builds a three-map sequence:
  `ENTRANCE -> MAIN -> EXIT`
- this rule applies to `first_order_tms` as well as to the active `tms`
- therefore an edge-aware family must always provide first-order entrance and
  exit hooks, even if its active tracking TM is a custom family-specific one

For slices built with `get_section_tms(...)`, the current contract is:

- `start_l == 0` includes `ENTRANCE`
- `start_l + delta_l == l` includes `EXIT`
- middle slices rebuild only the `MAIN` map for the requested `delta_l`
- copied edge maps are not rescaled
- `ignore_edges=True` suppresses entrance and exit maps explicitly

If a family declares support for a TM and also has `has_edge=True`, then that
TM must be able to build entrance, main, and exit maps from the atom hooks.
If it cannot, that is a bug in the wrapper contract, not a user mistake.

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
- for audited multi-method wrappers, `supported_tms` should match the TM families that the atom can really build for active tracking
- a family that intentionally exposes only one active TM should use `supported_tms = {default_tm}`, even if its atom can build more internal paths for optics or legacy helper code

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

| Family | Has edge | Linear optics / Twiss path | Declared / selectable active TMs | Extra internal or legacy-buildable paths | Notes |
| --- | --- | --- | --- | --- | --- |
| `Aperture` | No | `TransferMap` via generic `Element` first-order fallback | `TransferMap`, `SecondTM` | none | No dedicated aperture-specific transformation; mainly a holder for aperture metadata and offsets. |
| `Marker` | No | `TransferMap` via generic `Element` first-order fallback | `TransferMap`, `SecondTM` | none | Zero-length reference point. |
| `Monitor` | No | `TransferMap` via generic `Element` first-order fallback | `TransferMap`, `SecondTM` | none | Diagnostic state such as `x`, `y`, `x_ref`, `y_ref` lives on the atom. |
| `Matrix` | No | `TransferMap` from stored `R` / `B` data | `TransferMap`, `SecondTM` | none | Sliced sections currently fall back to a drift-like first-order map rather than slicing the stored matrix. |
| `UnknownElement` | No | `TransferMap` via inherited `Magnet` first-order hook | `TransferMap`, `SecondTM`, `KickTM` | none | Legacy placeholder family, but its generic `Magnet` hook surface is now declared explicitly. |
| `Pulse` | n/a | n/a | n/a | n/a | Separate helper object for time-dependent kicks; not part of the `OpticElement` wrapper/atom/TM stack. |

### Magnet-like families

| Family | Has edge | Linear optics / Twiss path | Declared / selectable active TMs | Extra internal or legacy-buildable paths | Notes |
| --- | --- | --- | --- | --- | --- |
| `Drift` | No | `TransferMap` via `DriftAtom` first-order hook | `TransferMap`, `SecondTM`, `KickTM`, `RungeKuttaTM`, `RungeKuttaTrTM` | none | Useful no-edge reference family. |
| `Quadrupole` | No | `TransferMap` via inherited `Magnet` first-order hook | `TransferMap`, `SecondTM`, `KickTM` | none | Wrapper also exposes convenience properties `k1l` and `k2l`. |
| `Sextupole` | No | `TransferMap` via inherited `Magnet` first-order hook | `TransferMap`, `SecondTM`, `KickTM` | none | `KickTM` is the important nonlinear tracking path here. |
| `Octupole` | No | `TransferMap` via inherited `Magnet` first-order hook | `TransferMap`, `SecondTM`, `KickTM` | none | `k3` only enters the dedicated kick path; `SecondTM` still comes from inherited generic `Magnet` second-order hooks. |
| `Solenoid` | No | `TransferMap` via solenoid-specific first-order hook | `TransferMap`, `SecondTM` | none | `SecondTM` currently comes from the generic `Element` second-order fallback, not from a solenoid-specific nonlinear model. |
| `Hcor` / `Vcor` | No | `TransferMap` via `CorAtom` first-order hook | `TransferMap`, `SecondTM` | none | These are not `Magnet` subclasses, so they do not inherit generic kick hooks. |
| `XYQuadrupole` | No | `TransferMap` via `XYQuadrupoleAtom` first-order hook | `TransferMap` only | the atom also has `SecondTM` and inherited `KickTM` hook surface, but the wrapper intentionally exposes only first-order active tracking | Important wrapper-level exception. |
| `Bend`, `SBend`, `RBend` | Yes | `TransferMap` via bend entrance/main/exit first-order hooks | `TransferMap`, `SecondTM`, `KickTM`, `RungeKuttaTM`, `RungeKuttaTrTM` | none | `KickTM` is declared here because the bend atoms do provide the full inherited `Magnet` kick hook family, including edges. |

### RF, field-integrated, and special families

| Family | Has edge | Linear optics / Twiss path | Declared / selectable active TMs | Extra internal or legacy-buildable paths | Notes |
| --- | --- | --- | --- | --- | --- |
| `Cavity` | Yes | `TransferMap` via cavity-specific first-order entrance/main/exit hooks | `CavityTM` only | first-order `TransferMap` exists for optics but is not exposed as an active tracking method | Most important complex reference family for edge handling and `delta_e`. |
| `TWCavity` | Yes | `TransferMap` via traveling-wave cavity first-order hooks | `TWCavityTM` only | first-order `TransferMap` exists for optics but is not exposed as an active tracking method | The atom currently warns that this family is unfinished. |
| `TDCavity` | No | `TransferMap` via `TDCavityAtom` first-order hook | `TransferMap`, `SecondTM` | none | The atom's `additional_tms = [SecondTM]` is informational only; `OpticElement` does not consume it. |
| `Undulator` | No | `TransferMap` via `UndulatorAtom` first-order hook | `TransferMap`, `SecondTM`, `RungeKuttaTM`, `RungeKuttaTrTM`, `UndulatorTestTM` | none | `MagneticLattice.update_transfer_maps()` also has special-case length handling when a field map is attached. |
| `Multipole` | No | `TransferMap` via `MultipoleAtom` first-order hook | `MultipoleTM` only | first-order `TransferMap` exists for optics but the wrapper intentionally exposes only `MultipoleTM` as an active tracking method | This first-order path is a linearized multipole optics map, not a pure drift. |

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

### Public attribute access rule

Code outside the element internals should use the public element API, not raw
storage.

- read parameters with `elem.k1` or `getattr(elem, "k1")`
- write parameters with `elem.k1 = value` or `setattr(elem, "k1", value)`
- do not use `elem.__dict__["k1"]` for wrappers

Why this matters:

- physics state may live on the atom, not in the wrapper `__dict__`
- wrapper forwarding makes `elem.k1` work even when storage lives on `elem.element`
- `setattr(...)` goes through wrapper cache invalidation, while direct `__dict__` writes bypass it

For generic serializer or adaptor code, the safe pattern is:

1. discover public constructor or exported fields
2. read them with `getattr(...)`
3. write them with `setattr(...)`

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

Here `*` is a placeholder for the map location:

- `main` for the body of the element
- `entrance` for the entrance edge
- `exit` for the exit edge

So, for example, `create_first_order_*_params(...)` means methods such as:

- `create_first_order_main_params(...)`
- `create_first_order_entrance_params(...)`
- `create_first_order_exit_params(...)`

Elements without edges usually implement only the `main` variant.

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

Before discussing further cleanup, it is important to preserve the advantages that already exist.

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

So cleanup cannot simply "remove the wrapper" unless that behavior is moved somewhere explicit and every caller still gets the same semantics.

### 2. TM selection policy must stay explicit

This part of the architecture used to be inconsistent. The current direction is
to keep it explicit and centralized in `OpticElement`.

Current rule:

- if `Quadrupole.supported_tms = {TransferMap, SecondTM, KickTM}` and the user requests `SecondTM`, the wrapper treats that as an explicitly supported active TM
- if the user requests an undeclared TM such as `RungeKuttaTM` directly through `quad.set_tm(RungeKuttaTM)` or `Quadrupole(..., tm=RungeKuttaTM)`, the wrapper raises because that is an explicit unsupported request
- if a lattice applies `method={"global": RungeKuttaTM}`, the same undeclared request is treated as permissive: the wrapper warns and falls back to `default_tm`
- if `Cavity.supported_tms = {CavityTM}` and the user requests `TransferMap` directly, the wrapper raises because `Cavity` exposes only one active tracking TM
- if a lattice applies `method={"global": SecondTM}` to a sequence containing a `Cavity`, the wrapper warns and falls back to `CavityTM`

Internally, `OpticElement.set_tm(...)` distinguishes these cases with
`request_source`:

- `Quadrupole(..., tm=...)` uses `request_source="explicit"`
- `quad.set_tm(...)` uses `request_source="explicit"`
- `MagneticLattice(method={Quadrupole: ...})` uses `request_source="explicit"`
- `MagneticLattice(method={"global": ...})` uses `request_source="global"`

That distinction is the reason the wrapper can stay strict for direct
family-specific requests while remaining permissive for broad global lattice
requests such as `method={"global": SecondTM}`.

Current outcomes are:

- declared request:
  keep the requested TM active
- explicit undeclared request:
  raise an error
- global undeclared request:
  warn and fall back to `default_tm`
- declared support that cannot actually be built from atom hooks:
  raise an error, because that is a broken wrapper declaration

So the intended user-facing rule is:

- specific requests should be strict
- broad `method={"global": ...}` requests should be permissive
- families with `supported_tms = {default_tm}` should keep that single active TM even under global reassignment

This is intentionally different from the old broad hook-based fallback path.
The remaining cleanup should continue by auditing families and either:

- declaring support explicitly in `supported_tms`
- or, for families with only one allowed active TM, keeping `supported_tms = {default_tm}`
- or, for rare unaudited families, leaving only the family `default_tm` as the safe path

### 3. The wrapper/atom split is still externally relevant

This is still a critical compatibility point.

Historically, parts of CPBD and related adaptors inspected `element.element` directly. The known serializer/adaptor cases have been migrated to the public wrapper API, so this is no longer the main blocker. The split is still not an implementation detail because:

- `OpticElement` itself stores physics on `self.element`
- wrapper forwarding and cache invalidation depend on that split
- some code paths still assume an `OpticElement` instance rather than a raw physics object

That means a direct one-class element design is not only a local change inside `elements/`. It still affects compatibility surfaces even after those serializer/adaptor cleanups.

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

## Important Invariants To Preserve

Any future cleanup should preserve the following behavior.

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
- the public wrapper API used by serialization and adaptor code remains stable

## Recommended Improvement Plan

The safest improvement path is staged.

### Phase 0: Freeze the current contract

Do this first:

- tighten this documentation
- add architecture-level unit tests
- list the current per-family behavior for supported transformations

Goal:

- remove guesswork about what the existing system guarantees

### Phase 1: Finish making capabilities explicit

This is now part of the main architecture, not just a future idea.

Current state:

- representative public wrappers already declare `default_tm`
- representative public wrappers already declare `supported_tms`
- architecture contract tests already freeze that behavior

What still needs to be completed:

- what its `TransferMap` / `first_order_tms` optics path means
- `has_edge`
- which edge hooks exist and how `ENTRANCE -> MAIN -> EXIT` is formed for that family

Current policy in the main architecture:

- `default_tm` is the family fallback for active tracking
- `supported_tms`, when declared, should list wrapper-selectable active tracking TMs
- wrappers treat constructor `tm=...`, direct `set_tm(...)`, and `MagneticLattice(method={Family: ...})` as explicit requests
- wrappers treat `MagneticLattice(method={"global": ...})` as a global request
- the `TransferMap` optics path used by `first_order_tms` is a separate concept and may exist even when `TransferMap` is not an allowed active TM for the wrapper
- declared TMs must actually build from the atom hooks; otherwise that is treated as a bug in the declaration
- explicit undeclared requests raise
- global undeclared requests warn and fall back to `default_tm`
- families with `supported_tms = {default_tm}` therefore stay on that single active TM unless they are asked for it explicitly

Goal:

- make supported behavior visible without relying on missing-method exceptions
- later remove constructor duplication by letting `tm=None` resolve to the class `default_tm` instead of repeating the same TM in both the `__init__` signature and the class metadata

### Phase 2: Make edge behavior more explicit

The next high-value clarification is edge behavior.

Current situation:

- `has_edge` already controls whether `ENTRANCE` and `EXIT` maps are built
- that behavior is important, but it is still easy to miss by just reading wrappers
- slicing rules in `get_section_tms(...)` depend on it directly

What to improve:

- make `has_edge` visually obvious in the family inventory
- document which families really implement entrance and exit hooks
- keep dedicated tests for `ENTRANCE -> MAIN -> EXIT` semantics
- avoid hidden assumptions that every magnetic family behaves like a bend

Goal:

- make edge-aware behavior readable without digging through hook names

### Phase 3: Keep TM-selection policy centralized

Current rule:

- declared support that cannot be built raises clearly
- explicit undeclared requests raise clearly
- global lattice requests warn and fall back to `default_tm`
- families that expose only `supported_tms = {default_tm}` therefore keep that TM active under global reassignment

This keeps user-facing behavior conservative while still making bad
declarations fail fast in tests and development.

There is still real code that relies on broad lattice-level TM requests such as
`MagneticLattice(method={"global": SecondTM})`. That is why the permissive
fallback now lives only on the explicit global lattice path, not on direct
family-specific or user-specific requests.

Goal:

- keep method selection readable, centralized, and testable

### Phase 4: Review `TMParams` containers carefully

`TMParams` already does something useful: it gives an explicit boundary between
atom physics code and transformation code.

Potential benefits of dataclasses:

- clearer field lists
- less boilerplate in simple containers
- easier debugging and repr output

Current risks / limits:

- some params classes already have behavior, not just fields
- inheritance is used (`SecondOrderParams` extends `FirstOrderParams`)
- refactoring them now gives less value than clarifying TM selection and edge behavior

Recommended approach:

- do not start with a broad dataclass conversion
- if this is revisited later, start with the simplest leaf containers first
- keep the current constructor contract stable until tests clearly freeze it

Goal:

- improve readability of params objects without creating churn in a stable boundary

### Phase 5: Extract the framework contract only if needed

Instead, identify the exact framework services currently provided by `OpticElement`:

- forwarding
- cache invalidation
- slicing
- TM creation
- compatibility with current public APIs

Then extract that contract into a reusable base or protocol.

Goal:

- make the current architecture easier to understand before considering deeper refactoring

## What To Do If You Add A New Element Today

Until the architecture is improved further, the safest extension path is still the current one:

1. start from the closest existing family
2. keep the wrapper-plus-atom split
3. reuse existing `TMParams` and transformations when the algorithm is unchanged
4. add new `TMParams` or a new transformation only when the tracking algorithm really changes
5. support first-order hooks even if the active method will usually be something else
6. add edge hooks only when `has_edge=True`
7. add tests before trying to optimize the structure

For new work today, that is lower risk than creating new public direct-style elements.

## Unit Tests That Protect Future Cleanup

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

These tests are high value because slicing is easy to break during cleanup.

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

These tests protect cleanup around cavities and other active elements.

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
- `Undulator` or another nontrivial custom-tracking family when cleanup reaches that area

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

3. Public wrapper based serialization / adaptor paths continue to work without reintroducing direct `wrapper.element` access.

## Problems In The Previous Version Of This Note

The earlier README was directionally useful, but it had several weaknesses.

### 1. It repeated itself

There were multiple sections saying similar things about architecture changes, explicit capabilities, and templates. That made the document longer without making the contract clearer.

### 2. It described the split more cleanly than the code actually behaves

The old note described wrappers mainly as facades, while current code shows that `OpticElement` still owns important framework behavior.

### 3. It moved toward architecture replacement too quickly

The earlier version talked too much about replacing the structure before finishing the basic clarity work on the current implementation.

### 4. It did not emphasize external compatibility enough

The previous note did not clearly separate two facts:

- external code still depends on public wrapper objects
- the known serializer/adaptor cases no longer need direct `wrapper.element` access

### 5. It mixed extension guidance with architecture replacement

Those are related topics, but they should not be confused. The current extension path and the future cleanup path are not the same thing.

## Recommendation

The right direction is evolutionary, not revolutionary.

The project should move toward:

- explicit capabilities
- consistent TM selection behavior
- stronger architecture tests
- less hidden wrapper magic
- maybe later, if clearly justified, a simpler internal model

But the immediate next step should be test-first hardening and clearer documentation of the current contract, not replacing the architecture by intuition.
