# QuantumSystems.jl API

QuantumSystems wraps named collections of [QuantumToolbox](https://github.com/qutip/QuantumToolbox.jl)
operators into *systems* that know their Hamiltonian, can be composed with `⊗`, and can label and
look up eigenstates.

```julia
using QuantumToolbox, QuantumSystems
```

## Types

### `AbstractQuantumSystem`

Supertype of all systems. Aliases: `Qsys`, `QuantumSystem`.

Other exported aliases:

- `StateIndex = Union{Int, AbstractVector{<:Int}, Symbol}`: a state label (see [State labels](#state-labels)).
- `Maybe{T} = Union{T, Nothing}`.

### `SimpleQuantumSystem(ops; params=nothing, name=nothing)`

A single (non-factorized) system.

- `ops::NamedTuple`: `QuantumObject` operators, all with the same dimensions, and containing the
  Hamiltonian under the key `:H`. Otherwise an `ArgumentError` is thrown.
- `params::Maybe{NamedTuple}`: arbitrary metadata, e.g. circuit parameters.
- `name::Maybe{Symbol}`: required if the system is used as a factor in `⊗`.

`Qsys(ops; params, name)` is equivalent. `Qsys(sys::SimpleQuantumSystem, name=nothing)` makes a copy
of `sys`, optionally renamed.

### `CompositeQuantumSystem(subsys, ops; params=nothing, U_basis=nothing, dressed_dim=nothing, name=nothing)`

A tensor product of named subsystems.

- `subsys::NamedTuple`: non-empty, with `AbstractQuantumSystem` values. Order defines the tensor
  ordering.
- `ops::NamedTuple`: operators on the full tensor-product space; must contain `:H`.
- `U_basis`: basis transform applied when lifting subsystem operators and bare states (set by
  `to_energy_basis`).
- `dressed_dim`: number of dressed states kept (set by `to_energy_basis`).

Normally built with `⊗` rather than directly.

## Accessing contents

| Expression | Returns |
|---|---|
| `sys.H`, `sys.ops`, `sys.params`, `sys.name` | Struct fields take priority, then operators in `ops` (`sys.H`, `sys.n`, ...) |
| `sys.I`, `I(sys)`, `qeye(sys)` | Identity on the Hilbert space of `sys.H` |
| `sys[:key]` | Operator `key` (composites: also the subsystem named `key`) |
| `sys[:sub, :op]` (composite) | Operator `op` of subsystem `sub`, lifted to the full space |
| `sys[:s1, :op1, :s2, :op2, ...]` (composite) | Product of lifted operators, one per subsystem (each subsystem at most once). If `U_basis` is set, returned as `U' * O * U` |
| `sys.sub` (composite) | The subsystem named `sub` |
| `tensor_dims(sys)` | Local dimensions of each subsystem (a vector for simple systems) |

## Composition and modification

### `A ⊗ B`

Tensor product of any combination of simple and composite systems. Simple systems must have a
`name`. The result is a `CompositeQuantumSystem` with:

- `subsys`: subsystems of `A` followed by those of `B` (name collisions throw `ArgumentError`),
- `H = H_A ⊗ I + I ⊗ H_B`,
- `params`: a `NamedTuple` of the factors' params (keyed by name for simple systems),
- `name = :A_B`.

### `sys + op`, `op + sys`, `sys - op`

`op` is a `QuantumObject` with the same dimensions as `sys.H`. Returns a new system with `op` added
to (or subtracted from) the Hamiltonian. Other operators, `params`, `name` (and `U_basis`,
`dressed_dim` for composites) are kept.

### `sys + ops`, `ops + sys`

`ops` is a `NamedTuple` of `QuantumObject`s with dimensions matching `sys.H`. Returns a new system
whose operators are `merge(sys.ops, ops)`: new keys are added and existing keys, including `:H`,
are replaced.

### `to_energy_basis(sys; set_energy_zero, etol)`

Returns `sys` rewritten in the eigenbasis of `sys.H`.

- `set_energy_zero`: `true` shifts the ground-state energy to zero, a number shifts by that value,
  `false` does nothing. Default is `false` for simple systems and `true` for composite systems.
- `etol`: `true` applies `tidyup` with the default tolerance, a `Float64` applies it with that
  tolerance, `false` (default) skips it.

For composites, only `:H` is kept in `ops`, and `U_basis` and `dressed_dim` are set.

### `truncate_sys(sys::SimpleQuantumSystem, N)`

Keeps the first `N` levels of every operator.

### `tidyup(sys, etol=nothing)`

Applies QuantumToolbox's `tidyup` to every operator (and, for composites, subsystems and
`U_basis`). `nothing` uses `QuantumToolbox.settings.tidyup_tol`.

## State labels

A `StateIndex` is one of:

- an `Int` `k`: the `k`-th eigenstate of `sys.H`, with `0` the ground state;
- a `Vector{Int}` with one 0-based level per subsystem (composite systems);
- a `Symbol` whose characters are in `g, e, f, h` (levels `0, 1, 2, 3`), one per subsystem, e.g.
  `:g`, `:ge`.

Labels that are out of range or have the wrong length throw an `ArgumentError`.

## Eigenstates

### `eigenstates(sys)`

Memoized (LRU, 100 entries) wrapper around `eigenstates(sys.H)`. A new system object, such as one
returned by `sys + op`, has an empty cache.

### `state(sys::SimpleQuantumSystem, lbl; bare_to_dressed_tracking=:overlap)`

Returns `(λ, ψ)`, the eigenenergy and eigenstate of `sys.H` for `lbl`. The tracking keyword is
validated but otherwise ignored.

### `state(sys::CompositeQuantumSystem, lbl, bare=false; bare_to_dressed_tracking=:overlap)`

Returns `(λ, ψ)`.

- `bare=true`: the product of the subsystem eigenstates for `lbl` (energy is the sum of the
  subsystem energies), mapped through `U_basis` if present. A single-integer label indexes the
  eigenstates of the bare Hamiltonian.
- `bare=false`: the eigenstate of `sys.H` that `lbl` is mapped to by `bare_to_dressed_tracking`:
  - `:overlap`: the eigenstate with maximal overlap with the bare product state. Results are cached
    in `sys.state_lookup`.
  - `:spectrum`: the eigenstate whose energy rank equals that of the bare product state. Ties among
    degenerate bare levels are broken by label order, so prefer `:overlap` there.

### `state_idx(sys::CompositeQuantumSystem, lbl, bare=false; bare_to_dressed_tracking=:overlap)`

The 1-based index into `eigenstates(sys.H)` for the state selected by `state`.

### `energy(sys, lbl; ...)`, `energy(sys::CompositeQuantumSystem, lbl, bare=false; ...)`

`state(...)[1]`; same arguments and keywords as `state`.

## Floquet states

### `state(sys, fb::FloquetBasis, lbl, t=0.0; bare_to_dressed_tracking=:overlap, min_overlap=0.5, cache=false, kwargs...)`

Returns `(ε, ψ)`: the folded quasienergy and Floquet state at time `t` that is adiabatically
connected to the dressed eigenstate of `sys` labelled by `lbl`.

- `fb` must be built with the tracked constructor `FloquetBasis(sys.H, H_dr, T)` so that its
  `k`-th Floquet state connects to the `k`-th eigenstate of `sys.H`. Mismatched dimensions throw
  `ArgumentError`.
- `min_overlap`: warns if the squared overlap between the static eigenstate and the Floquet state at
  `t=0` is below this value. Use `nothing` to skip the check.
- `cache=true` stores the micromotion propagator at `t` in `fb`.
- `kwargs` are passed to `propagator`.
- `ψ` is `U(t, 0) ψ(0)`, not the periodic Floquet mode, and its global phase is arbitrary.

### `energy(sys, fb::FloquetBasis, lbl; kwargs...)`

`state(sys, fb, lbl; kwargs...)[1]`.

### `fold_quasienergy(ε, T)`

Folds `ε` onto the interval ``[-\pi/T, \pi/T)``.
