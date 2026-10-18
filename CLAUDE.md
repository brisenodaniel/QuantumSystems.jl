# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

Julia package (`julia = "1.11"`). From the repo root:

- Run all tests: `julia --project -e 'using Pkg; Pkg.test()'`
- Run a single test file: `julia --project -e 'using Test, LinearAlgebra, QuantumToolbox, QuantumSystems; include("test/helpers.jl"); include("test/test_SimpleQuantumSystem.jl")'` (test files rely on `test/helpers.jl`, which defines the `fluxonium` builder)
- Install deps: `julia --project -e 'using Pkg; Pkg.instantiate()'`

`QuantumToolbox` is pulled from a fork (`brisenodaniel/QuantumToolbox.jl`, branch `main`) via `[sources]` in `Project.toml`, not the registered release; the code depends on fork-specific behavior (e.g. `qeye_like`, `FloquetBasis` support).

## Architecture

Built on QuantumToolbox's `QuantumObject`s. `src/QuantumSystems.jl` holds the module header (abstract type, aliases, exports) and `include`s, in order:

- `abstract_quantum_systems.jl`: code shared by both system types (validation helpers, state-label parsing, `getproperty`, `tensor_dims`, memoized `eigenstates`).
- `simple_quantum_systems.jl` / `composite_quantum_systems.jl`: the struct plus everything specific to it (`getindex`, `state`, `to_energy_basis`, `tidyup`, ...).
- `floquet.jl`: Floquet `state`/`energy` methods taking a `FloquetBasis`, and `fold_quasienergy`.
- `infix_operators.jl`: `⊗` and `+`/`-` on systems.

- `SimpleQuantumSystem`: a `NamedTuple` of operators (`ops`, must include `H`) plus optional `params` and `name`.
- `CompositeQuantumSystem`: a `NamedTuple` of named `subsys` (each must have a non-`nothing` `name`) plus operators on the tensor-product space. Optional `U_basis`/`dressed_dim` hold an energy-basis transform/truncation. `state_lookup` is a mutable array (bare index → dressed eigenstate index) filled lazily and used by `state`/`state_idx`.
- `Base.getproperty` is overloaded so `sys.H`, `sys.φ`, etc. resolve to entries in `ops` (and, for composites, to subsystems by name); `sys.I` gives the identity. `sys[:key]` indexes `ops`. Real struct fields take priority over operator names.
- `⊗` (from QuantumToolbox) composes Simple/Composite systems in all four combinations; `sys + Op` / `sys + NamedTuple` add operators. These lift local operators into the full tensor space via the `_lift_*` helpers.
- State labelling: `StateIndex` is an `Int`, a vector of `Int`s (one per subsystem for composites), or a symbol from `g,e,f,h` (`_STATE_SYMBOLS`). Indexing is 0-based in the user API and converted to Julia 1-based internally.
- `state`/`state_idx` take `bare_to_dressed_tracking` (`:overlap`, default, or `:spectrum`) to map bare labels to dressed eigenstates in composites. A `FloquetBasis` method of `state` returns Floquet states with quasienergies folded via `fold_quasienergy`.
- `eigenstates(sys)` is memoized with an LRU cache (size 100), keyed on the system object. Be mindful of this when mutating or reusing systems.
- `to_energy_basis`, `truncate_sys`, and `tidyup` (QuantumToolbox extension) transform/truncate systems; validation is done by `_check_ops`, `_check_subsys`, and `_check_idx`, which return exceptions that constructors throw.

Tests (`test/`) cover `SimpleQuantumSystem` and `CompositeQuantumSystem` using a test-local fluxonium builder, so there is no dependency on any circuit library.
