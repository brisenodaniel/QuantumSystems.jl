module QuantumSystems

using LinearAlgebra
using QuantumToolbox
using Memoize, LRUCache

abstract type AbstractQuantumSystem end
const Qsys = AbstractQuantumSystem
const QuantumSystem = AbstractQuantumSystem
const StateIndex = Union{Int, AbstractVector{<:Int}, Symbol}
const Maybe{T} = Union{T, Nothing} where T
const _STATE_SYMBOLS = (; g=1, e=2, f=3, h=4) # Julia indexing

export AbstractQuantumSystem, Qsys, QuantumSystem, StateIndex, Maybe
export SimpleQuantumSystem, CompositeQuantumSystem, tensor_dims
export state, energy, to_energy_basis, truncate_sys, I, state_idx, fold_quasienergy

include("abstract_quantum_systems.jl")
include("simple_quantum_systems.jl")
include("composite_quantum_systems.jl")
include("floquet.jl")
include("infix_operators.jl")

end # module QuantumSystems
