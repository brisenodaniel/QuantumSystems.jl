struct SimpleQuantumSystem{
    Ops<:NamedTuple,
    PT<:NamedTuple} <: AbstractQuantumSystem
    ops::Ops
    params::Maybe{PT}
    name::Maybe{Symbol}

    function SimpleQuantumSystem(
        ops::Ops;
        params::Maybe{PT} = nothing,
        name::Maybe{Symbol} = nothing
        ) where {
        Ops<:NamedTuple,
        PT<:NamedTuple
        }
        err = _check_ops(ops)
        isnothing(err) ? nothing : throw(err)
        return new{Ops, PT}(ops, params, name)
    end
end

function Qsys(
    ops::Ops;
    params::Maybe{PT}=nothing,
    name::Maybe{Symbol}=nothing,
    ) where {
        Ops <: NamedTuple,
        PT <: NamedTuple
    }
    return SimpleQuantumSystem(ops; params=params, name=name)
end

function Qsys(
    sys::SimpleQuantumSystem,
    name::Maybe{Symbol}=nothing
    )
    new_ops = (;sys.ops...)
    new_params = isnothing(sys.params) ? nothing : (;sys.params...)
    name = isnothing(name) ? sys.name : name
    return SimpleQuantumSystem(new_ops; params=new_params, name=name)
end

function Base.getindex(sys::SimpleQuantumSystem, key::Symbol)
    return Core.getfield(sys, :ops)[key]
end

function to_energy_basis(
    sys::SimpleQuantumSystem;
    set_energy_zero::Union{Bool, Float64}=false,
    etol::Union{Bool, Float64}=false)
    λ, _, U = eigenstates(sys.H)
    U = Qobj(U, dims=sys.H.dims)
    e0 = set_energy_zero === true ? λ[1] :
        set_energy_zero isa Float64 ? set_energy_zero : 0.0
    new_ops = (;(k=> U'*op*U for (k, op) in pairs(sys.ops))...)
    if e0 != 0.0
        H = new_ops[:H] - e0
        new_ops = (; new_ops..., H=H)
    end
    sys = SimpleQuantumSystem(new_ops; params=sys.params, name=sys.name)
    if etol isa Float64
        return tidyup(sys, etol)
    elseif etol
        return tidyup(sys)
    else
        return sys
    end
end

function truncate_sys(sys::SimpleQuantumSystem, N::Int)
    new_ops = (;(k => truncate_Qobj(op, N) for (k, op) in pairs(sys.ops))...)
    return SimpleQuantumSystem(new_ops; params=sys.params, name=sys.name)
end

function QuantumToolbox.tidyup(sys::SimpleQuantumSystem, etol::Maybe{Float64}=nothing)
    etol = isnothing(etol) ? QuantumToolbox.settings.tidyup_tol : etol
    new_ops = (;(k=> tidyup(op, etol) for (k, op) in pairs(sys.ops))...)
    return SimpleQuantumSystem(new_ops; params=sys.params, name=sys.name)
end

function truncate_Qobj(q::AbstractQuantumObject, N::Int)
    return QuantumObject(q.data[1:N, 1:N], q.type, (N,))
end

@doc raw"""
    `state(sys::SimpleQuantumSystem, lbl::StateIndex)`
Given a fock-space eigenstate label, returns the corresponding
eigenstate of the system Hamiltonian.

# Arguments
- `sys`: Simple QuantumSystem for which to return the Hamiltonian eigenstate
- `lbl`: The state label of the desired system eigenstate. Must be a positive integer, a symbol,
         or a vector of Ints of length 1.
- `bare_to_dressed_tracking`: Either `:overlap` (default) or `:spectrum`. Accepted for API symmetry with
         `CompositeQuantumSystem` and validated, but otherwise ignored since the eigenstates of a
         simple system are already energy-ordered.

# Notes:
- If `lbl` is a positive integer `n`, then the function will return the `nth` eigenstate of the
  system Hamiltonian where `n=0` corresponds to the ground state.
- If `lbl` is a Symbol, it must be a single character in the set ``\{g, e, f, h\}``, which will
  correspond to the state indices ``(0,1,2,3)``.
- If `lbl` is a symbol outside of the allowed set, more than a single character, or an integer not in
.  the interval `[0, size(sys.H[1]))`, the function will throw an `ArgumentError`

# Returns:
- The eigenstate of `sys.H` corresponding to the state label `lbl`.
 """
function state(sys::SimpleQuantumSystem, lbl::StateIndex; bare_to_dressed_tracking::Symbol=:overlap)
    _check_tracking(bare_to_dressed_tracking)
    λ, ψ, _ = eigenstates(sys.H)
    idx = _state_index(lbl, sys)[1]
    return λ[idx], ψ[idx]
end

energy(sys::SimpleQuantumSystem, lbl::StateIndex; kwargs...) = state(sys, lbl; kwargs...)[1]

# Index of the dressed eigenstate of `sys.H` labelled by `lbl`, following the rules of `state`.
function _dressed_idx(sys::SimpleQuantumSystem, lbl::StateIndex; bare_to_dressed_tracking::Symbol=:overlap)
    _check_tracking(bare_to_dressed_tracking)
    return _state_index(lbl, sys)[1]
end

