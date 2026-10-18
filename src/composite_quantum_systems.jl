struct CompositeQuantumSystem{
    Ops<:NamedTuple,
    Systms<:NamedTuple,
    PT<:NamedTuple,} <: AbstractQuantumSystem

    subsys::Systms
    ops::Ops
    params::Maybe{PT}
    U_basis :: Maybe{AbstractQuantumObject}
    dressed_dim :: Maybe{Int}
    name :: Maybe{Symbol}
    state_lookup::AbstractArray{Union{Int64, Missing}}

    function CompositeQuantumSystem(
        subsys::Systms,
        ops::Ops;
        params::Maybe{PT}=nothing,
        U_basis::Maybe{AbstractQuantumObject}=nothing,
        dressed_dim::Maybe{Int}=nothing,
        name::Maybe{Symbol}=nothing,
        ) where{
            Ops<:NamedTuple,
            Systms<:NamedTuple,
            PT<:NamedTuple,
        }
        err = _check_subsys(subsys)
        err = isnothing(err) ? _check_ops(ops) : err
        isnothing(err) ? nothing : throw(err)
        n_sys = length(subsys)
        tdims = (size(sys.H)[1] for (_, sys) in pairs(subsys)) |> collect
        state_lookup = Array{Union{Int64, Missing}, n_sys}(missing, tdims...)
        return new{Ops, Systms, PT}(subsys, ops, params, U_basis, dressed_dim, name, state_lookup)
    end
end

"""
    _tensor_identity(H)

Return identity operator on the Hilbert space of operator `H`.
"""
_tensor_identity(H::AbstractQuantumObject) = qeye_like(H)


"""
    _lift_left(op, Hright)

Return op ⊗ I
"""
_lift_left(op::AbstractQuantumObject, Hright::AbstractQuantumObject) =
    tensor(op, _tensor_identity(Hright))


"""
    _lift_right(Hleft, op)

Return I ⊗ op
"""
_lift_right(Hleft::AbstractQuantumObject, op::AbstractQuantumObject) =
    tensor(_tensor_identity(Hleft), op)

function _lift_local_operator(
    sys::CompositeQuantumSystem,
    subkey::Symbol,
    opkey::Symbol,
)
    # --- sanity checks ---
    hasproperty(sys.subsys, subkey) ||
        throw(KeyError("No subsystem named $subkey"))

    subsys = sys.subsys[subkey]

    hasproperty(subsys.ops, opkey) ||
        throw(KeyError("Subsystem $subkey has no operator $opkey"))

    # --- build tensor factors in correct order ---
    factors = map(keys(sys.subsys)) do name
        if name === subkey
            subsys.ops[opkey]
        else
            qeye(sys.subsys[name])
        end
    end

    # tensor in order
    O = reduce(tensor, factors)

    # --- apply dressed basis if present ---
    if sys.U_basis !== nothing
        U = sys.U_basis
        O = U' * O * U
    end

    return O
end

function _lift_multi_operator(
    sys::CompositeQuantumSystem,
    names::NTuple{N,Symbol},
) where {N}
    # --- must be pairs ---
    iseven(N) ||
        throw(ArgumentError("Expected (subsystem, operator) pairs"))

    # --- parse pairs ---
    pairs = Pair{Symbol,Symbol}[
        names[i] => names[i+1] for i in 1:2:N
    ]

    # --- forbid duplicate subsystems ---
    subs = first.(pairs)
    length(unique(subs)) == length(subs) ||
        throw(ArgumentError(
            "Each subsystem may appear at most once in operator lifting"
        ))

    # --- validate subsystems and operators ---
    for (subkey, opkey) in pairs
        hasproperty(sys.subsys, subkey) ||
            throw(ArgumentError("No subsystem named $subkey"))

        hasproperty(sys.subsys[subkey].ops, opkey) ||
            throw(ArgumentError(
                "Subsystem $subkey has no operator $opkey"
            ))
    end

    # --- map requested local operators ---
    local_ops = Dict{Symbol,AbstractQuantumObject}(
        subkey => sys.subsys[subkey].ops[opkey]
        for (subkey, opkey) in pairs
    )

    # --- build tensor factors in subsystem order ---
    factors = map(keys(sys.subsys)) do name
        if haskey(local_ops, name)
            local_ops[name]
        else
           qeye(sys.subsys[name])
        end
    end

    # --- tensor product ---
    O = reduce(tensor, factors)

    # --- apply dressed basis if present ---
    if sys.U_basis !== nothing
        U = sys.U_basis
        O = U' * O * U
    end

    return O
end

function _check_subsys(subsys::NamedTuple)
    msg = nothing
    if isempty(subsys) || !all(
        sys -> sys isa AbstractQuantumSystem, subsys
        )
        msg = "Parameter `subsys` must be a non-empty NamedTuple with " *
            "`AbstractQuantumSystem` values"
    end
    return isnothing(msg) ? msg : ArgumentError(msg)
end

function Base.getindex(
    sys::CompositeQuantumSystem,
    key::Symbol,
)
    if hasproperty(sys.ops, key)
        return sys.ops[key]
    elseif hasproperty(sys.subsys, key)
        return sys.subsys[key]
    else
        throw(KeyError("No operator or subsystem named $key"))
    end
end


function Base.getindex(
    sys::CompositeQuantumSystem,
    subkey::Symbol,
    opkey::Symbol,
)
    return _lift_local_operator(sys, subkey, opkey)
end


function Base.getindex(
    sys::CompositeQuantumSystem,
    keys::Vararg{Symbol},
)
    return _lift_multi_operator(sys, keys)
end

function to_energy_basis(
    sys::CompositeQuantumSystem;
    set_energy_zero::Union{Bool, Float64} = true,
    etol::Union{Bool, Float64} = false,
)
    # --- 1. Diagonalize composite Hamiltonian ---
    λ, _, U = eigenstates(sys.H)

    # Wrap U as a QuantumObject so algebra stays in QuantumToolbox
    Uq = QuantumObject(U, sys.H.type, sys.H.dims)

    # --- 2. Similarity transform of Hamiltonian ---
    H_diag = Uq' * sys.H * Uq

    # --- 3. Optional energy shift ---
    if set_energy_zero !== false
        e0 = set_energy_zero === true ? (minimum ∘ real)(λ) : set_energy_zero
        H_diag = H_diag - e0
    end

    # --- 4. Build new ops ---
    new_ops = (; H = H_diag)

    # --- 5. Construct dressed system ---
    dressed = CompositeQuantumSystem(
        sys.subsys,
        new_ops;
        params       = sys.params,
        U_basis      = Uq,
        dressed_dim  = length(λ),
        name         = sys.name,
    )

    # --- 6. Optional tidyup ---
    if etol isa Float64
        return tidyup(dressed, etol)
    elseif etol === true
        return tidyup(dressed)
    else
        return dressed
    end
end

function QuantumToolbox.tidyup(
    sys::CompositeQuantumSystem,
    etol::Maybe{Float64} = nothing,
)
    # Resolve tolerance
    tol = isnothing(etol) ? QuantumToolbox.settings.tidyup_tol : etol

    # --- 1. Tidy composite operators ---
    new_ops = (; (
        k => tidyup(op, tol)
        for (k, op) in pairs(sys.ops)
    )...)

    # --- 2. Tidy subsystems recursively ---
    new_subsys = (; (
        name => tidyup(subsys, tol)
        for (name, subsys) in pairs(sys.subsys)
    )...)

    # --- 3. Tidy basis transform if present ---
    new_U_basis = isnothing(sys.U_basis) ? nothing : tidyup(sys.U_basis, tol)

    # --- 4. Rebuild CompositeQuantumSystem ---
    return CompositeQuantumSystem(
        new_subsys,
        new_ops;
        params      = sys.params,
        U_basis     = new_U_basis,
        dressed_dim = sys.dressed_dim,
        name        = sys.name,
    )
end

@doc raw"""
    ``state(sys::CompositeQuantumSystem, lbl::StateIndex)

Given a fock-space eigenstate label, returns the corresponding eigenstate
of the system Hamiltonian.

# Arguments
- `sys`: CompositeQuantumSystem for which to return the Hamiltonian eigenstate.
- `lbl`: The state label of the desired system eigenstate.
- `bare`: Boolean flag which sets whether to return the eigenstate in the bare product
          basis or the dressed eigenbasis of `sys.H`.
- `bare_to_dressed_tracking`: Keyword, either `:overlap` (default) or `:spectrum`. Only used when
          `bare==false` and `lbl` is a multi-subsystem label.
        - `:overlap`: the dressed eigenstate with maximal overlap with the bare product state.
        - `:spectrum`: the dressed eigenstate whose energy rank equals the energy rank of the
          bare product state in the bare spectrum (e.g. if `:ge` is the 2nd excited bare state,
          the 2nd excited eigenstate of `sys.H` is returned). Degenerate bare levels are ordered
          by label order (first subsystem index varies fastest), so the result is somewhat
          arbitrary for them; prefer `:overlap` there.

# Notes
- For a `sys` with `N` subsystems, `lbl` must be a positive integer, a Symbol of `N` characters,
  or a `Vector{Int, N}`.
- If `lbl` is a single positive integer `k`, function will return `eigenstates(sys.H)[2][k]`
- If `lbl` is a `Vector{Int, N}` or a `Symbol` of `N` characters in the set ``\{g,e,f,h\}``:
      - When `bare==true`, will return the kronecker product
             ``math
             \hat{U}\left(\bigotimes_{n=0}^N |l[n]\rangle\right) \hat{U}^\dagger
             ``
             where ``l_n`` is the `nth` element or character in `lbl`, `|l_n\rangle` is the corresponding eigenstate of
             `nth` subsystem, and `U` is the basis transform to the basis defined by `sys.U_basis`. If `isnothing(sys.U_basis)`,
             then ``U=1``.
        - When `bare==false`, will return the eigenstate of `sys.H` with the maximal value of
             ``math
             |\langle l | k\rangle|
             ``
             where ``| l \langle`` is the output of `state(sys, lbl, true)` and `|k\rangle` is any eigenstate of
             sys.H`.

# Returns
- The dressed or bare eigenenergy and eigenstate pair ``(\lambda, \psi)`` of `sys` specified by `lbl`.
"""
function state(
    sys::CompositeQuantumSystem,
    lbl::StateIndex,
    bare::Bool=false;
    bare_to_dressed_tracking::Symbol=:overlap,
)
    _check_tracking(bare_to_dressed_tracking)
    if bare
        idx = _state_index(lbl, sys)
        return _bare_state(sys, idx)
    else
        idx = state_idx(sys, lbl; bare_to_dressed_tracking)
        λ, ψ, _ = eigenstates(sys.H)
        return λ[idx], ψ[idx]
    end
end

function _bare_state(sys::CompositeQuantumSystem, idx::AbstractVector{<:Int})

    if length(idx) == 1
        i = idx[1]
        H = tensor((sub.H for sub in sys.subsys)...)
        λ, ψ, _ = eigenstates(H)
        return λ[i], ψ[i]
    end
    ens_states = map(sys.subsys, idx) do sub, i
        λ_sub, ψ_sub, _ = eigenstates(sub.H)
        return λ_sub[i], ψ_sub[i]
    end
    λ, ψ = first.(ens_states) |> prod, last.(ens_states) |> x -> tensor(x...)
    U = sys.U_basis
    if isnothing(U)
        return λ, ψ
    else
        return λ, U'*ψ
    end
end

function state_idx(
    sys::CompositeQuantumSystem,
    lbl::StateIndex,
    bare::Bool=false;
    bare_to_dressed_tracking::Symbol=:overlap,
)
    _check_tracking(bare_to_dressed_tracking)
    idx = _state_index(lbl, sys)
    if length(idx) == 1
        return idx[1]
    elseif !bare && bare_to_dressed_tracking == :spectrum
        return _spectrum_idx(sys, idx)
    elseif !bare && !ismissing(sys.state_lookup[idx...])
        return sys.state_lookup[idx...]
    else
        H = bare ? tensor((sub.H for sub in sys.subsys)...) : sys.H
        _, ψ, _ = eigenstates(H)
        _, ψbare = _bare_state(sys, idx)
        best_idx = abs.(ψbare' * ψn for ψn in ψ) |> argmax
        if !bare
            sys.state_lookup[idx...] = best_idx
        end
        return best_idx
    end
end

# Energy rank of the bare product state `idx` among all bare product states.
function _spectrum_idx(sys::CompositeQuantumSystem, idx::AbstractVector{<:Int})
    sub_λ = [real.(first(eigenstates(sub.H))) for sub in sys.subsys]
    # first subsystem index varies fastest, matching the kron ordering of `tensor`
    inds = CartesianIndices(Tuple(length.(sub_λ)))
    energies = [sum(sub_λ[n][I[n]] for n in eachindex(sub_λ)) for I in inds]
    order = sortperm(vec(energies))  # stable: ties broken by label order
    target = LinearIndices(inds)[idx...]
    return findfirst(==(target), order)
end

energy(sys::CompositeQuantumSystem, lbl::StateIndex, bare::Bool=false; kwargs...) =
    state(sys, lbl, bare; kwargs...)[1]

_dressed_idx(sys::CompositeQuantumSystem, lbl::StateIndex; bare_to_dressed_tracking::Symbol=:overlap) =
    state_idx(sys, lbl; bare_to_dressed_tracking)

