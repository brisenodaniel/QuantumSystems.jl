function _check_ops(ops::NamedTuple)
    msg = nothing
    if !_all_qobj(ops)
        msg = "All elements of `ops` must be QuantumObjects"
    elseif !(values(ops) |> same_dims)
        msg = "All QuantumObjects in `ops` must have the same dimension"
    elseif !hasproperty(ops, :H)
        msg = "`ops` must contain Hamiltonian operator with key `:H`"
    end
    return isnothing(msg) ? msg : ArgumentError(msg)
end

function _state_index(lbl::StateIndex, sys::AbstractQuantumSystem)
    if lbl isa Int
        idx = [lbl + 1,]
    elseif lbl isa AbstractVector{<:Int}
        idx = lbl .+ 1
    else
        _lbl = (collect ∘ String)(lbl)
        idx =  map(_lbl) do l_n
            l_n = Symbol(l_n)
            return _STATE_SYMBOLS[l_n]
        end
    end
    _check_idx(idx, sys)
    return idx
end

function _check_idx(lbl::AbstractVector{<:Int}, sys::AbstractQuantumSystem)
    tdims = length(lbl) == 1 ? (size(sys.H)[1],) : tensor_dims(sys)
    length(tdims) == length(lbl) || throw(ArgumentError(
        "Number of characters or integers in `lbl` must match number of subsystems in `sys`"
    ))
    idx_checks = map(lbl, tdims) do l_n, t_n
        l_n ∉ 1:t_n
    end
    bad_idx = findall(idx_checks)
    return isempty(bad_idx) || throw(ArgumentError(
        "State indices at `lbl` indices $bad_idx lie outside of the range allowed  " *
            "by the Hilbert space dimension of their subsystem "
    ))
end

function same_dims(ops_vals)
    if isempty(ops_vals)
        return true
    end
    dim0 = first(ops_vals).dimensions
    return all(q -> q.dimensions == dim0, ops_vals)
end

function _all_qobj(ops::NamedTuple)
    if isempty(ops)
        return false
    end
    return all(q -> q isa AbstractQuantumObject, ops)
end



function tensor_dims(sys::AbstractQuantumSystem)
    if sys isa SimpleQuantumSystem
        return collect(sys.H.dims)
    else
        return ((tensor_dims(subsys)[1] for subsys in sys.subsys)...,)
    end
end

function QuantumToolbox.qeye(sys::QuantumSystem)
    return qeye_like(sys.H)
end

function I(sys::QuantumSystem)
    return qeye(sys)
end

function Base.getproperty(sys::AbstractQuantumSystem, key::Symbol)
    if hasfield(typeof(sys), key)
        return Core.getfield(sys, key)
    elseif hasproperty(Core.getfield(sys, :ops), key)
        return Core.getfield(sys, :ops)[key]
    elseif sys isa CompositeQuantumSystem &&
        hasproperty(Core.getfield(sys, :subsys), key)
        return Core.getfield(sys, :subsys)[key]
    elseif key==:I
        return qeye(sys)
    else
        throw(FieldError(typeof(sys), key))
    end
end

@memoize LRU{Any, Any}(maxsize=100) function QuantumToolbox.eigenstates(sys::AbstractQuantumSystem)
    eigenstates(sys.H)
end

function _check_tracking(mode::Symbol)
    mode in (:overlap, :spectrum) || throw(ArgumentError(
        "`bare_to_dressed_tracking` must be `:overlap` or `:spectrum`, got `:$mode`"
    ))
    return mode
end
