"""
    _require_name(sys)

Throw if `sys.name === nothing`, otherwise return `sys.name`.
"""
_require_name(sys::AbstractQuantumSystem) =
    isnothing(sys.name) ?
        throw(ArgumentError("All subsystems must have a non-nothing `name`")) :
        sys.name


function _check_ops_addition(
    sys::AbstractQuantumSystem,
    ops::NamedTuple,
)
    for (k, op) in pairs(ops)
        op isa AbstractQuantumObject ||
            throw(ArgumentError(
                "Value for key $k must be an AbstractQuantumObject"
            ))

        sys.H.dimensions == op.dimensions ||
            throw(DimensionMismatch(
                "Operator $k has incompatible dimensions with system Hamiltonian"
            ))
    end
end

## Composition ##########################

function QuantumToolbox.:⊗(A::SimpleQuantumSystem, B::SimpleQuantumSystem)
    # --- 1. Require names ---
    isnothing(A.name) && throw(ArgumentError(
        "Left-hand SimpleQuantumSystem must have a non-nothing `name`"
    ))
    isnothing(B.name) && throw(ArgumentError(
        "Right-hand SimpleQuantumSystem must have a non-nothing `name`"
    ))

    nameA = A.name
    nameB = B.name

    # --- 2. Build subsystem NamedTuple (order matters!) ---
    subsys = (; nameA => A, nameB => B)

    # --- 3. Build composite Hamiltonian ---
    H = tensor(A.H, qeye(B)) + tensor(qeye(A), B.H)

    ops = (; H = H)

    # --- 4. Build params NamedTuple ---
    params = (; nameA => A.params, nameB => B.params)

    # --- 5. Construct CompositeQuantumSystem ---
    return CompositeQuantumSystem(
        subsys,
        ops;
        params = params,
        name = Symbol(nameA, :_, nameB),
    )
end

function QuantumToolbox.:⊗(A::SimpleQuantumSystem, B::CompositeQuantumSystem)
    nameA = _require_name(A)

    # Subsystems (A comes first)
    subsys = (; nameA => A, B.subsys...)

    # Hamiltonian
    H = tensor(A.H, qeye(B)) + tensor(qeye(A), B.H)

    ops = (; H = H)

    # Params
    params = (; nameA => A.params, B.params...)

    return CompositeQuantumSystem(
        subsys,
        ops;
        params = params,
        name = Symbol(nameA, :_, B.name),
    )
end

function QuantumToolbox.:⊗(A::CompositeQuantumSystem, B::SimpleQuantumSystem)
    nameB = _require_name(B)

    subsys = (; A.subsys..., nameB => B)

    H = tensor(A.H, qeye(B)) + tensor(qeye(A), B.H)

    ops = (; H = H)

    params = (; A.params..., nameB => B.params)

    return CompositeQuantumSystem(
        subsys,
        ops;
        params = params,
        name = Symbol(A.name, :_, nameB),
    )
end

function QuantumToolbox.:⊗(A::CompositeQuantumSystem, B::CompositeQuantumSystem)
    # Ensure no name collisions
    overlap = intersect(keys(A.subsys), keys(B.subsys))
    isempty(overlap) || throw(ArgumentError(
        "Subsystem name collision in composite product: $(collect(overlap))"
    ))

    subsys = (; A.subsys..., B.subsys...)

    H = tensor(A.H, qeye(B)) + tensor(qeye(A), B.H)

    ops = (; H = H)

    params = (; A.params..., B.params...)

    return CompositeQuantumSystem(
        subsys,
        ops;
        params = params,
        name = Symbol(A.name, :_, B.name),
    )
end

function Base.:+(sys::SimpleQuantumSystem, Op::AbstractQuantumObject)
    # Dimension check (QuantumToolbox will throw if incompatible,
    # but we fail early with a clearer message)
    sys.H.dimensions == Op.dimensions ||
        throw(DimensionMismatch(
            "Operator dimensions do not match system Hamiltonian"
        ))

    new_H = sys.H + Op
    new_ops = (; sys.ops..., H = new_H)

    return SimpleQuantumSystem(
        new_ops;
        params = sys.params,
        name   = sys.name,
    )
end

function Base.:+(sys::CompositeQuantumSystem, Op::AbstractQuantumObject)
    sys.H.dimensions == Op.dimensions ||
        throw(DimensionMismatch(
            "Operator dimensions do not match composite Hamiltonian"
        ))

    new_H = sys.H + Op
    new_ops = (; sys.ops..., H = new_H)

    return CompositeQuantumSystem(
        sys.subsys,
        new_ops;
        params      = sys.params,
        U_basis     = sys.U_basis,
        dressed_dim = sys.dressed_dim,
        name        = sys.name,
    )
end

function Base.:+(sys::SimpleQuantumSystem, ops::NamedTuple)
    _check_ops_addition(sys, ops)

    new_ops = merge(sys.ops, ops)

    return SimpleQuantumSystem(
        new_ops;
        params = sys.params,
        name   = sys.name,
    )
end

function Base.:+(sys::CompositeQuantumSystem, ops::NamedTuple)
    _check_ops_addition(sys, ops)

    new_ops = merge(sys.ops, ops)

    return CompositeQuantumSystem(
        sys.subsys,
        new_ops;
        params      = sys.params,
        U_basis     = sys.U_basis,
        dressed_dim = sys.dressed_dim,
        name        = sys.name,
    )
end

Base.:+(ops::NamedTuple, sys::AbstractQuantumSystem) = sys + ops
Base.:+(Op::AbstractQuantumObject, sys::AbstractQuantumSystem) = sys + Op

Base.:-(sys::AbstractQuantumSystem, Op::AbstractQuantumObject) = sys +(-1*Op)

