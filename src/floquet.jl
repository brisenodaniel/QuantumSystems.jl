## Floquet states #######################

function _check_floquet_compat(sys::AbstractQuantumSystem, fb::FloquetBasis, idx::Int, min_overlap::Maybe{Real})
    (size(fb.U_T) == size(sys.H) && fb.U_T.dims == sys.H.dims) || throw(ArgumentError(
        "Dimensions of `fb` $(fb.U_T.dims) do not match those of `sys.H` $(sys.H.dims)"
    ))
    isnothing(min_overlap) && return nothing
    _, ψ, _ = eigenstates(sys)
    ovlp = abs2(ψ[idx]' * fb.ψfl0[idx])
    if ovlp < min_overlap
        @warn "Overlap $(round(ovlp, digits=3)) of Floquet state $idx with the static eigenstate is below " *
            "`min_overlap = $min_overlap`. Check that `fb` was built with the tracked constructor " *
            "`FloquetBasis(H0, H_dr, T)` with `H0 = sys.H`, or decrease its `δλ`."
    end
    return nothing
end

@doc raw"""
    state(sys::QuantumSystem, fb::FloquetBasis, lbl::StateIndex, t::Real=0.0; bare_to_dressed_tracking=:overlap, min_overlap=0.5, cache=false, kwargs...)

Returns the quasienergy and Floquet state at time `t` of `fb` that is adiabatically connected to the
dressed eigenstate of `sys` labelled by `lbl`.

# Arguments
- `sys`: QuantumSystem whose Hamiltonian `sys.H` is the static part of the Floquet Hamiltonian.
- `fb`: FloquetBasis built with the tracked constructor `FloquetBasis(sys.H, H_dr, T)`, so that its
        `k`-th quasienergy and Floquet state are connected to the `k`-th eigenstate of `sys.H`.
- `lbl`: State label, interpreted as in `state(sys, lbl)`.
- `t`: Time at which to return the Floquet state.
- `bare_to_dressed_tracking`: Either `:overlap` (default) or `:spectrum`, see `state(sys, lbl)`.
- `min_overlap`: A warning is raised if ``|\langle \lambda_k | \psi_k(0)\rangle|^2`` is below this value,
        which indicates that `fb` is untracked, built from a different static Hamiltonian, or strongly
        hybridized by the drive. Set to `nothing` to skip the check.
- `cache`: If `true`, stores the micromotion propagator at `t` in `fb` (see `propagator!`).
- `kwargs`: Passed to `propagator`.

# Notes
- The quasienergy is folded onto the first Brillouin zone ``[-\pi/T, \pi/T)`` (see `fold_quasienergy`),
  even though `fb.equasi` itself is unfolded onto the branch of the static energy.
- The Floquet state is ``|\psi_k(t)\rangle = \hat{U}(t, 0)|\psi_k(0)\rangle = e^{-i\epsilon_k t}|u_k(t)\rangle``,
  not the periodic mode ``|u_k(t)\rangle``. Its global phase is that of the eigenvector of ``\hat{U}(T, 0)``,
  and is arbitrary.

# Returns
- The pair ``(\epsilon_k, |\psi_k(t)\rangle)``.
"""
function state(
    sys::AbstractQuantumSystem,
    fb::FloquetBasis,
    lbl::StateIndex,
    t::Real=0.0;
    bare_to_dressed_tracking::Symbol=:overlap,
    min_overlap::Maybe{Real}=0.5,
    cache::Bool=false,
    kwargs...
)
    idx = _dressed_idx(sys, lbl; bare_to_dressed_tracking)
    _check_floquet_compat(sys, fb, idx, min_overlap)
    ψ0 = fb.ψfl0[idx]
    ε = fold_quasienergy(fb.equasi[idx], fb.T)
    if iszero(t)
        return ε, ψ0
    end
    U = cache ? propagator!(fb, t; kwargs...) : propagator(fb, t; kwargs...)
    return ε, U * ψ0
end

"""
    fold_quasienergy(ε::Real, T::Real)

Fold the quasienergy `ε` of a system driven with period `T` onto ``[-\\pi/T, \\pi/T)``.
"""
fold_quasienergy(ε::Real, T::Real) = mod(ε + π / T, 2π / T) - π / T

energy(sys::AbstractQuantumSystem, fb::FloquetBasis, lbl::StateIndex; kwargs...) =
    state(sys, fb, lbl; kwargs...)[1]

