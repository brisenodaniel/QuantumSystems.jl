# Test-only fluxonium builder (mirrors XoniumUnitCell.fluxonium) so the
# QuantumSystems tests do not depend on any circuit library.
function φn_operators(; EC::Real, EL::Real, nlev::Int)
    Z = EC / EL
    φ_zpf = (8.0 * Z)^(1/4)
    n_zpf = (8.0 * Z)^(-1/4)
    return φ_zpf * QuantumToolbox.position(nlev), n_zpf * QuantumToolbox.momentum(nlev)
end

function fluxonium(; EC::Real, EL::Real, EJ::Real, δϕ::Real = π, n_g::Real = 0.0,
                   nlev::Int = 40, name::Symbol = :Fluxonium)
    φ, n = φn_operators(; EC=EC, EL=EL, nlev=nlev)
    H = 4.0 * EC * (n - n_g)^2 + (1.0/2.0) * EL * (φ - δϕ)^2 - EJ * cos(φ)
    ops = (; φ=φ, n=n, H=H)
    params = (; EC=EC, EL=EL, EJ=EJ, δϕ=δϕ, Φext_incos=false)
    return SimpleQuantumSystem(ops; params=params, name=name)
end
