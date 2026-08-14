# 0PASchwarz inspiral forcing (InspiralModels/0PASchwarz.m): adiabatic (0PA) evolution of a
# circular Schwarzschild EMRI. The forcing `F0[r0]` is a closed-form prefactor times the total
# 1SF energy flux `F1(r0) = F_I(r0) + F_H(r0)` (infinity + horizon, order-8 interpolants). The
# inspiral evolves `(Ω, φ)` with `r0 = Ω^{-2/3}` (the constraint `Ω = √(1/r0³)`), `ν` and `M`
# constant, terminating at `r0 ≤ 6.01`.

"""
    Forcing0PASchwarz

0PASchwarz forcing model: the horizon and infinity 1SF energy-flux interpolants (`FEH`, `FEI`),
order-8 [`Interp1`](@ref). Build with [`forcing_0PASchwarz`](@ref); evaluate the forcing with
[`forcing_F0`](@ref).
"""
struct Forcing0PASchwarz
    FEH::Interp1{Float64,Float64}   # horizon 1SF flux
    FEI::Interp1{Float64,Float64}   # infinity 1SF flux
end

"""
    forcing_0PASchwarz(datadir) -> Forcing0PASchwarz

Load the `1SFCircShwarzDotE{Horizon,Inf}.m` 1SF flux tables from `datadir`
(`sf_data/1SF_Flux/Schwarz_Circ`) and build the order-8 flux interpolants.
"""
function forcing_0PASchwarz(datadir::AbstractString)
    FEH = interp1_order8(read_wl_table(joinpath(datadir, "1SFCircShwarzDotEHorizon.m")))
    FEI = interp1_order8(read_wl_table(joinpath(datadir, "1SFCircShwarzDotEInf.m")))
    return Forcing0PASchwarz(FEH, FEI)
end

"""
    flux_1SF(model, r0)

Total 1SF energy flux `F1(r0) = F_I(r0) + F_H(r0)` (infinity + horizon), 0PASchwarz.m.
"""
flux_1SF(m::Forcing0PASchwarz, r0::Real) = m.FEI(r0) + m.FEH(r0)

# Closed-form forcing prefactor `3(1-3M/r0)^{3/2}√(M/r0) / (M²(1-6M/r0))` (M = 1). Shared by
# the 0PA `F0` and 1PA `F1` forcing terms (0PASchwarz.m, 1PAT1.m).
_forcing_prefactor(r0::Real) = (3 * (1 - 3 / r0)^(3 / 2) * sqrt(1 / r0)) / (1 - 6 / r0)

"""
    forcing_F0(model, r0)

Adiabatic forcing `F0[r0] = 3(1-3M/r0)^{3/2}√(M/r0)/(M²(1-6M/r0)) · F1(r0)` (M = 1),
0PASchwarz.m.
"""
forcing_F0(m::Forcing0PASchwarz, r0::Real) = _forcing_prefactor(r0) * flux_1SF(m, r0)

"""
    omega_dot(model, r0, ν, M)

`dΩ/dt` of the inspiral. For 0PASchwarz this is `ν/M · F0(r0)` (0PASchwarz.m).
"""
omega_dot(m::Forcing0PASchwarz, r0::Real, ν::Real, M::Real) = (ν / M) * forcing_F0(m, r0)

"""
    stop_r0(model, ν)

Terminating radius of the inspiral: integration stops when `r0 ≤ stop_r0`. For 0PASchwarz this
is the constant `6.01` (0PASchwarz.m, independent of `ν`).
"""
stop_r0(::Forcing0PASchwarz, ν::Real) = 6.01

"""
    inspiral_rhs!(du, u, p, t)

In-place ODE right-hand side for the inspiral state `u = (Ω, φ)` with `p = (model, ν, M)`:
`dΩ/dt = omega_dot(model, r0, ν, M)`, `dφ/dt = Ω/M`, `r0 = Ω^{-2/3}`. Model-agnostic; the
per-model physics lives in [`omega_dot`](@ref).
"""
function inspiral_rhs!(du, u, p, t)
    m, ν, M = p
    Ω, φ = u
    r0 = Ω^(-2 / 3)
    du[1] = omega_dot(m, r0, ν, M)
    du[2] = Ω / M
    return nothing
end

"""
    inspiral_rhs(model, Ω, φ, ν, M) -> (dΩ, dφ)

Out-of-place form of [`inspiral_rhs!`](@ref): returns the `(dΩ/dt, dφ/dt)` rhs for `(Ω, φ)`.
"""
function inspiral_rhs(m, Ω::Real, φ::Real, ν::Real, M::Real)
    r0 = Ω^(-2 / 3)
    return (omega_dot(m, r0, ν, M), Ω / M)
end
