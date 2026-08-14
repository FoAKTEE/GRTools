# 1PAT1ea inspiral forcing (`InspiralModels/1PAT1ea.m`): the 1PAT1 Schwarzschild-circular
# forcing extended with the linear-in-spin and mass-loss channels of the `1PAT1ea` model —
# secondary-spin `Fσ = F2σI + F2σH`, primary-(slow-)spin `Fa = FaI + FaH` (both order-8
# interpolants over rational-abscissa tables), and the closed-form channel functions
# `F1s1` (primary spin χt1), `F1s2` (secondary spin χt2), `F1ExtraTerm` (horizon-flux χt1
# back-reaction) and `F1δM` (mass-loss δm). All take `M = 1` as upstream. The `F0`/`F1`
# channels are inherited unchanged from the 1PAT1 base ([`Forcing1PAT1`](@ref)).

"""
    Forcing1PAT1ea

1PAT1ea forcing model: the 1PAT1 base ([`Forcing1PAT1`](@ref): 1SF/2SF flux + redshift) plus
the secondary-spin (`F2σI`, `F2σH`) and slow-primary-spin (`FaI`, `FaH`) order-8 flux
interpolants. Build with [`forcing_1PAT1ea`](@ref). Evaluate channels with [`forcing_F0`](@ref),
[`forcing_F1`](@ref), [`forcing_F1s1`](@ref), [`forcing_F1s2`](@ref),
[`forcing_F1extra`](@ref), [`forcing_F1dm`](@ref).
"""
struct Forcing1PAT1ea
    base::Forcing1PAT1
    F2σI::Interp1{Float64,Float64}   # secondary-spin flux, infinity
    F2σH::Interp1{Float64,Float64}   # secondary-spin flux, horizon
    FaI::Interp1{Float64,Float64}    # slow-χ1 flux, infinity
    FaH::Interp1{Float64,Float64}    # slow-χ1 flux, horizon
end

"""
    forcing_1PAT1ea(sf_datadir) -> Forcing1PAT1ea

Load the 1PAT1 base tables plus the secondary-spin `Leading_S2_Flux/Schwarz_Circ/chi2_2SF*`
and slow-χ1 `1SF_Flux/Linear_Kerr_Circ/Slow_chi1_2SF*` flux tables under `sf_datadir`
(`InspiralModels/sf_data`) and build the order-8 spin-flux interpolants.
"""
function forcing_1PAT1ea(sf_datadir::AbstractString)
    base = forcing_1PAT1(sf_datadir)
    dsig = joinpath(sf_datadir, "Leading_S2_Flux", "Schwarz_Circ")
    da = joinpath(sf_datadir, "1SF_Flux", "Linear_Kerr_Circ")
    F2σI = interp1_order8(read_wl_table(joinpath(dsig, "chi2_2SFCircShwarzDotEInf.m")))
    F2σH = interp1_order8(read_wl_table(joinpath(dsig, "chi2_2SFCircShwarzDotEHorizon.m")))
    FaI = interp1_order8(read_wl_table(joinpath(da, "Slow_chi1_2SFCircShwarzDotEInf.m")))
    FaH = interp1_order8(read_wl_table(joinpath(da, "Slow_chi1_2SFCircShwarzDotEHorizon.m")))
    return Forcing1PAT1ea(base, F2σI, F2σH, FaI, FaH)
end

# Base-flux delegates (the 1PAT1 base owns the 1SF/2SF flux + redshift).
flux_1SF(m::Forcing1PAT1ea, r0::Real) = flux_1SF(m.base, r0)
flux_2SF(m::Forcing1PAT1ea, r0::Real) = flux_2SF(m.base, r0)

"""
    flux_LH(model::Forcing1PAT1ea, r0) -> Float64

Horizon angular-momentum flux `FLH(r0) = r0^{3/2}·F_H(r0)` (`1PAT1ea.m`/`1PAT1R.m`), driving
the primary-spin evolution `χt1'`.
"""
flux_LH(m::Forcing1PAT1ea, r0::Real) = r0^(3 / 2) * m.base.FEH(r0)

"""
    flux_sigma(model::Forcing1PAT1ea, r0) -> Float64

Secondary-spin energy-flux sum `Fσ(r0) = F2σI(r0) + F2σH(r0)` (`1PAT1ea.m`).
"""
flux_sigma(m::Forcing1PAT1ea, r0::Real) = m.F2σI(r0) + m.F2σH(r0)

"""
    flux_a(model::Forcing1PAT1ea, r0) -> Float64

Slow-primary-spin energy-flux sum `Fa(r0) = FaI(r0) + FaH(r0)` (`1PAT1ea.m`).
"""
flux_a(m::Forcing1PAT1ea, r0::Real) = m.FaI(r0) + m.FaH(r0)

"""
    forcing_F0(model::Forcing1PAT1ea, r0) -> Float64

0PA forcing `F0[r0]` (inherited from the 1PAT1 base; `1PAT1ea.m`).
"""
forcing_F0(m::Forcing1PAT1ea, r0::Real) = forcing_F0(m.base, r0)

"""
    forcing_F1(model::Forcing1PAT1ea, r0) -> Float64

1PA forcing `F1[r0]` (inherited from the 1PAT1 base; `1PAT1ea.m`).
"""
forcing_F1(m::Forcing1PAT1ea, r0::Real) = forcing_F1(m.base, r0)

"""
    forcing_F1s1(model::Forcing1PAT1ea, r0) -> Float64

Primary-spin (χt1) forcing channel `F1s1[r0] = pref·Fa(r0)
- 2√(1-3/r0)(36-33r0+10r0²)/(r0²(r0-6)²)·F1flux(r0)` (M = 1; `1PAT1ea.m`).
"""
function forcing_F1s1(m::Forcing1PAT1ea, r0::Real)
    pref = _forcing_prefactor(r0)
    back = (2 * sqrt(1 - 3 / r0) * (36 - 33r0 + 10r0^2)) / (r0^2 * (r0 - 6)^2)
    return pref * flux_a(m, r0) - back * flux_1SF(m, r0)
end

"""
    forcing_F1s2(model::Forcing1PAT1ea, r0) -> Float64

Secondary-spin (χt2) forcing channel `F1s2[r0] = pref·Fσ(r0)
- (3/r0²)·(1-3/r0)^{3/2}(5-12/r0)/(1-6/r0)²·F1flux(r0)` (M = 1; `1PAT1ea.m`).
"""
function forcing_F1s2(m::Forcing1PAT1ea, r0::Real)
    pref = _forcing_prefactor(r0)
    back = (3 / r0^2) * ((1 - 3 / r0)^(3 / 2) * (5 - 12 / r0)) / (1 - 6 / r0)^2
    return pref * flux_sigma(m, r0) - back * flux_1SF(m, r0)
end

"""
    forcing_F1extra(model::Forcing1PAT1ea, r0) -> Float64

Horizon-flux χt1 back-reaction channel `F1ExtraTerm[r0]
= -(1/r0)^{3/2}·(5-12/r0)/(1-6/r0)·F_H(r0)` (M = 1; `1PAT1ea.m`).
"""
forcing_F1extra(m::Forcing1PAT1ea, r0::Real) =
    -(1 / r0)^(3 / 2) * ((5 - 12 / r0) / (1 - 6 / r0)) * m.base.FEH(r0)

"""
    forcing_F1dm(model::Forcing1PAT1ea, r0) -> Float64

Mass-loss (δm) forcing channel `F1δM[r0] = pref·(√(1/r0³)(-2r0^{5/2}/3)·F1flux'(r0))
+ √(1-3/r0)√(1/r0)(3(7-6/r0)/r0-2)/(1-6/r0)²·F1flux(r0)` (M = 1; `1PAT1ea.m`).
Reuses [`flux_1SF_deriv`](@ref) for `F1flux'`.
"""
function forcing_F1dm(m::Forcing1PAT1ea, r0::Real)
    pref = _forcing_prefactor(r0)
    t1 = pref * (sqrt(1 / r0^3) * (-2 * r0^(5 / 2) / 3) * flux_1SF_deriv(m.base, r0))
    t2 = (sqrt(1 - 3 / r0) * sqrt(1 / r0) * (3 * (7 - 6 / r0) / r0 - 2)) /
         (1 - 6 / r0)^2 * flux_1SF(m, r0)
    return t1 + t2
end

# --- 7-var inspiral models (1PAT1ea additive / 1PAT1R ratio form) -----------

# The frequency floor shared by both stop conditions: Ω ≥ Min[1/6.06^{3/2}, 1/6.26^{3/2}] =
# 1/6.26^{3/2} (the larger radius → smaller frequency; `1PAT1ea.m`/`1PAT1R.m` stopcondition).
const _OMEGA_FLOOR_1PA = 1 / 6.26^(3 / 2)

"""
    Model1PAT1ea

7-variable `1PAT1ea` inspiral model: the [`Forcing1PAT1ea`](@ref) channels plus the 3D
`Ωcritea(ν, χt1, χt2)` stop interpolant ([`Interp3Lagrange`](@ref) over `Stop1PAT1ea.m`).
State `(Ω, φ, ν, M, χt1, χt2, δm)`; `Ω'` is the additive `F0 + F1 + spin/δm` forcing
(`1PAT1ea.m`). Build with [`model_1PAT1ea`](@ref).
"""
struct Model1PAT1ea
    forcing::Forcing1PAT1ea
    ωcrit::Interp3Lagrange{Float64}
end

"""
    Model1PAT1R

7-variable `1PAT1R` inspiral model: same data as [`Model1PAT1ea`](@ref) but `Ω'` uses the ratio
form `dΩdt·dEdΩ / dEdΩ` (`1PAT1R.m`) and the `ΩcritR` stop interpolant (`Stop1PAT1R.m`). Build
with [`model_1PAT1R`](@ref).
"""
struct Model1PAT1R
    forcing::Forcing1PAT1ea
    ωcrit::Interp3Lagrange{Float64}
end

# Stop-condition directory shared by both models (`InspiralModels/sf_data/stop_conditions`).
_stop_interp(sf_datadir, fname) =
    interp3_lagrange(read_wl_grid3d(joinpath(sf_datadir, "stop_conditions", fname)))

"""
    model_1PAT1ea(sf_datadir) -> Model1PAT1ea

Build the `1PAT1ea` inspiral model from `sf_datadir` (`InspiralModels/sf_data`): the spin/δm
forcing channels plus the `Stop1PAT1ea.m` 3D Ωcrit interpolant.
"""
model_1PAT1ea(sf_datadir::AbstractString) =
    Model1PAT1ea(forcing_1PAT1ea(sf_datadir), _stop_interp(sf_datadir, "Stop1PAT1ea.m"))

"""
    model_1PAT1R(sf_datadir) -> Model1PAT1R

Build the `1PAT1R` inspiral model from `sf_datadir` (`InspiralModels/sf_data`): the spin/δm
forcing channels plus the `Stop1PAT1R.m` 3D Ωcrit interpolant.
"""
model_1PAT1R(sf_datadir::AbstractString) =
    Model1PAT1R(forcing_1PAT1ea(sf_datadir), _stop_interp(sf_datadir, "Stop1PAT1R.m"))

"""
    omega_crit(model::Union{Model1PAT1ea,Model1PAT1R}, ν, χt1, χt2) -> Float64

Effective stop frequency `Ωstop = min(1/6.26^{3/2}, Ωcrit(ν, χt1, χt2))` (`1PAT1ea.m`/
`1PAT1R.m`); evaluated once from the *initial* `(ν, χt1, χt2)` as upstream substitutes the ICs
into the stop condition.
"""
omega_crit(model::Union{Model1PAT1ea,Model1PAT1R}, ν::Real, χt1::Real, χt2::Real) =
    min(_OMEGA_FLOOR_1PA, model.ωcrit(ν, χt1, χt2))

# Numerator of the 1PAT1R Ω' ratio: dΩ/dt · dEdΩ (`1PAT1R.m` `dΩdttimesdEdΩ`, M = 1, x = 1/r0).
function _domega_times_dEdΩ(f::Forcing1PAT1ea, r0::Real, ν::Real, at1::Real, at2::Real, δm1::Real)
    x = 1 / r0
    dm_term = δm1 * (sqrt(1 / r0^3) * (-2 * r0^(5 / 2) / 3) * flux_1SF_deriv(f.base, r0))
    return -ν * flux_1SF(f, r0) +
           ν^2 * (-flux_2SF(f, r0) - dm_term - at1 * flux_a(f, r0) / ν -
                  at2 * flux_sigma(f, r0) / ν -
                  (2 * x^(5 / 2) * (-2 + 3x) * flux_LH(f, r0)) / (3 * (1 - 3x)^(3 / 2)) -
                  E_dm(x) * f.base.FEH(r0))
end

"""
    _spin_eom!(du, u, model::Union{Model1PAT1ea,Model1PAT1R})

In-place 7-variable spinning-secondary equation of motion for the state
`u = (Ω, φ, ν, M, χt1, χt2, δm)`; `ν`, `M`, `χt2` are constants of motion while the
primary-spin correction `χt1` and the secondary-mass correction `δm` evolve alongside
`(Ω, φ)`. The `Model1PAT1ea` method drives `dΩ/dt` with the 1PA forcing-channel sum
`ν F₀ + ν²(F₁ + F₁ₑₓₜᵣₐ + δm F₁δm) + ν χt1 F₁s1 + ν χt2 F₁s2`; the `Model1PAT1R` method
instead uses the resummed ratio `(dΩ/dt · dE/dΩ)/(dE/dΩ)` (`1PAT1R.m`). Wrapped by the
`Hybrid.Drivers` SciML rhs (`_spin_rhs!`) for the `SpinInspiralSolution` evolution.
"""
function _spin_eom!(du, u, model::Model1PAT1ea)
    Ω, _, ν, M, χt1, χt2, δm = u
    f = model.forcing
    r0 = Ω^(-2 / 3)
    du[1] = (1 / M) * (ν * forcing_F0(f, r0) +
                       ν^2 * (forcing_F1(f, r0) + forcing_F1extra(f, r0) + δm * forcing_F1dm(f, r0)) +
                       ν * χt1 * forcing_F1s1(f, r0) + ν * χt2 * forcing_F1s2(f, r0))
    du[2] = Ω / M
    du[3] = 0.0
    du[4] = 0.0
    du[5] = (ν^2 / M) * flux_LH(f, r0)
    du[6] = 0.0
    du[7] = (ν / M) * f.base.FEH(r0)
    return nothing
end

function _spin_eom!(du, u, model::Model1PAT1R)
    Ω, _, ν, M, χt1, χt2, δm = u
    f = model.forcing
    r0 = Ω^(-2 / 3)
    num = _domega_times_dEdΩ(f, r0, ν, χt1, χt2, δm)
    den = dEdΩ(f.base.z, r0, ν, χt1, χt2, δm)
    du[1] = (1 / M) * num / den
    du[2] = Ω / M
    du[3] = 0.0
    du[4] = 0.0
    du[5] = (1 / M) * ν^2 * flux_LH(f, r0)
    du[6] = 0.0
    du[7] = (ν / M) * f.base.FEH(r0)
    return nothing
end
