# WaSABI-C SF-PN hybrid: 9-variable forcing assembly + evolution equations + stop condition +
# δx-coupled trajectory (`InspiralModels/WaSABI-C.m`, §WF frequency forcing term … §Evolution
# equations). The nine forcing channels Fω/Fφ/Fν/Fm/Fχ1/Fχ2 are MACHINE-GENERATED from the WL
# InputForm (wasabic_codegen.jl); the three shift channels Fδm/Fδν/Fδχ are the upstream
# closed-form combinations of Fm/Fν/Fχ1.
#
# State u = (ω, φ, ν, M, χ1, χ2, δm, δν, δχ). Evolution (`WaSABI-C.m` L1511-1520, forcing
# evaluated at the literal mass m=1):
#   ω' = Fω/M, φ' = ω/M, ν' = Fν/M, M' = Fm, χ1' = Fχ1/M, χ2' = Fχ2/M,
#   δm' = Fδm/M, δν' = Fδν/M, δχ' = Fδχ/M.
# Stop: ω ≥ min(1/((1.05 rISCO[χ1])^{3/2}+χ1), 1/(6.26^{3/2}+χ1), ωcrit(ν,χ1,χ2)), the stop
# frequency fixed by the initial (ν,χ1,χ2) as upstream substitutes the ICs.

# --- nine forcing channels -------------------------------------------------

"""
    wc_Fφ(ω) = ω

WaSABI-C phase-evolution forcing `Fφ = ω` (`WaSABI-C.m` §WF phase evolution).
"""
wc_Fφ(ω::Real) = ω

"""
    wc_forcing(mdl, ω, φ, ν, m, χ1, χ2, δm, δν, δχ) -> NTuple{9}

The full WaSABI-C forcing channel vector `(Fω, Fφ, Fν, Fm, Fχ1, Fχ2, Fδm, Fδν, Fδχ)`
(`WaSABI-C.m`, the `ForcingTerms["WaSABI-C"]` order). `Fδm/Fδν/Fδχ` are the upstream shift
combinations of `Fm/Fν/Fχ1` (computed once and reused).
"""
function wc_forcing(mdl::WaSABICModel, ω, φ, ν, m, χ1, χ2, δm, δν, δχ)
    Fω = _wc_Fω(mdl, ω, φ, ν, m, χ1, χ2, δm, δν, δχ)
    Fν = _wc_Fν(mdl, ω, φ, ν, m, χ1, χ2, δm, δν, δχ)
    Fm = _wc_Fm(mdl, ω, φ, ν, m, χ1, χ2, δm, δν, δχ)
    Fχ1 = _wc_Fχ1(mdl, ω, φ, ν, m, χ1, χ2, δm, δν, δχ)
    Fχ2 = _wc_Fχ2(mdl, ω, φ, ν, m, χ1, χ2, δm, δν, δχ)
    Fδm = (Fm - δm * Fν) / ν                       # §Total mass shift
    Fδν = -((-1 + 2 * δν * ν) * Fν) / ν^2          # §Symmetric mass ratio shift
    Fδχ = (-δχ * Fν + Fχ1) / ν                     # §Primary spin shift
    return (Fω, ω, Fν, Fm, Fχ1, Fχ2, Fδm, Fδν, Fδχ)
end

# --- stop frequency (fixed by the initial conditions) ----------------------

"""
    omega_crit(mdl::WaSABICModel, ν, χ1, χ2) -> Float64

WaSABI-C effective stop frequency
`min(1/((1.05 rISCO[χ1])^{3/2}+χ1), 1/(6.26^{3/2}+χ1), ωcrit(ν,χ1,χ2))` (`WaSABI-C.m`
§Stop Condition). Integration terminates when `ω` reaches this value. Extends the `Models`
`omega_crit` generic to the WaSABI-C model.
"""
function omega_crit(mdl::WaSABICModel, ν::Real, χ1::Real, χ2::Real)
    ri = kerr_isco(χ1)
    return min(1 / ((1.05 * ri)^(3 / 2) + χ1), 1 / (6.26^(3 / 2) + χ1), mdl.ωcrit(ν, χ1, χ2))
end

# --- 9-var equations of motion ---------------------------------------------

# In-place EOM for u = (ω, φ, ν, M, χ1, χ2, δm, δν, δχ); forcing at the literal mass m=1.
function _wasabic_eom!(du, u, mdl::WaSABICModel)
    ω, φ, ν, M, χ1, χ2, δm, δν, δχ = u
    Fω, Fφ, Fν, Fm, Fχ1, Fχ2, Fδm, Fδν, Fδχ =
        wc_forcing(mdl, ω, φ, ν, 1, χ1, χ2, δm, δν, δχ)
    du[1] = Fω / M
    du[2] = Fφ / M
    du[3] = Fν / M
    du[4] = Fm
    du[5] = Fχ1 / M
    du[6] = Fχ2 / M
    du[7] = Fδm / M
    du[8] = Fδν / M
    du[9] = Fδχ / M
    return nothing
end

_wasabic_rhs!(du, u, p, t) = _wasabic_eom!(du, u, p)

# --- driver + solution -----------------------------------------------------

"""
    WaSABICSolution

Result of the 9-var WaSABI-C [`int_inspiral`](@ref): the integrated ODE solution for
`u = (ω, φ, ν, M, χ1, χ2, δm, δν, δχ)`. Query with `omega_at`, `phi_at`, `duration`, and
[`wc_state_at`](@ref).
"""
struct WaSABICSolution{S}
    sol::S
    mdl::WaSABICModel
end

"""
    int_inspiral(mdl::WaSABICModel, ics::NamedTuple; reltol=1e-10, abstol=1e-10, tmax=1e12)
        -> WaSABICSolution

Integrate the 9-var WaSABI-C inspiral from `ics` (`(; ω, ν, χ1, χ2, M, φ, δm, δν, δχ)`; `φ` and
the three shifts default to 0, `M` to 1), terminating when `ω` reaches
[`omega_crit`](@ref)`(mdl, ν, χ1, χ2)`. Adds a `WaSABICModel` method to the `Drivers.int_inspiral`
generic. `reltol`/`abstol` map to upstream `PrecisionGoal`/`AccuracyGoal`.
"""
function int_inspiral(mdl::WaSABICModel, ics::NamedTuple;
                      reltol::Real = 1e-10, abstol::Real = 1e-10, tmax::Real = 1e12)
    g(k, d) = haskey(ics, k) ? float(getfield(ics, k)) : d
    u0 = [float(ics.ω), g(:φ, 0.0), float(ics.ν), g(:M, 1.0),
          float(ics.χ1), float(ics.χ2), g(:δm, 0.0), g(:δν, 0.0), g(:δχ, 0.0)]
    ωstop = omega_crit(mdl, ics.ν, ics.χ1, ics.χ2)
    cond(u, t, integrator) = ωstop - u[1]          # ω increases; crosses ωstop → terminate.
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(cond, affect!)
    prob = ODEProblem(_wasabic_rhs!, u0, (0.0, float(tmax)), mdl)
    sol = solve(prob, Vern9(); reltol = reltol, abstol = abstol, callback = cb)
    return WaSABICSolution(sol, mdl)
end

"""
    duration(s::WaSABICSolution) -> Float64

Inspiral duration (stop time) of the WaSABI-C trajectory.
"""
duration(s::WaSABICSolution) = s.sol.t[end]

"""
    omega_at(s::WaSABICSolution, t) -> Float64
    phi_at(s::WaSABICSolution, t) -> Float64

Frequency `ω(t)` / phase `φ(t)` along the WaSABI-C trajectory.
"""
omega_at(s::WaSABICSolution, t::Real) = s.sol(t)[1]
phi_at(s::WaSABICSolution, t::Real) = s.sol(t)[2]

"""
    wc_state_at(s::WaSABICSolution, t) -> NamedTuple

Full 9-var state `(ω, φ, ν, M, χ1, χ2, δm, δν, δχ)` along the trajectory at time `t`.
"""
function wc_state_at(s::WaSABICSolution, t::Real)
    u = s.sol(t)
    return (ω = u[1], φ = u[2], ν = u[3], M = u[4], χ1 = u[5],
            χ2 = u[6], δm = u[7], δν = u[8], δχ = u[9])
end
