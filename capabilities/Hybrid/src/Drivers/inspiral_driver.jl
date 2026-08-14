# Inspiral ODE driver (`Kernel/Inspiral.m` `IntInspiral`): integrate the forcing rhs to a
# terminating event. The upstream uses `NDSolveValue[..., PrecisionGoal->prec,
# AccuracyGoal->acc]` with a `WhenEvent[StopCondition, "StopIntegration"]`; here the modern
# `OrdinaryDiffEqVerner` `solve(..., Vern9(); reltol, abstol)` with a `ContinuousCallback`
# zero-crossing (swapped off the monolithic `OrdinaryDiffEq` of the legacy package — Vern9 is the
# only integrator WaSABI ever used). Three state shapes share this driver: the 2-var (Ω, φ)
# Schwarzschild inspiral, the 4-arg (ω, φ; χ1) Kerr inspiral, and the 7-var
# (Ω, φ, ν, M, χt1, χt2, δm) spin/δm inspiral.

"""
    InspiralSolution

Result of [`int_inspiral`](@ref): the integrated ODE solution for the state `u = (Ω, φ)`, the
forcing `model`, and the mass parameters `ν`, `M` (`χ1` for the Kerr `(ω, φ)` state, else `NaN`).
Query with [`omega_at`](@ref), [`r0_at`](@ref), [`phi_at`](@ref) and [`duration`](@ref).
"""
struct InspiralSolution{S,Mdl}
    sol::S
    model::Mdl
    ν::Float64
    M::Float64
    χ1::Float64
end

# Schwarzschild models carry no spin; default χ1 = NaN keeps the 4-arg constructor calls.
InspiralSolution(sol, model, ν::Float64, M::Float64) = InspiralSolution(sol, model, ν, M, NaN)

"""
    int_inspiral(model, r0_init, ν, M; φ0=0.0, reltol=1e-10, abstol=1e-10, tmax=1e12)
        -> InspiralSolution

Integrate the 2-var (Ω, φ) inspiral of `model` (0PASchwarz / 1PAT1) from initial radius
`r0_init` (so `Ω0 = r0_init^{-3/2}`) and phase `φ0`, terminating when `r0 = Ω^{-2/3}` reaches
`stop_r0(model, ν)`. `reltol`/`abstol` map to the upstream `PrecisionGoal`/`AccuracyGoal`.
"""
function int_inspiral(model, r0_init::Real, ν::Real, M::Real;
                      φ0::Real = 0.0, reltol::Real = 1e-10, abstol::Real = 1e-10,
                      tmax::Real = 1e12)
    Ω0 = float(r0_init)^(-3 / 2)
    sc = stop_r0(model, ν)
    # r0 - stop_r0 crosses zero from above (r0 decreasing) → terminate.
    cond(u, t, integrator) = u[1]^(-2 / 3) - sc
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(cond, affect!)
    prob = ODEProblem(inspiral_rhs!, [Ω0, float(φ0)], (0.0, float(tmax)),
                      (model, float(ν), float(M)))
    sol = solve(prob, Vern9(); reltol = reltol, abstol = abstol, callback = cb)
    return InspiralSolution(sol, model, float(ν), float(M))
end

"""
    int_inspiral(model::Forcing0PAKerr, ω0, ν, M, χ1; φ0=0.0, reltol=1e-10, abstol=1e-10, tmax=1e12)
        -> InspiralSolution

Integrate the Kerr (ω, φ) inspiral from initial frequency `ω0` and phase `φ0`, terminating when
`ω` reaches [`omega_crit`](@ref). Same driver pattern as the Schwarzschild variant, but with the
Kerr state and the `ω ≥ Ωcrit` stop (`0PAKerr.m`).
"""
function int_inspiral(model::Forcing0PAKerr, ω0::Real, ν::Real, M::Real, χ1::Real;
                      φ0::Real = 0.0, reltol::Real = 1e-10, abstol::Real = 1e-10,
                      tmax::Real = 1e12)
    ωc = omega_crit(model, χ1)
    cond(u, t, integrator) = u[1] - ωc           # ω increases; crosses ωc → terminate.
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(cond, affect!)
    prob = ODEProblem(kerr_inspiral_rhs!, [float(ω0), float(φ0)], (0.0, float(tmax)),
                      (model, float(ν), float(M), float(χ1)))
    sol = solve(prob, Vern9(); reltol = reltol, abstol = abstol, callback = cb)
    return InspiralSolution(sol, model, float(ν), float(M), float(χ1))
end

"""
    duration(s::InspiralSolution) -> Float64

Inspiral duration: the terminating time (`Max[..."Domain"]` upstream).
"""
duration(s::InspiralSolution) = s.sol.t[end]

"""
    omega_at(s::InspiralSolution, t) -> Float64

Orbital frequency `Ω(t)` (Kerr: `ω(t)`) along the trajectory.
"""
omega_at(s::InspiralSolution, t::Real) = s.sol(t)[1]

"""
    r0_at(s::InspiralSolution, t) -> Float64

Orbital radius `r0(t) = Ω(t)^{-2/3}` along the trajectory.
"""
r0_at(s::InspiralSolution, t::Real) = s.sol(t)[1]^(-2 / 3)

"""
    phi_at(s::InspiralSolution, t) -> Float64

Orbital phase `φ(t)` along the trajectory.
"""
phi_at(s::InspiralSolution, t::Real) = s.sol(t)[2]

"""
    int_inspiral_batch(model, ics; threaded=false, kwargs...) -> Vector{InspiralSolution}

Integrate a batch of 2-var inspirals, one per `(r0_init, ν, M)` tuple in `ics`. Each solve is
independent (the model interpolants are read-only shared state), so the `threaded` path is
numerically identical to the serial path — assert that parity before claiming any speedup.
"""
function int_inspiral_batch(model, ics; threaded::Bool = false, kwargs...)
    out = Vector{InspiralSolution}(undef, length(ics))
    if threaded
        # §8 BLAS guard: pin OpenBLAS to 1 thread inside the @threads region (deadlock-safe,
        # bit-deterministic). Launch policy: keep Julia nthreads ≤ 20 CPUs.
        blas_single_thread() do
            Threads.@threads for i in eachindex(ics)
                r0i, νi, Mi = ics[i]
                out[i] = int_inspiral(model, r0i, νi, Mi; kwargs...)
            end
        end
    else
        for i in eachindex(ics)
            r0i, νi, Mi = ics[i]
            out[i] = int_inspiral(model, r0i, νi, Mi; kwargs...)
        end
    end
    return out
end

# --- 7-var spin/δm driver (1PAT1ea / 1PAT1R) -------------------------------
#
# The spin models integrate the full upstream 7-var state (Ω, φ, ν, M, χt1, χt2, δm) with the
# additive (1PAT1ea) or ratio (1PAT1R) `Ω'` and the dual frequency stop
# `Ω ≥ min(1/6.26^{3/2}, Ωcrit(ν,χt1,χt2))`, where Ωcrit is fixed by the *initial*
# (ν, χt1, χt2) (upstream substitutes the ICs).

const _SpinModel = Union{Model1PAT1ea,Model1PAT1R}

"""
    SpinInspiralSolution

Result of the 7-var [`int_inspiral`](@ref) for a [`Model1PAT1ea`](@ref)/[`Model1PAT1R`](@ref):
the integrated ODE solution for `u = (Ω, φ, ν, M, χt1, χt2, δm)`. Query with [`omega_at`](@ref),
[`r0_at`](@ref), [`phi_at`](@ref), [`chit1_at`](@ref), [`dm_at`](@ref), [`duration`](@ref).
"""
struct SpinInspiralSolution{S,Mdl<:_SpinModel}
    sol::S
    model::Mdl
end

# The spin ODE rhs; `p` is the model, dispatched to the per-variant EOM in Models.
_spin_rhs!(du, u, p, t) = _spin_eom!(du, u, p)

"""
    int_inspiral(model::Union{Model1PAT1ea,Model1PAT1R}, ics::NamedTuple;
                 reltol=1e-10, abstol=1e-10, tmax=1e12) -> SpinInspiralSolution

Integrate the 7-var spin/δm inspiral from initial conditions `ics`
(`(; Ω, φ, ν, M, χt1, χt2, δm)`; `φ` defaults to 0), terminating when `Ω` reaches
`min(1/6.26^{3/2}, Ωcrit(ν, χt1, χt2))`.
"""
function int_inspiral(model::_SpinModel, ics::NamedTuple;
                      reltol::Real = 1e-10, abstol::Real = 1e-10, tmax::Real = 1e12)
    φ0 = haskey(ics, :φ) ? float(ics.φ) : 0.0
    u0 = [float(ics.Ω), φ0, float(ics.ν), float(ics.M),
          float(ics.χt1), float(ics.χt2), float(ics.δm)]
    Ωstop = omega_crit(model, ics.ν, ics.χt1, ics.χt2)
    # Ω increases; (Ωstop - Ω) crosses zero from above → terminate.
    cond(u, t, integrator) = Ωstop - u[1]
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(cond, affect!)
    prob = ODEProblem(_spin_rhs!, u0, (0.0, float(tmax)), model)
    sol = solve(prob, Vern9(); reltol = reltol, abstol = abstol, callback = cb)
    return SpinInspiralSolution(sol, model)
end

duration(s::SpinInspiralSolution) = s.sol.t[end]

"""
    omega_at(s::SpinInspiralSolution, t) -> Float64

Orbital frequency `Ω(t)` along the 7-var trajectory.
"""
omega_at(s::SpinInspiralSolution, t::Real) = s.sol(t)[1]

"""
    r0_at(s::SpinInspiralSolution, t) -> Float64

Orbital radius `r0(t) = Ω(t)^{-2/3}` along the 7-var trajectory.
"""
r0_at(s::SpinInspiralSolution, t::Real) = s.sol(t)[1]^(-2 / 3)

"""
    phi_at(s::SpinInspiralSolution, t) -> Float64

Orbital phase `φ(t)` along the 7-var trajectory.
"""
phi_at(s::SpinInspiralSolution, t::Real) = s.sol(t)[2]

"""
    chit1_at(s::SpinInspiralSolution, t) -> Float64

Primary reduced spin `χt1(t)` along the 7-var trajectory.
"""
chit1_at(s::SpinInspiralSolution, t::Real) = s.sol(t)[5]

"""
    dm_at(s::SpinInspiralSolution, t) -> Float64

Dimensionless secondary-mass change `δm(t)` along the 7-var trajectory.
"""
dm_at(s::SpinInspiralSolution, t::Real) = s.sol(t)[7]

"""
    int_inspiral_batch(model::Union{Model1PAT1ea,Model1PAT1R}, ics; threaded=false, kwargs...)
        -> Vector{SpinInspiralSolution}

Integrate a batch of 7-var spin inspirals, one per `NamedTuple` IC in `ics`. Independent solves
⇒ the `threaded` path is bit-identical to the serial path (assert parity before any speedup).
"""
function int_inspiral_batch(model::_SpinModel, ics; threaded::Bool = false, kwargs...)
    out = Vector{SpinInspiralSolution}(undef, length(ics))
    if threaded
        # §8 BLAS guard: pin OpenBLAS to 1 thread inside the @threads region (deadlock-safe,
        # bit-deterministic). Launch policy: keep Julia nthreads ≤ 20 CPUs.
        blas_single_thread() do
            Threads.@threads for i in eachindex(ics)
                out[i] = int_inspiral(model, ics[i]; kwargs...)
            end
        end
    else
        for i in eachindex(ics)
            out[i] = int_inspiral(model, ics[i]; kwargs...)
        end
    end
    return out
end
