"""
    EMRIWaveforms.Trajectory

Flux-driven inspiral ODE evolution (design §6.8: FEW `trajectory/{emri_inspiral,ode_base,pn5,pn5_rhs}.jl`;
KerrBackwards `trajectory.jl`, `pn5fluxes.jl`). The 5PN-e10 flux coefficient table
(`Pn5Fluxes.dEdt8H_5PNe10` / `dLdt8H_5PNe10`) feeds `pn5_rhs`, integrated by the vendored
`NumericsCore` DOPR853 over the shared `ode_rhs` with separatrix-refined termination
(`Spacetime.get_separatrix`).

State `[ṗ, ė, Ẏ, Ω_φ, Ω_θ, Ω_r]`. The equatorial (`Y = 1`) `pn5_rhs` branch maps the `(Ė, L̇z)` 5PN
fluxes through `Mappings.ELdot_to_PEdot_Jacobian`; the generic-inclination (`Y ≠ 1`) branch drives the
direct `(ṗ, ė, Ẏ)` 5PN-e10 fluxes (`Pn5Fluxes.dpdt/dedt/dYdt8H_5PNe10`) with the FEW `Ω_θ` sign
convention. Kerr-only: `a = 0` (Schwarzschild) uses a separate FEW ODE upstream, so `pn5_rhs(0, …)` is
out of domain (the flux table's removable 0/0). [`run_inspiral`](@ref) takes a `backend::Symbol` selector
(`:few` default, `:kb`): both canonicalize on the single adaptive DOPR853 driver of the identical 5PN ODE
(the AK/KB Verner integrator and the FEW DOPR853 are the SAME adaptive-integration decision).
"""
module Trajectory

include("pn5_fluxes.jl")
using .Pn5Fluxes: dEdt8H_5PNe10, dLdt8H_5PNe10, dpdt8H_5PNe10, dedt8H_5PNe10, dYdt8H_5PNe10

import NumericsCore
using NumericsCore: MTSUN_SI, YRSID_SI
using Spacetime: get_separatrix, kerr_frequencies
using ..Mappings: Y_to_xI, ELdot_to_PEdot_Jacobian
import Roots

export dEdt8H_5PNe10, dLdt8H_5PNe10, dpdt8H_5PNe10, dedt8H_5PNe10, dYdt8H_5PNe10,
       pn5_rhs, ODEModel, PN5FluxODE, EMRIInspiral, run_inspiral

function _symmetric_mass_ratio(m1::Real, m2::Real)
    m1f = Float64(m1); m2f = Float64(m2)
    (isfinite(m1f) && m1f > 0 && isfinite(m2f) && m2f > 0) ||
        throw(ArgumentError("m1 and m2 must be finite, representable, and positive, got ($m1, $m2)"))
    q = min(m1f, m2f) / max(m1f, m2f)
    η = q / (1 + q)^2
    (isfinite(η) && η > 0) ||
        throw(ArgumentError("the symmetric mass ratio is not representable for masses ($m1, $m2)"))
    return η
end

# =====================================================================================
# pn5_rhs — the analytic 5PN RHS [ṗ, ė, Ẏ, Ω_φ, Ω_θ, Ω_r] (FEW trajectory/pn5_rhs.jl).
# =====================================================================================

"""
    pn5_rhs(a, p, e, Y) -> Vector{Float64}

`[ṗ, ė, Ẏ, Ω_φ, Ω_θ, Ω_r]` for the analytic 5PN (PN5) EMRI trajectory. `Y == 1` maps the equatorial
`(Ė, L̇z)` 5PN fluxes (orders (10,10)) through `ELdot_to_PEdot_Jacobian(|a|, p, e, sign(a), …)` with
`Ẏ = 0`. The generic-inclination `Y ≠ 1` branch drives the KerrBackwards direct 5PN fluxes
`dpdt`/`dedt`/`dYdt8H_5PNe10` at orders (10,10)/(10,8)/(7,10) (`Ẏ = 0` recovered at `Y = ±1`). The
frequencies are `kerr_frequencies` at `xI = Y_to_xI(a, p, e, Y)`; the generic branch restores the FEW
`Ω_θ` sign via `copysign(|Ω_θ|, sign(xI))` (BHPToolkit returns the magnitude, but `Υ_θ ∝ L_z`), so
`Ω_θ < 0` for retrograde orbits. (`a = 0` is out of domain.)
"""
function pn5_rhs(a, p, e, Y)
    iszero(a) && throw(DomainError(a, "pn5_rhs: a=0 requires the Schwarzschild flux model"))
    xI = Y_to_xI(a, p, e, Y)
    f = kerr_frequencies(a, p, e, xI)
    if Y == 1
        Edot = dEdt8H_5PNe10(a, p, e, Y, 10, 10)
        # The flux table uses the original signed-spin frame (Lz > 0 for Y = 1).
        # The Jacobian is evaluated in the equivalent a >= 0 frame, where a retrograde
        # orbit has Lz < 0, so reflect Lzdot together with Lz when a < 0.
        Lzdot = sign(a) * dLdt8H_5PNe10(a, p, e, Y, 10, 10)
        pdot, edot = ELdot_to_PEdot_Jacobian(abs(a), p, e, sign(a), Edot, Lzdot)
        Ydot = 0.0
        Omega_theta = f.Omegath
    else
        pdot = dpdt8H_5PNe10(a, p, e, Y, 10, 10)
        edot = dedt8H_5PNe10(a, p, e, Y, 10, 8)
        Ydot = dYdt8H_5PNe10(a, p, e, Y, 7, 10)
        Omega_theta = copysign(f.Omegath, sign(xI))
    end
    return [pdot, edot, Ydot, f.Omegaph, Omega_theta, f.Omegar]
end

# =====================================================================================
# ODE model interface + shared RHS driver (FEW trajectory/ode_base.jl).
# =====================================================================================

"Abstract trajectory ODE model (`few.trajectory.ode.base.ODEBase`)."
abstract type ODEModel end

ode_spin(model::ODEModel) = model.a
converts_Y(::ODEModel) = false
function evaluate_rhs end
separatrix_buffer_dist(::ODEModel) = 0.1
_separatrix_frame(::ODEModel, a, x) = (a, x)

"""
    ode_rhs(model, y) -> Vector{Float64}

Trajectory RHS at state `y = [c1, c2, c3, Φ_φ, Φ_θ, Φ_r]` (FEW `ODEBase.__call__`): out of bounds
(`e < 0` or `p ≤ p_sep`) ⇒ all-NaN; in bounds ⇒ `evaluate_rhs(model, p, e, x_or_Y, p_sep)`. PN5's
`evaluate_rhs` already returns pex fluxes, so no ELQ Jacobian here.
"""
function ode_rhs(model::ODEModel, y)
    a = ode_spin(model)
    p, e, x_or_Y = y[1], y[2], y[3]
    in_bounds = e >= 0
    p_sep = NaN
    x = x_or_Y
    if in_bounds
        converts_Y(model) && (x = Y_to_xI(a, p, e, x_or_Y))
        asep, xsep = _separatrix_frame(model, a, x)
        p_sep = get_separatrix(asep, e, xsep)
        in_bounds = p > p_sep
    end
    out = zeros(Float64, 6)
    if in_bounds
        out .= evaluate_rhs(model, p, e, x_or_Y, p_sep)
    else
        fill!(out, NaN)
    end
    return out
end

"""
    PN5FluxODE(a)

Analytic 5PN-e10 Kerr inspiral model (`few.trajectory.ode.pn5.PN5`) at spin `a`. Integrates `Y = cos ι`
(`converts_Y = true`); `evaluate_rhs` delegates to [`pn5_rhs`].
"""
struct PN5FluxODE <: ODEModel
    a::Float64
end
ode_spin(m::PN5FluxODE) = m.a
converts_Y(::PN5FluxODE) = true
_separatrix_frame(::PN5FluxODE, a, x) = (abs(a), sign(a) * x)
evaluate_rhs(model::PN5FluxODE, p, e, Y, p_sep) = pn5_rhs(model.a, p, e, Y)

# =====================================================================================
# Adaptive forward inspiral driver — canonical DOPR853 (FEW trajectory/emri_inspiral.jl).
# =====================================================================================

"""
    EMRIInspiral(model; err=1e-11, enforce_schwarz_sep=false)

Adaptive forward inspiral integrator wrapping an [`ODEModel`] (`few.trajectory.inspiral.EMRIInspiral`).
`err` is the DOPR853 absolute tolerance; `enforce_schwarz_sep` forces `p_sep = 6 + 2e`. (The ELQ
integration path is not carried in v0.1.)
"""
struct EMRIInspiral{M<:ODEModel}
    model::M
    err::Float64
    enforce_schwarz_sep::Bool
    function EMRIInspiral(model::M, err::Float64,
                          enforce_schwarz_sep::Bool) where {M<:ODEModel}
        (isfinite(err) && err > 0) ||
            throw(ArgumentError("err must be finite and positive, got $err"))
        new{M}(model, err, enforce_schwarz_sep)
    end
end
EMRIInspiral(model::ODEModel; err::Real = 1e-11, enforce_schwarz_sep::Bool = false) =
    EMRIInspiral(model, Float64(err), enforce_schwarz_sep)

function _insp_psep(insp::EMRIInspiral, p, e, x_or_Y)
    a = ode_spin(insp.model)
    (insp.enforce_schwarz_sep || a == 0.0) && return 6.0 + 2.0 * e
    x = converts_Y(insp.model) ? Y_to_xI(a, p, e, x_or_Y) : x_or_Y
    asep, xsep = _separatrix_frame(insp.model, a, x)
    return get_separatrix(asep, e, xsep)
end

"""
    run_inspiral(insp::EMRIInspiral, m1, m2, a, y0; backend=:few, T=1.0, dt=10.0, max_iter=100_000) -> NamedTuple

Integrate the EMRI inspiral from `y0 = [p, e, x, 0, 0, 0]` for masses `m1`, `m2` (M_⊙) and spin `a` to
the separatrix buffer or to time `T` (years). Adimensionalised exactly as `few.trajectory.integrate`
(`massratio = m1 m2/(m1+m2)²`, `tmax = T·YRSID/(M·MTSUN)·massratio`, phases ÷ massratio); the forward
stop is placed on the buffer by a dense-output Brent root-find. Returns dimensionful `t` (s), `p`, `e`,
`x`, `Phi_phi`, `Phi_theta`, `Phi_r`.

The call-site `a` must equal the spin carried by `insp.model`, and `max_iter` must be positive. If the
step budget is exhausted (or the integrator produces a non-finite rejected step) before reaching the
requested time or a terminal dense-output segment, the function throws instead of returning a silent
partial trajectory.

`backend::Symbol` (`:few` default, `:kb`) selects the adaptive integrator: both map to the single
canonical DOPR853 driver of the identical 5PN ODE (the AK/KB Verner integrator and the FEW DOPR853 are
the SAME adaptive-integration decision), so the two backends are bit-identical.
"""
function run_inspiral(insp::EMRIInspiral, m1::Real, m2::Real, a::Real, y0;
                      backend::Symbol = :few, T::Real = 1.0, dt::Real = 10.0,
                      max_iter::Integer = 100_000)
    backend in (:few, :kb) ||
        throw(ArgumentError("backend must be :few or :kb, got $backend"))
    model = insp.model
    model_spin = Float64(ode_spin(model))
    call_spin = Float64(a)
    (isfinite(call_spin) && isfinite(model_spin) && call_spin == model_spin) ||
        throw(ArgumentError("call-site spin a=$a does not match model spin $(ode_spin(model))"))
    max_iter > 0 || throw(ArgumentError("max_iter must be positive, got $max_iter"))
    massratio = _symmetric_mass_ratio(m1, m2)
    (isfinite(T) && T >= 0) ||
        throw(ArgumentError("T must be finite and nonnegative, got $T"))
    (isfinite(dt) && dt > 0) ||
        throw(ArgumentError("dt must be finite and positive, got $dt"))
    length(y0) == 6 || throw(DimensionMismatch("y0 must contain [p,e,x/Y,Φφ,Φθ,Φr] (6 values)"))
    all(isfinite, y0) || throw(ArgumentError("all y0 values must be finite"))
    buffer = separatrix_buffer_dist(model)
    Mtot = Float64(m1) + Float64(m2)
    (isfinite(Mtot) && Mtot > 0) ||
        throw(ArgumentError("total mass must be finite and representable, got m1+m2=$Mtot"))
    tmax = T * YRSID_SI / (Mtot * MTSUN_SI) * massratio
    Msec = MTSUN_SI * Mtot / massratio
    (isfinite(tmax) && tmax >= 0 && isfinite(Msec) && Msec > 0) ||
        throw(ArgumentError("derived time/mass scales are not finite for the supplied masses and T"))

    ode! = function (dydx, x, yv, args)
        dydx .= ode_rhs(model, yv)
        return nothing
    end
    solver = NumericsCore.DOPR853(ode!; abstol = insp.err)
    solver.errOld = 1e-4
    solver.previousReject = false

    y = collect(float.(y0))
    p0, e0, x0 = y[1], y[2], y[3]
    e0 >= 0 || throw(DomainError(e0, "initial eccentricity must be nonnegative"))
    psep0 = _insp_psep(insp, p0, e0, x0)
    p0 > psep0 + buffer ||
        throw(DomainError(p0, "initial p must exceed the separatrix buffer p_sep+$buffer = $(psep0 + buffer)"))
    n = length(y)
    ts = Float64[0.0]
    ys = Vector{Vector{Float64}}([copy(y)])

    h = min(dt / (Mtot * MTSUN_SI) * massratio, tmax / 10)
    x = 0.0
    niter = 0
    bad = 0
    have_seg = false
    seg_xa = 0.0
    seg_xb = 0.0
    seg_cont = Matrix{Float64}(undef, 8, n)
    failure_reason = nothing

    while x < tmax && niter < max_iter
        success, x1, y1, hnew, cont = NumericsCore.take_step(solver, x, y, h)
        if !success
            if !isfinite(hnew) || any(!isfinite, y1)
                failure_reason = "the adaptive stepper produced a non-finite rejected step"
                break
            end
            h = hnew
            bad += 1
            bad > 1000 && error("run_inspiral: exceeded 1000 consecutive rejected steps")
            continue
        end
        if any(!isfinite, y1)
            error("run_inspiral: adaptive stepper accepted a non-finite state")
        end
        bad = 0

        p, e, x_or_Y = y1[1], y1[2], y1[3]
        p_sep = _insp_psep(insp, p, e, x_or_Y)
        if (p - p_sep) < buffer
            seg_xa = x; seg_xb = x1; seg_cont = cont; have_seg = true
            break
        end
        if x1 > tmax
            seg_xa = x; seg_xb = x1; seg_cont = cont; have_seg = true
            break
        end

        push!(ts, x1)
        push!(ys, copy(y1))
        x = x1
        y .= y1
        h = hnew
        niter += 1
    end

    if !have_seg && x < tmax
        reason = failure_reason === nothing ?
            "max_iter=$max_iter was exhausted after $niter accepted steps" : failure_reason
        error("run_inspiral: integration terminated before the requested endpoint: $reason")
    end

    if have_seg
        seg_ts = [seg_xa, seg_xb]
        seg_co = [seg_cont]
        fdist(tt) = begin
            yv = NumericsCore.dense_eval(seg_ts, seg_co, tt)
            p, e, x_or_Y = yv[1], yv[2], yv[3]
            p_sep = _insp_psep(insp, p, e, x_or_Y)
            return p - (p_sep + buffer)
        end
        if seg_xb > tmax && fdist(tmax) > 0
            push!(ts, tmax)
            push!(ys, NumericsCore.dense_eval(seg_ts, seg_co, tmax))
        else
            t_star = Roots.find_zero(fdist, (seg_xa, seg_xb), Roots.Brent(); atol = 1e-12)
            push!(ts, t_star)
            push!(ys, NumericsCore.dense_eval(seg_ts, seg_co, t_star))
        end
    end

    npts = length(ts)
    a_spin = ode_spin(model)
    p_out = Vector{Float64}(undef, npts)
    e_out = Vector{Float64}(undef, npts)
    x_out = Vector{Float64}(undef, npts)
    Phi_phi = Vector{Float64}(undef, npts)
    Phi_theta = Vector{Float64}(undef, npts)
    Phi_r = Vector{Float64}(undef, npts)
    for i in 1:npts
        p, e, x_or_Y = ys[i][1], ys[i][2], ys[i][3]
        p_out[i] = p
        e_out[i] = e
        x_out[i] = converts_Y(model) ? Y_to_xI(a_spin, p, e, x_or_Y) : x_or_Y
        Phi_phi[i] = ys[i][4] / massratio
        Phi_theta[i] = ys[i][5] / massratio
        Phi_r[i] = ys[i][6] / massratio
    end

    return (t = ts .* Msec, p = p_out, e = e_out, x = x_out,
            Phi_phi = Phi_phi, Phi_theta = Phi_theta, Phi_r = Phi_r)
end

end # module Trajectory
