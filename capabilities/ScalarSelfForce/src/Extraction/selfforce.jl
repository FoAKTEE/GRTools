# Extraction/selfforce.jl — F_r^m (and F_t^m, F_φ^m, F_θ^m) from a converged m-mode solution:
# the number arXiv:2606.04998 is about (ESSF node `mode-self-force-extraction`, GLOBAL_DAG
# row 44; DESIGN §2.2 `Extraction/selfforce.jl`).
#
# WHAT THIS FILE OWNS — and nothing else:
#   1. `mode_self_force`            — Osburn–Nishimura Eq. (3.10) on raw ingredients, on a
#                                     `ModeSolution`, or on any level of an `AmrHistory`.
#   2. `mode_self_force_history`    — the per-AMR-level sequence (rows 50 / 53's input).
#   3. `mode_force_convergence`     — its successive-difference / ratio table and the
#                                     resolution floor (main.tex:705-724).
#   4. `regularized_mode_force`     — `F_r^m − F^m_{r[order]}`, the per-m term row 46 SUMS.
# The field and its gradient come from the SEALED row-42 read-out surface (`field_at`,
# `field_grad_at`) or from the SEALED row-43 per-level diagnostics (`AmrLevel.psi/dpsi_r/
# dpsi_w`); no stencil, no interpolation and no basis is re-derived here. THE MODE SUM
# (row 46), the process fan-out (row 47) and every figure (rows 48-53) are other nodes.
#
# SOURCES (read-only mirrors):
#   * Osburn & Nishimura, arXiv:2206.07031 `src/main.tex:178-184` — Eq. (3.10), the formula;
#     `:89-94` — their mode ansatz and `Δφ(r) = a/(r₊−r₋) ln((r−r₊)/(r−r₋))`.
#   * CRP arXiv:2606.04998 `main.tex:428-430` (`F^self_α = q ∇_α Φ^R|_{x_p}`),
#     `main.tex:701-703` ("we extract the radial self-force contributions F_r^m ... following
#     Eq. (3.10) in [Osburn_2022]"), `main.tex:719-724` (the resolution floor is the
#     difference between the highest and second-highest AMR level).
#   * spectre `src/Elliptic/Systems/SelfForce/Scalar/Events/ObserveSelfForce.cpp:52-77` — the
#     reference implementation; its regression rows
#     `tests/InputFiles/SelfForce/ScalarCircularOrbit.yaml:19-23` are gated in
#     test/test_selfforce.jl.

export mode_self_force, mode_self_force_history, mode_force_convergence,
       regularized_mode_force

# ============================== the formula ================================================

"""
    mode_self_force(orbit::CircularOrbit, m, psi, dpsi_r, dpsi_z = zero(psi); phi_p = 0)
    mode_self_force(sol::ModeSolution, m = sol.m; element = nothing, phi_p = sol.phi_p)
    mode_self_force(hist::AmrHistory, m = hist.m; level = lastindex(hist.levels), phi_p = …)

The `m`-mode contribution to the scalar self-force at the particle, Osburn–Nishimura
Eq. (3.10) (`arXiv:2206.07031` `main.tex:178-184`) in this package's locked conventions.
Returns a `NamedTuple`; `Fr` is `(M²/q²) F_r^m`, the number the Warburton–Barack tables print
and the one `fig:ssf_error_budget` sums (`q = 1` and `M = 1` package-wide).

# THE NORMATIVE CONVENTION STATEMENT for extraction (the ansatz side is locked in
# [`mode_term`](@ref), the subtraction side in [`mmode_reg_factor_over_pi`](@ref))

The self-force is `F_α = q ∇_α Φ^ℛ|_{x_p}` (CRP `main.tex:428-430`) and the residual field is
the locked mode sum of [`mode_term`](@ref),

    Φ^ℛ = Σ_{m∈ℤ} Ψ^ℛ_m(r, z) · sin^{|m|}θ / r · e^{+i m (φ̄ − Ω t)},   φ̄ = φ + χ(r),

with `χ = kerr_varphi_offset(a, r)` and `dχ/dr = a/Δ` (`kerr_varphi_offset_dr`). Folding the
conjugate pairs `Ψ^ℛ_{−m} = conj(Ψ^ℛ_m)` gives Eq. (3.10) exactly:

    F^{m=0}_α = q ∂_α (Ψ^ℛ_0 / r)                              — UNFOLDED, no ×2, no phase
    F^{m>0}_α = 2q Re[ e^{+i m (φ̄_p − Ω t_p)} ∂_α (Ψ^ℛ_m / r) ]

evaluated at `x_p = (t_p, r₀, π/2, φ_p)`. Written out at the equator, where `sin^{|m|}θ = 1`
and its θ-derivative `|m| sin^{|m|−1}θ cosθ` vanishes:

  * **radial** — `∂_r` at fixed BOYER–LINDQUIST `φ` hits three things, `Ψ^ℛ_m`, the `1/r`,
    AND the phase, because `φ̄` depends on `r`:

        F_r^m / (2q Re[e^{imφ̄_p} · … ]) :   ∂_rΨ^ℛ_m/r₀ − Ψ^ℛ_m/r₀² + i m (a/Δ) Ψ^ℛ_m/r₀

    The last term is the DEPHASING DERIVATIVE. Dropping it is an `O(a)` error (MEASURED 9.6 %
    at `a = 0.5`, `m = 2`); it is the term spectre adds by hand at
    `ObserveSelfForce.cpp:70-72`, with the sign fixed by the `+im φ̄` of the ansatz — the
    conjugate of the `−i a m/Δ` the sealed puncture applies to `Ψ^𝒫_m`
    ([`puncture_psi_chart`](@ref), `CircularOrbit.cpp:337-339`).
  * **temporal / azimuthal** — the phase gives `∂_t → −imΩ` and `∂_φ → +im`, so
    `F_t^m = −Ω F_φ^m` exactly, mode by mode (the circular-orbit identity), and both vanish
    at `m = 0`.
  * **polar** — `F_θ^m = −2q Re[e^{imφ̄_p}] ∂_zΨ^ℛ_m / r₀` through `∂_θ = −sinθ ∂_z`. In the
    production `z²` chart `∂_z = 2z ∂_{z²}` and `z = 0` at the equator, so `F_θ^m` is a
    BITWISE zero — the `= 0` spectre hard-codes at `ObserveSelfForce.cpp:49-51`. In the full
    `z` chart it is zero only up to truncation error and is returned as measured.

`m ≥ 0` always (retrograde orbits are `a < 0`, never `m < 0`). `phi_p` is the particle's
Boyer–Lindquist azimuth and `t_p = 0`; the CRP and spectre both use `φ_p = 0`, which makes the
phase `e^{+imχ(r₀)}` — spectre's `rotation` at `ObserveSelfForce.cpp:73-74`.

# FIELDS OF THE RESULT

`Fr`, `Ft`, `Ftheta`, `Fphi` — the real force components; `Fr_complex`, `Ft_complex`,
`Ftheta_complex`, `Fphi_complex` — the complex quantities BEFORE the real part is taken
(spectre's regression rows print `Re`/`Im` of `Fr_complex`, so this is what a comparison
against them uses); `m`, `phase` (`e^{im(φ_p+χ(r₀))}`), `psi`, `dpsi_r`, `dpsi_z` — the
ingredients; `element` — which one-sided trace produced them (`0` for the raw-ingredient
method). The `AmrHistory` method adds `step` and `ndofs`.

# WHICH FIELD, AND WHICH TRACE

`sol.u` IS `Ψ^ℛ_m` inside the worldtube (row 42's contract) and the particle is inside it, so
the extraction reads `sol.u` directly — never [`reconstruct_full_field`](@ref), which would
add `Ψ^𝒫_m` back. The puncture sits on an element CORNER, so the value is a one-sided DG
trace and `element` selects which; omitted, the first containing element in mesh order is used
— the same choice [`field_at`](@ref) and the sealed [`AmrLevel`](@ref) diagnostics make.
spectre instead AVERAGES over its `NumContributingElements` (`ObserveSelfForce.hpp:97-107`);
the two agree to the trace spread, MEASURED `~1e-11` in `F_r` at the paper's 8 AMR levels,
far below the resolution floor.

Type-generic: `BigFloat` orbits and ingredients give `BigFloat` components at the ambient
`precision`.
"""
function mode_self_force(orbit::CircularOrbit, m::Integer, psi::Number, dpsi_r::Number,
                         dpsi_z::Number = zero(psi); phi_p::Real = 0)
    m >= 0 || throw(ArgumentError(
        "mode_self_force: modes are folded (m ≥ 0; retrograde orbits are a < 0), got m = $m"))
    RT = float(promote_type(typeof(orbit.r0), real(typeof(psi)), real(typeof(dpsi_r)),
                            real(typeof(dpsi_z)), typeof(phi_p)))
    CT = Complex{RT}
    r0 = RT(orbit.r0)
    p, dr, dz = CT(psi), CT(dpsi_r), CT(dpsi_z)
    brad = dr / r0 - p / r0^2                    # ∂_r(Ψ^ℛ_m / r), phase excluded
    bpol = -dz / r0                              # ∂_θ(Ψ^ℛ_m / r) = −∂_z(...) at the equator
    if m == 0
        z = zero(CT)
        return _sf_result(m, one(CT), brad, z, z, bpol, p, dr, dz, 0)
    end
    ph = CT(cis(m * kerr_varphi(RT(orbit.a), r0, RT(phi_p))))
    chi1 = RT(kerr_varphi_offset_dr(RT(orbit.a), r0))              # χ′ = a/Δ (sealed)
    Frc = 2 * ph * (brad + im * m * chi1 * p / r0)
    Fphc = 2 * ph * (im * m) * p / r0
    Ftc = -RT(orbit.Omega) * Fphc                                  # ∂_t = −Ω ∂_φ on the phase
    Fthc = 2 * ph * bpol
    return _sf_result(m, ph, Frc, Ftc, Fphc, Fthc, p, dr, dz, 0)
end

_sf_result(m, phase, Frc, Ftc, Fphc, Fthc, psi, dpsi_r, dpsi_z, element) =
    (Fr = real(Frc), Ft = real(Ftc), Ftheta = real(Fthc), Fphi = real(Fphc),
     Fr_complex = Frc, Ft_complex = Ftc, Ftheta_complex = Fthc, Fphi_complex = Fphc,
     m = Int(m), phase = phase, psi = psi, dpsi_r = dpsi_r, dpsi_z = dpsi_z,
     element = Int(element))

# `∂_z Ψ` from the domain's own angular chart: `w = z²` under equatorial symmetry (so
# `∂_z = 2z ∂_w`, an EXACT zero at the equator) and `z = cosθ` without it (row 39's charts).
function _sf_dpsi_z(chart::Symbol, x2::Real, dpsi_x2::Number)
    chart === :z2 && return 2 * sqrt(x2) * dpsi_x2
    chart === :z && return dpsi_x2
    throw(ArgumentError("mode_self_force: unknown angular chart :$chart " *
                        "(expected :z2 or :z — see `block_lattice_bounds`)"))
end

function mode_self_force(sol::ModeSolution, m::Integer = sol.m; element = nothing,
                         phi_p::Real = sol.phi_p)
    Int(m) == sol.m || throw(ArgumentError(
        "mode_self_force: this ModeSolution solves m = $(sol.m), got m = $m — a different " *
        "mode is a different problem (a different RunConfig), never a keyword"))
    xp = puncture_chart_position(sol.dom, sol.orbit)
    e = element === nothing ? _dr_locate(sol.mesh, xp[1], xp[2]) : Int(element)
    psi = field_at(sol, xp[1], xp[2]; element = e)
    g = field_grad_at(sol, xp[1], xp[2]; element = e)
    r = mode_self_force(sol.orbit, sol.m, psi, g[1],
                        _sf_dpsi_z(sol.dom.chart, xp[2], g[2]); phi_p = phi_p)
    return merge(r, (element = e,))
end

function mode_self_force(hist::AmrHistory, m::Integer = hist.m;
                         level::Integer = lastindex(hist.levels),
                         phi_p::Real = hist.solution.phi_p)
    Int(m) == hist.m || throw(ArgumentError(
        "mode_self_force: this AmrHistory is the m = $(hist.m) sweep, got m = $m"))
    l = hist.levels[level]
    r = mode_self_force(hist.orbit, hist.m, l.psi, l.dpsi_r,
                        _sf_dpsi_z(hist.solution.dom.chart, hist.x_p[2], l.dpsi_w);
                        phi_p = phi_p)
    return merge(r, (element = l.element,), (step = l.step, ndofs = l.mesh.ndofs))
end

# ========================= the per-level sequence and its convergence ======================

"""
    mode_self_force_history(hist::AmrHistory; phi_p = hist.solution.phi_p) -> Vector

[`mode_self_force`](@ref) at EVERY level of an hp-AMR sweep, level 0 first — the sequence
`fig:ssf_convergence` plots and the error budget's resolution floor is measured from
(CRP `main.tex:719-724`). Each entry carries the extra `step` and `ndofs` fields.

The per-level ingredients are the SEALED [`AmrLevel`](@ref) diagnostics recorded by
[`solve_mode_amr`](@ref), so this costs nothing beyond the sweep itself and needs no stored
per-level solution vector.
"""
mode_self_force_history(hist::AmrHistory; phi_p::Real = hist.solution.phi_p) =
    [mode_self_force(hist; level = k, phi_p = phi_p) for k in eachindex(hist.levels)]

"""
    mode_force_convergence(hist::AmrHistory; component = :Fr, phi_p = …) -> NamedTuple

The successive-difference convergence table of one force component over an hp-AMR sweep:
`steps`, `ndofs`, `values`, `diffs` (`|Δ|` between consecutive levels), `ratios`
(consecutive `diffs`), `two_step_ratios` (`diffs[k+2]/diffs[k]`), and the two headline
numbers `resolution_error = diffs[end]` and `relative_resolution_error`.

`resolution_error` is EXACTLY the CRP's own "error due to numerical resolution ... the
difference between the highest and second-highest AMR level" (`main.tex:723`), which is what
`fig:ssf_error_budget` plots — this function is that measurement, on `F_r^m` rather than on a
raw field diagnostic. It is a DISCRETIZATION/round-off floor, NOT the GMRES tolerance: the
sealed row-43 node re-ran the sweep at `abstol = 1e-14` and the floor did not move, so a
caller must not attribute it to the linear solve or try to buy it down there. The DEFINITION
is the row-43 one — the last successive difference, relative denominator `|values[end]|` — so
this is the `F_r^m` analogue of that node's `resolution_error(hist; field)` helper on a raw
field diagnostic, and the two numbers are directly comparable.

THE PER-LEVEL VALUES ARE MESH-DEPENDENT, so every number quoted here is a MEASUREMENT of one
AMR sequence, not a contract: an `hp_refine` order-map change moves them slightly. Gate on
tolerances and on the shape of the sequence, never on a golden per-level value.

TWO-STEP, NOT ONE-STEP: the single-step ratio is NOT universally below 0.5 — MEASURED 1.313
at level 5 of the `a = 0.5, r₀ = 10, m = 1` sweep (the sealed row-43 node saw the same on its
raw diagnostics). The geometric statement that does hold is the two-step ratio (MEASURED max
0.238 there) and the envelope over the whole sweep (mean ratio 0.168). Iteration counts and
convergence rates are CONFIGURATION-SCOPED — quote them per `(a, r₀, m)`, never in general.
`component` must name a real field of [`mode_self_force`](@ref)'s result
(`:Fr`, `:Ft`, `:Ftheta`, `:Fphi`).
"""
function mode_force_convergence(hist::AmrHistory; component::Symbol = :Fr,
                                phi_p::Real = hist.solution.phi_p)
    component in (:Fr, :Ft, :Ftheta, :Fphi) || throw(ArgumentError(
        "mode_force_convergence: component must be one of (:Fr, :Ft, :Ftheta, :Fphi), " *
        "got :$component"))
    fs = mode_self_force_history(hist; phi_p = phi_p)
    v = [getproperty(f, component) for f in fs]
    d = [abs(v[k + 1] - v[k]) for k in 1:length(v)-1]
    rat = [iszero(d[k]) ? oftype(d[k], NaN) : d[k + 1] / d[k] for k in 1:length(d)-1]
    two = [iszero(d[k]) ? oftype(d[k], NaN) : d[k + 2] / d[k] for k in 1:length(d)-2]
    resolution = isempty(d) ? zero(eltype(v)) : d[end]
    return (steps = [f.step for f in fs], ndofs = [f.ndofs for f in fs], values = v,
            diffs = d, ratios = rat, two_step_ratios = two,
            resolution_error = resolution,
            relative_resolution_error = iszero(v[end]) ? oftype(resolution, NaN) :
                                        resolution / abs(v[end]))
end

# ============================ the regularized per-mode term ================================

# The term to subtract from an Osburn–Nishimura Eq. (3.10) force mode: the Heffernan parameter
# carries the ×2 real-part folding INSIDE it at every m (the sealed rule of
# `mmode_reg_factor_over_pi`), while Eq. (3.10)'s m = 0 branch is UNFOLDED — hence the halving.
_sf_reg_term(order::Integer, orbit::CircularOrbit, m::Integer) =
    m == 0 ? mmode_reg_Fr(order, orbit, 0) / 2 : mmode_reg_Fr(order, orbit, m)

"""
    regularized_mode_force(orbit::CircularOrbit, m, Fr::Real; order = 4) -> T
    regularized_mode_force(sol::ModeSolution, m = sol.m; order = 4, kwargs...) -> T
    regularized_mode_force(hist::AmrHistory, m = hist.m; order = 4, kwargs...) -> T

`F_r^m − F^m_{r[order]}`: the extracted radial force mode with the Heffernan large-m
regularization parameter subtracted ([`mmode_reg_Fr`](@ref); `order ∈ (2, 4, 6, 8)`, the CRP
uses 4, `main.tex:723`).

THE m = 0 HALVING RULE (the sealed contract of [`mmode_reg_factor_over_pi`](@ref)): the ×2
real-part folding is inside `F^m_{r[order]}` at every `m`, `m = 0` included, whereas
Osburn–Nishimura Eq. (3.10)'s `m = 0` branch is UNFOLDED — so the subtracted term is
`mmode_reg_Fr(order, orbit, 0)/2` at `m = 0` and `mmode_reg_Fr(order, orbit, m)` above it.
`regularized_mode_force` applies that rule; a caller must never re-apply it.

THIS IS THE PER-m TERM THE MODE SUM ADDS UP. Because the doubly-infinite parameter sum
vanishes (`F^0/2 + Σ_{m≥1} F^m = 0`, gated in `test_regularization.jl`), the accelerated
mode sum of `main.tex:723` is literally

    F_r ≈ Σ_{m=0}^{m_max} regularized_mode_force(…, m) ,

identical to `Σ_{m≤m_max} F_r^m + `[`mmode_reg_tail_Fr`](@ref)`(order, orbit, m_max+1)`.
Performing that sum is node `mode-sum-total-force` (row 46), not this file.

MEASURED at `a = 0.5`, `r₀ = 10`, 6 AMR levels: the residual falls like `m^{-6}` — log-log
slope `−6.13` over `m ∈ {2, 5, 10}` — against the raw `F_r^m`'s `m^{-4}`
(`main.tex:520-522`).
"""
regularized_mode_force(orbit::CircularOrbit, m::Integer, Fr::Real; order::Integer = 4) =
    Fr - _sf_reg_term(order, orbit, m)

regularized_mode_force(sol::ModeSolution, m::Integer = sol.m; order::Integer = 4, kwargs...) =
    regularized_mode_force(sol.orbit, sol.m, mode_self_force(sol, m; kwargs...).Fr;
                           order = order)

regularized_mode_force(hist::AmrHistory, m::Integer = hist.m; order::Integer = 4, kwargs...) =
    regularized_mode_force(hist.orbit, hist.m, mode_self_force(hist, m; kwargs...).Fr;
                           order = order)
