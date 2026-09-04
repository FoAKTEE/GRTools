# Extraction/selfforce.jl — F_r^m (and F_t^m, F_φ^m, F_θ^m) from a converged m-mode solution:
# the number arXiv:2606.04998 is about (ESSF node `mode-self-force-extraction`, GLOBAL_DAG
# row 44; DESIGN §2.2 `Extraction/selfforce.jl`).
#
# WHAT THIS FILE OWNS — and nothing else:
#   1. `mode_self_force`            — Osburn–Nishimura Eq. (3.10) on raw ingredients, on a
#                                     `ModeSolution`, on any level of an `AmrHistory`, or on a
#                                     bare `AmrLevel` + `(orbit, config, m, x_p)`.
#   2. `mode_self_force_history`    — the per-AMR-level sequence (rows 50 / 53's input).
#   3. `mode_force_convergence`     — its successive-difference / ratio table and the
#                                     resolution floor (main.tex:719-724).
#      2 and 3 come in two forms: the `AmrHistory` one and the LEVEL one, which needs only
#      what the row-47 fan-out ships back (`ModeResult.levels`, `.x_p`) because the live
#      `ModeSolution` is process-local and deliberately does not travel. Both go through the
#      single `_sf_level_force` / `_sf_convergence` pair, so they cannot drift.
#   4. `regularized_mode_force`     — `F_r^m − F^m_{r[order]}`, the per-m term row 46 SUMS;
#                                     the folded/unfolded subtraction term itself is
#                                     [`mmode_reg_term_Fr`](@ref), owned by
#                                     `Extraction/regularization.jl` (row 45).
# The field and its gradient come from the SEALED row-42 read-out surface (`field_at`,
# `field_grad_at`) or from the SEALED row-43 per-level diagnostics (`AmrLevel.psi/dpsi_r/
# dpsi_w`); no stencil, interpolation or basis is re-derived here. THE MODE SUM (row 46), the
# process fan-out (row 47) and every figure (rows 48-53) are other nodes.
#
# MEASURED numbers do NOT live in this file's docstrings — they are configuration-scoped and
# mesh-dependent. They are in `docs/src/capabilities/scalarselfforce.md` (with their
# configuration and status) and in the node evidence.
#
# SOURCES (read-only mirrors):
#   * Osburn & Nishimura, arXiv:2206.07031 `src/main.tex:178-184` — Eq. (3.10), the formula;
#     `:89-94` — their mode ansatz and `Δφ(r) = a/(r₊−r₋) ln((r−r₊)/(r−r₋))`.
#   * CRP arXiv:2606.04998 `main.tex:428-430` (the field normalization and
#     `F^self_α = q ∇_α Φ^R|_{x_p}`), `main.tex:701-703`, `main.tex:719-724`.
#   * spectre `src/Elliptic/Systems/SelfForce/Scalar/Events/ObserveSelfForce.cpp:52-77`
#     (the formula) and `.../ObserveSelfForce.hpp:89-107` (the reduction: `funcl::Plus`
#     then `funcl::Divides` by `NumContributingElements`); regression rows
#     `tests/InputFiles/SelfForce/ScalarCircularOrbit.yaml:19-23`, column legend
#     `tests/Unit/Elliptic/Systems/SelfForce/Scalar/Events/Test_ObserveSelfForce.cpp:195-199`.

export mode_self_force, mode_self_force_history, mode_force_convergence,
       regularized_mode_force

# `∂_z Ψ` from the domain's own angular chart: `w = z²` under equatorial symmetry (so
# `∂_z = 2z ∂_w`, an EXACT zero at the equator) and `z = cosθ` without it (row 39's charts).
# `x2` is clamped at 0 first: a chart coordinate that is a round-off-negative `w` is still the
# equator, and `sqrt` of it must not throw.
function _sf_dpsi_z(chart::Symbol, x2::Real, dpsi_x2::Number)
    chart === :z2 && return 2 * sqrt(max(x2, zero(x2))) * dpsi_x2
    chart === :z && return dpsi_x2
    throw(ArgumentError("mode_self_force: unknown angular chart :$chart " *
                        "(expected :z2 or :z — see `block_lattice_bounds`)"))
end

# The radial chart contract, checked rather than assumed (see the `mode_self_force` docstring):
# `x1` must be the horizon-penetrating Boyer–Lindquist `r`.
function _sf_check_chart(config::RunConfig)
    config.orbit.penetrating_horizon || throw(ArgumentError(
        "mode_self_force: this configuration declares PenetratingHorizon = false, so its " *
        "radial chart is not Boyer–Lindquist r and ∂_{x1} is not ∂_r. This node deliberately " *
        "implements only the penetrating chart (see the docstring's RADIAL CHART note); a " *
        "tortoise-chart extraction needs the α = Δ/(r²+a²) factor of spectre " *
        "ObserveSelfForce.cpp:61-65 and a test that exercises it."))
    return nothing
end

# ============================== the formula ================================================

"""
    mode_self_force(orbit::CircularOrbit, m, psi, dpsi_r, dpsi_z = zero(psi); phi_p = 0)
    mode_self_force(sol::ModeSolution, m = sol.m; element = nothing, phi_p = sol.phi_p)
    mode_self_force(hist::AmrHistory, m = hist.m; level = lastindex(hist.levels),
                    reduction = :trace, phi_p = …)

The `m`-mode contribution to the scalar self-force at the particle, Osburn–Nishimura
Eq. (3.10) (`arXiv:2206.07031` `main.tex:178-184`) in this package's locked conventions.
Returns a `NamedTuple` (fields below).

# UNITS

`q = 1` and `M = 1` package-wide, so every returned component is the NORMALIZED
`(M²/q²) F_α^m` — the product the Warburton–Barack tables print and the quantity
`fig:ssf_error_budget` sums. The field normalization is the CRP's own,
`F^self_α = q ∇_α Φ^ℛ|_{x_p}` (`main.tex:428-430`).

# THE NORMATIVE CONVENTION STATEMENT for extraction

(The ansatz side is locked in [`mode_term`](@ref), the subtraction side in
[`mmode_reg_factor_over_pi`](@ref).) The residual field is the locked mode sum

    Φ^ℛ = Σ_{m∈ℤ} Ψ^ℛ_m(r, z) · sin^{|m|}θ / r · e^{+i m (φ̄ − Ω t)},   φ̄ = φ + χ(r),

with `χ = kerr_varphi_offset(a, r)` and `dχ/dr = a/Δ` (`kerr_varphi_offset_dr`). Folding the
conjugate pairs `Ψ^ℛ_{−m} = conj(Ψ^ℛ_m)` gives Eq. (3.10):

    F^{m=0}_α = q ∂_α (Ψ^ℛ_0 / r)                              — UNFOLDED, no ×2, no phase
    F^{m>0}_α = 2q Re[ e^{+i m (φ̄_p − Ω t_p)} ∂_α (Ψ^ℛ_m / r) ]

evaluated at `x_p = (t_p, r₀, π/2, φ_p)`. At the equator `sin^{|m|}θ = 1` and its
θ-derivative `|m| sin^{|m|−1}θ cosθ` vanishes, so only `Ψ^ℛ_m` and the explicit `r`- and
phase-dependence survive.

## The radial component — the spectre-compatible contract

    F_r^m ∝ ∂_rΨ^ℛ_m/r₀ − Ψ^ℛ_m/r₀² + i m (a/Δ) Ψ^ℛ_m/r₀

`∂_r` is taken at fixed BOYER–LINDQUIST `φ`, so it hits `Ψ^ℛ_m`, the `1/r` AND the phase; the
last term is the DEPHASING DERIVATIVE. Dropping it is an `O(a)` error. It is the term spectre
adds by hand at `ObserveSelfForce.cpp:70-72`, with the sign fixed by the `+im φ̄` of the
ansatz — the conjugate of the `−i a m/Δ` the sealed puncture applies to `Ψ^𝒫_m`
([`puncture_psi_chart`](@ref), `CircularOrbit.cpp:337-339`). With `φ_p = 0` this is exactly
spectre's `extract_self_force`, whose complex result is [`Fr_complex`](@ref
`mode_self_force`).

RADIAL CHART: `x1` is the horizon-penetrating Boyer–Lindquist `r`, so `∂_{x1} ≡ ∂_r` and
spectre's `/= alpha` at `ObserveSelfForce.cpp:61-65` — the `α = Δ/(r²+a²)` that converts
`∂_{r_*}` to `∂_r` in its NON-penetrating chart — is deliberately absent here. The
`ModeSolution` / `AmrHistory` methods THROW if a configuration ever declares
`PenetratingHorizon = false` rather than silently returning a `∂_{r_*}`-scaled number.

## The other components — this package's extension beyond spectre

spectre observes only the radial force (and hard-codes `F_θ = 0`); the rest follows from the
same ansatz phase and is provided here, gated against a finite-difference oracle of the full
folded field:

  * `∂_t → −imΩ` and `∂_φ → +im`, so `F_t^m = −Ω F_φ^m` exactly, mode by mode (the
    circular-orbit relation of Warburton–Barack `Warburton.tex:134`), and both vanish at
    `m = 0`;
  * `F_θ^m = −2q Re[e^{imφ̄_p} ∂_zΨ^ℛ_m] / r₀` through `∂_θ = −sinθ ∂_z`. In the production
    `z²` chart `∂_z = 2z ∂_{z²}` and `z = 0` at the equator, so `F_θ^m` is a BITWISE zero —
    what spectre hard-codes at `ObserveSelfForce.cpp:49-51`. In the full `z` chart it is zero
    only up to truncation error and is returned as measured;
  * `phi_p` is the particle's Boyer–Lindquist azimuth (`t_p = 0`). The CRP and spectre both
    use `φ_p = 0`, which makes the phase `e^{+imχ(r₀)}` — spectre's `rotation`.

`m ≥ 0` always (retrograde orbits are `a < 0`, never `m < 0`).

# WHICH FIELD, AND WHICH REDUCTION

`sol.u` IS `Ψ^ℛ_m` inside the worldtube (row 42's contract) and the particle is inside it, so
the extraction reads `sol.u` directly — never [`reconstruct_full_field`](@ref), which would
add `Ψ^𝒫_m` back.

The puncture sits on an element CORNER, so a DG field has one value per touching element.
`mode_self_force(sol)` DEFAULTS TO THE MEAN over every containing element — spectre parity
(`ObserveSelfForce.hpp:89-107` sums the per-element forces and divides by
`NumContributingElements`). Because Eq. (3.10) is LINEAR in `(Ψ^ℛ_m, ∂_rΨ^ℛ_m, ∂_zΨ^ℛ_m)`,
averaging the ingredients and applying the formula once is the same reduction as averaging the
forces, and it is what this file does (gated equal to the mean-of-forces in the node tests).
`element = e` selects one element's one-sided trace instead — the diagnostic form.

The `AmrHistory` method reads the SEALED per-level `AmrLevel` diagnostics, which record ONE
one trace per containing element plus their mean, so its default `reduction = :trace` is a
single-element trace at EVERY level and the per-level sequence is homogeneous, while
`reduction = :mean` is spectre's reduction and is now available at
the final level, where `hist.solution` exists; it throws elsewhere rather than silently mixing
EVERY level — neither needs a `ModeSolution`. Do not mix reductions inside one convergence
sequence; `AmrLevel.trace_spread` bounds the difference, for `dpsi_w` as well as `psi`/`dpsi_r`.

# FIELDS OF THE RESULT

`Fr`, `Ft`, `Ftheta`, `Fphi` — the real (normalized) force components; `Fr_complex`,
`Ft_complex`, `Ftheta_complex`, `Fphi_complex` — the complex quantities BEFORE the real part
is taken (spectre's regression rows print `Re`/`Im` of `Fr_complex`); `m`, `phase`
(`e^{im(φ_p+χ(r₀))}`), `psi`, `dpsi_r`, `dpsi_z` — the ingredients actually used;
`reduction` ∈ `(:direct, :mean, :trace)`; `element` (the trace's element, `0` for a mean or
for the raw-ingredient method) and `elements` (all containing elements). The `AmrHistory`
method adds `step` and `ndofs`.

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
        return _sf_result(m, one(CT), brad, z, z, bpol, p, dr, dz, :direct, 0, Int[])
    end
    ph = CT(cis(m * kerr_varphi(RT(orbit.a), r0, RT(phi_p))))
    chi1 = RT(kerr_varphi_offset_dr(RT(orbit.a), r0))              # χ′ = a/Δ (sealed)
    Frc = 2 * ph * (brad + im * m * chi1 * p / r0)
    Fphc = 2 * ph * (im * m) * p / r0
    Ftc = -RT(orbit.Omega) * Fphc                                  # ∂_t = −Ω ∂_φ on the phase
    Fthc = 2 * ph * bpol
    return _sf_result(m, ph, Frc, Ftc, Fphc, Fthc, p, dr, dz, :direct, 0, Int[])
end

_sf_result(m, phase, Frc, Ftc, Fphc, Fthc, psi, dpsi_r, dpsi_z, reduction, element, elements) =
    (Fr = real(Frc), Ft = real(Ftc), Ftheta = real(Fthc), Fphi = real(Fphc),
     Fr_complex = Frc, Ft_complex = Ftc, Ftheta_complex = Fthc, Fphi_complex = Fphc,
     m = Int(m), phase = phase, psi = psi, dpsi_r = dpsi_r, dpsi_z = dpsi_z,
     reduction = reduction, element = Int(element), elements = elements)

function mode_self_force(sol::ModeSolution, m::Integer = sol.m; element = nothing,
                         phi_p::Real = sol.phi_p)
    Int(m) == sol.m || throw(ArgumentError(
        "mode_self_force: this ModeSolution solves m = $(sol.m), got m = $m — a different " *
        "mode is a different problem (a different RunConfig), never a keyword"))
    _sf_check_chart(sol.config)
    xp = puncture_chart_position(sol.dom, sol.orbit)
    # `Solve/driver.jl`'s `_dr_containing` is the ONE containment locator of the package
    # (`_dr_locate` is its first entry, which is what `field_at` uses by default); it returns
    # an empty list rather than throwing, so the emptiness check belongs to the caller.
    es = _dr_containing(sol.mesh, xp[1], xp[2])
    isempty(es) && throw(ArgumentError(
        "mode_self_force: the puncture $(xp) lies in no element of this mesh"))
    if element === nothing
        n = length(es)
        psi = sum(field_at(sol, xp[1], xp[2]; element = e) for e in es) / n
        g1 = sum(field_grad_at(sol, xp[1], xp[2]; element = e)[1] for e in es) / n
        g2 = sum(field_grad_at(sol, xp[1], xp[2]; element = e)[2] for e in es) / n
        r = mode_self_force(sol.orbit, sol.m, psi, g1,
                            _sf_dpsi_z(sol.dom.chart, xp[2], g2); phi_p = phi_p)
        return merge(r, (reduction = :mean, element = 0, elements = es))
    end
    e = Int(element)
    psi = field_at(sol, xp[1], xp[2]; element = e)
    g = field_grad_at(sol, xp[1], xp[2]; element = e)
    r = mode_self_force(sol.orbit, sol.m, psi, g[1],
                        _sf_dpsi_z(sol.dom.chart, xp[2], g[2]); phi_p = phi_p)
    return merge(r, (reduction = :trace, element = e, elements = es))
end

function mode_self_force(hist::AmrHistory, m::Integer = hist.m;
                         level::Integer = lastindex(hist.levels),
                         reduction::Symbol = :trace,
                         phi_p::Real = hist.solution.phi_p)
    Int(m) == hist.m || throw(ArgumentError(
        "mode_self_force: this AmrHistory is the m = $(hist.m) sweep, got m = $m"))
    _sf_check_chart(hist.config)
    l = hist.levels[level]
    # BOTH reductions are available at EVERY level: `AmrLevel` stores one trace per containing
    # element and their mean, so neither needs the (process-local, and for every level but the
    # last discarded) solution vector. At the final level the stored mean is BITWISE
    # `mode_self_force(hist.solution)` (gated).
    return _sf_level_force(hist.orbit, hist.config, hist.m, l, hist.x_p, phi_p;
                           reduction = reduction)
end

# ===================== the level route: extraction WITHOUT a live ModeSolution ==============
# The `AmrHistory` methods above and the `(orbit, config, m, levels, x_p)` methods below are
# the SAME code through this one function, so the two routes cannot drift.

# BOTH reductions from a level alone: `:mean` is spectre's (what `AmrLevel.psi/dpsi_r/dpsi_w`
# record), `:trace` is one containing element's one-sided trace — `element`, or the level's
# representative `l.element` when none is named. `:trace` is the DEFAULT because a per-level
# SEQUENCE must be homogeneous and because it is what every recorded table was measured with;
# ask for `:mean` explicitly to get the reduction spectre reports.
function _sf_level_force(orbit::CircularOrbit, config::RunConfig, m::Integer, l::AmrLevel,
                         x_p, phi_p::Real; reduction::Symbol = :trace,
                         element::Union{Nothing,Integer} = nothing)
    _sf_check_chart(config)
    reduction in (:trace, :mean) || throw(ArgumentError(
        "mode_self_force: reduction must be :trace or :mean, got :$(reduction)"))
    chart = block_lattice_bounds(config)[3]          # the config IS the chart's source
    if reduction === :mean
        element === nothing || throw(ArgumentError(
            "mode_self_force: reduction = :mean is the reduction over ALL containing " *
            "elements, so it takes no `element`"))
        r = mode_self_force(orbit, m, l.psi, l.dpsi_r, _sf_dpsi_z(chart, x_p[2], l.dpsi_w);
                            phi_p = phi_p)
        return merge(r, (reduction = :mean, element = 0, elements = l.elements),
                     (step = l.step, ndofs = l.mesh.ndofs))
    end
    e = element === nothing ? l.element : Int(element)
    k = findfirst(==(e), l.elements)
    k === nothing && throw(ArgumentError(
        "mode_self_force: element $(e) does not contain the puncture at this level; the " *
        "containing elements are $(l.elements)"))
    t = l.traces[k]
    r = mode_self_force(orbit, m, t[1], t[2], _sf_dpsi_z(chart, x_p[2], t[3]); phi_p = phi_p)
    return merge(r, (reduction = :trace, element = e, elements = l.elements),
                 (step = l.step, ndofs = l.mesh.ndofs))
end

# the mode number is the one thing an `AmrLevel` does NOT record, so the level route checks it
# against the configuration that produced the level (`mode_config(config, m)` in the fan-out,
# `crp_paper_config(; m_mode = m)` in a driver). A silently wrong `m` would rotate the ansatz
# phase by the wrong amount and still return a plausible number.
function _sf_check_mode(config::RunConfig, m::Integer)
    Int(config.orbit.m_mode_number) == Int(m) || throw(ArgumentError(
        "mode_self_force: this configuration solves m = $(config.orbit.m_mode_number) but the " *
        "levels were passed as m = $m. `AmrLevel`s carry no mode number, so the config must " *
        "be the one that produced them — `mode_config(config, m)` for a fan-out result."))
    return nothing
end

"""
    mode_self_force(orbit, config, m, level::AmrLevel, x_p; reduction = :trace, phi_p = 0)

`F^m_α` from ONE hp-AMR level's sealed diagnostics, addressed WITHOUT a live
[`ModeSolution`](@ref) — the form the process fan-out's payload supports.

A worker ships back a [`ModeResult`](@ref) carrying `levels::Vector{AmrLevel}` and `x_p`; the
`ModeSolution` is process-local by design and does not travel. This method takes exactly those
two pieces plus the `orbit` and the per-mode `config` (`mode_config(config, m)`), which is
where the angular chart comes from ([`block_lattice_bounds`](@ref) — `:z2` under equatorial
symmetry, so `∂_zΨ = 2z ∂_{z²}Ψ`) and where the radial-chart and mode-number contracts are
checked. The result is identical, field by field, to `mode_self_force(hist; level = k)` — both
go through the same internal function.
"""
function mode_self_force(orbit::CircularOrbit, config::RunConfig, m::Integer, l::AmrLevel,
                         x_p; reduction::Symbol = :trace, phi_p::Real = 0)
    _sf_check_mode(config, m)
    return _sf_level_force(orbit, config, m, l, x_p, phi_p; reduction = reduction)
end

# ========================= the per-level sequence and its convergence ======================

"""
    mode_self_force_history(hist::AmrHistory; phi_p = hist.solution.phi_p) -> Vector

[`mode_self_force`](@ref) at EVERY level of an hp-AMR sweep, level 0 first — the sequence
`fig:ssf_convergence` plots and the error budget's resolution floor is measured from
(CRP `main.tex:719-724`). Each entry carries the extra `step` and `ndofs` fields.

`reduction` selects the per-level reduction and BOTH are available at EVERY level, from the
traces the sealed [`AmrLevel`](@ref) stores — no level needs a `ModeSolution`, which matters
because only the last one keeps one and a fanned-out sweep ships none.

  - `:trace` (the DEFAULT) — the level's representative one-sided trace, `AmrLevel.element`.
    Homogeneous across the sequence and what every recorded table in this package was measured
    with, so it is what a convergence sequence keeps using.
  - `:mean` — SPECTRE'S OWN REDUCTION over the level's `NumContributingElements`
    (`ObserveSelfForce.hpp:89-107`), i.e. `AmrLevel.psi/dpsi_r/dpsi_w`. At the final level it is
    BITWISE `mode_self_force(hist.solution)` (gated). It is also the better resolution measure —
    the face jump it cancels is first order in `h`: MEASURED, the sweep's
    `relative_resolution_error` is 6.62e-6 under `:mean` against 3.39e-5 under `:trace` on the
    production half domain, and 1.60e-7 against 6.66e-6 on the full-`z` lattice.

The two differ by the trace spread, which `AmrLevel.trace_spread` bounds — now for `dpsi_w` as
well as `psi` and `dpsi_r`, because in the full-`z` chart `dpsi_w` is what carries `F_θ`.
"""
mode_self_force_history(hist::AmrHistory; reduction::Symbol = :trace,
                       phi_p::Real = hist.solution.phi_p) =
    [_sf_level_force(hist.orbit, hist.config, hist.m, l, hist.x_p, phi_p;
                     reduction = reduction) for l in hist.levels]

"""
    mode_self_force_history(orbit, config, m, levels::AbstractVector{<:AmrLevel};
                            x_p, phi_p = 0) -> Vector

The per-level sequence from a fanned-out sweep's `levels` — the same list
`mode_self_force_history(hist)` returns, addressed through the pieces a
[`ModeResult`](@ref) carries (`res.levels`, `res.x_p`) instead of through a live
[`ModeSolution`](@ref). `config` must be the per-mode configuration that solved `m`
(`mode_config(config, m)`); `x_p` is the puncture in the discretization chart.
"""
function mode_self_force_history(orbit::CircularOrbit, config::RunConfig, m::Integer,
                                 levels::AbstractVector{<:AmrLevel};
                                 x_p = puncture_chart_position(block_domain(orbit, config),
                                                               orbit),
                                 reduction::Symbol = :trace, phi_p::Real = 0)
    _sf_check_mode(config, m)
    return [_sf_level_force(orbit, config, m, l, x_p, phi_p; reduction = reduction)
            for l in levels]
end

"""
    mode_force_convergence(hist::AmrHistory; component = :Fr, phi_p = …) -> NamedTuple

The successive-difference convergence table of one force component over an hp-AMR sweep:
`steps`, `ndofs`, `values`, `diffs` (`|Δ|` between consecutive levels), `ratios`
(consecutive `diffs`), `two_step_ratios` (`diffs[k+2]/diffs[k]`), and the two headline
numbers `resolution_error = diffs[end]` and `relative_resolution_error`.

`resolution_error` is the CRP's own "error due to numerical resolution ... the difference
between the highest and second-highest AMR level" (`main.tex:722`), applied to `F_r^m` rather
than to a raw field diagnostic. The DEFINITION is the row-43 one — the last successive
difference, relative denominator `|values[end]|` — so it is the `F_r^m` analogue of that
node's `resolution_error(hist; field)` helper and the two are directly comparable. It is a
DISCRETIZATION/round-off floor, not the GMRES tolerance (row 43 re-ran its sweep at
`abstol = 1e-14` and the floor did not move), so a caller must not try to buy it down there.

TWO-LEVEL CONTRACTION, NOT ONE: the single-step ratio of an exponentially converging sequence
with plateaus is not bounded by a constant below 1 — `ratios` is reported, not contracted.
`two_step_ratios` is the quantity that contracts. Requires at least 4 levels (three
differences), so that a two-step ratio exists at all; fewer throws rather than returning a
vacuously empty table. Per-level values are MESH-DEPENDENT: gate on tolerances and on the
shape of the sequence, never on a golden per-level value. `component` must name a real field
of [`mode_self_force`](@ref)'s result (`:Fr`, `:Ft`, `:Ftheta`, `:Fphi`).
"""
mode_force_convergence(hist::AmrHistory; component::Symbol = :Fr,
                       phi_p::Real = hist.solution.phi_p) =
    _sf_convergence(mode_self_force_history(hist; phi_p = phi_p), component,
                    length(hist.levels))

"""
    mode_force_convergence(orbit, config, m, levels::AbstractVector{<:AmrLevel};
                           x_p, component = :Fr, phi_p = 0) -> NamedTuple

`mode_force_convergence` on a fanned-out sweep's `levels` — the table `fig:ssf_convergence`
plots, reachable from a [`ModeResult`](@ref) (`res.levels`, `res.x_p`) rather than only from a
live [`ModeSolution`](@ref). Identical, field by field, to the [`AmrHistory`](@ref) method on
the same sweep. `config` must be the per-mode configuration that solved `m`.
"""
function mode_force_convergence(orbit::CircularOrbit, config::RunConfig, m::Integer,
                                levels::AbstractVector{<:AmrLevel};
                                x_p = puncture_chart_position(block_domain(orbit, config),
                                                              orbit),
                                component::Symbol = :Fr, phi_p::Real = 0)
    _sf_check_mode(config, m)
    fs = mode_self_force_history(orbit, config, m, levels; x_p = x_p, phi_p = phi_p)
    return _sf_convergence(fs, component, length(levels))
end

function _sf_convergence(fs, component::Symbol, nlev::Integer)
    component in (:Fr, :Ft, :Ftheta, :Fphi) || throw(ArgumentError(
        "mode_force_convergence: component must be one of (:Fr, :Ft, :Ftheta, :Fphi), " *
        "got :$component"))
    nlev >= 4 || throw(ArgumentError(
        "mode_force_convergence: needs at least 4 AMR levels for a two-step ratio to " *
        "exist, got $(nlev) — an emptier table would be a vacuous gate"))
    v = [getproperty(f, component) for f in fs]
    d = [abs(v[k + 1] - v[k]) for k in 1:length(v)-1]
    rat = [iszero(d[k]) ? oftype(d[k], NaN) : d[k + 1] / d[k] for k in 1:length(d)-1]
    two = [iszero(d[k]) ? oftype(d[k], NaN) : d[k + 2] / d[k] for k in 1:length(d)-2]
    resolution = d[end]
    return (steps = [f.step for f in fs], ndofs = [f.ndofs for f in fs], values = v,
            diffs = d, ratios = rat, two_step_ratios = two,
            resolution_error = resolution,
            relative_resolution_error = iszero(v[end]) ? oftype(resolution, NaN) :
                                        resolution / abs(v[end]))
end

# ============================ the regularized per-mode term ================================

"""
    regularized_mode_force(orbit::CircularOrbit, m, Fr::Real; order = 4) -> T
    regularized_mode_force(sol::ModeSolution, m = sol.m; order = 4, kwargs...) -> T
    regularized_mode_force(hist::AmrHistory, m = hist.m; order = 4, kwargs...) -> T

`F_r^m − F^m_{r[order]}`: the extracted radial force mode with the Heffernan large-m
regularization parameter subtracted. The subtracted term is [`mmode_reg_term_Fr`](@ref),
which owns the `m = 0` halving rule (`Extraction/regularization.jl`, row 45) — a caller must
never re-apply it. `order ∈ HEFFERNAN_MMODE_ORDERS`; the CRP uses 4 (`main.tex:722`).

THIS IS THE PER-m TERM THE MODE SUM ADDS UP. Because the doubly-infinite parameter sum
vanishes (`F^0/2 + Σ_{m≥1} F^m = 0`, gated in `test_regularization.jl`), the accelerated mode
sum of `main.tex:722` is literally

    F_r ≈ Σ_{m=0}^{m_max} regularized_mode_force(…, m) ,

identical to `Σ_{m≤m_max} F_r^m + `[`mmode_reg_tail_Fr`](@ref)`(order, orbit, m_max+1)`.
Performing that sum is node `mode-sum-total-force` (row 46), not this file.
"""
regularized_mode_force(orbit::CircularOrbit, m::Integer, Fr::Real; order::Integer = 4) =
    Fr - mmode_reg_term_Fr(order, orbit, m)

regularized_mode_force(sol::ModeSolution, m::Integer = sol.m; order::Integer = 4, kwargs...) =
    regularized_mode_force(sol.orbit, sol.m, mode_self_force(sol, m; kwargs...).Fr;
                           order = order)

regularized_mode_force(hist::AmrHistory, m::Integer = hist.m; order::Integer = 4, kwargs...) =
    regularized_mode_force(hist.orbit, hist.m, mode_self_force(hist, m; kwargs...).Fr;
                           order = order)
