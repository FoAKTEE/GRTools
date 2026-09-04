# Extraction/modesum.jl — the TOTAL radial scalar self-force: the accelerated sum over
# m-modes plus the large-m tail that the truncation leaves behind (ESSF node
# `mode-sum-total-force`, GLOBAL_DAG row 46; DESIGN §2.2 `Extraction/modesum.jl`).
#
# WHAT THIS FILE OWNS — and nothing else:
#   1. `total_self_force`  — `F_r = Σ_{m=0}^{m_max} [F_r^m − F^m_{r[order]}] + tail`, from raw
#                            per-mode numbers, from a `ModeSet` (row 47's fan-out), from a
#                            vector of `AmrHistory` (row 43), or solved end-to-end from an
#                            `(orbit, config)` pair.
#   2. `mode_sum_tail`     — the estimate of the REMAINDER `Σ_{m>m_max} [F_r^m − F^m_{r[order]}]`
#                            with its uncertainty, in closed form (never a numerical sum).
#   3. `mode_decay_fit`    — the log-log decay exponent of a mode sequence with the uncertainty
#                            the per-mode resolution errors propagate into it, and the
#                            headroom rule that decides which modes may enter a fit at all.
#
# The per-mode force `F_r^m` comes from the SEALED row-44 surface (`mode_self_force`,
# `regularized_mode_force`) and the subtraction term from the SEALED row-45 helper
# (`mmode_reg_term_Fr`, which owns the `m = 0` halving — it is NEVER re-applied here).
# Nothing about the DG solve, the AMR loop or the fan-out is re-derived.
#
# THE TWO TAILS — the distinction that this file exists to keep straight:
#
#   (T1) THE ANALYTIC PARAMETER TAIL, `mmode_reg_tail_Fr(order, orbit, m_max+1)`. Because the
#        doubly-infinite parameter sum vanishes (`F^0/2 + Σ_{m≥1}F^m = 0`, gated in
#        `test_regularization.jl`), summing the REGULARIZED per-mode terms already contains it:
#            Σ_{m≤M} [F_r^m − F^m_{r[4]}]  ==  Σ_{m≤M} F_r^m + mmode_reg_tail_Fr(4, orbit, M+1).
#        This IS the acceleration of CRP `main.tex:722`. It is not optional and not a model —
#        it is an identity, gated to round-off in `test_modesum.jl` (a).
#
#   (T2) THE REMAINDER OF THE ACCELERATED SUM, `Σ_{m>M} [F_r^m − F^m_{r[4]}]`. With a
#        fourth-order puncture the `[2]` and `[4]` parameters are gone, so the leading
#        unsubtracted content is `F^m_{r[6]} ∼ m^{-6}` ([H14]:1675-1680) — against the
#        `m^{-4}` of the UNREGULARIZED `F_r^m`, which asymptotes to `F^m_{r[4]}` itself.
#        `F_{r[6]}` has NO closed form in [H14] (`:1795-1797` fits it numerically), so its
#        amplitude is FITTED here on the resolved large-m modes while its m-dependence and its
#        infinite sum are taken EXACTLY from `mmode_reg_factor(6, ·)` /
#        `mmode_reg_factor_tail(6, ·)`. That is `method = :heffernan`. `method = :fit` frees
#        the exponent too and sums `C m^{-p}` by the Hurwitz zeta. `method = :none` sets the
#        remainder to zero and reports the CRP's own truncation estimate `|T_{m_max}|`
#        (`main.tex:722`) as its uncertainty.
#
#        Note the size relation this makes explicit: for `T_m ≈ C m^{-p}`,
#        `Σ_{m>M} T_m ≈ C M^{1-p}/(p-1) = |T_M| · M/(p-1)`, so the paper's
#        omit-the-last-mode estimator is smaller than the true remainder by ≈ `M/(p−1)`
#        (≈ 4 at `M = 20`, `p = 6`) — the "constant offset" of `main.tex:724`.
#
# SOURCES (read-only mirrors):
#   * CRP arXiv:2606.04998 `main.tex:712` (extract `F_r^m` per mode, sum over m),
#     `:722-723` (the `F_{r[4]}` acceleration), `:726-731` (the error budget: resolution,
#     truncation, true error).
#   * A. Heffernan, A. Ottewill, B. Wardell, arXiv:1211.6446 `KerrRegularization.tex:1675-1680`
#     (the m-dependence at orders 2, 4, 6, 8), `:1784-1798` (the `m^{-6}` residual).
#   * N. Warburton and L. Barack, arXiv:1003.1860 `Warburton.tex:620-645` Table III — the
#     `(M²/q²)F_r` reference values this node is gated against (baked, string-equal to the
#     print, in `test/fixtures/wb_reference.jl`).
#
# MEASURED numbers do NOT live in this file's docstrings — they are configuration-scoped.
# They are in `docs/src/capabilities/scalarselfforce.md` and in the node evidence.

export total_self_force, TotalSelfForce, mode_sum_tail, ModeSumTail,
       mode_decay_fit, ModeDecayFit

# ---------------------------------------------------------------- Hurwitz zeta (internal)
# `ζ(s, a) = Σ_{k≥0} (a+k)^{-s}` for real `s > 1`, `a > 0`, by Euler–Maclaurin: `n` explicit
# terms, then the integral plus the first three Bernoulli corrections. Needed ONLY by
# `method = :fit`, whose exponent is not an integer; the `:heffernan` shape sums in closed form
# through `mmode_reg_factor_tail`. Internal: a general special function is not this package's
# business, and the name would collide with every special-function library.
# ACCURACY: `n = 24` explicit terms plus three Bernoulli corrections give ~1e-16 RELATIVE for
# `s ≥ 1.5`, `a ≥ 1` — Float64-grade. The result inherits the argument TYPE (so a `BigFloat`
# call returns a `BigFloat`) but NOT arbitrary precision: measured 2.2e-18 relative at
# `BigFloat(256)`, i.e. the Euler–Maclaurin truncation, not the arithmetic, is the limit.
# That is far below the ~1e-2 uncertainty of the tail model it serves.
function _ms_hurwitz_zeta(s::Real, a::Real; n::Integer = 24)
    T = float(promote_type(typeof(s), typeof(a)))
    sT, aT = T(s), T(a)
    (sT > 1) || throw(DomainError(s, "_ms_hurwitz_zeta: needs s > 1 for convergence"))
    (aT > 0) || throw(DomainError(a, "_ms_hurwitz_zeta: needs a > 0"))
    acc = zero(T)
    for k in 0:(Int(n) - 1)
        acc += (aT + k)^(-sT)
    end
    b = aT + T(n)
    acc += b^(1 - sT) / (sT - 1)                       # ∫_b^∞ x^{-s} dx
    acc += b^(-sT) / 2                                 # B_2 term of Euler–Maclaurin
    acc += sT * b^(-sT - 1) / 12                       # B_2/2! · f'(b)
    acc -= sT * (sT + 1) * (sT + 2) * b^(-sT - 3) / 720
    acc += sT * (sT + 1) * (sT + 2) * (sT + 3) * (sT + 4) * b^(-sT - 5) / 30240
    return acc
end

# ============================== the decay fit =============================================

"""
    ModeDecayFit{T}

The log-log decay fit of a mode sequence, `|x_m| ≈ amplitude · m^{slope}`.

| field | what it is |
|---|---|
| `ms` | the modes that ENTERED the fit (in order) |
| `slope` | the least-squares slope of `log|x_m|` against `log m` (NEGATIVE for a decay) |
| `slope_sigma` | the uncertainty the per-mode errors propagate into `slope` (`NaN` when no errors were supplied) |
| `amplitude` | `exp(intercept)`, so `amplitude · m^slope` is the fitted magnitude |
| `excluded` | modes dropped for insufficient headroom, and `m ≤ 0` |
| `headroom` | `|x_m| / error_m` for the modes that entered (empty when no errors) |

See [`mode_decay_fit`](@ref).
"""
struct ModeDecayFit{T<:AbstractFloat}
    ms::Vector{Int}
    slope::T
    slope_sigma::T
    amplitude::T
    excluded::Vector{Int}
    headroom::Vector{T}
end

function Base.show(io::IO, f::ModeDecayFit)
    print(io, "ModeDecayFit(m ∈ ", isempty(f.ms) ? "∅" : "$(first(f.ms)):$(last(f.ms))",
          ", slope = ", f.slope, " ± ", f.slope_sigma, ")")
end

"""
    mode_decay_fit(ms, values; errors = nothing, headroom = 3, window = nothing) -> ModeDecayFit

Least-squares fit of `log|values[k]|` against `log(ms[k])`: the log-log DECAY EXPONENT of a
mode sequence, with the uncertainty its per-mode errors propagate into the slope.

* `errors[k]` is the uncertainty of `values[k]` — for a solved sweep, the per-mode resolution
  error (the last successive AMR difference, `mode_force_convergence(...).resolution_error`,
  CRP `main.tex:722`). It enters twice: as `σ(log|x_k|) = errors[k]/|values[k]|` in the
  propagated `slope_sigma`, and through the HEADROOM RULE.
* THE HEADROOM RULE (the discipline row 44 established and this node inherits): a mode may
  enter a fit only when `|values[k]| ≥ headroom · errors[k]`. A mode whose residual has sunk
  into its own resolution floor carries no information about the decay, and a fit that
  includes it is resolution-limited, not physical. Excluded modes are listed in `excluded`,
  never silently dropped. Pass `errors = nothing` to disable the rule (then `slope_sigma` is
  `NaN`, and it is the caller's job not to quote a bare number as a measurement).
* `window` keeps only the LAST `window` admissible modes (the asymptotic end of the sequence);
  `nothing` keeps all of them. `m ≤ 0` is always excluded — `log m` does not exist — and a
  value of exactly zero is excluded as well.

At least two admissible modes are required; fewer throws rather than returning a fabricated
slope. `slope_sigma` is the quadrature propagation
`σ² = Σ_k (X_k − X̄)² σ_{Y,k}² / (Σ_k (X_k − X̄)²)²` with `X = log m`, `Y = log|x|`, exact for
uncorrelated errors and reducing to `(σ_1 ⊕ σ_2)/Δ log m` on a two-point fit.
"""
function mode_decay_fit(ms, values::AbstractVector; errors = nothing,
                        headroom::Real = 3, window = nothing)
    mv = collect(Int, ms)
    length(mv) == length(values) || throw(ArgumentError(
        "mode_decay_fit: got $(length(mv)) modes and $(length(values)) values"))
    errors === nothing || length(errors) == length(values) || throw(ArgumentError(
        "mode_decay_fit: got $(length(errors)) errors and $(length(values)) values"))
    T = float(promote_type(eltype(values) <: Number ? real(eltype(values)) : Float64,
                           errors === nothing ? Float64 : real(eltype(errors)), Float64))
    keep, drop, head = Int[], Int[], T[]
    for k in eachindex(mv)
        m, v = mv[k], abs(T(values[k]))
        if m <= 0 || iszero(v)
            push!(drop, m); continue
        end
        if errors !== nothing
            e = abs(T(errors[k]))
            h = iszero(e) ? T(Inf) : v / e
            if h < headroom
                push!(drop, m); continue
            end
            push!(keep, k); push!(head, h)
        else
            push!(keep, k)
        end
    end
    if window !== nothing && length(keep) > Int(window)
        n = Int(window)
        head = isempty(head) ? head : head[end-n+1:end]
        keep = keep[end-n+1:end]
    end
    length(keep) >= 2 || throw(ArgumentError(
        "mode_decay_fit: only $(length(keep)) mode(s) pass the headroom rule " *
        "(headroom = $headroom, excluded $(drop)) — a slope needs at least two. Solve more " *
        "AMR levels, lower `headroom` deliberately, or do not quote a slope."))
    X = T[log(T(mv[k])) for k in keep]
    Y = T[log(abs(T(values[k]))) for k in keep]
    xb, yb = sum(X) / length(X), sum(Y) / length(Y)
    Sxx = sum((X .- xb) .^ 2)
    slope = sum((X .- xb) .* (Y .- yb)) / Sxx
    inter = yb - slope * xb
    sig = T(NaN)
    if errors !== nothing
        sy = T[abs(T(errors[k])) / abs(T(values[k])) for k in keep]
        sig = sqrt(sum(((X .- xb) .^ 2) .* (sy .^ 2))) / Sxx
    end
    return ModeDecayFit{T}(mv[keep], slope, sig, exp(inter), drop, head)
end

# ============================== the large-m tail ==========================================

"""
    ModeSumTail{T}

The estimate of the mode sum's REMAINDER `Σ_{m>m_max} [F_r^m − F^m_{r[order]}]` — tail (T2) of
this file's header — with everything needed to judge it.

| field | what it is |
|---|---|
| `value` | the estimate itself, in the same units as the sum (`(M²/q²)F_r`) |
| `uncertainty` | statistical (amplitude spread) and model (`:heffernan` vs `:fit`) added in quadrature |
| `method` | `:heffernan`, `:fit` or `:none` — which estimate `value` is |
| `alt_value` | the OTHER model's estimate (`:heffernan` ↔ `:fit`); `NaN` when only one is available |
| `exponent` | the decay exponent used: `6` exactly for `:heffernan`, the fitted `p` for `:fit` |
| `exponent_sigma` | the propagated uncertainty of the FITTED exponent (`NaN` for `:heffernan`) |
| `amplitude` | the fitted amplitude: `F_{r[6]}` for `:heffernan`, `C` in `C m^{-p}` for `:fit` |
| `amplitude_sigma` | its uncertainty: the per-mode SHAPE MISMATCH and the per-mode RESOLUTION errors, propagated separately and added in quadrature |
| `fit` | the [`ModeDecayFit`](@ref) of the regularized sequence, `nothing` for `:none` |
| `truncation_error` | `|T_{m_max}|`, the CRP's own omit-the-last-mode estimator (`main.tex:722`) |
| `m_start` | the first mode the tail covers, `m_max + 1` |

See [`mode_sum_tail`](@ref).
"""
struct ModeSumTail{T<:AbstractFloat}
    value::T
    uncertainty::T
    method::Symbol
    alt_value::T
    exponent::T
    exponent_sigma::T
    amplitude::T
    amplitude_sigma::T
    fit::Union{Nothing,ModeDecayFit{T}}
    truncation_error::T
    m_start::Int
end

function Base.show(io::IO, t::ModeSumTail)
    print(io, "ModeSumTail(", t.method, ": ", t.value, " ± ", t.uncertainty,
          ", m ≥ ", t.m_start, ")")
end

"""
    mode_sum_tail(orbit, ms, terms; method = :heffernan, order = 4, window = 6,
                  errors = nothing, headroom = 3) -> ModeSumTail

The remainder `Σ_{m>m_max} [F_r^m − F^m_{r[order]}]` of an accelerated mode sum whose terms
`terms[k]` (= [`regularized_mode_force`](@ref) at `ms[k]`) are known up to `m_max = last(ms)`.

`ms` must be a contiguous `0:m_max`-style run of folded modes. THE ESTIMATE IS ALWAYS A CLOSED
FORM — no numerical sum is ever run to a cut-off.

# THE THREE METHODS

* `:heffernan` (default). The remainder is modelled as the NEXT Heffernan parameter,
  `T_m ≈ F_{r[6]} · F^m_{r[6]}/F_{r[6]}`, i.e. the exact printed m-dependence
  [`mmode_reg_factor`](@ref)`(6, m)` ([H14] `KerrRegularization.tex:1675-1680`) with the one
  unknown amplitude `F_{r[6]}` — which [H14] never puts in closed form (`:1795-1797` fits it)
  — obtained here as the mean of `terms[k] / mmode_reg_factor(6, ms[k])` over the fit window.
  The sum is then EXACT: `amplitude · `[`mmode_reg_factor_tail`](@ref)`(6, m_max+1)`.
  This is the physically-motivated estimate: it imposes the analytically known shape and fits
  only what is genuinely unknown.
* `:fit`. Both the exponent and the amplitude are free: `|T_m| ≈ C m^{-p}` from
  [`mode_decay_fit`](@ref), summed by the Hurwitz zeta, `± C ζ(p, m_max+1)`, with the sign of
  the fitted modes. Weaker prior, and it is the CONTROL on `:heffernan`: the two agreeing is
  evidence that the `m^{-6}` model is the right one.
* `:none`. `value = 0`; the sum stops where it stops. `uncertainty` is then the CRP's own
  truncation estimator `|T_{m_max}|` (`main.tex:722`), which for a `C m^{-p}` tail is smaller
  than the true remainder by a factor `≈ m_max/(p−1)`.

# THE WINDOW, THE HEADROOM AND THE UNCERTAINTY

`window` (default 6) is the number of trailing modes the fit uses; `errors` (per-mode
resolution errors) enforce the headroom rule of [`mode_decay_fit`](@ref) — a mode that has sunk
into its own resolution floor never enters. `uncertainty` is the quadrature sum of

  * STATISTICAL: `|value| · amplitude_sigma/|amplitude|`, where `amplitude_sigma` combines the
    per-mode SHAPE MISMATCH (the spread of the per-mode amplitudes about their mean) with the
    per-mode RESOLUTION errors PROPAGATED through the same mean, `√(Σ (err_k/shape_k)²)/n`.
    The headroom rule only decides which modes may enter; it does not propagate anything, so
    the propagation is done explicitly here; and
  * MODEL: `|value − alt_value|`, the disagreement between the two independent models.

Both models are always computed, so `alt_value` is populated whenever the window admits a fit.

`order ∈ HEFFERNAN_MMODE_ORDERS` is the order ALREADY subtracted from `terms` (4 for the CRP's
fourth-order puncture); the tail model uses `order + 2`, the leading parameter that survives.
"""
function mode_sum_tail(orbit::CircularOrbit, ms, terms::AbstractVector;
                       method::Symbol = :heffernan, order::Integer = 4,
                       window::Integer = 6, errors = nothing, headroom::Real = 3)
    method in (:heffernan, :fit, :none) || throw(ArgumentError(
        "mode_sum_tail: method must be :heffernan, :fit or :none, got :$method"))
    mv = collect(Int, ms)
    length(mv) == length(terms) || throw(ArgumentError(
        "mode_sum_tail: got $(length(mv)) modes and $(length(terms)) terms"))
    isempty(mv) && throw(ArgumentError("mode_sum_tail: no modes"))
    mv == collect(minimum(mv):maximum(mv)) || throw(ArgumentError(
        "mode_sum_tail: the modes must be a contiguous ascending run, got $mv"))
    next = Int(order) + 2
    next in HEFFERNAN_MMODE_ORDERS || throw(ArgumentError(
        "mode_sum_tail: the tail model needs the printed m-dependence of order $(next) " *
        "(one of $HEFFERNAN_MMODE_ORDERS); order = $order has none"))
    T = float(promote_type(typeof(orbit.r0), real(eltype(terms))))
    mmax = maximum(mv)
    mstart = mmax + 1
    trunc_err = abs(T(terms[end]))

    if method === :none
        return ModeSumTail{T}(zero(T), trunc_err, :none, T(NaN), T(NaN), T(NaN), T(NaN),
                              T(NaN), nothing, trunc_err, mstart)
    end

    # --- the free fit (also the control on the Heffernan shape) ----------------------------
    fit = mode_decay_fit(mv, terms; errors = errors, headroom = headroom, window = window)
    sgn = T(sign(T(terms[findfirst(==(last(fit.ms)), mv)])))
    p = -fit.slope
    method === :fit && p <= 1 && throw(ArgumentError(
        "mode_sum_tail: the free fit found an exponent p = $(p) ≤ 1 over m ∈ $(fit.ms), so " *
        "Σ_{m>m_max} C m^{-p} does not converge — the sequence is not in its asymptotic " *
        "regime (or is resolution-limited). Use method = :heffernan, or solve more modes."))
    v_fit = p > 1 ? sgn * T(fit.amplitude) * _ms_hurwitz_zeta(p, T(mstart)) : T(NaN)

    # --- the Heffernan shape: exact m-dependence, one fitted amplitude ----------------------
    shape = T[mmode_reg_factor(T, next, m) for m in fit.ms]
    amps = T[T(terms[findfirst(==(m), mv)]) / s for (m, s) in zip(fit.ms, shape)]
    A = sum(amps) / length(amps)
    # TWO INDEPENDENT contributions to the amplitude's uncertainty, kept separate and then
    # added in quadrature:
    #   SHAPE MISMATCH — how far the per-mode amplitudes spread around their mean, i.e. how
    #     badly a single order-(order+2) parameter describes the window (it also carries the
    #     next unsubtracted order); and
    #   RESOLUTION — the per-mode numerical-resolution errors, PROPAGATED THROUGH THE SAME
    #     MEAN: `σ_A = √(Σ (err_k/shape_k)²)/n` for uncorrelated per-mode errors. Gating mode
    #     ADMISSION on the headroom rule does NOT propagate them, so this is done here.
    spread = maximum(abs.(amps .- A))
    sig_res = zero(T)
    if errors !== nothing
        idx = [findfirst(==(m), mv) for m in fit.ms]
        sig_res = sqrt(sum((T(abs(errors[i])) / abs(sh))^2 for (i, sh) in zip(idx, shape))) /
                  length(shape)
    end
    amp_sigma = sqrt(spread^2 + sig_res^2)
    v_heff = A * mmode_reg_factor_tail(T, next, mstart)

    value, alt, expo, esig, amp = method === :heffernan ?
        (v_heff, v_fit, T(next), T(NaN), A) :
        (v_fit, v_heff, p, T(fit.slope_sigma), T(fit.amplitude))
    # the amplitude uncertainty is RELATIVE to the fitted amplitude, so it transfers to either
    # model's tail value unchanged (both are linear in their amplitude)
    stat = iszero(A) ? zero(T) : abs(value) * amp_sigma / abs(A)
    model = isfinite(alt) ? abs(value - alt) : zero(T)
    sigma_amp = method === :heffernan ? amp_sigma :
                (iszero(A) ? T(NaN) : abs(T(fit.amplitude)) * amp_sigma / abs(A))
    return ModeSumTail{T}(value, sqrt(stat^2 + model^2), method, alt, expo, esig, amp,
                          sigma_amp, fit, trunc_err, mstart)
end

# ============================== the total force ===========================================

"""
    TotalSelfForce{T}

The total radial scalar self-force and the whole budget behind it — the product of ESSF node
`mode-sum-total-force` and the input of `fig:ssf_error_budget` (CRP `main.tex:722-724`).

| field | what it is |
|---|---|
| `orbit` | the orbit (`a`, `r0`, `M`) |
| `ms` | the modes summed, ascending and contiguous |
| `order` | the regularization order subtracted per mode (4 for the CRP) |
| `Fr_modes` | the raw `F_r^m` (the UNREGULARIZED sequence, `∼ m^{-order}`) |
| `reg_terms` | the subtracted `F^m_{r[order]}`, `m = 0` ALREADY halved ([`mmode_reg_term_Fr`](@ref)) |
| `terms` | `Fr_modes − reg_terms` — the regularized per-m contributions, `∼ m^{-(order+2)}` |
| `partial` | `Σ terms` — the accelerated sum truncated at `m_max` |
| `tail` | the [`ModeSumTail`](@ref) remainder estimate |
| `Fr` | `partial + tail.value` — **THE NUMBER** |
| `truncation_error` | `|terms[end]|`, the CRP's omit-the-last-mode estimator (`main.tex:722`) |
| `resolution_errors` | per-mode resolution errors `|Δ_m|` (empty when the input carried none) |
| `resolution_error` | the CRP's own definition (`main.tex:722`) applied to the COMPLETE TOTAL: `|F_r(L) − F_r(L−1)|`, i.e. the signed sum of the per-mode differences PLUS the change of the fitted tail, since the tail is re-fitted on the penultimate level's own terms |
| `resolution_error_tail` | `|tail(L) − tail(L−1)|` — how much of that is the tail model moving |
| `resolution_error_bound` | `Σ_m |Δ_m| + resolution_error_tail` — the incoherent bound, always ≥ `resolution_error` |
| `resolution_includes_tail` | `false` only when the penultimate level admits no tail fit (then `resolution_error` degrades to the partial-sum proxy `|Σ_m Δ_m|` and says so) |
| `reduction` | which per-element reduction produced the modes (`:trace`, `:mean`, `:direct`) |
| `levels` | how many AMR levels each mode had (`0` when unknown) |

THE DISSIPATIVE COMPANION, filled when the caller supplies `Ft_modes` (both solved routes do):

| field | what it is |
|---|---|
| `Ft_modes` | the raw `F_t^m` (empty when not supplied) |
| `Ft` | `Σ_m F_t^m` — NO regularization and NO tail: the dissipative mode sum converges GEOMETRICALLY (measured ratio ≈ 0.19 at `a = 0`, `r₀ = 10M`), so truncation at `m_max = 20` is far below round-off |
| `Fphi` | `−Ft/Ω`, from the circular-orbit identity `F_t + Ω F_φ = 0` (`Warburton.tex:134`) |
| `Ft_truncation_error` | `|Ft_modes[end]|`, the same omit-the-largest-`m` estimator |
| `Ft_resolution_errors` / `Ft_resolution_error` | per mode, and the coherent `|Σ_m Δ_m^t|` |

`Ft` is a SECOND, INDEPENDENT external gate on the same pipeline: Warburton–Barack print
`(M²/q²)F_t` to nine significant figures (Table II, `Warburton.tex:580`) and the sum needs none
of the Heffernan machinery, so a disagreement there localizes to the extraction/solve while a
disagreement in `Fr` alone localizes to the regularization and the tail.

`Fr` is the normalized `(M²/q²)F_r` of the Warburton–Barack tables directly (`q = M = 1`
package-wide); there is no sign map and no unit conversion anywhere downstream.
"""
struct TotalSelfForce{T<:AbstractFloat}
    orbit::CircularOrbit{T}
    ms::Vector{Int}
    order::Int
    Fr_modes::Vector{T}
    reg_terms::Vector{T}
    terms::Vector{T}
    partial::T
    tail::ModeSumTail{T}
    Fr::T
    truncation_error::T
    resolution_errors::Vector{T}
    resolution_error::T
    resolution_error_tail::T
    resolution_error_bound::T
    resolution_includes_tail::Bool
    reduction::Symbol
    levels::Int
    Ft_modes::Vector{T}
    Ft::T
    Fphi::T
    Ft_truncation_error::T
    Ft_resolution_errors::Vector{T}
    Ft_resolution_error::T
end

function Base.show(io::IO, ::MIME"text/plain", r::TotalSelfForce)
    println(io, "TotalSelfForce  a = ", r.orbit.a, "  r0 = ", r.orbit.r0,
            "  m = ", first(r.ms), ":", last(r.ms), "  order = ", r.order)
    println(io, "  (M²/q²)F_r        = ", r.Fr)
    println(io, "  Σ_{m≤m_max} terms = ", r.partial)
    println(io, "  tail (", r.tail.method, ")   = ", r.tail.value, " ± ", r.tail.uncertainty)
    println(io, "  truncation error  = ", r.truncation_error)
    isempty(r.Ft_modes) || println(io, "  (M²/q²)F_t        = ", r.Ft,
                                   "   F_φ = ", r.Fphi)
    print(io,   "  resolution error  = ", r.resolution_error,
          r.resolution_includes_tail ? " (tail part " : " (PARTIAL-SUM PROXY, tail part ",
          r.resolution_error_tail, ", bound ", r.resolution_error_bound,
          ", reduction :", r.reduction, ")")
end

Base.show(io::IO, r::TotalSelfForce) =
    print(io, "TotalSelfForce(a = ", r.orbit.a, ", r0 = ", r.orbit.r0, ", Fr = ", r.Fr, ")")

"""
    total_self_force(orbit, Fr_modes::AbstractVector{<:Real}; ms = 0:length(Fr_modes)-1,
                     order = 4, tail = :heffernan, resolution_deltas = nothing, window = 6,
                     headroom = 3, reduction = :direct, levels = 0) -> TotalSelfForce
    total_self_force(set::ModeSet; ...)                              # row 47's fan-out
    total_self_force(hists::AbstractVector{<:AmrHistory}; ...)       # row 43, one per mode
    total_self_force(orbit, config::RunConfig; m_max = 20, pool = nothing,
                     parallel = true, ...)                           # solve, then sum

`F_r`, the total radial scalar self-force, as the CRP computes it (`main.tex:712`,
`:722-723`):

    F_r  =  Σ_{m=0}^{m_max} [ F_r^m − F^m_{r[order]} ]  +  tail
         =  Σ_{m=0}^{m_max}   F_r^m  +  mmode_reg_tail_Fr(order, orbit, m_max+1)  +  tail

The two lines are the SAME NUMBER — an identity, because the doubly-infinite parameter sum
vanishes (see this file's header, tail (T1)) — and both are computed and gated equal.

The `tail` keyword selects the model for the REMAINDER of the accelerated sum,
`Σ_{m>m_max}[F_r^m − F^m_{r[order]}]` (tail (T2)); see [`mode_sum_tail`](@ref) for `:heffernan`
(default), `:fit` and `:none`.

`resolution_deltas[k]` is the SIGNED difference `F_r^{m_k}(L) − F_r^{m_k}(L−1)` between the two
highest AMR levels (the `ModeSet` / `AmrHistory` methods supply it automatically). The sign
matters: its MAGNITUDE gates mode admission in the tail fit and propagates into the tail
amplitude, while the SIGNED values reconstruct the penultimate level's terms — and therefore
the penultimate level's own tail — so that `resolution_error` is the difference of two COMPLETE
totals rather than of two partial sums.

# UNITS AND CONVENTIONS

`q = M = 1`, so `Fr` IS the printed `(M²/q²)F_r` of Warburton–Barack. Modes are FOLDED `m ≥ 0`
with the ×2 real part inside them; retrograde orbits are `a < 0`, never `m < 0`. The `m = 0`
halving of the subtraction term is owned by [`mmode_reg_term_Fr`](@ref) and is applied exactly
once, inside [`regularized_mode_force`](@ref) — **never re-apply it in a caller**.

# THE SOLVING METHOD

`total_self_force(orbit, config)` runs the row-47 fan-out ([`solve_modes_parallel`](@ref),
or [`solve_modes_serial`](@ref) with `parallel = false`) over `0:m_max` and reduces the result.
Pass a `pool` from [`mode_pool`](@ref) and HOLD IT across a whole sweep of orbits — a fresh
worker pays its entire JIT before its first solve. Any mode that fails is an error, not a
silently shorter sum. Extra keywords (`iterations`, `timeout_s`, `log`, …) go to the fan-out.

The per-mode reduction available from a `ModeSet` is the one-sided `:trace` of the sealed
per-level diagnostics: a `ModeResult` deliberately does not ship a live `ModeSolution`, so the
spectre-parity `:mean` is not reachable from the fan-out. `reduction` records which one the
returned numbers are; how far apart the two reductions sit is a measurement of one mesh
sequence and lives in `docs/src/capabilities/scalarselfforce.md`, not here.
"""
function total_self_force(orbit::CircularOrbit, Fr_modes::AbstractVector;
                          ms = 0:(length(Fr_modes) - 1), order::Integer = 4,
                          tail::Symbol = :heffernan, resolution_deltas = nothing,
                          Ft_modes = nothing, Ft_resolution_deltas = nothing,
                          window::Integer = 6, headroom::Real = 3,
                          reduction::Symbol = :direct, levels::Integer = 0)
    mv = collect(Int, ms)
    length(mv) == length(Fr_modes) || throw(ArgumentError(
        "total_self_force: got $(length(mv)) modes and $(length(Fr_modes)) values"))
    isempty(mv) && throw(ArgumentError("total_self_force: no modes to sum"))
    mv == collect(minimum(mv):maximum(mv)) || throw(ArgumentError(
        "total_self_force: the modes must be a contiguous ascending run 0:m_max " *
        "(the mode sum has no gaps), got $mv"))
    first(mv) == 0 || throw(ArgumentError(
        "total_self_force: the sum starts at m = 0 (the unfolded Eq. (3.10) branch), " *
        "got first mode $(first(mv))"))
    resolution_deltas === nothing || length(resolution_deltas) == length(Fr_modes) ||
        throw(ArgumentError("total_self_force: got $(length(resolution_deltas)) " *
                            "resolution deltas and $(length(Fr_modes)) values"))
    T = float(promote_type(typeof(orbit.r0), real(eltype(Fr_modes))))
    orb = CircularOrbit(T(orbit.a), T(orbit.r0))
    fr = T[T(x) for x in Fr_modes]
    # row 45 owns the m = 0 halving; row 44 owns the subtraction. Both are called, never redone.
    rt = T[mmode_reg_term_Fr(order, orb, m) for m in mv]
    tm = T[regularized_mode_force(orb, mv[k], fr[k]; order = order) for k in eachindex(mv)]
    partial = sum(tm)
    dl = resolution_deltas === nothing ? T[] : T[T(d) for d in resolution_deltas]
    t = mode_sum_tail(orb, mv, tm; method = tail, order = order, window = window,
                      errors = isempty(dl) ? nothing : abs.(dl), headroom = headroom)
    # ---- the CRP's resolution measure applied to the COMPLETE total (main.tex:722) --------
    # "the difference between the highest and second-highest AMR level" is a statement about
    # THE NUMBER, and the number is `partial + tail`. The penultimate level's per-mode values
    # are `Fr_modes − Δ_m` (that is what the signed deltas ARE), so its terms, its own tail fit
    # and its own total are all reconstructible here without a further solve — and they must
    # be, because the tail is FITTED on the terms and therefore moves with them. Reporting
    # `|Σ Δ_m|` alone would be a partial-sum proxy.
    #
    # CORRELATION ASSUMPTION, stated rather than hidden: `resolution_error` is the COHERENT
    # difference — the per-mode errors are treated as one correlated realisation of the same
    # discretization, which is what they are (one mesh sequence), so they are allowed to
    # cancel. `resolution_error_bound` makes the opposite, worst-case assumption (every mode's
    # error aligned) and adds magnitudes. The truth is bracketed by the two; the CRP's own
    # definition is the coherent one.
    res, res_tail, bnd, with_tail = T(NaN), T(NaN), T(NaN), false
    if !isempty(dl)
        prev_fr = fr .- dl
        prev_tm = T[regularized_mode_force(orb, mv[k], prev_fr[k]; order = order)
                    for k in eachindex(mv)]
        prev_partial = sum(prev_tm)
        prev_tail = try
            mode_sum_tail(orb, mv, prev_tm; method = tail, order = order, window = window,
                          errors = abs.(dl), headroom = headroom).value
        catch err
            err isa ArgumentError || rethrow()
            nothing                       # the penultimate level admits no fit of its own
        end
        if prev_tail === nothing
            res, res_tail, bnd = abs(sum(dl)), T(NaN), sum(abs, dl)
        else
            res_tail = abs(t.value - prev_tail)
            res = abs((partial + t.value) - (prev_partial + prev_tail))
            bnd = sum(abs, dl) + res_tail
            with_tail = true
        end
    end
    # ---- the dissipative companion: a mode sum with NO regularization and NO tail ---------
    ft = Ft_modes === nothing ? T[] : T[T(x) for x in Ft_modes]
    isempty(ft) || length(ft) == length(mv) || throw(ArgumentError(
        "total_self_force: got $(length(ft)) F_t modes and $(length(mv)) modes"))
    ftd = Ft_resolution_deltas === nothing ? T[] : T[T(x) for x in Ft_resolution_deltas]
    isempty(ftd) || length(ftd) == length(mv) || throw(ArgumentError(
        "total_self_force: got $(length(ftd)) F_t deltas and $(length(mv)) modes"))
    Ft = isempty(ft) ? T(NaN) : sum(ft)
    Fphi = isempty(ft) ? T(NaN) : -Ft / T(orb.Omega)
    Ft_tr = isempty(ft) ? T(NaN) : abs(ft[end])
    Ft_res = isempty(ftd) ? T(NaN) : abs(sum(ftd))
    return TotalSelfForce{T}(orb, mv, Int(order), fr, rt, tm, partial, t, partial + t.value,
                             abs(tm[end]), abs.(dl), res, res_tail, bnd, with_tail,
                             reduction, Int(levels),
                             ft, Ft, Fphi, Ft_tr, abs.(ftd), Ft_res)
end

# ------------------------------------------------------- from solved sweeps ---------------
# The final value and the SIGNED last successive difference of one mode's `F_r^m` sequence.
# The sequence itself comes from the SEALED row-44 LEVEL ROUTE
# (`mode_self_force_history(orbit, config, m, levels; x_p)`), which is the same internal
# function the `AmrHistory` methods use — so nothing here can drift from row 44, and the
# angular chart, the radial-chart contract and the mode-number check are all its business.
# `config` MUST be the per-mode configuration that solved `m` (`mode_config(config, m)`).
function _ms_mode_value(orbit::CircularOrbit, config::RunConfig, m::Integer, levels, x_p;
                        phi_p::Real = 0)
    length(levels) >= 2 || throw(ArgumentError(
        "total_self_force: mode $m has $(length(levels)) AMR level(s); the resolution " *
        "error IS the difference between the highest and second-highest (main.tex:722)"))
    h = mode_self_force_history(orbit, config, m, levels; x_p = x_p, phi_p = phi_p)
    vr = [f.Fr for f in h]
    vt = [f.Ft for f in h]
    # the per-level reduction is the SEALED row-44 route's business, not this node's: report
    # what it actually used rather than asserting a constant, so a change of its default is a
    # re-measurement here and never a silently mislabelled table.
    return (Fr = vr[end], dFr = vr[end] - vr[end - 1],
            Ft = vt[end], dFt = vt[end] - vt[end - 1], reduction = h[end].reduction)
end

"""
    total_self_force(set::ModeSet; kwargs...) -> TotalSelfForce

The mode sum of a completed row-47 fan-out. Every requested mode must be present — a
`ModeSet` with failures throws rather than summing a shorter series — and the modes must be a
contiguous `0:m_max`. The per-mode resolution deltas come from the sealed per-level
diagnostics of each `ModeResult`, so the total's `resolution_error` needs no extra solve.
`reduction` is `:trace` (see the main docstring).
"""
function total_self_force(set::ModeSet; phi_p::Real = 0, kwargs...)
    isempty(set.failures) || throw(ArgumentError(
        "total_self_force: this ModeSet has failed modes $( [f.m for f in set.failures] ) — " *
        "a mode sum with a hole is not a shorter mode sum"))
    ms = sort(collect(Int, set.ms))
    # `mode_config(config, m)` is the per-mode configuration the fan-out itself built; the
    # sealed level route checks the mode number and the radial chart against it.
    vd = [_ms_mode_value(set.orbit, mode_config(set.config, m), m, set[m].levels, set[m].x_p;
                         phi_p = phi_p) for m in ms]
    vals, dels = [x.Fr for x in vd], [x.dFr for x in vd]
    tvals, tdels = [x.Ft for x in vd], [x.dFt for x in vd]
    red = unique(x.reduction for x in vd)
    length(red) == 1 || throw(ArgumentError(
        "total_self_force: the modes were reduced differently ($(red)); a mode sum must be " *
        "homogeneous"))
    lv = maximum(length(set[m].levels) for m in ms)
    return total_self_force(set.orbit, vals; ms = ms, resolution_deltas = dels,
                            Ft_modes = tvals, Ft_resolution_deltas = tdels,
                            reduction = only(red), levels = lv, kwargs...)
end

"""
    total_self_force(hists::AbstractVector{<:AmrHistory}; kwargs...) -> TotalSelfForce

The mode sum of a set of row-43 AMR sweeps, one per mode (the serial route; the fan-out's
`ModeSet` is the parallel one). The histories may come in any order; their `m` values must
form a contiguous `0:m_max`. Per-mode values and resolution deltas are read through
[`mode_self_force_history`](@ref), so this route uses whatever per-level reduction that sealed
route reports — the same one the `ModeSet` route gets — and the two are directly comparable.
`reduction` records which one it was; a set of modes that were not reduced alike throws.
"""
function total_self_force(hists::AbstractVector{<:AmrHistory}; kwargs...)
    isempty(hists) && throw(ArgumentError("total_self_force: no histories"))
    orb = first(hists).orbit
    p = sortperm([h.m for h in hists])
    ms = [hists[k].m for k in p]
    seqs = map(p) do k
        h = mode_self_force_history(hists[k])
        length(h) >= 2 || throw(ArgumentError(
            "total_self_force: mode $(hists[k].m) has $(length(h)) AMR level(s); the " *
            "resolution error IS the difference between the highest and second-highest " *
            "(main.tex:722)"))
        ([f.Fr for f in h], [f.Ft for f in h], h[end].reduction)
    end
    vals = [v[1][end] for v in seqs]
    dels = [v[1][end] - v[1][end - 1] for v in seqs]
    tvals = [v[2][end] for v in seqs]
    tdels = [v[2][end] - v[2][end - 1] for v in seqs]
    lv = maximum(length(v[1]) for v in seqs)
    red = unique(v[3] for v in seqs)
    length(red) == 1 || throw(ArgumentError(
        "total_self_force: the modes were reduced differently ($(red)); a mode sum must be " *
        "homogeneous"))
    return total_self_force(orb, vals; ms = ms, resolution_deltas = dels,
                            Ft_modes = tvals, Ft_resolution_deltas = tdels,
                            reduction = only(red), levels = lv, kwargs...)
end

"""
    total_self_force(orbit::CircularOrbit, config::RunConfig; m_max = last(CRP_PAPER_MMODES),
                     pool = nothing, parallel = true, phi_p = 0, kwargs...) -> TotalSelfForce

Solve `m = 0:m_max` and sum them: the end-to-end number. `config` is the m-INDEPENDENT
production configuration (`crp_paper_config`); each mode's own configuration is derived by the
sealed [`mode_config`](@ref), so the `m = 0` outer-boundary branch is handled once, in one
place. Keywords not listed here (`order`, `tail`, `window`, `headroom` for the sum;
`iterations`, `timeout_s`, `log`, `nworkers`, … for the fan-out) are routed to the sum and to
the solver respectively.

HOLD ONE POOL ACROSS A SWEEP OF ORBITS — see [`mode_pool`](@ref); building a pool per point
pays every worker's JIT again.
"""
function total_self_force(orbit::CircularOrbit, config::RunConfig;
                          m_max::Integer = last(CRP_PAPER_MMODES), pool = nothing,
                          parallel::Bool = true, phi_p::Real = 0,
                          order::Integer = 4, tail::Symbol = :heffernan,
                          window::Integer = 6, headroom::Real = 3, kwargs...)
    m_max >= 1 || throw(ArgumentError(
        "total_self_force: m_max must be ≥ 1 to have a mode sum at all, got $m_max"))
    # the fan-out solves at the package-wide φ_p = 0 (row 47 has no φ_p knob); `phi_p` here is
    # therefore the EXTRACTION phase and must stay 0 for a spectre/CRP-comparable number.
    set = parallel ?
        solve_modes_parallel(orbit, config, 0:Int(m_max); pool = pool, kwargs...) :
        solve_modes_serial(orbit, config, 0:Int(m_max); kwargs...)
    return total_self_force(set; phi_p = phi_p, order = order, tail = tail, window = window,
                            headroom = headroom)
end
