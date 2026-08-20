# Tail-rate estimators and the RT-chart excision policy for the generic-spin sector.
#
# Promoted from the legacy arXiv:1301.1591 reproduction drivers (emri-disk
# `results/dg-teukolsky-spin/scripts/{driver_common,fit_rates}.jl`, campaigns A/B/C): each function
# here PRODUCES a published number — the measured late-time tail exponents of the Kerr s = −2 tails
# and mode-scattering matrix, and the self-convergence orders of the DG panels — so it belongs in the
# package with the estimator it generalises (`fit_tail_window`, common.jl) rather than in a script.
#
# One least-squares log-log kernel only: every estimator below reduces to `fit_tail_window`.

# ---- oscillation-robust envelope fit ------------------------------------------------------------
"""
    fit_tail_envelope(τs, ψ, win; floor = 1e-14) -> (n, npeaks)

Late-time exponent fitted THROUGH THE LOCAL MAXIMA of `|ψ|` inside the window
`win = (τlo, τhi)`. Intended for a tail riding a numerical oscillation — the diagnosed
κ-mediated |s| = 2 oscillation of the RT-chart sector, whose deep nulls pull a plain log-log
fit off the envelope: on an exact `A τⁿ (½+½cos)²` record this returns `n` while
[`fit_tail_window`](@ref) does not.

A MONOTONE record has no interior maxima, so `(NaN, 0)` is returned — as it is whenever fewer
than four maxima clear `floor`. Never a fabricated number. `npeaks` is the number of maxima used.
"""
function fit_tail_envelope(τs, ψ, win; floor = 1e-14)
    t = Float64[]; a = Float64[]
    for (τ, v) in zip(τs, ψ)
        win[1] <= τ <= win[2] && abs(v) > floor && (push!(t, τ); push!(a, abs(v)))
    end
    length(a) < 6 && return (NaN, 0)
    pk = [i for i in 2:length(a)-1 if a[i] >= a[i-1] && a[i] >= a[i+1]]
    length(pk) < 4 && return (NaN, length(pk))
    return (fit_tail_window(t[pk], a[pk], (-Inf, Inf)), length(pk))
end

# ---- adaptive window: fit until the record turns up into its floor -------------------------------
"""
    fit_tail_adaptive(τs, ψ; t0 = 800.0, floor = 1e-13) -> n

Late-time exponent from `t0` until the series stops decreasing (turn-up into the numerical
floor, detected as `|ψ| > 1.5 |ψ_prev|`) or drops below `floor`. This is the least
floor-biased uniform window choice when the floor time is configuration-dependent — a fixed
window either clips the tail or fits through the floor. `NaN` with fewer than four samples.
"""
function fit_tail_adaptive(τs, ψ; t0 = 800.0, floor = 1e-13)
    t = Float64[]; a = Float64[]
    last = Inf
    for (τ, v) in zip(τs, ψ)
        τ < t0 && continue
        av = abs(v)
        (av < floor || av > 1.5 * last) && break
        push!(t, τ); push!(a, av); last = av
    end
    length(t) < 4 && return NaN
    return fit_tail_window(t, a, (-Inf, Inf))
end

"Window starts of the sliding-window scan used by [`fit_tail_sliding`](@ref) (legacy campaign choice)."
const TAIL_SLIDING_STARTS = (600.0, 750.0, 900.0, 1100.0, 1300.0, 1500.0, 1700.0)

"""
    fit_tail_sliding(τs, ψ; tf = maximum(τs), floor = 1e-13,
                     starts = TAIL_SLIDING_STARTS, span = 1.8, minwidth = 250.0) -> (plateau, nmin, nmax, nwin)

Scan of geometric windows `(w₀, min(span·w₀, tf))` over `starts`, keeping only windows whose
whole record stays above `floor`. `plateau` is the fit of the window pair with the smallest
consecutive difference (the flattest part of the scan — a window choice that is NOT
law-guided, so it cannot be tuned toward an expected exponent); `nmin`/`nmax` bracket the
window-to-window scatter and `nwin` counts the windows used. All `NaN`/`0` when no window
qualifies.
"""
function fit_tail_sliding(τs, ψ; tf = maximum(τs), floor = 1e-13,
                          starts = TAIL_SLIDING_STARTS, span = 1.8, minwidth = 250.0)
    fits = Tuple{Float64,Float64}[]
    for w0 in starts
        w1 = min(span * w0, tf)
        w1 - w0 < minwidth && continue
        sel = [j for j in eachindex(τs) if w0 <= τs[j] <= w1]
        isempty(sel) && continue
        minimum(abs.(ψ[sel])) < floor && continue
        n = fit_tail_window(τs, ψ, (w0, w1))
        isnan(n) || push!(fits, (w0, n))
    end
    isempty(fits) && return (NaN, NaN, NaN, 0)
    length(fits) == 1 && return (fits[1][2], fits[1][2], fits[1][2], 1)
    dmin, best = Inf, fits[end][2]
    for i in 1:length(fits)-1
        d = abs(fits[i+1][2] - fits[i][2])
        d < dmin && (dmin = d; best = fits[i+1][2])
    end
    ns = [f[2] for f in fits]
    return (best, minimum(ns), maximum(ns), length(fits))
end

# ---- self-convergence order ----------------------------------------------------------------------
"""
    fit_convergence_order(Ks, errs; lo = 1e-13, hi = 0.5) -> p

Observed order `p` of an error ladder `errs ∝ Ks^(-p)`: the negated least-squares log-log slope
over the points whose error lies in the ordinate window `(lo, hi)` — `lo` removes points already
at the round-off/self-convergence floor, `hi` the un-converged coarse end. `NaN` with fewer than
three usable points. Same kernel as [`fit_tail_window`](@ref) (abscissa `K`, ordinate `E`).
"""
function fit_convergence_order(Ks, errs; lo = 1e-13, hi = 0.5)
    x = Float64[]; y = Float64[]
    for (K, E) in zip(Ks, errs)
        lo < E < hi && (push!(x, float(K)); push!(y, float(E)))
    end
    length(x) < 3 && return NaN
    return -fit_tail_window(x, y, (-Inf, Inf))
end

# ---- RT-chart excision policy ---------------------------------------------------------------------
"""
    spin_excise_policy(M, a; frac = 0.7, cap = 0.05) -> (; excise, Rplus, Rminus)

Width of the excised strip inside the outer horizon for a [`SpinModesProblem`](@ref) at spin `a`:
`excise = min(cap, frac·(R₊ − R₋))` in RT coordinates.

The inner boundary must stay strictly BETWEEN the two horizons: inside `R₊` the trapped region
makes every characteristic speed ≤ 0 (the constructor asserts this), while past `R₋` (where Δ > 0
again) an outgoing characteristic reappears and the build fails. At moderate spin the horizon gap
is wide and the constant `cap` applies; near extremality the gap closes and `frac` keeps the inner
edge inside it. `Rminus = 0` at `a = 0` (no inner horizon).

The excised edge must also SURVIVE the mesh build: `build_mesh` snaps a forced coordinate to the
nearest uniform edge, so `excise < h/2` silently deletes the excised strip and re-pins the inner
boundary at `R₊`. Callers assert `p.mesh.edges[1] == R₊ − excise` after construction.
"""
function spin_excise_policy(M, a; frac = 0.7, cap = 0.05)
    Rp = rt_Rplus(M, a)
    rm = M - sqrt(M^2 - a^2)
    Rm = rm == 0 ? zero(Rp) : rt_R_of_r(rm)      # a = 0: no inner horizon, R₋ = 0
    return (excise = min(cap, frac * (Rp - Rm)), Rplus = Rp, Rminus = Rm)
end
