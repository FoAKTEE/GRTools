# QNM/qnm_frequency.jl — Kerr QNM frequency ω_{s,l,m,n}(a): the angular spheroidal eigenvalue
# λ(s,l,m,aω) is inserted into the radial Leaver CF's c₃ and the single CF condition is
# root-found for ω (Leaver 1985 eqs 12-14; the n-th overtone via the r=n inversion).
# Coefficients per Paper I (2401.15516) app:Leaver tex:1156-1174 with σ₊ = (ωr₊ − am/2)/√(1−a²)
# [tex:329] — with TWO constraints the printed paper does not show:
#  • c₄: the paper's −(2s+1)iω [tex:1172] is a typo; the correct Cook–Zalutskiy form is
#    −4(s+1)iω (1410.7698 Eq.31; ref-code/qnm radial.py D[4]) — the wrong form misses the a=0
#    anchor by 1.2e-2. The test gates c₀..c₄ elementwise against CZ D[0..4].
#  • λ in c₃ is the Harmonics PUBLIC (Teukolsky) eigenvalue, NOT the raw CZ separation
#    constant A = λ + 2mγ − γ²: they coincide at a=0 only, and the raw-A insertion errs by
#    5e-2 at a≠0 (the s=±2 independence gate pins the choice).
# {T}-generic: Float64 λ via the spectral path, BigFloat λ via the Harmonics Leaver backend.

"""
    _qnm_leaver_c(s, m, a, ω, λ) -> (c0, c1, c2, c3, c4)

The five Leaver constants (Paper I eq:leaver_alpha_beta_gamma; double-entry transcription in the
file header). `λ` is the BHPT/Teukolsky spheroidal eigenvalue (Harmonics convention). {T}-generic:
`a` real, `ω`/`λ` complex; returns `Complex{T}` constants.
"""
@inline function _qnm_leaver_c(s::Integer, m::Integer, a::Real, ω::Number, λ::Number)
    T = float(real(promote_type(typeof(a), real(typeof(ω)), real(typeof(λ)))))
    aT = T(a); mT = T(m); sT = T(s); i = im
    b = sqrt(one(T) - aT^2)                         # b = √(1−a²)          [tex:1174]
    rp = one(T) + b                                 # r₊ = 1 + √(1−a²)     [tex:186]
    σp = (ω * rp - aT * mT / 2) / b                 # σ₊ = (ω r₊ − am/2)/b  [tex:329]
    c0 = 1 - sT - 2i * σp                           #                       [tex:1167]
    c1 = -4 + 4i * σp + 4i * rp * ω                 #                       [tex:1168]
    c2 = sT + 3 - 4i * ω - 2i * σp                  #                       [tex:1169]
    c3 = 4 * ω^2 * (2 + 2b) - λ - 4 * aT * mT * ω - sT - 1 +
         2 * (1 + b) * i * ω + (8 * ω + 2i) * σp    #                       [tex:1170-1171]
    # c₄ = CZ D[4] (typo-corrected; header [TYPO-CORRECTION]). Paper tex:1172 writes −(2s+1)iω;
    # the correct Cook-Zalutskiy coefficient is −4(s+1)iω (qnm/radial.py D[4]=α(α−γ+1)).
    c4 = sT + 1 - 4 * (sT + 1) * i * ω - (8 * ω + 2i) * σp
    return (c0, c1, c2, c3, c4)
end

"""
    qnm_cf_residual(s, l, m, a, ω, λ; ninv, rtol) -> Complex

Leaver radial continued-fraction residual at the `ninv`-th inversion, with the angular spheroidal
eigenvalue `λ` inserted into `c₃`. A QNM ω is a root: `qnm_cf_residual(...) == 0`.

Assembly = `β_r − γ_r·conv₁ − α_r·conv₂`, the `r`-th-inversion of Leaver's CF (Leaver 1985 eqs
12-14; identical to `qnm.radial.leaver_cf_trunc_inversion`, ref-code/qnm/qnm/radial.py). The FINITE
"down" leg `conv₁` (indices r-1…0) reuses `NumericsCore.cf_finite` (it fits exactly — a finite
K-CF). The "up" tail `conv₂` (indices r+1…∞) uses the ADAPTIVE truncated backward recurrence
`_qnm_uptail` rather than `NumericsCore.cf_convergent`: at the Schwarzschild algebraically-special
overtone (l=2 n=8, ω≈−2i) the minimal-solution tail converges too slowly for `cf_convergent`'s
iterate-to-`rtol`-or-error contract, whereas the truncated inversion (the h5 oracle's OWN method)
always returns a value. For well-behaved modes `conv₂ == −cf_convergent(...)` to precision (test
cross-check). {T}/BigFloat-inherited throughout.
"""
function qnm_cf_residual(s::Integer, l::Integer, m::Integer, a::Real, ω::Number, λ::Number;
                         ninv::Integer=0, rtol::Real=_qnm_default_rtol(ω), Nmax::Integer=100_000)
    c0, c1, c2, c3, c4 = _qnm_leaver_c(s, m, a, ω, λ)
    αn(n) = n^2 + (c0 + 1) * n + c0
    βn(n) = -2 * n^2 + (2 + c1) * n + c3
    γn(n) = n^2 + (c2 - 3) * n + c4 - c2 + 2
    r = ninv
    # finite "down" leg (indices r-1 … 0), reindexed as in the angular Leaver front door.
    # cf_finite here equals −γ_r·conv₁ (qnm's conv₁ = the forward down recurrence).
    down = r >= 1 ? cf_finite(k -> -αn(r - k) * γn(r - k + 1), k -> βn(r - k), 1, r - 1) : zero(c3)
    # convergent "up" tail (indices r+1 … ∞) — the minimal-solution ratio CF, truncated-adaptive
    conv2 = _qnm_uptail(αn, βn, γn, r, typeof(complex(one(c3))); rtol=rtol, Nmax=Nmax)
    return βn(r) + down - αn(r) * conv2
end

"""
    _qnm_uptail(αn, βn, γn, r, ::Type{CT}; rtol, Nmax) -> conv₂::CT

The Leaver "up" tail ratio `conv₂ = γ_{r+1}/(β_{r+1} − α_{r+1}γ_{r+2}/(β_{r+2} − ⋯))` by backward
recurrence from a truncation depth `N` with seed `r_N = 0` (qnm's `leaver_cf_trunc_inversion`). `N`
is doubled (from `max(2r+20, 60)`) until `conv₂` stabilises to `rtol` or hits `Nmax` — so a
BigFloat run automatically takes enough terms to reach its floor. Never errors.
"""
function _qnm_uptail(αn, βn, γn, r::Integer, ::Type{CT}; rtol::Real, Nmax::Integer) where {CT}
    RT = real(CT)
    N = max(2r + 20, 60)
    prev = CT(Inf)
    conv2 = zero(CT)
    while true
        c2 = zero(CT)                                   # r_N seed = 0
        @inbounds for i in N:-1:(r + 1)
            c2 = γn(i) / (βn(i) - αn(i) * c2)
        end
        conv2 = c2
        if isfinite(abs(conv2)) && abs(conv2 - prev) <= rtol * max(abs(conv2), one(RT))
            break
        end
        prev = conv2
        N >= Nmax && break
        N = min(2N, Nmax)
    end
    return conv2
end

@inline _qnm_default_rtol(ω::Number) = eps(float(real(typeof(ω)))) * 8

"""
    _qnm_secant(f, x0, x1; rtol, maxiter) -> Complex

Derivative-free complex secant (kept local: the CF residual is complex-valued, so the tolerance
test compares the real step `|Δω|` to a real bound — the generic `NumericsCore.secant`'s complex
`eps(T)` cannot express this, exactly the reason Harmonics keeps `_swsh_secant` local).
"""
function _qnm_secant(f, x0::Number, x1::Number; rtol::Real, maxiter::Integer=200)
    f0 = f(x0); f1 = f(x1)
    x0c, x1c = x0, x1
    for _ in 1:maxiter
        denom = f1 - f0
        iszero(denom) && return x1c
        x2 = x1c - f1 * (x1c - x0c) / denom
        abs(x2 - x1c) <= rtol * abs(x2) && return x2
        x0c, f0, x1c, f1 = x1c, f1, x2, f(x2)
    end
    return x1c
end

"""
    _qnm_angular_lambda(s, l, m, a, ω) -> Complex

The spheroidal eigenvalue λ(s,l,m,γ=aω) via Harmonics (complex γ, `method=:Auto`: ComplexF64 →
spectral, Complex{BigFloat} → Leaver CF). Handles γ=0 exactly (returns l(l+1)−s(s+1)).
"""
@inline function _qnm_angular_lambda(s::Integer, l::Integer, m::Integer, a::Real, ω::Number)
    γ = complex(a * ω)
    return SpinWeightedSpheroidalEigenvalue(s, l, m, γ)
end

"""
    _qnm_eikonal_seed(l, n, ::Type{T}) -> Complex{T}

Geometric-optics (eikonal) QNM seed for a=0: ω ≈ ((l+½) − i(n+½))/(3√3) (M=1). A DERIVED starting
guess (not a table lookup), used only when no explicit `seed`/continuation value is supplied; the
returned ω is always refined by the CF secant and gated by the CF residual.
"""
@inline function _qnm_eikonal_seed(l::Integer, n::Integer, ::Type{T}) where {T}
    base = (T(l) + T(1) / 2 - im * (T(n) + T(1) / 2)) / (3 * sqrt(T(3)))
    return base
end

"""
    qnm_frequency(s, l, m, n, a; seed=nothing, seed_A=nothing, ninv=n, rtol, maxiter=200,
                  seed_step=..., return_residual=false, backend=:leaver) -> ω  (or (ω, |residual|))

Kerr quasinormal-mode frequency `ω_{s,l,m,n}(a)` (M=1, prograde/physical branch; the mirror mode is
`−conj(ω_{s,l,−m,n})`). Root-finds the Leaver radial CF condition (`qnm_cf_residual`) at the `ninv`-th
inversion (default = overtone `n`), with the complex-γ spheroidal eigenvalue coupled in.

`seed` — complex initial guess (χ-continuation value or a published value used ONLY as a starting
guess, always refined+gated). Defaults to the eikonal estimate `_qnm_eikonal_seed`.

`return_residual=true` returns `(ω, |qnm_cf_residual(…,ω,λ(ω))|)` — the CF defect at the root, the
BF-256 precision-convergence gate diagnostic.

`backend` — the angular-λ engine (all root-find the same `qnm_cf_residual`):
  * `:leaver` (default, the reference) — λ from `Harmonics.SpinWeightedSpheroidalEigenvalue`.
  * `:fast` — λ from shift-invert inverse iteration on the pentadiagonal, warm-started from
    `seed_A`. A Float64 speed path by design: a `BigFloat` `a` routes to `:leaver`.
  * `:hyperboloidal` — the collocation QEP (an independent formulation; see qnm_hyperboloidal.jl).

{T}/BigFloat-inherited: pass `a::BigFloat` and a `Complex{BigFloat}` `seed` for the BF-256 path.
"""
function qnm_frequency(s::Integer, l::Integer, m::Integer, n::Integer, a::Real;
                       seed=nothing, seed_A=nothing, ninv::Integer=n,
                       rtol::Real=eps(float(real(typeof(a)))) * 64,
                       maxiter::Integer=200, seed_step=nothing, n_coll::Integer=48,
                       return_residual::Bool=false, backend::Symbol=:leaver)
    (l < abs(s) || abs(m) > l) &&
        throw(ArgumentError("invalid (s,l,m)=($s,$l,$m): need l>=|s| and |m|<=l"))
    T = float(real(typeof(a)))
    # n_coll=48 reaches ≤1e-13 for |a|≲0.9 prograde; the retrograde near-extremal corner
    # (e.g. (2,-2,a≳0.95)) needs n_coll≳56.
    if backend === :hyperboloidal
        return _qnm_frequency_hyperboloidal(s, l, m, n, a; seed=seed, rtol=rtol, n_coll=n_coll,
                                            return_residual=return_residual)
    end
    # :fast is Float64-only; BigFloat routes to :leaver (tier split).
    if backend === :fast && !(T <: BigFloat)
        return _qnm_frequency_fast(s, l, m, n, a; seed=seed, seed_A=seed_A, ninv=ninv, rtol=rtol,
                                   maxiter=maxiter, seed_step=seed_step, return_residual=return_residual)
    end
    backend in (:leaver, :fast, :hyperboloidal) ||
        throw(ArgumentError("unknown backend $backend (expected :leaver, :fast, or :hyperboloidal)"))
    ω0::Complex{T} = seed === nothing ? _qnm_eikonal_seed(l, n, T) : Complex{T}(seed)
    G(ω) = qnm_cf_residual(s, l, m, a, ω, _qnm_angular_lambda(s, l, m, a, ω);
                           ninv=ninv, rtol=min(rtol, _qnm_default_rtol(ω)))
    # second secant seed: a small complex perturbation of ω0 (scaled to ω0's magnitude)
    δ = seed_step === nothing ? Complex{T}(T(1) / 10^4 + im * T(1) / 10^4) : Complex{T}(seed_step)
    ω1 = ω0 * (1 + δ)
    ω = _qnm_secant(G, ω0, ω1; rtol=rtol, maxiter=maxiter)
    return return_residual ? (ω, abs(G(ω))) : ω
end
