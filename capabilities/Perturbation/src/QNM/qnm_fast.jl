# QNM/qnm_fast.jl — the `:fast` Kerr-QNM backend: a native-Julia port of the `qnm` package's
# speed engineering (ref-code/qnm @f3abd18). Default stays `:leaver`; `:fast` is the speed option.
#
#  • Angular eigenpair by SHIFT-INVERT inverse iteration on the pentadiagonal matrix — only the
#    one eigenpair near the tracked seed is needed (the reference computes all N by dense geev).
#    The matrix comes from `Harmonics._swsh_spectral_setup` (band elements single-sourced in
#    Harmonics; ONE non-exported-symbol coupling — a rename breaks at load, test-caught).
#    Convention bridge: A_sep = −μ_sel − s(s+1), λ = A_sep − 2mγ + γ²; a seed A_seed maps to the
#    shift μ₀ = −A_seed − s(s+1).
#  • Spin-sequence continuation (spinsequence.py:228-301 semantics): quadratic extrapolation of
#    BOTH ω(a) and A(a) through the last 3 nodes, curvature-adaptive Δa = 0.05/√(max d²) clamped
#    [1e-5, 4e-3]; a non-converging polish throws (never a silent skip).
#  • The radial CF is REUSED from `qnm_cf_residual` (not forked); the angular solve dominates the
#    per-root cost, so no inlined Lentz variant. Adaptive l_max: grow the band until the tail
#    eigenvector component < `ltol`.
# `:fast` is a Float64 speed path by design (the banded-LU/RQI is itself {T}-generic); BigFloat
# callers are routed to `:leaver` at the front door.

# ─────────────────────────────────────────────────────────────────────────────────────────────
# Banded pentadiagonal LU (lower/upper bandwidth p=q=2), no pivoting, O(N).
# The matrix is stored dense (N≲25, ~5 KB) but the elimination touches ONLY the 5 bands, so the
# flop count is O(N)·O(pq) = O(N) — the genuine banded win over dense O(N³) `lu`/geev. Fill-in stays
# strictly inside the band (Golub & Van Loan Alg. 4.3.1). Returns the smallest |pivot| seen, so the
# caller can nudge an (intentionally near-singular, by shift-invert design) shift off an EXACT-zero
# pivot; inverse-iteration accuracy is otherwise insensitive to the ill-conditioning (Peters-
# Wilkinson). {T}-generic (Float64/Complex/BigFloat). Validated == dense `\\` in the tests.
# ─────────────────────────────────────────────────────────────────────────────────────────────
"""
    _qnm_penta_lu!(A) -> min|pivot|

In-place pentadiagonal LU (unit-lower, bandwidth 2 each side, no pivot). Overwrites `A` with its
`L\\U` band factors. Returns the minimum |diagonal pivot| encountered.
"""
function _qnm_penta_lu!(A::AbstractMatrix{T}) where {T}
    n = size(A, 1)
    mp = typemax(real(T))
    @inbounds for k in 1:(n - 1)
        piv = A[k, k]
        mp = min(mp, abs(piv))
        for i in (k + 1):min(k + 2, n)
            A[i, k] /= piv
            for j in (k + 1):min(k + 2, n)
                A[i, j] -= A[i, k] * A[k, j]
            end
        end
    end
    return min(mp, abs(A[n, n]))
end

"""
    _qnm_penta_solve!(A, b) -> b

Solve `A x = b` in place (`b` overwritten by `x`) given the band factors from `_qnm_penta_lu!`.
Forward (unit-lower, 2 subdiagonals) then back (upper, 2 superdiagonals) substitution, O(N).
"""
function _qnm_penta_solve!(A::AbstractMatrix{T}, b::AbstractVector{T}) where {T}
    n = length(b)
    @inbounds for i in 1:n                       # L y = b   (unit diagonal)
        s = b[i]
        for k in max(1, i - 2):(i - 1)
            s -= A[i, k] * b[k]
        end
        b[i] = s
    end
    @inbounds for i in n:-1:1                     # U x = y
        s = b[i]
        for k in (i + 1):min(i + 2, n)
            s -= A[i, k] * b[k]
        end
        b[i] = s / A[i, i]
    end
    return b
end

# ─────────────────────────────────────────────────────────────────────────────────────────────
# Shift-invert inverse iteration (Rayleigh-quotient iteration) for the ONE eigenpair
# of the complex-SYMMETRIC pentadiagonal `M` nearest a shift `μ₀`. RQI updates the shift to the
# Rayleigh quotient each step → cubic convergence; each step is one banded LU + one banded solve
# (O(N)). The bilinear form is the complex-SYMMETRIC `vᵀ·v` (transpose, NOT the Hermitian `v'·v`),
# matching `Harmonics._swsh_eigenvalue_dgamma`'s complex-symmetric left eigenvector. A near-singular
# (M−ρI) is expected AT convergence (ρ→eigenvalue); the direction stays accurate. Falls back to a
# dense pivoted `lu` solve only if the band solve produces a non-finite iterate.
# ─────────────────────────────────────────────────────────────────────────────────────────────
"""
    _qnm_shift_invert(M, μ0, v0; tol, maxiter) -> (μ, v, iters)

Rayleigh-quotient inverse iteration on the pentadiagonal `M` (dense-stored) from the initial shift
`μ0` and unit-ish vector `v0`. Returns the converged eigenvalue `μ`, complex-symmetric-unit
eigenvector `v` (`vᵀv=1`), and the iteration count. `M` is NOT modified (work copies are made).
"""
function _qnm_shift_invert(M::AbstractMatrix{CT}, μ0::Number, v0::AbstractVector;
                           tol::Real=1e-14, maxiter::Integer=50) where {CT}
    n = size(M, 1)
    v = CT.(v0)
    v ./= sqrt(sum(x -> x^2, v))                  # complex-symmetric normalisation (vᵀv=1)
    ρ = CT(μ0)
    A = Matrix{CT}(undef, n, n)
    y = Vector{CT}(undef, n)
    iters = 0
    for it in 1:maxiter
        iters = it
        @inbounds for j in 1:n, i in 1:n          # A = M − ρI
            A[i, j] = M[i, j] - (i == j ? ρ : zero(CT))
        end
        mp = _qnm_penta_lu!(A)
        copyto!(y, v)
        if mp < eps(real(CT))^2                   # shift landed on an exact pivot zero → nudge
            ρn = ρ * (1 + real(CT)(1e-13)) + real(CT)(1e-30)
            @inbounds for j in 1:n, i in 1:n
                A[i, j] = M[i, j] - (i == j ? ρn : zero(CT))
            end
            _qnm_penta_lu!(A)
        end
        _qnm_penta_solve!(A, y)
        if any(!isfinite, y)                       # rare: dense pivoted fallback
            y .= (M - ρ * I) \ v
            any(!isfinite, y) && break
        end
        nrm = sqrt(sum(x -> x^2, y))
        v .= y ./ nrm
        ρnew = (transpose(v) * M * v) / (transpose(v) * v)
        if abs(ρnew - ρ) <= tol * max(abs(ρnew), one(real(CT)))
            ρ = ρnew
            break
        end
        ρ = ρnew
    end
    return ρ, v, iters
end

# ─────────────────────────────────────────────────────────────────────────────────────────────
# The angular spheroidal eigenvalue via shift-invert, adaptive l_max.
# ─────────────────────────────────────────────────────────────────────────────────────────────
"""
    _qnm_angular_fast(s, l, m, γ, A_seed; ltol, growmax, tol, maxiter)
        -> (λ, A_sep, tailfrac, lmax_used, iters)

Spheroidal eigenvalue λ(s,l,m,γ=aω) and the raw separation constant `A_sep = −μ_sel − s(s+1)` via
shift-invert around the seed `A_seed` (from continuation / the previous ω-iterate).
`λ = A_sep − 2mγ + γ²` (the Harmonics public / Teukolsky convention `qnm_cf_residual` expects).

Adaptive l_max: starts from Harmonics' `:auto` band width and GROWS `numterms` (+4 per
step, up to `growmax`) until the tail eigenvector component `tailfrac = max(boundary |vᵢ|)/‖v‖∞`
drops below `ltol`. Returns `lmax_used = l + nUp` (the highest l' retained) and the RQI count.

γ=0 (a=0) is returned exactly (`λ = l(l+1)−s(s+1)`, `A_sep = λ`), no eigensolve. {T}-generic, but
the SPEED target is ComplexF64 (see the file header tier split).
"""
function _qnm_angular_fast(s::Integer, l::Integer, m::Integer, γ::Number, A_seed::Number;
                           ltol::Real=1e-9, growmax::Integer=60, tol::Real=1e-14,
                           maxiter::Integer=50)
    if iszero(γ)
        λ0 = l * (l + 1) - s * (s + 1)
        CT = complex(float(typeof(γ)))
        return (CT(λ0), CT(λ0), zero(real(CT)), l, 0)
    end
    CT = complex(float(typeof(γ)))
    RT = real(CT)
    numterms = :auto
    local λ, A_sep, tailfrac, lmax_used, iters, nt_used
    while true
        M0, nUp, nDown, _ = Harmonics._swsh_spectral_setup(s, l, m, γ, numterms)
        N = size(M0, 1)
        M = Matrix{CT}(M0)
        v0 = zeros(CT, N)
        v0[nDown + 1] = one(CT)                    # seed the l-th spherical component (physical branch)
        μ0 = -A_seed - s * (s + 1)
        μ, v, iters = _qnm_shift_invert(M, μ0, v0; tol=tol, maxiter=maxiter)
        E = -μ - s * (s + 1)                        # = A_sep (qnm/CZ raw separation constant)
        A_sep = E
        λ = E - 2 * m * γ + γ^2
        peak = maximum(abs, v)
        tailfrac = max(abs(v[N]), nDown > 0 ? abs(v[1]) : zero(RT)) / peak
        nt_used = nUp
        lmax_used = l + nUp
        if tailfrac <= ltol || nUp >= growmax
            break
        end
        numterms = nUp + 4
    end
    return (λ, A_sep, tailfrac, lmax_used, iters)
end

# ─────────────────────────────────────────────────────────────────────────────────────────────
# PRIORITY-2a — single-root polish. Complex secant on the SAME `qnm_cf_residual` (reused, not
# forked), with the angular λ per residual supplied by shift-invert WARM-STARTED from the previous
# ω-iterate's `A_sep` (the continuation idea at the inner scale). `NumericsCore.secant` is the same complex
# secant `:leaver` uses.
# ─────────────────────────────────────────────────────────────────────────────────────────────
"""
    _qnm_polish_fast(s, l, m, a, ω_seed, A_seed; ninv, rtol, maxiter, seed_step, angtol, ltol, growmax)
        -> (ω, λ, A, residual, lmax, angiters, nfev)

Polish one QNM root at spin `a` from the seeds `(ω_seed, A_seed)`, via the `:fast` angular shift-
invert coupled into the reused radial CF residual. Returns the root, its λ / separation constant,
the CF residual at the root, the adaptive l_max, the total RQI iterations, and the secant fevals.
"""
function _qnm_polish_fast(s::Integer, l::Integer, m::Integer, a::Real, ω_seed::Number, A_seed::Number;
                          ninv::Integer=0, rtol::Real=eps(Float64) * 64, maxiter::Integer=200,
                          seed_step=nothing, angtol::Real=1e-14, ltol::Real=1e-9,
                          growmax::Integer=60)
    T = float(real(typeof(a)))
    Aref = Ref(complex(A_seed))
    lmax_used = Ref(l)
    angiters = Ref(0)
    nfev = Ref(0)
    function G(ω)
        γ = complex(a * ω)
        λ, A_sep, _, lmx, it = _qnm_angular_fast(s, l, m, γ, Aref[]; ltol=ltol, growmax=growmax,
                                                  tol=angtol)
        Aref[] = A_sep
        lmax_used[] = lmx
        angiters[] += it
        nfev[] += 1
        return qnm_cf_residual(s, l, m, a, ω, λ; ninv=ninv, rtol=min(rtol, _qnm_default_rtol(ω)))
    end
    ω0 = Complex{T}(ω_seed)
    δ = seed_step === nothing ? Complex{T}(T(1) / 10^4 + im * T(1) / 10^4) : Complex{T}(seed_step)
    ω1 = ω0 * (1 + δ)
    ω = secant(G, ω0, ω1; rtol=rtol, maxiter=maxiter, throw_on_fail=false)
    # final consistent recompute (λ, A_sep, residual) at the root
    γ = complex(a * ω)
    λ, A_sep, tailfrac, lmx, _ = _qnm_angular_fast(s, l, m, γ, Aref[]; ltol=ltol, growmax=growmax,
                                                   tol=angtol)
    res = abs(qnm_cf_residual(s, l, m, a, ω, λ; ninv=ninv))
    return (ω=ω, λ=λ, A=A_sep, residual=res, lmax=lmx, tailfrac=tailfrac,
            angiters=angiters[], nfev=nfev[])
end

# ─────────────────────────────────────────────────────────────────────────────────────────────
# PRIORITY-2b — the spin-sequence continuation (KerrSpinSeq semantics).
# ─────────────────────────────────────────────────────────────────────────────────────────────
"quadratic (Newton divided-difference) fit through 3 (x,y) nodes; returns (eval::Function, d2)."
function _qnm_quadfit(x1, x2, x3, y1, y2, y3)
    d01 = (y2 - y1) / (x2 - x1)
    d12 = (y3 - y2) / (x3 - x2)
    c2 = (d12 - d01) / (x3 - x1)
    evalf(x) = y1 + d01 * (x - x1) + c2 * (x - x1) * (x - x2)
    return evalf, 2 * c2                          # p''(x) ≡ 2·c2 (curvature)
end

"""
    qnm_sequence_fast(s, l, m, n; a_max=0.99, seed=nothing, A0=nothing, init_step=1e-3,
                      dmin=1e-5, dmax=4e-3, curv=0.05, rtol, maxiter, angtol, ltol, growmax,
                      fail_tol=1e-6, build_interp=true)
        -> NamedTuple

Port of `qnm.spinsequence.KerrSpinSeq.do_find_sequence`: continue ω_{s,l,m,n}(a) from a=0 to
`a_max`. The a=0 root is polished from a DERIVED eikonal seed (or the supplied `seed`) with
`A0 = l(l+1)−s(s+1)`; each next node is seeded by the quadratic extrapolation of BOTH ω(a) and A(a)
through the last 3 nodes (linear for 2, fixed for 1), with the curvature-adaptive step
`Δa = curv/√(max d²)` clamped to `[dmin, dmax]`. A polish whose CF residual exceeds `fail_tol`
(or is non-finite) raises a `[BLOCKING]`-per-mode error (never a silent skip).

Returns `(; a, omega, A, lmax, nfev, iters, nnodes, omega_of_a, A_of_a)` — the node vectors, the
per-node diagnostics, and (if `build_interp`) cubic-spline interpolants (Re/Im splined jointly) for
ω(a) and A(a) over the built grid. Float64 SPEED path.
"""
function qnm_sequence_fast(s::Integer, l::Integer, m::Integer, n::Integer; a_max::Real=0.99,
                           seed=nothing, A0=nothing, init_step::Real=1e-3, dmin::Real=1e-5,
                           dmax::Real=4e-3, curv::Real=0.05, rtol::Real=eps(Float64) * 64,
                           maxiter::Integer=200, angtol::Real=1e-14, ltol::Real=1e-9,
                           growmax::Integer=60, fail_tol::Real=1e-6, build_interp::Bool=true)
    (l < abs(s) || abs(m) > l) &&
        throw(ArgumentError("invalid (s,l,m)=($s,$l,$m): need l>=|s| and |m|<=l"))
    T = Float64
    ω0 = seed === nothing ? _qnm_eikonal_seed(l, n, T) : Complex{T}(seed)
    A0v = A0 === nothing ? complex(T(l * (l + 1) - s * (s + 1))) : complex(A0)
    r0 = _qnm_polish_fast(s, l, m, zero(T), ω0, A0v; ninv=n, rtol=rtol, maxiter=maxiter,
                          angtol=angtol, ltol=ltol, growmax=growmax)
    (!isfinite(r0.residual) || r0.residual > fail_tol) &&
        error("[BLOCKING] qnm_sequence_fast (s,l,m,n)=($s,$l,$m,$n): a=0 root did not converge " *
              "(residual=$(r0.residual)); the enclosing continuation cannot start.")
    as = T[0.0]; ωs = Complex{T}[r0.ω]; As = Complex{T}[r0.A]
    lmaxs = Int[r0.lmax]; nfevs = Int[r0.nfev]; iters = Int[r0.angiters]
    a = zero(T)
    while a < a_max
        k = length(as)
        if k == 1
            Δa = T(init_step); ωseed = ωs[end]; Aseed = As[end]
            anext = min(a + Δa, T(a_max))
        elseif k == 2
            Δa = T(init_step)
            anext = min(a + Δa, T(a_max))
            slω = (ωs[end] - ωs[end-1]) / (as[end] - as[end-1])   # linear extrapolation seed
            slA = (As[end] - As[end-1]) / (as[end] - as[end-1])
            ωseed = ωs[end] + slω * (anext - as[end])
            Aseed = As[end] + slA * (anext - as[end])
        else
            fω, d2ω = _qnm_quadfit(as[end-2], as[end-1], as[end], ωs[end-2], ωs[end-1], ωs[end])
            fA, d2A = _qnm_quadfit(as[end-2], as[end-1], as[end], As[end-2], As[end-1], As[end])
            d2 = max(abs(d2ω), abs(d2A))
            Δa = clamp(T(curv) / sqrt(max(d2, eps(T))), T(dmin), T(dmax))
            anext = min(a + Δa, T(a_max))
            ωseed = fω(anext); Aseed = fA(anext)
        end
        anext <= a && break                         # reached a_max
        r = _qnm_polish_fast(s, l, m, anext, ωseed, Aseed; ninv=n, rtol=rtol, maxiter=maxiter,
                             angtol=angtol, ltol=ltol, growmax=growmax)
        (!isfinite(r.residual) || r.residual > fail_tol) &&
            error("[BLOCKING] qnm_sequence_fast (s,l,m,n)=($s,$l,$m,$n): polish failed at a=$anext " *
                  "(residual=$(r.residual), |ω−seed|=$(abs(r.ω-ωseed))). Owner must widen the step " *
                  "clamp or reseed; the sequence cannot advance past this node.")
        push!(as, anext); push!(ωs, r.ω); push!(As, r.A)
        push!(lmaxs, r.lmax); push!(nfevs, r.nfev); push!(iters, r.angiters)
        a = anext
    end
    ω_of_a = nothing; A_of_a = nothing
    if build_interp && length(as) >= 3
        ωm = permutedims(hcat(real.(ωs), imag.(ωs)))     # 2×N  [Re; Im]
        Am = permutedims(hcat(real.(As), imag.(As)))
        sω = CubicSplineInterpolant(as, ωm); sA = CubicSplineInterpolant(as, Am)
        ω_of_a = x -> (v = sω(x); complex(v[1], v[2]))
        A_of_a = x -> (v = sA(x); complex(v[1], v[2]))
    end
    return (a=as, omega=ωs, A=As, lmax=lmaxs, nfev=nfevs, iters=iters, nnodes=length(as),
            omega_of_a=ω_of_a, A_of_a=A_of_a)
end

# ─────────────────────────────────────────────────────────────────────────────────────────────
# The single-frequency `:fast` front door — called by `qnm_frequency(...; backend=:fast)`.
# Same signature semantics as `:leaver`, plus `seed_A` (the tracked separation-constant seed).
# ─────────────────────────────────────────────────────────────────────────────────────────────
function _qnm_frequency_fast(s::Integer, l::Integer, m::Integer, n::Integer, a::Real;
                             seed=nothing, seed_A=nothing, ninv::Integer=n,
                             rtol::Real=eps(Float64) * 64, maxiter::Integer=200, seed_step=nothing,
                             return_residual::Bool=false, angtol::Real=1e-14, ltol::Real=1e-9,
                             growmax::Integer=60)
    T = Float64
    ω0 = seed === nothing ? _qnm_eikonal_seed(l, n, T) : Complex{T}(seed)
    A0 = seed_A === nothing ? complex(T(l * (l + 1) - s * (s + 1))) : complex(seed_A)
    r = _qnm_polish_fast(s, l, m, T(a), ω0, A0; ninv=ninv, rtol=rtol, maxiter=maxiter,
                         seed_step=seed_step, angtol=angtol, ltol=ltol, growmax=growmax)
    return return_residual ? (r.ω, r.residual) : r.ω
end
