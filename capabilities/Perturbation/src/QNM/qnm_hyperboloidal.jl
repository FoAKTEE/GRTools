# QNM/qnm_hyperboloidal.jl — the hyperboloidal Kerr-QNM backend: Chebyshev collocation of the
# compactified horizon-penetrating radial operator (2410.14529 supp:42-55 == Ripley 2202.03837),
# solved as the QEP L(ω) = C₀ + ωC₁ + ω²C₂ (F64 companion seed → BigFloat RQI). Regularity at
# scri (R=0) and horizon (R=R₊) is implicit in the compactified analytic basis — NO boundary rows.
#
#   α ψ° + β ∂_R ψ° − ½Δ̃R² ∂²_R ψ° = −(Λ/2 + 1) ψ°,   ψ_s = R^{1+2s}ψ°_s,   R = 1/r.
#
# Two constraints the printed sources do not show:
#  • Δ̃ = 1 − 2MR + a²R² — the printed supp:44 `1−2MR²` is a typo (its principal part vanishes at
#    the WRONG radius, breaking the horizon structure); the test asserts Δ̃(1/r₊) ≡ 0 exactly.
#  • the ½-form: the printed unit principal with half-scaled α,β is internally inconsistent; the
#    consistent operator is R22/2 (global scale — same spectrum and eigenfunctions).
# λ bookkeeping: the raw constant A_raw = λ + 2maω − a²ω² (λ = the Harmonics public eigenvalue)
# enters the potential; in the QEP split λ is FROZEN at the trial ω (outer self-consistency
# iterated) while the polynomial 2maω − a²ω² folds exactly into C₁, C₂.

"""
    _hphc_alpha(s, m, a, ω, R; M=1) -> Complex
    _hphc_beta(s, m, a, ω, R; M=1)  -> Complex

The closed-form α, β of the hyperboloidal radial operator (2410 supp:50-55). Kept as scalar
functions of `R` so the tests can gate the ω-power collected assembly against them.
{T}-generic; `s,m` integers, `a` real, `ω` complex.
"""
@inline function _hphc_alpha(s::Integer, m::Integer, a::Real, ω::Number, R::Number; M::Real=1)
    i = im; h = one(real(float(ω))) / 2
    return a^2 * R^2 * (8 * M^2 * ω^2 + 6i * M * ω - 1) -
           R * (4 * M * ω + i) * (-a * (a * ω + m) + 4 * M^2 * ω + i * M * (s + 1)) +
           h * ω^2 * (a^2 - 16 * M^2) + a * m * ω + s * (-h + 2i * M * ω)
end
@inline function _hphc_beta(s::Integer, m::Integer, a::Real, ω::Number, R::Number; M::Real=1)
    i = im
    return 2 * a^2 * R^3 * (-1 + 2i * M * ω) + R * (-s - 1) + i * ω +
           R^2 * (i * a * (a * ω + m) - 8i * M^2 * ω + M * (s + 3))
end

"""
    hphc_radial_operator(s, m, a, ω, λ; n=48, T=..., M=1) -> NamedTuple

Assemble the `(n+1)×(n+1)` Chebyshev collocation operator of the hyperboloidal Teukolsky radial QEP
at frequency `ω` with the FROZEN angular eigenvalue `λ` (Harmonics public convention;
`A_raw = λ + 2maω − a²ω²`). Returns `(; C0, C1, C2, L, R, Rp, rp, iscri, ihor)` where
`L = C0 + ω·C1 + ω²·C2` (the operator at `ω`), `(C0,C1,C2)` the exposed QEP split for the eigen
path, `R` the CGL nodes (`R[iscri]=0` scri, `R[ihor]=R₊` horizon). Assembled in the ½-form
(`−½Δ̃R² ∂²_R`, `Δ̃=1−2MR+a²R²`) from the verbatim α,β. {T}/BigFloat-inherited (pass `T=BigFloat`).
"""
function hphc_radial_operator(s::Integer, m::Integer, a::Real, ω::Number, λ::Number;
                              n::Integer=48,
                              T::Type=float(real(promote_type(typeof(a), real(typeof(ω)),
                                                              real(typeof(λ))))),
                              M::Real=1)
    RT = T; CT = Complex{RT}
    Mv = RT(M); av = RT(a); sT = RT(s); mT = RT(m)
    rp = Mv + sqrt(Mv^2 - av^2)            # outer horizon r₊
    Rp = one(RT) / rp                       # compactified horizon R₊ = 1/r₊
    x, D1r, D2r = cheb_collocation(RT, n; lo=zero(RT), hi=Rp)   # x[1]=R₊ (hor), x[end]=0 (scri)
    D1 = CT.(D1r); D2 = CT.(D2r); R = x; i = im
    # ω-power coefficient vectors (collected from the verbatim closed forms — see header)
    aw0 = @. -av^2 * R^2 + R * (i * av * mT + Mv * (sT + 1)) - sT / 2
    aw1 = @. 6i * Mv * av^2 * R^2 + (av * mT + 2i * Mv * sT) -
             R * (-4 * Mv * av * mT + 4i * Mv^2 * (sT + 1) + i * (4 * Mv^2 - av^2))
    aw2 = @. 8 * Mv^2 * av^2 * R^2 + (av^2 - 16 * Mv^2) / 2 - 4 * Mv * (4 * Mv^2 - av^2) * R
    bw0 = @. -2 * av^2 * R^3 + R * (-sT - 1) + R^2 * (i * av * mT + Mv * (sT + 3))
    bw1 = @. 4i * Mv * av^2 * R^3 + i + R^2 * (i * av^2 - 8i * Mv^2)
    Δ̃  = @. 1 - 2 * Mv * R + av^2 * R^2
    P   = @. -(Δ̃ * R^2) / 2                # ½-form principal part (−½Δ̃R² == R22/2)
    λc  = CT(λ)
    # QEP split, λ frozen: A_raw/2 = λ/2 + maω − ½a²ω² ⇒ +ma → C1, −½a² → C2; +s/2 → C0.
    C0 = Diagonal(CT.(P)) * D2 + Diagonal(CT.(bw0)) * D1 +
         Diagonal(CT.(aw0 .+ λc / 2 .+ sT / 2))
    C1 = Diagonal(CT.(bw1)) * D1 + Diagonal(CT.(aw1 .+ mT * av))
    C2 = Diagonal(CT.(aw2 .- av^2 / 2))
    ωc = CT(ω)
    L  = Matrix(C0) .+ ωc .* Matrix(C1) .+ ωc^2 .* Matrix(C2)
    return (; C0=Matrix(C0), C1=Matrix(C1), C2=Matrix(C2), L=L,
            R=R, Rp=Rp, rp=rp, iscri=n + 1, ihor=1)
end

"""
    _hphc_companion(C0, C1, C2) -> Matrix

Standard first-companion linearization of the QEP `(C0 + ωC1 + ω²C2)v = 0`: the `2n×2n` matrix whose
eigenvalues are the QEP frequencies ω and whose eigenvectors are `[v; ωv]`. {T}/BigFloat-clean
(diagonal `C2` inverse). Mirrors `cond_bf.jl::companion`.
"""
function _hphc_companion(C0::AbstractMatrix, C1::AbstractMatrix, C2::AbstractMatrix)
    CT = eltype(C0); n = size(C0, 1)
    Z = zeros(CT, n, n); Id = Matrix{CT}(I, n, n)
    C2i = inv(C2)
    return [Z Id; -(C2i * C0) -(C2i * C1)]
end

"""
    _hphc_rqi(L, μ0, v0; maxit=60, tol=...) -> (μ, v, res)

Rayleigh-quotient iteration for the eigenpair of `L` nearest the shift `μ0`, seeded by `v0`. Each
step is one generic-LU shift-invert solve (`cheb_solve`, BigFloat-clean) + a Rayleigh-quotient
update — the PRODUCTIONIZED `cond_bf.jl::rqi`. Returns the eigenvalue `μ`, unit eigenvector `v`, and
the residual `‖Lv − μv‖`. {T}/BigFloat-inherited.
"""
function _hphc_rqi(L::AbstractMatrix, μ0::Number, v0::AbstractVector;
                   maxit::Integer=60, tol::Real=eps(real(eltype(L))) * 64)
    CT = eltype(L); RT = real(CT); N = size(L, 1)
    v = CT.(v0); v ./= norm(v); μ = CT(μ0)
    Id = Matrix{CT}(I, N, N); res = RT(Inf)
    for _ in 1:maxit
        w = cheb_solve(L - μ * Id, v)
        any(!isfinite, w) && break
        v = w ./ norm(w)
        μ = dot(v, L * v)                   # Rayleigh quotient (v Hermitian-unit)
        res = norm(L * v - μ * v)
        res <= tol && break
    end
    return μ, v, res
end

"""
    _hphc_nullvector(L; maxit=40, tol=...) -> (v, defect)

The null eigenvector of the (near-singular at a QNM ω) collocation matrix `L(ω)` by inverse
iteration at zero shift (`L` amplifies its smallest-singular direction). Returns the unit vector `v`
and the annihilation defect `‖Lv‖`. {T}/BigFloat-clean — the eigenfunction path when generic SVD is
unavailable at BigFloat.
"""
function _hphc_nullvector(L::AbstractMatrix; maxit::Integer=40, tol::Real=eps(real(eltype(L))) * 64)
    CT = eltype(L); RT = real(CT); N = size(L, 1)
    v = ones(CT, N); v ./= norm(v)
    defect = RT(Inf); prev = RT(Inf)
    for _ in 1:maxit
        w = cheb_solve(L, v)
        (any(!isfinite, w) || iszero(norm(w))) && break
        v = w ./ norm(w)
        defect = norm(L * v)
        (defect <= tol || abs(defect - prev) <= tol * max(defect, one(RT))) && break
        prev = defect
    end
    return v, defect
end

# The :hyperboloidal frequency front door: F64 self-consistent companion eigensolve → BigFloat
# RQI polish. The branch is SELECTED by an independent Leaver/:fast seed, but the VALUE is the
# collocation operator's own residual-gated root (null-defect ‖Lv‖ returned as the residual).
function _qnm_frequency_hyperboloidal(s::Integer, l::Integer, m::Integer, n_ovt::Integer, a::Real;
                                      seed=nothing, n_coll::Integer=48, sc_iters::Integer=12,
                                      rtol::Real=eps(float(real(typeof(a)))) * 64,
                                      return_residual::Bool=false)
    T = float(real(typeof(a)))
    a64 = Float64(a)
    # 1. branch-selection seed (independent Leaver reference; only localizes the collocation root)
    ωsel::ComplexF64 = seed === nothing ?
        ComplexF64(qnm_frequency(s, l, m, n_ovt, a64; backend=:leaver)) : ComplexF64(seed)
    # 2. Float64 self-consistent companion eigensolve (freeze λ, eigensolve, pick nearest, iterate λ)
    ωf = ωsel; vf = ComplexF64[]
    for _ in 1:sc_iters
        λ = _qnm_angular_lambda(s, l, m, a64, ωf)
        op = hphc_radial_operator(s, m, a64, ωf, λ; n=n_coll, T=Float64)
        F = eigen(_hphc_companion(op.C0, op.C1, op.C2))
        k = argmin(abs.(F.values .- ωf))
        ωnew = F.values[k]
        vf = ComplexF64.(F.vectors[1:(n_coll + 1), k])
        (abs(ωnew - ωf) <= 1e-14 * max(abs(ωnew), 1.0)) && (ωf = ωnew; break)
        ωf = ωnew
    end
    if !(T <: BigFloat)
        ω = Complex{T}(ωf)
        if return_residual
            λ = _qnm_angular_lambda(s, l, m, a64, ωf)
            op = hphc_radial_operator(s, m, a64, ωf, λ; n=n_coll, T=Float64)
            _, defect = _hphc_nullvector(op.L)
            return (ω, T(defect))
        end
        return ω
    end
    # 3. BigFloat polish — BRANCH-LOCKED Rayleigh-quotient iteration on the QEP companion with
    #    λ(ω)-continuation. ONE shift-invert per iteration (not a full inner RQI), carrying the companion vector forward
    #    so the iteration stays on the PHYSICAL branch seeded by the F64 eigenvector — essential for
    #    s=+2, whose collocation spectrum carries a spurious near-zero mode that a zero-shift inverse
    #    iteration on `L(ω)` would drift to. λ is refreshed each step (folded polynomial ω-parts stay
    #    exact). {T}/BigFloat-inherited; converges cubically once λ stabilises.
    # λ at full BigFloat (Harmonics Leaver path) — the conditioning-sensitive radial LU/RQI dominates
    # the per-mode cost, so a full-precision λ (down to the retrograde near-extremal corner where the
    # F64 spheroidal eigenvalue would lose ~1e-12) is retained at negligible marginal cost.
    aT = T(a); ω = Complex{T}(ωf)
    vc = Complex{T}.(vcat(vf, ComplexF64(ωf) .* vf)); vc ./= norm(vc)   # companion seed [v; ωv]
    N2 = 2 * (n_coll + 1)
    for _ in 1:sc_iters
        λ = _qnm_angular_lambda(s, l, m, aT, ω)
        op = hphc_radial_operator(s, m, aT, ω, λ; n=n_coll, T=T)
        Lc = _hphc_companion(op.C0, op.C1, op.C2)
        w = cheb_solve(Lc - ω * Matrix{Complex{T}}(I, N2, N2), vc)
        any(!isfinite, w) && break
        vc = w ./ norm(w)
        μ = dot(vc, Lc * vc)                               # Rayleigh quotient (→ ω)
        Δ = abs(μ - ω); ω = μ
        (Δ <= rtol * max(abs(ω), one(T))) && break
    end
    if return_residual
        λ = _qnm_angular_lambda(s, l, m, aT, ω)
        op = hphc_radial_operator(s, m, aT, ω, λ; n=n_coll, T=T)
        _, defect = _hphc_nullvector(op.L)
        return (ω, defect)
    end
    return ω
end

"""
    qnm_radial_hphc(s, l, m, n_ovt, a; n=48, T=..., seed=nothing) -> NamedTuple

The hyperboloidal QNM RADIAL eigenfunction ℛ(R) for `ω_{s,l,m,n_ovt}(a)`. Solves the frequency via
`backend=:hyperboloidal`, then extracts the null eigenfunction of `L(ω)` (inverse iteration) and
pins the in-house normalization `ℛ(scri) = ℛ(R=0) = 1` (STATED convention — the overall complex
constant is arbitrary; the scri/horizon R-RATIO is normalization-invariant, PREP-B L4/OBS-2).

Returns `(; omega, lambda, R, Rvals, R_of, coeffs, ratio_scri_hor, residual, n, rp, Rp)` where
`Rvals = ℛ` at the CGL nodes (scri-normalized), `R_of(Rq)` a barycentric evaluator, `coeffs` the
Chebyshev `T_k(2R/R₊ − 1)` coefficients, and `ratio_scri_hor = ℛ(scri)/ℛ(horizon)` (the invariant
compared to the MMA fixture's `I`). {T}/BigFloat-inherited.
"""
function qnm_radial_hphc(s::Integer, l::Integer, m::Integer, n_ovt::Integer, a::Real;
                         n::Integer=48, T::Type=float(real(typeof(a))), seed=nothing)
    aT = T(a)
    ω = _qnm_frequency_hyperboloidal(s, l, m, n_ovt, aT; n_coll=n, seed=seed)
    λ = _qnm_angular_lambda(s, l, m, aT, ω)
    op = hphc_radial_operator(s, m, aT, ω, λ; n=n, T=T)
    v, defect = _hphc_nullvector(op.L)
    scri = v[op.iscri]; hor = v[op.ihor]
    Rvals = v ./ scri                         # ℛ(scri) = 1 (in-house pin)
    R = op.R
    w = bary_weights(R)
    R_of(Rq) = bary_eval(R, w, Rvals, T(Rq))
    coeffs = _hphc_cheb_coeffs(R, Rvals, op.Rp)
    return (; omega=ω, lambda=λ, R=R, Rvals=Rvals, R_of=R_of, coeffs=coeffs,
            ratio_scri_hor=scri / hor, residual=defect, n=n, rp=op.rp, Rp=op.Rp,
            iscri=op.iscri, ihor=op.ihor)
end

"Chebyshev coefficients c_k of ℛ ≈ Σ c_k T_k(ξ), ξ = 2R/R₊ − 1 (scri ξ=−1, horizon ξ=+1)."
function _hphc_cheb_coeffs(R::AbstractVector{RT}, Rvals::AbstractVector, Rp::RT) where {RT<:Real}
    N = length(R); ξ = @. 2 * R / Rp - 1
    V = Matrix{RT}(undef, N, N)
    @inbounds for r in 1:N
        V[r, 1] = one(RT); N >= 2 && (V[r, 2] = ξ[r])
        for k in 2:(N - 1)
            V[r, k + 1] = 2 * ξ[r] * V[r, k] - V[r, k - 1]
        end
    end
    return cheb_solve(Complex{RT}.(V), Complex{RT}.(Rvals))
end

"""
    hphc_mirror_carry(ω, m) -> (ωΦ, mΦ) = (−conj(ω), −m)

The Hertz mirror carry (PREP-B ledger L3; eq:qnm_kerr_symmetry `ω^p_{ℓmn} = −ω̄^{−p}_{ℓ,−m,n}`). To
place the physical QNM strain amplitude of mode `(ω, m)` (read from `Φ̄`), the s=+2 Hertz potential
`Φ` (the `ζ⁻⁴Φ` variable, the Q3 source input) must carry frequency `−ω̄` and azimuthal number `−m`.
Documented helper — the out/in boundary conditions are NOT mixed by the conjugation.
"""
hphc_mirror_carry(ω::Number, m::Integer) = (-conj(ω), -m)
