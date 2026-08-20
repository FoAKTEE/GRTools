# SpheroidalEigenvalues — λ(s,l,m,γ) + ∂_γλ: pentadiagonal spectral solve (Float64/geev),
# Leaver continued fraction (the BigFloat home), and the small-/large-γ asymptotic series.
# Real AND complex γ. LAPACK eigensolves run under blas_single_thread (BigFloat suffers
# GC+BLAS contention under threads).

# --- pentadiagonal spectral matrix elements: real γ ⇒ symmetric, complex γ ⇒ complex-SYMMETRIC
#     (the elements are real-integer-coefficient polynomials in γ; bit-identical on real input) ---

"diagonal spectral matrix element (`l==0&&m==0` special case first)."
function _kHat(s::Integer, l::Integer, m::Integer, g::Number)
    (l == 0 && m == 0) && return g^2 / 3
    return -l * (1 + l) + (2 * m * s^2 * g) / (l + l^2) +
           (1 / 3) * (1 + (2 * (l + l^2 - 3 * m^2) * (l + l^2 - 3 * s^2)) /
                          (l * (-3 + l + 8 * l^2 + 4 * l^3))) * g^2
end

"`±2` band spectral matrix element."
function _k2(s::Integer, l::Integer, m::Integer, g::Number)
    num = (1 + l - m) * (2 + l - m) * (1 + l + m) * (2 + l + m) *
          (1 + l - s) * (2 + l - s) * (1 + l + s) * (2 + l + s)
    return sqrt(num / ((1 + 2 * l) * (5 + 2 * l))) * g^2 / ((1 + l) * (2 + l) * (3 + 2 * l))
end

"`±1` band spectral matrix element (`l==0&&m==0` special case first)."
function _kTilde2(s::Integer, l::Integer, m::Integer, g::Number)
    (l == 0 && m == 0) && return -2 * s * sqrt((1 - s^2) / 3) * g
    return -(2 * s * sqrt(((1 + 2 * l + l^2 - m^2) * (1 + 2 * l + l^2 - s^2)) /
                          (3 + 8 * l + 4 * l^2)) * g * (2 * l + l^2 + m * g)) /
            (l * (2 + 3 * l + l^2))
end

"""
    _swsh_spectral_setup(s, l, m, γ, numterms) -> (M, nUp, nDown, lmin)

Build the pentadiagonal symmetric (real γ) spectral matrix, shared by the eigenvalue and the
harmonic (eigenvector) paths.
"""
function _swsh_spectral_setup(s::Integer, l::Integer, m::Integer, g::Number, numterms)
    nUp = numterms === :auto ? ceil(Int, abs(3 * g / 2)) + 5 : Int(numterms)
    lmin = max(abs(s), abs(m))
    nDown = min(nUp, l - lmin)
    n = nUp + nDown + 1
    M = zeros(float(typeof(g)), n, n)
    for i in 1:n
        Ld = l - nDown - 1 + i                 # l-index of row i (diagonal)
        M[i, i] = _kHat(s, Ld, m, g)
        i - 2 >= 1 && (M[i, i-2] = _k2(s, Ld - 2, m, g))
        i - 1 >= 1 && (M[i, i-1] = _kTilde2(s, Ld - 1, m, g))
        i + 1 <= n && (M[i, i+1] = _kTilde2(s, Ld, m, g))
        i + 2 <= n && (M[i, i+2] = _k2(s, Ld, m, g))
    end
    return M, nUp, nDown, lmin
end

# --- γ-derivatives of the spectral matrix elements (analytic; Hellmann-Feynman dλ/dγ) ---

"`d/dγ` of `_kHat` (diagonal)."
function _dkHat(s::Integer, l::Integer, m::Integer, g::Number)
    (l == 0 && m == 0) && return 2 * g / 3
    a1 = (2 * m * s^2) / (l + l^2)
    a2 = (1 / 3) * (1 + (2 * (l + l^2 - 3 * m^2) * (l + l^2 - 3 * s^2)) /
                        (l * (-3 + l + 8 * l^2 + 4 * l^3)))
    return a1 + 2 * a2 * g                                   # _kHat = a0 + a1 γ + a2 γ²
end

"`d/dγ` of `_k2` (`±2` band; `_k2 ∝ γ²`)."
_dk2(s::Integer, l::Integer, m::Integer, g::Number) =
    iszero(g) ? zero(g) : 2 * _k2(s, l, m, g) / g

"`d/dγ` of `_kTilde2` (`±1` band; `_kTilde2 = -A((2l+l²)γ + m γ²)`)."
function _dkTilde2(s::Integer, l::Integer, m::Integer, g::Number)
    (l == 0 && m == 0) && return -2 * s * sqrt((1 - s^2) / 3)
    A = (2 * s * sqrt(((1 + 2 * l + l^2 - m^2) * (1 + 2 * l + l^2 - s^2)) /
                      (3 + 8 * l + 4 * l^2))) / (l * (2 + 3 * l + l^2))
    return -A * ((2 * l + l^2) + 2 * m * g)
end

"build dM/dγ over the same index layout as `_swsh_spectral_setup`."
function _swsh_dM_dgamma(s::Integer, l::Integer, m::Integer, g::Number, nUp::Integer, nDown::Integer)
    n = nUp + nDown + 1
    dM = zeros(float(typeof(g)), n, n)
    for i in 1:n
        Ld = l - nDown - 1 + i
        dM[i, i] = _dkHat(s, Ld, m, g)
        i - 2 >= 1 && (dM[i, i-2] = _dk2(s, Ld - 2, m, g))
        i - 1 >= 1 && (dM[i, i-1] = _dkTilde2(s, Ld - 1, m, g))
        i + 1 <= n && (dM[i, i+1] = _dkTilde2(s, Ld, m, g))
        i + 2 <= n && (dM[i, i+2] = _dk2(s, Ld, m, g))
    end
    return dM
end

# --- spectral eigenvalue / eigenvector (BLAS-guarded LAPACK) ---

"""
    _swsh_eigenvalue_spectral(s, l, m, γ; numterms=:auto) -> Real

Spectral eigenvalue `E` (the `-s(s+1)`-shifted spheroidal eigenvalue): the `(nDown+1)`-th-largest
eigenvalue of the pentadiagonal matrix, `E = -Sort[eigvals][[-(nDown+1)]] - s(s+1)`.
"""
function _swsh_eigenvalue_spectral(s::Integer, l::Integer, m::Integer, g::Real; numterms=:auto)
    M, _, nDown, _ = _swsh_spectral_setup(s, l, m, g, numterms)
    ev = blas_single_thread(() -> eigvals(Symmetric(M)))   # ascending real (symmetric for real γ)
    return -ev[end-nDown] - s * (s + 1)
end

"""
    _swsh_spectral_coeffs(s, l, m, γ; numterms=:auto) -> (λ, coeffs, lmin', lmax')

Public eigenvalue `λ = E - 2mγ + γ²` and the sign-fixed unit eigenvector `coeffs`
(spherical-expansion coefficients, ordered `l'=lmin'..lmax'`). Sign convention: pin a central
component positive.
"""
function _swsh_spectral_coeffs(s::Integer, l::Integer, m::Integer, g::Real; numterms=:auto)
    M, nUp, nDown, lmin = _swsh_spectral_setup(s, l, m, g, numterms)
    F = blas_single_thread(() -> eigen(Symmetric(M)))      # values ascending, vectors aligned
    k = lastindex(F.values) - nDown            # the (nDown+1)-th-largest eigenvalue
    evalraw = -F.values[k]
    evec = F.vectors[:, k]
    idx = min(l - lmin + 1, div(nUp + nDown, 2) + 1)
    sgn = sign(evec[idx]); sgn == 0 && (sgn = one(sgn))
    lambda = evalraw - s * (s + 1) - 2 * m * g + g^2
    return lambda, sgn .* evec, l - nDown, l + nUp
end

"""
    _swsh_eigenvalue_dgamma(s, l, m, γ; numterms=:auto) -> Real

`dλ/dγ` of the spectral eigenvalue via Hellmann-Feynman: `dE_raw/dγ = vᵀ(dM/dγ)v` for the
selected unit eigenvector `v`. With `λ = -E_sel - s(s+1) - 2mγ + γ²`,
`dλ/dγ = -vᵀ(dM/dγ)v - 2m + 2γ`.
"""
function _swsh_eigenvalue_dgamma(s::Integer, l::Integer, m::Integer, g::Real; numterms=:auto)
    M, nUp, nDown, _ = _swsh_spectral_setup(s, l, m, g, numterms)
    F = blas_single_thread(() -> eigen(Symmetric(M)))
    v = F.vectors[:, lastindex(F.values)-nDown]
    dM = _swsh_dM_dgamma(s, l, m, g, nUp, nDown)
    return -dot(v, Symmetric(dM), v) - 2 * m + 2 * g
end

"""
    _swsh_coeffs_dgamma(s, l, m, γ; numterms=:auto) -> (dcoeffs, lmin', lmax')

`d/dγ` of the sign-fixed unit eigenvector via first-order eigenvector perturbation:
`dv/dγ = Σ_{j≠sel} ⟨v_j,(dM/dγ)v_sel⟩/(μ_sel-μ_j) v_j`, sign-fixed consistently with
`_swsh_spectral_coeffs`.
"""
function _swsh_coeffs_dgamma(s::Integer, l::Integer, m::Integer, g::Real; numterms=:auto)
    M, nUp, nDown, lmin = _swsh_spectral_setup(s, l, m, g, numterms)
    F = blas_single_thread(() -> eigen(Symmetric(M)))
    n = nUp + nDown + 1
    k = lastindex(F.values) - nDown
    vsel = F.vectors[:, k]
    dMv = Symmetric(_swsh_dM_dgamma(s, l, m, g, nUp, nDown)) * vsel
    dv = zeros(float(typeof(g)), n)
    @inbounds for j in 1:n
        j == k && continue
        dv .+= (dot(view(F.vectors, :, j), dMv) / (F.values[k] - F.values[j])) .* view(F.vectors, :, j)
    end
    idx = min(l - lmin + 1, div(nUp + nDown, 2) + 1)
    sgn = sign(vsel[idx]); sgn == 0 && (sgn = one(sgn))
    return sgn .* dv, l - nDown, l + nUp
end

# --- complex-γ spectral tier (Q1a) ---
# For COMPLEX γ=aω the pentadiagonal matrix is complex-SYMMETRIC (M=Mᵀ, not Hermitian), so the
# LAPACK symmetric/Hermitian eigensolvers do NOT apply. We use the dense NONsymmetric eigensolver
# (LAPACK geev) — exact for these small pentadiagonal sizes (Q0 dossier: "dense nonsymmetric eigen
# of the pentadiagonal is fine"). geev returns UNSORTED complex eigenpairs, so the physical branch
# (the one that continuously deforms from ₛY_lm as γ→0) is selected by EIGENVECTOR TRACKING: its
# spherical expansion stays dominated by the l-th component (matrix row `nDown+1`, whose diagonal
# is `_kHat(s,l,m,γ)`) for the moderate |γ|≲1.2 regime this backend targets. This is a Float64
# LAPACK path (ComplexF64); a Complex{BigFloat} γ must use the Leaver CF backend (which seeds off
# this at ComplexF64). The Hellmann-Feynman derivative uses the complex-symmetric left eigenvector
# `= transpose(v)` (NOT the conjugate `v'`).

"select the geev eigenpair dominated by the l-th spherical component (matrix row `selrow`)."
function _select_spheroidal_index(vectors::AbstractMatrix, selrow::Integer)
    n = size(vectors, 2)
    kbest = 1
    wbest = abs(vectors[selrow, 1])
    @inbounds for k in 2:n
        w = abs(vectors[selrow, k])
        if w > wbest
            wbest = w
            kbest = k
        end
    end
    return kbest
end

"""
    _swsh_eigenvalue_spectral(s, l, m, γ::Complex; numterms=:auto) -> Complex

Complex-γ spectral eigenvalue `E` (the `-s(s+1)`-shifted spheroidal eigenvalue): the geev
eigenvalue of the branch tracked to the l-th spherical component, `E = -μ_sel - s(s+1)`.
"""
function _swsh_eigenvalue_spectral(s::Integer, l::Integer, m::Integer, g::Complex; numterms=:auto)
    M, _, nDown, _ = _swsh_spectral_setup(s, l, m, g, numterms)
    F = blas_single_thread(() -> eigen(M))                 # general (geev): complex eigenpairs
    k = _select_spheroidal_index(F.vectors, nDown + 1)
    return -F.values[k] - s * (s + 1)
end

"""
    _swsh_spectral_coeffs(s, l, m, γ::Complex; numterms=:auto) -> (λ, coeffs, lmin', lmax')

Complex-γ public eigenvalue `λ = E - 2mγ + γ²` and the phase-fixed unit eigenvector `coeffs`
(spherical-expansion coefficients, `l'=lmin'..lmax'`, Euclidean-unit-norm — same normalisation as
the real path). Phase convention: the central component `evec[idx]` is pinned real-positive (the
complex generalisation of the real path's `sign` fix, `sign(x)=x/|x|`), so `coeffs → e_central` as
`γ→0`.
"""
function _swsh_spectral_coeffs(s::Integer, l::Integer, m::Integer, g::Complex; numterms=:auto)
    M, nUp, nDown, lmin = _swsh_spectral_setup(s, l, m, g, numterms)
    F = blas_single_thread(() -> eigen(M))
    k = _select_spheroidal_index(F.vectors, nDown + 1)
    evec = F.vectors[:, k]
    idx = min(l - lmin + 1, div(nUp + nDown, 2) + 1)
    ph = evec[idx]
    ph = iszero(ph) ? one(ph) : ph / abs(ph)               # pin central component real-positive
    lambda = -F.values[k] - s * (s + 1) - 2 * m * g + g^2
    return lambda, evec ./ ph, l - nDown, l + nUp
end

"""
    _swsh_eigenvalue_dgamma(s, l, m, γ::Complex; numterms=:auto) -> Complex

`dλ/dγ` at complex γ via Hellmann-Feynman for the complex-SYMMETRIC matrix: the left eigenvector
is `transpose(v)` (not `v'`), so `dμ/dγ = (vᵀ dM v)/(vᵀ v)`. With `λ = -μ_sel - s(s+1) - 2mγ + γ²`,
`dλ/dγ = -dμ/dγ - 2m + 2γ`.
"""
function _swsh_eigenvalue_dgamma(s::Integer, l::Integer, m::Integer, g::Complex; numterms=:auto)
    M, nUp, nDown, _ = _swsh_spectral_setup(s, l, m, g, numterms)
    F = blas_single_thread(() -> eigen(M))
    k = _select_spheroidal_index(F.vectors, nDown + 1)
    v = F.vectors[:, k]
    dM = _swsh_dM_dgamma(s, l, m, g, nUp, nDown)
    dmu = (transpose(v) * dM * v) / (transpose(v) * v)     # complex-symmetric Hellmann-Feynman
    return -dmu - 2 * m + 2 * g
end

# --- Leaver continued-fraction backend (type-generic; consumes NumericsCore CF) ---

"""
    _swsh_eigenvalue_leaver_residual(s, l, m, γ) -> (λ, residual)

Leaver continued-fraction eigenvalue (Leaver 1985) WITH the converged CF residual: secant-solve
`Eq(A)=LHS(A)-RHS(A)=0` for the separation constant `A` from the spectral initial guess, then
return `λ = A - 2mγ + γ²` (BHPT spheroidal convention) AND `residual = |Eq(A)|` at the root.

`residual` is the diagnostic the BigFloat radial precision-convergence gate consumes (≈ machine
floor of the active precision — it is the Leaver-equation defect at the secant root). This is the
single CF body; `_swsh_eigenvalue_leaver` is the bare-λ wrapper over it (AHA — one CF home).

Type-generic: `T` INHERITS the input precision (Float64 OR BigFloat), so a BigFloat γ runs the
whole CF/secant at full precision. The CF legs are `NumericsCore.cf_finite` (LHS finite K-CF) and
`NumericsCore.cf_convergent` (RHS tail). Folds KerrLorenzCirc `precision.jl/swsh_eigenvalue_big`
into this canonical Harmonics path — ONE BigFloat eigenvalue home.
"""
function _swsh_eigenvalue_leaver_residual(s::Integer, l::Integer, m::Integer, γ::Number)
    RT = typeof(float(real(γ)))                               # REAL underlying float (Float64/BigFloat)
    g = complex(float(γ)); CT = typeof(g)                     # complex spheroidicity, Complex{RT}
    # cheap Float64/ComplexF64 spectral seed, promoted to CT (mirrors the real path's Float64 seed).
    seedγ = γ isa Complex ? ComplexF64(g) : Float64(real(γ))
    Aini = CT(_swsh_eigenvalue_spectral(s, l, m, seedγ))
    nInv = l - max(abs(m), abs(s))                             # NumInversions
    α = abs(m + s)
    β = abs(m - s)
    ab = (α + β) / 2
    αn(n) = (-4g * (n + α + 1) * (n + β + 1) * (n + ab + 1 + s)) /
            ((2n + α + β + 2) * (2n + α + β + 3))
    βn(n, A) = A + s * (s + 1) + g^2 - (n + ab) * (n + ab + 1) +
               (s != 0 ? (8 * m * s^2 * g) / ((2n + α + β) * (2n + α + β + 2)) : zero(g))
    γn(n) = (4g * n * (n + α + β) * (n + ab - s)) /
            ((2n + α + β - 1) * (2n + α + β))
    LHS(A) = βn(nInv, A) +
             (nInv >= 1 ? cf_finite(n -> -αn(nInv - n) * γn(nInv - n + 1),
                                    n -> βn(nInv - n, A), 1, nInv - 1) : zero(A))
    rt = eps(RT) * 8                                           # precision-INHERITED secant/CF tolerance
    RHS(A) = -cf_convergent(n -> -αn(n - 1) * γn(n), n -> βn(n, A), nInv + 1; rtol=rt)
    Eq(A) = LHS(A) - RHS(A)
    # complex-seeded secant (`NumericsCore.secant`, real tolerances on a complex iterate);
    # `throw_on_fail=false` keeps the .m `FindRoot[..., Method->"Secant"]` behaviour of returning the
    # last iterate — the caller judges the CF residual returned alongside λ.
    A = secant(Eq, Aini, Aini * (1 + RT(1)/RT(10)^6) + RT(1)/RT(10)^8;
               rtol=rt, maxiter=400, throw_on_fail=false)
    λ = A - 2 * m * γ + g^2
    return (iszero(imag(λ)) ? real(λ) : λ), abs(Eq(A))
end

"""
    _swsh_eigenvalue_leaver(s, l, m, γ) -> spheroidal λ

Bare-λ wrapper over `_swsh_eigenvalue_leaver_residual` (single CF body): the Leaver
continued-fraction eigenvalue `λ = A - 2mγ + γ²`, discarding the CF residual.
"""
_swsh_eigenvalue_leaver(s::Integer, l::Integer, m::Integer, γ::Number) =
    first(_swsh_eigenvalue_leaver_residual(s, l, m, γ))

# --- public dispatcher ---

"""
    SpinWeightedSpheroidalEigenvalue(s, l, m, γ; method=:Auto, numterms=:auto, return_residual=false) -> λ

Spin-weighted spheroidal eigenvalue `λ` (BHPT convention; `γ=0` limit `l(l+1)-s(s+1)`), with
`γ=aω` real and integer `s`. `λ = E - 2mγ + γ²`.

`method=:Auto` routes a `BigFloat` γ to the type-generic Leaver continued-fraction backend
(INHERITS precision) and a `Float64` γ to the LAPACK spectral backend (fast). `:SphericalExpansion`
forces spectral; `:Leaver`/`:Automatic` force the CF backend.

The asymptotic-series backends (design §4.2 THIRD backend; SpheroidalEigenvalueSeries.jl) are
reachable here too: `:SmallGamma` evaluates the small-γ Taylor series `Σ c_i γ^i`
(coefficients `SpinWeightedSpheroidalEigenvalueSeries`; validity window |γ|≲1e-2 at machine
precision, measured) and `:LargeGamma` the large-γ Laurent series `evalpoly(γ, cl)/γ^order`
(coefficients `SpinWeightedSpheroidalEigenvalueSeriesLargeGamma`; asymptotic error decays
monotonically, measured ≲2e-5 rel at γ=20, ≲2e-7 at γ=40 for the default order). For these two
routes `numterms` is the series `order` (`:auto` → the series defaults, 6 small-γ / 4 large-γ).
The coefficient recurrences are Float64, so these routes carry a Float64 accuracy floor even
for BigFloat γ — the BigFloat home stays the Leaver backend.

`return_residual=true` routes the type-generic Leaver path unconditionally and returns the tuple
`(λ, |Eq(A)|)` — the eigenvalue plus the converged CF residual the BigFloat radial
precision-convergence gate consumes (folds KerrLorenzCirc `swsh_eigenvalue_big`; ONE BigFloat
eigenvalue home). The default scalar-λ return is unchanged.
"""
function SpinWeightedSpheroidalEigenvalue(s::Integer, l::Integer, m::Integer, gamma::Real;
                                          method::Symbol=:Auto, numterms=:auto,
                                          return_residual::Bool=false)
    (l < abs(s) || abs(m) > l) &&
        throw(ArgumentError("invalid (s,l,m)=($s,$l,$m): need l>=|s| and |m|<=l"))
    if return_residual
        iszero(gamma) && return (oftype(float(gamma), l * (l + 1) - s * (s + 1)),
                                 abs(zero(float(gamma))))
        return _swsh_eigenvalue_leaver_residual(s, l, m, gamma)
    end
    if method === :Auto
        method = (gamma isa BigFloat) ? :Leaver : :SphericalExpansion
    end
    iszero(gamma) && return oftype(float(gamma), l * (l + 1) - s * (s + 1))
    if method === :SphericalExpansion
        return _swsh_eigenvalue_spectral(s, l, m, gamma; numterms=numterms) - 2 * m * gamma + gamma^2
    elseif method === :Leaver || method === :Automatic
        return _swsh_eigenvalue_leaver(s, l, m, gamma)
    elseif method === :SmallGamma
        order = numterms === :auto ? 6 : Int(numterms)
        return evalpoly(gamma, SpinWeightedSpheroidalEigenvalueSeries(s, l, m; order=order))
    elseif method === :LargeGamma
        order = numterms === :auto ? 4 : Int(numterms)
        return evalpoly(gamma, SpinWeightedSpheroidalEigenvalueSeriesLargeGamma(s, l, m; order=order)) / gamma^order
    else
        throw(ArgumentError("unknown method $method (expected :SphericalExpansion, :Leaver, :SmallGamma, or :LargeGamma)"))
    end
end

"""
    SpinWeightedSpheroidalEigenvalueDeriv(s, l, m, γ; numterms=:auto) -> Real

`dλ/dγ` of the spheroidal eigenvalue w.r.t. the spheroidicity `γ=aω` (so `dλ/dω = a·dλ/dγ`),
via analytic Hellmann-Feynman on the spectral matrix. At `γ=0`: `dλ/dγ = -2m s²/(l+l²) - 2m`.
"""
function SpinWeightedSpheroidalEigenvalueDeriv(s::Integer, l::Integer, m::Integer, gamma::Real;
                                               numterms=:auto)
    (l < abs(s) || abs(m) > l) &&
        throw(ArgumentError("invalid (s,l,m)=($s,$l,$m): need l>=|s| and |m|<=l"))
    if iszero(gamma)
        # The only valid l=0 mode is (s,m)=(0,0); its derivative is exactly zero. Avoid the
        # unsimplified analytic expression's removable 0/0 and preserve gamma's floating type.
        iszero(l) && return zero(float(gamma))
        return float(-2 * m * s^2 / (l + l^2) - 2 * m)
    end
    return _swsh_eigenvalue_dgamma(s, l, m, gamma; numterms=numterms)
end

# --- complex-γ front doors (Q1a; the module header's declared [FUTURE], now BUILT) ---

"""
    SpinWeightedSpheroidalEigenvalue(s, l, m, γ::Complex; method=:Auto, numterms=:auto,
                                     return_residual=false) -> Complex

Complex-γ (`γ=aω` with complex QNM frequency ω) spin-weighted spheroidal eigenvalue. `method=:Auto`
routes a `Complex{BigFloat}` γ to the type-generic Leaver CF backend (INHERITS precision, the
BF-256 QQNM home) and a `ComplexF64` γ to the geev complex-symmetric spectral backend (moderate
precision). `:SphericalExpansion` forces the spectral geev path (ComplexF64 only — LAPACK has no
complex-BigFloat eigensolver); `:Leaver`/`:Automatic` force the CF backend (BigFloat-capable);
`:SmallGamma`/`:LargeGamma` evaluate the (real-coefficient) asymptotic series at complex γ.
`return_residual=true` routes the Leaver path and returns `(λ, |Eq(A)|)`.

Symmetry (verified): `λ(s,l,m,γ̄) == conj(λ(s,l,m,γ))` (the CF coefficients carry real integer
structure). `γ→0` gives the exact spherical value `l(l+1)-s(s+1) = (l-s)(l+s+1)`.
"""
function SpinWeightedSpheroidalEigenvalue(s::Integer, l::Integer, m::Integer, gamma::Complex;
                                          method::Symbol=:Auto, numterms=:auto,
                                          return_residual::Bool=false)
    (l < abs(s) || abs(m) > l) &&
        throw(ArgumentError("invalid (s,l,m)=($s,$l,$m): need l>=|s| and |m|<=l"))
    λ0 = oftype(complex(float(gamma)), l * (l + 1) - s * (s + 1))
    if return_residual
        iszero(gamma) && return (λ0, abs(zero(real(float(gamma)))))
        return _swsh_eigenvalue_leaver_residual(s, l, m, gamma)
    end
    if method === :Auto
        method = (real(gamma) isa BigFloat) ? :Leaver : :SphericalExpansion
    end
    iszero(gamma) && return λ0
    if method === :SphericalExpansion
        return _swsh_eigenvalue_spectral(s, l, m, gamma; numterms=numterms) - 2 * m * gamma + gamma^2
    elseif method === :Leaver || method === :Automatic
        return _swsh_eigenvalue_leaver(s, l, m, gamma)
    elseif method === :SmallGamma
        order = numterms === :auto ? 6 : Int(numterms)
        return evalpoly(gamma, SpinWeightedSpheroidalEigenvalueSeries(s, l, m; order=order))
    elseif method === :LargeGamma
        order = numterms === :auto ? 4 : Int(numterms)
        return evalpoly(gamma, SpinWeightedSpheroidalEigenvalueSeriesLargeGamma(s, l, m; order=order)) / gamma^order
    else
        throw(ArgumentError("unknown method $method (expected :SphericalExpansion, :Leaver, :SmallGamma, or :LargeGamma)"))
    end
end

"""
    SpinWeightedSpheroidalEigenvalueDeriv(s, l, m, γ::Complex; numterms=:auto) -> Complex

`dλ/dγ` at complex γ via complex-symmetric Hellmann-Feynman (geev; the left eigenvector is
`transpose(v)`). ComplexF64 path (the spectral tier). At `γ=0`: `dλ/dγ = -2m s²/(l+l²) - 2m`.
"""
function SpinWeightedSpheroidalEigenvalueDeriv(s::Integer, l::Integer, m::Integer, gamma::Complex;
                                               numterms=:auto)
    (l < abs(s) || abs(m) > l) &&
        throw(ArgumentError("invalid (s,l,m)=($s,$l,$m): need l>=|s| and |m|<=l"))
    if iszero(gamma)
        iszero(l) && return zero(float(gamma))
        return oftype(complex(float(gamma)), -2 * m * s^2 / (l + l^2) - 2 * m)
    end
    return _swsh_eigenvalue_dgamma(s, l, m, gamma; numterms=numterms)
end
