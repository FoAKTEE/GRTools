# SpheroidalHarmonics — the callable Float64-VALUE tier (design §4.2): `{}_sS_{lm}(θ,φ)` carriers
# in the SphericalExpansion and Leaver representations, with mixed ∂_θ/∂_φ/∂_γ. Depends up the
# in-batch DAG on the Y seed (SphericalHarmonics), the sign-fixed eigenvector coeffs + γ-deriv
# (SpheroidalEigenvalues), and the small-γ d-tables (SpheroidalEigenvalueSeries).
#
# Precision asymmetry (design §4.2, load-bearing): the eigenvalue/Leaver paths are BigFloat-
# inheriting, but the harmonic VALUE tier is Float64-capped — `SWSHSFunction`/`SWSHSLeaverFunction`
# store Float64 coeffs. The Leaver coeff NORM consumes `NumericsCore.kummer_1F1` (no local Kummer
# copy, AHA); the half-angle term-list eval/derivative machinery is the promoted
# `NumericsCore._eval_terms`/`_dtheta_terms`/`_y_dtheta` seed.

# ============================== SphericalExpansion representation ==============================

"""
    SWSHSFunction

Callable carrier for `{}_sS_{lm}(·;γ)` in the SphericalExpansion representation: `coeffs` is the
sign-/phase-fixed unit eigenvector over `l'=lmin..lmax`. Evaluate as `S(θ,φ)`.

Type-parametrised (Q1a): for real γ, `Tg=Tλ=Tc=Float64` (bit-identical to the original concrete
struct); for a complex QNM spheroidicity `γ=aω`, `Tg=Tλ=Tc=ComplexF64` (the geev complex-symmetric
spectral tier). The value tier stays Float64-capped by design (§4.2 precision asymmetry): the Y
seed is Float64, so `S(θ)` carries the Float64 ceiling even for a `Complex{BigFloat}` request.
"""
struct SWSHSFunction{Tg<:Number,Tλ<:Number,Tc<:Number}
    s::Int
    l::Int
    m::Int
    gamma::Tg
    eigenvalue::Tλ            # public λ (real for real γ, complex for complex γ)
    method::Symbol
    coeffs::Vector{Tc}
    lmin::Int
    lmax::Int
end

"""
    SpinWeightedSpheroidalHarmonicS(s, l, m, γ; method=:SphericalExpansion, numterms=:auto)

Construct the spheroidal-harmonic carrier `{}_sS_{lm}(·;γ)`. `method=:SphericalExpansion`
(default) returns a [`SWSHSFunction`](@ref); `:Leaver` returns a [`SWSHSLeaverFunction`](@ref).
"""
function SpinWeightedSpheroidalHarmonicS(s::Integer, l::Integer, m::Integer, gamma::Real;
                                         method::Symbol=:SphericalExpansion, numterms=:auto)
    (l < abs(s) || abs(m) > l) &&
        throw(ArgumentError("invalid (s,l,m)=($s,$l,$m): need l>=|s| and |m|<=l"))
    if method === :SphericalExpansion || method === :Automatic
        lambda, coeffs, lmn, lmx = _swsh_spectral_coeffs(s, l, m, gamma; numterms=numterms)
        return SWSHSFunction(s, l, m, float(gamma), lambda, :SphericalExpansion, coeffs, lmn, lmx)
    elseif method === :Leaver
        lambda, coeffs, twok1, twok2 = _swsh_leaver_coeffs(s, l, m, gamma; numterms=numterms)
        return SWSHSLeaverFunction(s, l, m, float(gamma), lambda, :Leaver, coeffs, twok1, twok2)
    else
        throw(ArgumentError("unknown method $method (expected :SphericalExpansion or :Leaver)"))
    end
end

"""
    SpinWeightedSpheroidalHarmonicS(s, l, m, γ::Complex; method=:SphericalExpansion, numterms=:auto)

Complex-γ (`γ=aω`, complex QNM spheroidicity) value-tier carrier (Q1a). Only the
`:SphericalExpansion` representation is available at complex γ (the geev complex-symmetric spectral
coeffs); the Leaver VALUE representation at complex γ is [FUTURE] (needs a complex Kummer
normalisation). The eigensolve runs at `ComplexF64` — the value tier is Float64-capped by design
(§4.2), so a `Complex{BigFloat}` request is downcast; the BigFloat-inheriting surface is the
eigenvalue `SpinWeightedSpheroidalEigenvalue`.
"""
function SpinWeightedSpheroidalHarmonicS(s::Integer, l::Integer, m::Integer, gamma::Complex;
                                         method::Symbol=:SphericalExpansion, numterms=:auto)
    (l < abs(s) || abs(m) > l) &&
        throw(ArgumentError("invalid (s,l,m)=($s,$l,$m): need l>=|s| and |m|<=l"))
    if method === :SphericalExpansion || method === :Automatic
        g = ComplexF64(gamma)                         # value tier Float64-capped (geev is ComplexF64)
        lambda, coeffs, lmn, lmx = _swsh_spectral_coeffs(s, l, m, g; numterms=numterms)
        return SWSHSFunction(s, l, m, g, lambda, :SphericalExpansion, coeffs, lmn, lmx)
    elseif method === :Leaver
        throw(ArgumentError("SpinWeightedSpheroidalHarmonicS: the Leaver VALUE representation at " *
                            "complex γ is [FUTURE] (needs a complex Kummer/normalisation); use " *
                            ":SphericalExpansion"))
    else
        throw(ArgumentError("unknown method $method (expected :SphericalExpansion or :Leaver)"))
    end
end

"Evaluate `S(θ,φ) = (Σ_{l'} coeffs[l'] {}_sY_{l'm}(θ,0)) e^{imφ}`."
function (S::SWSHSFunction)(theta::Real, phi::Real)
    acc = zero(ComplexF64)
    for lp in S.lmin:S.lmax
        acc += S.coeffs[lp-S.lmin+1] * SpinWeightedSphericalHarmonicY(S.s, lp, S.m, theta, 0.0)
    end
    return acc * cis(S.m * phi)
end

"SphericalExpansion coefficients as `l' => coefficient` (mirrors `ExpansionCoefficients`)."
expansion_coefficients(S::SWSHSFunction) = Dict(lp => S.coeffs[lp-S.lmin+1] for lp in S.lmin:S.lmax)

"""
    swsh_derivative(S::SWSHSFunction, d1, d2, θ, φ) -> Complex

Mixed derivative `∂_θ^{d1} ∂_φ^{d2} {}_sS_{lm}(θ,φ)` (SphericalExpansion):
`(Σ_{l'} c_{l'} ∂_θ^{d1} {}_sY_{l'm}(θ,0)) · (im)^{d2} · e^{imφ}`.
"""
function swsh_derivative(S::SWSHSFunction, d1::Integer, d2::Integer, theta::Real, phi::Real)
    S.method === :SphericalExpansion ||
        throw(ArgumentError("swsh_derivative: SWSHSFunction must be :SphericalExpansion"))
    re = 0.0
    for lp in S.lmin:S.lmax
        re += S.coeffs[lp-S.lmin+1] * NumericsCore._y_dtheta(S.s, lp, S.m, theta, d1)
    end
    return re * (im * S.m)^d2 * cis(S.m * phi)
end

"""
    swsh_gamma_derivative(s, l, m, γ, d1, d2, θ, φ; numterms=:auto) -> Complex

Mixed spheroidicity/angle derivative `∂_γ ∂_θ^{d1} ∂_φ^{d2} {}_sS_{lm}(θ,φ)` (the `angparNew`
angular correction). Only the expansion coefficients depend on γ, so this is
`(Σ_{l'} (∂_γ c_{l'}) ∂_θ^{d1} {}_sY_{l'm}(θ,0)) · (im)^{d2} · e^{imφ}` with `∂_γ c` from
eigenvector perturbation. `∂_ω = a · ∂_γ`. SphericalExpansion only.
"""
function swsh_gamma_derivative(s::Integer, l::Integer, m::Integer, gamma::Real,
                               d1::Integer, d2::Integer, theta::Real, phi::Real; numterms=:auto)
    (l < abs(s) || abs(m) > l) &&
        throw(ArgumentError("invalid (s,l,m)=($s,$l,$m): need l>=|s| and |m|<=l"))
    dc, lmn, lmx = _swsh_coeffs_dgamma(s, l, m, gamma; numterms=numterms)
    re = 0.0
    for lp in lmn:lmx
        re += dc[lp-lmn+1] * NumericsCore._y_dtheta(s, lp, m, theta, d1)
    end
    return re * (im * m)^d2 * cis(m * phi)
end

"""
    swsh_gamma_derivative(s, l, m, γ::Complex, d1, d2, θ, φ; numterms=:auto)

Complex-γ spheroidicity derivative is [FUTURE]: the complex-symmetric eigenvector γ-derivative
carries a phase-normalisation drift term (the real path's sign fix is piecewise-constant, so its
drift vanishes — the complex phase-pin's does not). The complex VALUE + θ-derivative surface
(`SpinWeightedSpheroidalHarmonicS(::Complex)` + `swsh_derivative`) IS available.
"""
function swsh_gamma_derivative(s::Integer, l::Integer, m::Integer, gamma::Complex,
                               d1::Integer, d2::Integer, theta::Real, phi::Real; numterms=:auto)
    throw(ArgumentError("swsh_gamma_derivative at complex γ is [FUTURE] (complex eigenvector " *
                        "γ-derivative with phase-normalisation drift); the complex value + " *
                        "θ-derivative surface is available via SpinWeightedSpheroidalHarmonicS + " *
                        "swsh_derivative"))
end

# ============================== Leaver representation ==============================

"""
    SWSHSLeaverFunction

Callable carrier for `{}_sS_{lm}(·;γ)` via Leaver's series (Leaver 1985 Eq.18). `coeffs[n+1]` are
the sign-fixed, normalised series coefficients (n=0..nmax). Evaluate as `S(θ,φ)`; `swsh_derivative`
gives `∂_θ^{d1} ∂_φ^{d2} S`.
"""
struct SWSHSLeaverFunction
    s::Int
    l::Int
    m::Int
    gamma::Float64
    eigenvalue::Float64
    method::Symbol
    coeffs::Vector{Float64}        # n = 0..nmax
    twok1::Int                     # 2k1 = |m-s|
    twok2::Int                     # 2k2 = |m+s|
end

"""
    _leaver_coeff_precision() -> Int

Working `BigFloat` precision of the Leaver coefficient recurrence: `max(256, precision(BigFloat))`.

256 bits is a FLOOR, not a cap — the recurrence excites the dominant solution and the `normterm`
sum cancels over ~10⁹, so it needs far more than the caller may have active; but a caller who has
raised `setprecision` (BigFloat-512/1024 chains) must not have it silently truncated back to 256.
This is the same `max(256, precision(BigFloat))` policy used by `Perturbation.renorm_angmom` and
`MetricReconstruction.GHZ`.
"""
_leaver_coeff_precision() = max(256, precision(BigFloat))

"""
    _swsh_leaver_coeffs(s, l, m, γ; numterms=:auto, maxterms=200) -> (λ, coeffs, 2k1, 2k2)

Leaver-series coefficients & normalisation (Leaver 1985 Eq.20). `λ` from
[`SpinWeightedSpheroidalEigenvalue`](@ref); the `₁F₁` norm term sum (Pochhammer/Γ reduced to
integer factorials) is consumed from `NumericsCore.kummer_1F1`. Computed in `BigFloat` at
[`_leaver_coeff_precision`](@ref) (the forward recurrence excites the dominant solution and the
alternating `normterm` sum cancels over ~10⁹ — Float64 would lose the overall scale), then returned
as sign-fixed `Float64 coeffs/norm` (§4.2 value-tier ceiling).
"""
function _swsh_leaver_coeffs(s::Integer, l::Integer, m::Integer, γ::Real;
                             numterms=:auto, maxterms::Integer=200)
    λf = SpinWeightedSpheroidalEigenvalue(s, l, m, γ)
    twok1 = abs(m - s)
    twok2 = abs(m + s)
    coeffsbig = setprecision(BigFloat, _leaver_coeff_precision()) do
        g = BigFloat(γ); λ = BigFloat(λf)
        k1 = big(twok1) / 2; k2 = big(twok2) / 2
        αn(n) = -2 * (n + 1) * (n + twok1 + 1)
        βn(n) = n * (n - 1) + 2 * n * (k1 + k2 + 1 - 2g) -
                (2g * (2k1 + s + 1) - (k1 + k2) * (k1 + k2 + 1)) -
                (s * (s + 1) + λ + 2 * m * g)
        γn(n) = 2g * (n + k1 + k2 + s)

        a = BigFloat[one(BigFloat)]
        push!(a, -βn(0) / αn(0))
        grow!(nmax) = while length(a) <= nmax
            n = length(a)
            push!(a, (-βn(n - 1) * a[n] - γn(n - 1) * a[n - 1]) / αn(n - 1))
        end

        function normterm(i::Integer)
            grow!(i)
            poch = BigFloat(factorial(big(i + twok1)) // factorial(big(i + 1 + twok1 + twok2)))
            m1f1 = kummer_1F1(big(1 + i + twok1), big(i + 2 + twok1 + twok2), 4 * g)
            conv = sum(a[j + 1] * a[i - j + 1] for j in 0:i)
            return big(2)^i * poch * m1f1 * conv
        end

        nmax = 1
        normsum = normterm(0)
        if numterms === :auto
            while nmax < maxterms
                t = normterm(nmax)
                normsum += t
                abs(t) <= eps(Float64) * abs(normsum) && break
                nmax += 1
            end
        else
            nmax = Int(numterms) - 1
            normsum = sum(normterm(i) for i in 0:nmax)
        end
        grow!(nmax)

        norm = sqrt(2 * BigFloat(π)) *
               sqrt(big(2)^(1 + twok1 + twok2) * exp(-2g) * BigFloat(factorial(big(twok2))) * normsum)

        sgn = ((isodd(l) && iseven(m + s)) ||
               (isodd(l) && isodd(m + s) && m >= s) ||
               (iseven(l) && isodd(m - s) && m <= s)) ? -one(BigFloat) : one(BigFloat)

        (sgn / norm) .* a[1:(nmax + 1)]
    end
    return float(λf), Float64.(coeffsbig), twok1, twok2
end

"""
    _leaver_F_terms(coeffs, 2k1, 2k2) -> Vector{(coeff, a, b)}

Half-angle term list for the non-exponential Leaver factor
`F(θ) = Σ_n c_n 2^{(2k1+2k2)/2+n} sin(θ/2)^{2k2} cos(θ/2)^{2k1+2n}` as `coeff·sin(θ/2)^a cos(θ/2)^b`.
Then `S = e^{γcosθ}·F·e^{imφ}` (Leaver Eq.18).
"""
function _leaver_F_terms(coeffs::Vector{Float64}, twok1::Integer, twok2::Integer)
    terms = Tuple{Float64,Int,Int}[]
    for n in 0:(length(coeffs) - 1)
        c = coeffs[n + 1] * 2.0^((twok1 + twok2) / 2 + n)
        push!(terms, (c, twok2, twok1 + 2n))
    end
    return terms
end

"Evaluate `{}_sS_{lm}(θ,φ)` via Leaver's series."
function (S::SWSHSLeaverFunction)(theta::Real, phi::Real)
    terms = _leaver_F_terms(S.coeffs, S.twok1, S.twok2)
    F = NumericsCore._eval_terms(1.0, terms, theta)
    return exp(S.gamma * cos(theta)) * F * cis(S.m * phi)
end

"""
    swsh_derivative(S::SWSHSLeaverFunction, d1, d2, θ, φ) -> Complex

Mixed derivative `∂_θ^{d1} ∂_φ^{d2} {}_sS_{lm}(θ,φ)` for the Leaver series. With
`S = e^{γcosθ}·F(θ)·e^{imφ}` and `F` an exact half-angle polynomial: product-rule the `e^{γcosθ}`
factor (`E'=-γsinθ·E`, `E''=(γ²sin²θ-γcosθ)·E`) against the exact `F^{(d1)}`. Supports `d1∈{0,1,2}`.
"""
function swsh_derivative(S::SWSHSLeaverFunction, d1::Integer, d2::Integer, theta::Real, phi::Real)
    terms = _leaver_F_terms(S.coeffs, S.twok1, S.twok2)
    γ = S.gamma
    E = exp(γ * cos(theta))
    st = sin(theta); ct = cos(theta)
    if d1 == 0
        re = E * NumericsCore._eval_terms(1.0, terms, theta)
    elseif d1 == 1
        F = NumericsCore._eval_terms(1.0, terms, theta)
        dterms = NumericsCore._dtheta_terms(terms)
        Fp = NumericsCore._eval_terms(1.0, dterms, theta)
        re = E * (-γ * st * F + Fp)
    elseif d1 == 2
        F = NumericsCore._eval_terms(1.0, terms, theta)
        dterms = NumericsCore._dtheta_terms(terms)
        Fp = NumericsCore._eval_terms(1.0, dterms, theta)
        Fpp = NumericsCore._eval_terms(1.0, NumericsCore._dtheta_terms(dterms), theta)
        re = E * ((γ^2 * st^2 - γ * ct) * F - 2 * γ * st * Fp + Fpp)
    else
        throw(ArgumentError("Leaver swsh_derivative supports d1∈{0,1,2}, got d1=$d1"))
    end
    return re * (im * S.m)^d2 * cis(S.m * phi)
end

# ============================== small-γ harmonic Series ==============================

"""
    SpinWeightedSpheroidalHarmonicSSeries(s, l, m, θ, φ; order=6) -> Vector{ComplexF64}

Small-γ Taylor coefficients `[a_0,…,a_order]` of `{}_sS_{lm}(·;γ)(θ,φ)` about `γ=0`, i.e.
`S = Σ_i a_i γ^i + O(γ^{order+1})`. Each `a_i(θ,φ) = Σ_{j=-i}^{i} d[s,l,m][i,j] · {}_sY_{l+j,m}(θ,φ)`
(zero-ing `l+j < max(|s|,|m|)`), from the module-3 `_swsh_series_tables` d-tables and the module-1
Y seed. `a_0 = {}_sY_{lm}` (the γ→0 seed).
"""
function SpinWeightedSpheroidalHarmonicSSeries(s::Integer, l::Integer, m::Integer,
        theta::Real, phi::Real; order::Integer=6)
    order >= 0 || throw(ArgumentError("order must be ≥ 0; got $order"))
    _, Dd = _swsh_series_tables(s, l, m, order)
    coeffs = Vector{ComplexF64}(undef, order + 1)
    for i in 0:order
        a = 0.0 + 0.0im
        for j in -i:i
            d = get(Dd, (i, j), 0.0 + 0.0im)
            (l + j < abs(s) || l + j < abs(m)) && continue
            a += d * SpinWeightedSphericalHarmonicY(s, l + j, m, theta, phi)
        end
        coeffs[i + 1] = a
    end
    return coeffs
end
