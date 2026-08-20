# Complex{BigFloat}-correct Γ / logΓ + factorial. SpecialFunctions' Γ/logΓ are reliable only on
# the right half-plane Re z ≥ 1/2 (a real non-positive argument falls through to the real branch
# that raises DomainError on older versions, and the analytic continuation loses the validated
# precision on newer ones). The Euler reflection Γ(z)=π/(sin(πz)·Γ(1−z)) routes every Re z<1/2
# argument to Γ(1−z) (Re ≥ 1/2), which is where upstream is trustworthy. This consolidates
# Teukolsky `_gamma_big`/`_g`/`_lg`/`_needs_refl` + Punctures/H1Lorenz `_fact` into one canonical
# type-generic pair, so BigFloat precision is inherited (never hard-Float64'd).

"""
    needs_refl(z) -> Bool

True iff `z` is a `Complex{BigFloat}` on the real non-positive axis — the precise argument set
on which `SpecialFunctions.{gamma,loggamma}` fall through to the real branch and fail. A
documentary predicate; `gamma_big`/`loggamma_big` reflect on the broader `Re z < 1/2` to stay
on the trustworthy half-plane for all complex arguments.
"""
needs_refl(z) = z isa Complex{BigFloat} && imag(z) == 0 && real(z) <= 0

@inline _is_gamma_pole(z::Number) =
    iszero(imag(z)) && isfinite(real(z)) && real(z) <= zero(real(z)) && isinteger(real(z))
@inline _complex_pole(z::Number) = begin
    T = typeof(float(real(z)))
    Complex{T}(T(Inf), zero(T))
end
@inline _gamma_pole(z::Real) = oftype(float(z), Inf)
@inline _gamma_pole(z::Complex) = _complex_pole(z)

"""
    gamma_big(z) -> typeof(z)

Γ(z) valid for ALL `z` (incl. negative real part / real-negative half-integers). Exact
nonpositive integers are poles; because the signed infinity depends on approach direction, this API
uses the direction-independent convention `+Inf` (or `+Inf + 0im` for complex input). For
`Re z < 1/2` away from a pole the Euler reflection `Γ(z)=π/(sin(πz)·Γ(1−z))` is used (the inner Γ
then sees `Re ≥ 1/2`); otherwise the direct `SpecialFunctions.gamma`. Type-generic, so
`Complex{BigFloat}` precision is inherited.
"""
function gamma_big(z::Number)
    _is_gamma_pole(z) && return _gamma_pole(z)
    if real(z) < oftype(real(z), 1) / 2
        p = oftype(real(z), pi)
        return p / (sin(p * z) * SpecialFunctions.gamma(one(z) - z))
    end
    return SpecialFunctions.gamma(z)
end

"""
    loggamma_big(z) -> Complex

logΓ(z) for `Complex{BigFloat}` via the log-reflection on `Re z < 1/2`:
`logΓ(z)=log π − log sin(πz) − logΓ(1−z)`. The branch may differ from the principal logΓ by a
`2πi` multiple, but `exp(loggamma_big(z)) == gamma_big(z)` away from poles. At exact
nonpositive integers it follows the same direction-independent convention as `gamma_big`, returning
`+Inf + 0im`.
"""
function loggamma_big(z::Number)
    _is_gamma_pole(z) && return _complex_pole(z)
    zc = complex(z)
    if real(zc) < oftype(real(zc), 1) / 2
        p = complex(oftype(real(zc), pi))
        return log(p) - log(sin(p * zc)) - SpecialFunctions.loggamma(one(zc) - zc)
    end
    return SpecialFunctions.loggamma(zc)
end

"""
    fact_big(n) -> typeof(gamma_big(n+1))

Factorial as `Γ(n+1)` (accepts non-integer / complex `n`); the Γ-based factorial Punctures uses.
"""
fact_big(n) = gamma_big(n + one(n))

"""
    factorial_big(::Type{T}, n::Integer) -> T

Exact integer factorial `n!` accumulated in `T` (e.g. `BigFloat`) — overflow-safe at high `l`
where an `Int64` product would wrap. `factorial_big(T, 0) = factorial_big(T, 1) = one(T)`.
The BigSpecial name of the single body `NumericsCore.factorial_int` (forwarded, not re-implemented).
"""
factorial_big(::Type{T}, n::Integer) where {T} = NumericsCore.factorial_int(T, n)
