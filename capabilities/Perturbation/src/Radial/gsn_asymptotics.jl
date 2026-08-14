# :GSN backend — asymptotic-expansion coefficients of the GSN solutions at infinity and at the
# horizon (GSN paper arXiv:2306.16469 App. of Eqs. (33)-(41): X ~ e^{±iωr*}·Σ aᵢ/(ωr)^i at ∞,
# X ~ e^{±ipr*}·Σ aᵢ(r−r₊)^i at the horizon). Two REDESIGNS vs upstream
# src/Homogeneous/AsymptoticExpansionCoefficients.jl (v0.8.5 commit 982b087), both per
# GSN_PLAN §1/§3c and validated row-by-row against the captured oracle fixtures:
#   1. RECURRENCE-ONLY: upstream duplicated hand-coded closed forms for ∞-orders 1..3 (with
#      Float64 rational literals 1/6, 1/8, 1/48 — an audited precision leak); the general
#      recurrences below generate every order ≥ 1 from a₀ = 1, so the closed forms are dropped.
#      The fixture net FIX_ASYMP_COEFFS spans orders 0..8, i.e. BOTH upstream regimes.
#   2. NO TaylorSeries.jl (dep-hygiene ruling OQ-3): the P/Q coefficient functions are evaluated
#      on a minimal in-house truncated-Taylor arithmetic (_GSNTaylor below). P/Q are rational in
#      the expansion variable (verified: no sqrt/exp of the series variable in any branch), so
#      +,−,×,÷,^Int is a complete operation set.
#   3. TYPE-PARAMETRIC CACHES (audited hard blocker GSN_PLAN §3c-2): upstream's module-global
#      caches were annotated Vector{ComplexF64} and silently downcast BigFloat coefficients on
#      store; here the cache is keyed by (kind, s, m, a, ω, λ, T, precision) under a lock.
#
# The P/Q generator bodies at the end of this file are MECHANICALLY EXTRACTED verbatim
# (rename I→im only; scratchpad/gsn-port/extract_gsn.py).

# ----------------------------------------------------------------------------------------------
# Minimal truncated-Taylor arithmetic (fixed truncation length; c[k+1] = coeff of t^k)
# ----------------------------------------------------------------------------------------------

"""
    _GSNTaylor{T}

Truncated Taylor series Σₖ c[k+1]·tᵏ over Complex{T}, fixed length. Supports exactly the
operation set the GSN P/Q rational bodies need: ±, ×, ÷ (series division), ^Int≥0, and mixed
scalar ops. In-house replacement for TaylorSeries.taylor_expand (OQ-3 ruling).
"""
struct _GSNTaylor{T<:Real}
    c::Vector{Complex{T}}
end

_gsn_taylor_var(::Type{T}, N::Int) where {T<:Real} = begin
    c = zeros(Complex{T}, N + 1)
    N >= 1 && (c[2] = one(Complex{T}))
    _GSNTaylor{T}(c)
end
_gsn_taylor_const(x::Number, N::Int, ::Type{T}) where {T<:Real} = begin
    c = zeros(Complex{T}, N + 1); c[1] = Complex{T}(x); _GSNTaylor{T}(c)
end
_gsn_getcoeff(t::_GSNTaylor, k::Int) = t.c[k + 1]

Base.:+(x::_GSNTaylor{T}, y::_GSNTaylor{T}) where {T} = _GSNTaylor{T}(x.c .+ y.c)
Base.:-(x::_GSNTaylor{T}, y::_GSNTaylor{T}) where {T} = _GSNTaylor{T}(x.c .- y.c)
Base.:-(x::_GSNTaylor{T}) where {T} = _GSNTaylor{T}(-x.c)
Base.:+(x::_GSNTaylor{T}, y::Number) where {T} = begin c = copy(x.c); c[1] += y; _GSNTaylor{T}(c) end
Base.:+(y::Number, x::_GSNTaylor) = x + y
Base.:-(x::_GSNTaylor, y::Number) = x + (-y)
Base.:-(y::Number, x::_GSNTaylor) = (-x) + y
Base.:*(x::_GSNTaylor{T}, y::Number) where {T} = _GSNTaylor{T}(x.c .* Complex{T}(y))
Base.:*(y::Number, x::_GSNTaylor) = x * y
Base.:/(x::_GSNTaylor{T}, y::Number) where {T} = _GSNTaylor{T}(x.c ./ Complex{T}(y))

function Base.:*(x::_GSNTaylor{T}, y::_GSNTaylor{T}) where {T}
    n = length(x.c)
    c = zeros(Complex{T}, n)
    @inbounds for i in 1:n, j in 1:(n - i + 1)
        c[i + j - 1] += x.c[i] * y.c[j]
    end
    _GSNTaylor{T}(c)
end

"Series division (denominator must have a nonzero constant term — true for every GSN P/Q branch)."
function Base.:/(x::_GSNTaylor{T}, y::_GSNTaylor{T}) where {T}
    n = length(x.c)
    iszero(y.c[1]) && throw(DomainError("series division by a series with zero constant term"))
    c = zeros(Complex{T}, n)
    @inbounds for k in 0:(n - 1)
        acc = x.c[k + 1]
        for j in 1:k
            acc -= y.c[j + 1] * c[k - j + 1]
        end
        c[k + 1] = acc / y.c[1]
    end
    _GSNTaylor{T}(c)
end
Base.:/(y::Number, x::_GSNTaylor{T}) where {T} = _gsn_taylor_const(y, length(x.c) - 1, T) / x

function Base.:^(x::_GSNTaylor{T}, p::Integer) where {T}
    p >= 0 || throw(DomainError(p, "negative Taylor power unsupported (and unused by GSN P/Q)"))
    r = _gsn_taylor_const(one(T), length(x.c) - 1, T)
    b = x
    while p > 0
        isodd(p) && (r = r * b)
        b = b * b
        p >>= 1
    end
    r
end

# ----------------------------------------------------------------------------------------------
# Coefficient recurrences (∞: GSN upstream l.361-412/763-814; horizon Frobenius: l.1205-1260/
# 1651-1706 — transcribed with typed accumulators, GSN_PLAN §3c-6/7)
# ----------------------------------------------------------------------------------------------

"""
    _gsn_coeff_series(kind, s, m, a, omega, lambda, order, ::Type{T}) -> Vector{Complex{T}}

Full coefficient list a₀..a_order for one boundary/direction. `kind ∈ (:in_inf, :out_inf,
:in_hor, :out_hor)`. ∞ recurrence: aᵢ (stored ω^i-scaled, as upstream) from all previous terms,
    aᵢ = ω^i·[ i(i−1)·aᵢ₋₁/ω^{i−1} + Σₖ₌₁^i (Q_{k+1} − (i−k)P_k)·a_{i−k}/ω^{i−k} ]/(P₀·i).
Horizon (Frobenius): aᵢ = −Σᵣ₌₀^{i−1} aᵣ(r·P_{i−r} + Q_{i−r}) / [i(i−1) + P₀·i + Q₀].
[GSN AsymptoticExpansionCoefficients.jl recurrences, verbatim algebra]
"""
function _gsn_coeff_series(kind::Symbol, s::Int, m::Int, a, omega, lambda, order::Int, ::Type{T}) where {T<:Real}
    coeffs = Vector{Complex{T}}(undef, order + 1)
    coeffs[1] = one(Complex{T})
    order == 0 && return coeffs
    aT = T(a); omegaT = T(omega); lamT = Complex{T}(lambda)
    if kind === :in_inf || kind === :out_inf
        zP = _gsn_taylor_var(T, order)
        zQ = _gsn_taylor_var(T, order + 1)
        P = kind === :in_inf ? _gsn_P_in_inf(s, m, aT, omegaT, lamT, zP) :
                               _gsn_P_out_inf(s, m, aT, omegaT, lamT, zP)
        Q = kind === :in_inf ? _gsn_Q_in_inf(s, m, aT, omegaT, lamT, zQ) :
                               _gsn_Q_out_inf(s, m, aT, omegaT, lamT, zQ)
        P0 = _gsn_getcoeff(P, 0)
        for i in 1:order
            acc = zero(Complex{T})
            for k in 1:i
                acc += (_gsn_getcoeff(Q, k + 1) - (i - k)*_gsn_getcoeff(P, k)) * (coeffs[i - k + 1]/omegaT^(i - k))
            end
            coeffs[i + 1] = omegaT^i * ((i*(i - 1)*(coeffs[i]/omegaT^(i - 1)) + acc)/(P0*i))
        end
    elseif kind === :in_hor || kind === :out_hor
        xP = _gsn_taylor_var(T, order)
        P = kind === :in_hor ? _gsn_P_in_hor(s, m, aT, omegaT, lamT, xP) :
                               _gsn_P_out_hor(s, m, aT, omegaT, lamT, xP)
        Q = kind === :in_hor ? _gsn_Q_in_hor(s, m, aT, omegaT, lamT, xP) :
                               _gsn_Q_out_hor(s, m, aT, omegaT, lamT, xP)
        P0 = _gsn_getcoeff(P, 0); Q0 = _gsn_getcoeff(Q, 0)
        indicial(nu) = nu*(nu - 1) + P0*nu + Q0
        for i in 1:order
            acc = zero(Complex{T})
            for r in 0:(i - 1)
                acc += coeffs[r + 1]*(r*_gsn_getcoeff(P, i - r) + _gsn_getcoeff(Q, i - r))
            end
            coeffs[i + 1] = -acc/indicial(i)
        end
    else
        throw(ArgumentError("unknown asymptotic-coefficient kind :$kind"))
    end
    return coeffs
end

# Per-(mode, type, precision) cache — the {T}-generic replacement for the upstream ComplexF64-
# annotated module globals (audited BigFloat hard blocker, GSN_PLAN §3c-2). Lock-guarded
# (same pattern as radial_sn's _SN_LOCK).
const _GSN_COEFF_LOCK = ReentrantLock()
const _GSN_COEFF_CACHE = Dict{Tuple{Symbol,Int,Int,Any,Any,Any,DataType,Int},Vector}()

function _gsn_asym_coeff(kind::Symbol, s::Int, m::Int, a, omega, lambda, order::Int, ::Type{T}) where {T<:Real}
    order >= 0 || throw(DomainError(order, "Only positive expansion order is supported"))
    prec = T <: BigFloat ? precision(BigFloat) : 0
    key = (kind, s, m, a, omega, lambda, T, prec)
    lock(_GSN_COEFF_LOCK) do
        v = get(_GSN_COEFF_CACHE, key, nothing)
        if v === nothing || length(v) < order + 1
            v = _gsn_coeff_series(kind, s, m, a, omega, lambda, order, T)
            _GSN_COEFF_CACHE[key] = v
        end
        return (v::Vector{Complex{T}})[order + 1]
    end
end

"Ingoing ∞-expansion coefficient aᵢ (ω^i-scaled). [GSN ingoing_coefficient_at_inf]"
_gsn_coeff_in_inf(s::Int, m::Int, a, omega, lambda, order::Int, ::Type{T}) where {T} =
    _gsn_asym_coeff(:in_inf, s, m, a, omega, lambda, order, T)
"Outgoing ∞-expansion coefficient. [GSN outgoing_coefficient_at_inf]"
_gsn_coeff_out_inf(s::Int, m::Int, a, omega, lambda, order::Int, ::Type{T}) where {T} =
    _gsn_asym_coeff(:out_inf, s, m, a, omega, lambda, order, T)
"Ingoing horizon-expansion coefficient. [GSN ingoing_coefficient_at_hor]"
_gsn_coeff_in_hor(s::Int, m::Int, a, omega, lambda, order::Int, ::Type{T}) where {T} =
    _gsn_asym_coeff(:in_hor, s, m, a, omega, lambda, order, T)
"Outgoing horizon-expansion coefficient. [GSN outgoing_coefficient_at_hor]"
_gsn_coeff_out_hor(s::Int, m::Int, a, omega, lambda, order::Int, ::Type{T}) where {T} =
    _gsn_asym_coeff(:out_hor, s, m, a, omega, lambda, order, T)

# ----------------------------------------------------------------------------------------------
# P/Q generator functions (mechanically extracted; see file header)
# ----------------------------------------------------------------------------------------------

function _gsn_P_in_inf(s::Int, m::Int, a, omega, lambda, z)
    if s == 0
        return begin
            -((2*im*(omega + a^2*z^4*(-im + a^2*omega) + z^2*(im + 2*a^2*omega)))/
            ((1 + a^2*z^2)*(1 - 2*z + a^2*z^2)))
        end
    elseif s == 1
        return begin
            ((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - (2*z*(1 - 2*z + a^2*z^2))/
            (1 + a^2*z^2)^2 + (2*a*z^2*(1 - 2*z + a^2*z^2)*((-im)*m*(1 + 3*a^2*z^2) + 
            a*z*(3 + 3*z + 2*lambda + 2*a^2*z^2*(1 + lambda))))/((1 + a^2*z^2)*
            (2 - 2*im*a*m*z - 2*im*a^3*m*z^3 + lambda + a^4*z^4*(1 + lambda) + 
            a^2*z^2*(3 + 2*z + 2*lambda))) - 2*im*omega))/(1 - 2*z + a^2*z^2)
        end
    elseif s == -1
        return begin
            (1/(1 - 2*z + a^2*z^2))*((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - 
            (2*z*(1 - 2*z + a^2*z^2))/(1 + a^2*z^2)^2 + 
            (2*a*z^2*(1 - 2*z + a^2*z^2)*(im*m*(1 + 3*a^2*z^2) + 
            a*z*(-1 + 3*z + 2*a^2*z^2*(-1 + lambda) + 2*lambda)))/
            ((1 + a^2*z^2)*(2*im*a*m*z + 2*im*a^3*m*z^3 + a^4*z^4*(-1 + lambda) + lambda + 
            a^2*z^2*(-1 + 2*z + 2*lambda))) - 2*im*omega))
        end
    elseif s == 2
        return begin
            (1/(1 - 2*z + a^2*z^2))*((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - 
            (2*z*(1 - 2*z + a^2*z^2))/(1 + a^2*z^2)^2 - 2*im*omega + 
            (8*a*z^2*(1 - 2*z + a^2*z^2)*(-6*a*m^2*z + im*m*(-4 + 9*a^2*z^2 - lambda + 
            z*(-6 - 12*im*a^2*omega)) + a*(6*a^2*z^3 + im*(1 + lambda)*omega + z^2*(-9 - 9*im*a^2*omega) + 
            z*(3 + 6*im*omega - 6*a^2*omega^2))))/((1 + a^2*z^2)*(24 + 10*lambda + lambda^2 - 
            4*im*a*m*(6*z^2 + 2*z*(4 + lambda) + 3*im*omega) + 12*im*omega + 24*a^3*m*z^2*(im*z + 2*omega) + 
            12*a^4*z^2*(z^2 - 2*im*z*omega - 2*omega^2) - 4*a^2*(6*z^3 + z^2*(-3 + 6*m^2 - 6*im*omega) - 
            2*im*z*(1 + lambda)*omega + 3*omega^2)))))
        end
    elseif s == -2
        return begin
            (1/(1 - 2*z + a^2*z^2))*((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - (2*z*(1 - 2*z + a^2*z^2))/(1 + a^2*z^2)^2 - 2*im*omega + 
            (8*a*z^2*(1 - 2*z + a^2*z^2)*(im*m*(6*z + lambda) + 3*a^2*m*z*(-3*im*z + 4*omega) - 
            a*(9*z^2 + z*(-3 + 6*m^2 + 6*im*omega) + im*(-3 + lambda)*omega) + 3*a^3*z*(2*z^2 + 3*im*z*omega - 2*omega^2)))/
            ((1 + a^2*z^2)*(2*lambda + lambda^2 - 12*im*omega + 24*a^3*m*z^2*((-im)*z + 2*omega) + 4*a*m*(6*im*z^2 + 2*im*z*lambda + 3*omega) + 
            12*a^4*z^2*(z^2 + 2*im*z*omega - 2*omega^2) - 4*a^2*(6*z^3 + z^2*(-3 + 6*m^2 + 6*im*omega) + 2*im*z*(-3 + lambda)*omega + 3*omega^2)))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_Q_in_inf(s::Int, m::Int, a, omega, lambda, z)
    if s == 0
        return begin
            -((z^2*(-2*z*(-1 + lambda) + lambda - 4*a^2*z^3*(1 + lambda) - 2*a^4*z^5*(3 + lambda) + 2*a*m*omega + 
            a^6*z^6*(1 - m^2 + lambda + 2*a*m*omega) + a^2*z^4*(8 + a^2*(2 - 2*m^2 + 3*lambda) + 
            6*a^3*m*omega) + z^2*(-4 + a^2*(1 - m^2 + 3*lambda) + 6*a^3*m*omega)))/
            ((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2))
        end
    elseif s == 1
        return begin
            (z^2*(-4 - lambda^2 + 2*a^12*z^12*(1 + lambda) - 8*a*m*omega - 4*im*a^10*z^11*(1 + lambda)*
            (-3*im + a^2*omega) - 2*lambda*(2 + a*m*omega) - a^4*z^6*(52 + 32*lambda + 2*im*a*m*(31 + 6*lambda) + 
            a^2*(71 + 10*lambda^2 - 3*m^2*(3 + 2*lambda) + 8*lambda*(7 - 6*im*omega)) + 
            4*a^3*m*(17 + 5*lambda)*omega) + 2*a^6*z^9*(8 - 6*im*a*m + 
            a^2*(-3 + 3*m^2 - 4*lambda + lambda^2 + 4*im*omega) + 2*im*a^4*(m^2 - 5*(1 + lambda))*omega + 
            a^3*m*(4*im - im*m^2 + im*lambda + 2*omega)) + 
            2*z*((2 + lambda)^2 + 2*im*a^2*(-1 + m^2 - lambda)*omega + a*m*(2*im + im*lambda + 4*omega)) + 
            2*a^3*z^5*(14*im*m + a*(69 - m^2 + 48*lambda + 6*lambda^2 + 28*im*omega) + 
            4*im*a^3*(3*m^2 - 5*(1 + lambda))*omega + 3*a^2*m*(6*im - im*m^2 + 2*im*lambda + 6*omega)) + 
            2*a^4*z^7*(-4 + 8*im*a*m + 4*a^2*(10 + m^2 + 7*lambda + lambda^2 + 5*im*omega) + 
            4*im*a^4*(2*m^2 - 5*(1 + lambda))*omega + a^3*m*(14*im - 3*im*m^2 + 4*im*lambda + 10*omega)) + 
            2*a^2*z^3*(36 - 2*m^2 + 26*lambda + 4*lambda^2 + 12*im*omega + 2*im*a^2*(4*m^2 - 5*(1 + lambda))*omega + 
            a*m*(10*im - im*m^2 + 4*im*lambda + 14*omega)) - 
            a^8*z^10*(-6*im*a*m - 8*(1 + 2*lambda) + 2*a^3*m*(3 + lambda)*omega + 
            a^2*(1 - m^2*(-1 + lambda) + lambda^2 - 4*im*omega - 8*im*lambda*omega)) - 
            a*z^2*(4*im*m*(2 + lambda) + 2*a^2*m*(19 + 5*lambda)*omega + 
            a*(28 + 24*lambda + 5*lambda^2 - m^2*(4 + lambda) + 4*im*omega - 8*im*lambda*omega)) - 
            a^2*z^4*(24*(2 + lambda) + 2*im*a*m*(23 + 6*lambda) + 4*a^3*m*(18 + 5*lambda)*omega + 
            a^2*(67 + 54*lambda + 10*lambda^2 - m^2*(11 + 4*lambda) + 8*im*omega - 32*im*lambda*omega)) + 
            a^6*z^8*(4 + 8*lambda - 2*im*a*m*(9 + 2*lambda) - 2*a^3*m*(16 + 5*lambda)*omega + 
            a^2*(-31 - 24*lambda - 5*lambda^2 + m^2*(1 + 4*lambda) + 8*im*omega + 32*im*lambda*omega))))/
            ((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*(2 - 2*im*a*m*z - 2*im*a^3*m*z^3 + lambda + 
            a^4*z^4*(1 + lambda) + a^2*z^2*(3 + 2*z + 2*lambda)))
        end
    elseif s == -1
        return begin
            (z^2*((-1 + 2*z)*lambda^2 + a^2*z*(m^2*(-4*z^2 + z*(2 + lambda) - 4*im*omega) + 
            lambda*(-24*z^3 + 4*z^2*(5 + 2*lambda) + z*(-4 - 5*lambda + 8*im*omega) - 4*im*omega)) + 
            2*im*a*m*lambda*(-z + 2*z^2 + im*omega) + 2*a^12*z^11*(-1 + lambda)*(z - 2*im*omega) - 
            2*a^11*m*z^10*(-5 + lambda)*omega + 2*a^9*m*z^8*(-3*im*z^2 + im*z*(-2 + m^2 - lambda + 10*im*omega) - 
            5*(-4 + lambda)*omega) + 2*im*a^7*m*z^6*(6*z^3 + z^2*(5 + 2*lambda) + 
            z*(-6 + 3*m^2 - 4*lambda + 30*im*omega) + 10*im*(-3 + lambda)*omega) + 
            2*im*a^5*m*z^4*(-8*z^3 + z^2*(19 + 6*lambda) + 3*z*(-2 + m^2 - 2*lambda + 10*im*omega) + 
            10*im*(-2 + lambda)*omega) + 2*im*a^3*m*z^2*(-14*z^3 + z^2*(11 + 6*lambda) + 
            z*(-2 + m^2 - 4*lambda + 10*im*omega) + 5*im*(-1 + lambda)*omega) - 
            a^10*z^9*(12*z^2*(-1 + lambda) + 4*im*(-4 + m^2 + 5*lambda)*omega + 
            z*(5 - m^2*(-3 + lambda) - 4*lambda + lambda^2 + 16*im*omega - 8*im*lambda*omega)) + 
            a^4*z^3*(-8*z^4 - 4*z^3*(-3 + 8*lambda) - 2*z^2*(3 + m^2 - 24*lambda - 6*lambda^2 - 8*im*omega) - 
            4*im*(-1 + 4*m^2 + 5*lambda)*omega + z*(1 - 14*lambda - 10*lambda^2 + m^2*(3 + 4*lambda) - 16*im*omega + 
            32*im*lambda*omega)) + a^8*z^7*(8*z^3*(-3 + 2*lambda) + 
            2*z^2*(9 + 3*m^2 - 8*lambda + lambda^2 + 8*im*omega) - 8*im*(-3 + 2*m^2 + 5*lambda)*omega + 
            z*(-3 - 4*lambda - 5*lambda^2 + m^2*(-7 + 4*lambda) - 48*im*omega + 32*im*lambda*omega)) + 
            a^6*z^5*(16*z^4 + 4*z^3*(-3 + 2*lambda) + 8*z^2*(m^2 + 3*lambda + lambda^2 + 4*im*omega) - 
            8*im*(-2 + 3*m^2 + 5*lambda)*omega + z*(1 - 16*lambda - 10*lambda^2 + m^2*(-3 + 6*lambda) - 48*im*omega + 
            48*im*lambda*omega))))/((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*
            (2*im*a*m*z + 2*im*a^3*m*z^3 + a^4*z^4*(-1 + lambda) + lambda + a^2*z^2*(-1 + 2*z + 2*lambda)))
        end
    elseif s == 2
        return begin
            (z^2*(24*a^12*z^12 - lambda^3 - 48*im*a^10*z^11*(-3*im + a^2*omega) - 2*lambda^2*(8 + a*m*omega) + 
            4*lambda*(-21 - 3*(im + 4*a*m)*omega + 7*a^2*omega^2) - 
            12*a^8*z^10*(-24 + 6*im*a*m - a^2*(5 + 3*m^2 + 20*im*omega) + 4*a^4*omega^2) + 
            8*im*a^6*z^9*(24*im + 18*a*m + 3*im*a^2*(9 + m^2 + 12*im*omega) + 
            a^3*m*(19 - 3*m^2 + 4*lambda + 18*im*omega) + a^4*(-31 + 15*m^2 - 4*lambda - 24*im*omega)*omega - 
            24*a^5*m*omega^2 + 12*a^6*omega^3) + 8*(-18 - (9*im + 23*a*m)*omega + 
            a*(-3*im*m + a*(11 - 3*m^2))*omega^2 + 3*a^3*m*omega^3) + 
            a^3*z^6*(384*im*m + 8*im*a^2*m*(-177 + 18*m^2 - 53*lambda - 2*lambda^2 + 24*im*omega) + 
            16*a*(-33 + 12*m^2 - 10*lambda - lambda^2 - 36*im*omega) + 
            a^3*(-756 - 48*m^4 - lambda^3 + lambda^2*(-41 + 16*im*omega) + lambda*(-334 + 364*im*omega) + 
            m^2*(564 + 106*lambda + lambda^2 - 756*im*omega) + 780*im*omega) + 
            2*a^4*m*(-506 + 126*m^2 - 132*lambda - lambda^2 + 588*im*omega)*omega + 
            4*a^5*(49 - 93*m^2 + 37*lambda - 144*im*omega)*omega^2 + 168*a^6*m*omega^3) + 
            a*z^4*(96*im*m + 8*im*a^2*m*(-202 + 9*m^2 - 65*lambda - 4*lambda^2) - 
            8*a*(120 + 50*lambda + 5*lambda^2 + 72*im*omega) + a^3*(-1068 - 24*m^4 - 3*lambda^3 + 
            lambda^2*(-74 + 32*im*omega) + 2*m^2*(246 + 58*lambda + lambda^2 - 288*im*omega) + 
            lambda*(-512 + 388*im*omega) + 488*im*omega) + 6*a^4*m*(-224 + 36*m^2 - 60*lambda - lambda^2 + 
            172*im*omega)*omega + 12*a^5*(38 - 34*m^2 + 17*lambda - 48*im*omega)*omega^2 + 216*a^6*m*omega^3) + 
            4*im*a*z^3*(12*m*(3 + lambda) + im*a*(-408 - 194*lambda - 27*lambda^2 - lambda^3 + 16*m^2*(4 + lambda) - 
            168*im*omega) - 12*a^4*m*(17 + lambda)*omega^2 + 96*a^5*omega^3 + 
            2*a^3*omega*(-53 - 26*lambda - 3*lambda^2 + m^2*(79 + 7*lambda) + 60*im*omega + 24*im*lambda*omega) - 
            2*a^2*m*(-98 - 35*lambda - 3*lambda^2 + m^2*(10 + lambda) + 128*im*omega + 32*im*lambda*omega)) + 
            2*im*a^3*z^5*(24*m*(9 + 4*lambda) - 24*a^4*m*(25 + lambda)*omega^2 + 288*a^5*omega^3 + 
            4*a^3*omega*(-89 - 35*lambda - 3*lambda^2 + m^2*(113 + 8*lambda) + 51*im*omega + 27*im*lambda*omega) - 
            4*a^2*m*(-143 - 44*lambda - 3*lambda^2 + m^2*(23 + 2*lambda) + 181*im*omega + 43*im*lambda*omega) + 
            im*a*(-924 - 394*lambda - 47*lambda^2 - lambda^3 + 4*m^2*(91 + 16*lambda) - 300*im*omega + 
            84*im*lambda*omega)) + 4*a^6*z^8*(36 - 24*m^2 + 6*a^3*m*(4*m^2 - 3*(4 + lambda - 6*im*omega))*
            omega - 2*a^4*(7 + 15*m^2 - 5*lambda + 24*im*omega)*omega^2 + 12*a^5*m*omega^3 + 
            2*a*m*(-44*im + 9*im*m^2 - 11*im*lambda + 36*omega) - 
            a^2*(39 + 6*m^4 + 20*lambda + 2*lambda^2 + m^2*(-65 - 8*lambda + 78*im*omega) - 142*im*omega - 
            22*im*lambda*omega + 48*omega^2)) + 8*im*a^4*z^7*(-12*im + 6*a*m*(-5 + lambda) + 
            im*a^2*(m^2*(62 + 8*lambda) - 3*(24 + lambda^2 + 2*lambda*(5 - im*omega) + 10*im*omega)) + 
            a^3*m*(88 + lambda^2 - m^2*(16 + lambda) + lambda*(23 - 18*im*omega) - 66*im*omega) - 
            2*a^5*m*(49 + lambda)*omega^2 + 48*a^6*omega^3 + a^4*omega*(-73 - lambda^2 + 3*m^2*(23 + lambda) - 
            14*im*omega + 10*im*lambda*(2*im + omega))) + 2*z*((9 + lambda)*(24 + 10*lambda + lambda^2 + 12*im*omega) - 
            8*im*a^3*m*(13 + lambda)*omega^2 + 48*im*a^4*omega^3 - 
            4*a^2*omega*(12*im + im*lambda^2 - 2*im*m^2*(10 + lambda) + 19*omega + 7*lambda*(im + omega)) + 
            4*a*m*(24*im + im*lambda^2 + 31*omega + lambda*(10*im + 7*omega))) + 
            z^2*(-8*im*a*m*(60 + 23*lambda + 2*lambda^2 + 6*im*omega) - 12*(24 + 10*lambda + lambda^2 + 12*im*omega) + 
            6*a^3*m*(-134 + 10*m^2 - 36*lambda - lambda^2 + 44*im*omega)*omega + 
            4*a^4*(85 - 45*m^2 + 31*lambda - 48*im*omega)*omega^2 + 120*a^5*m*omega^3 + 
            a^2*(-3*lambda^3 + lambda^2*(-57 + 16*im*omega) + lambda*(-342 + 100*im*omega) + 
            m^2*(152 + 42*lambda + lambda^2 - 132*im*omega) - 12*(54 + 3*im*omega + 4*omega^2)))))/
            ((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*(24 + 10*lambda + lambda^2 - 
            4*im*a*m*(6*z^2 + 2*z*(4 + lambda) + 3*im*omega) + 12*im*omega + 24*a^3*m*z^2*(im*z + 2*omega) + 
            12*a^4*z^2*(z^2 - 2*im*z*omega - 2*omega^2) - 4*a^2*(6*z^3 + z^2*(-3 + 6*m^2 - 6*im*omega) - 
            2*im*z*(1 + lambda)*omega + 3*omega^2)))
        end
    elseif s == -2
        return begin
            (z^2*(-((-1 + 2*z)*(-2 + 6*z - lambda)*(2*lambda + lambda^2 - 12*im*omega)) - 48*a^11*m*z^8*omega*(3*z^2 - omega^2) + 
            24*a^12*z^10*(z^2 - 2*im*z*omega + 4*omega^2) - 4*a^10*z^8*(36*z^3 - 3*z^2*(5 + 3*m^2 + 8*im*omega) + 6*(-17 + 5*m^2 - lambda)*omega^2 + 
            2*z*omega*(39*im + 3*im*m^2 - 4*im*lambda + 24*omega)) + 8*a^9*m*z^6*(9*im*z^4 + z^2*(-72 + 12*m^2 - 7*lambda - 6*im*omega)*omega + 
            2*im*z*(-3 + lambda)*omega^2 + 21*omega^3 + z^3*(-3*im + 3*im*m^2 - 4*im*lambda + 30*omega)) - 
            2*im*a^7*m*z^4*(72*z^5 + 4*z^4*(9*m^2 - 11*lambda + 12*im*omega) + im*z^2*(-426 + 126*m^2 - 92*lambda - lambda^2 - 12*im*omega)*omega - 
            24*z*(-3 + lambda)*omega^2 + 108*im*omega^3 + 4*z^3*(12 + 15*lambda + lambda^2 - m^2*(12 + lambda) + 150*im*omega + 14*im*lambda*omega)) - 
            2*im*a^3*m*(192*z^6 + 24*z^5*(-7 + 4*lambda) + 4*z^4*(-6 + 9*m^2 - 33*lambda - 4*lambda^2 - 72*im*omega) - 8*z*(-3 + lambda)*omega^2 + 12*im*omega^3 - 
            z^2*omega*(90*im - 30*im*m^2 + 52*im*lambda + 3*im*lambda^2 + 60*omega) + 4*z^3*(6 + 11*lambda + 3*lambda^2 - m^2*(6 + lambda) + 84*im*omega + 
            20*im*lambda*omega)) - 2*im*a^5*m*z^2*(24*z^5*(-9 + lambda) + 4*z^4*(3 + 18*m^2 - 37*lambda - 2*lambda^2 - 96*im*omega) + 
            3*im*z^2*(-96 + 36*m^2 - 36*lambda - lambda^2 + 20*im*omega)*omega - 24*z*(-3 + lambda)*omega^2 + 60*im*omega^3 + 
            4*z^3*(15 + 20*lambda + 3*lambda^2 - m^2*(15 + 2*lambda) + 189*im*omega + 31*im*lambda*omega)) + 
            4*a^8*z^6*(-6*(9 + m^2)*z^3 + 72*z^4 + z^2*(9 - 6*m^4 - 2*lambda^2 + lambda*(-4 - 22*im*omega) + m^2*(33 + 8*lambda + 30*im*omega) + 258*im*omega) + 
            3*(55 - 31*m^2 + 7*lambda)*omega^2 - 2*z*omega*(87*im - im*lambda^2 + 3*im*m^2*(3 + lambda) + 114*omega + 6*lambda*(-2*im + omega))) + 
            2*a*m*(-48*im*z^4 - 24*im*z^3*(-1 + lambda) + 4*im*z^2*(7*lambda + 2*lambda^2 + 18*im*omega) - (12 + 8*lambda + lambda^2 - 12*im*omega)*omega + 
            4*z*((-im)*lambda^2 + 15*omega + lambda*(-2*im + 3*omega))) - a^2*(8*z^4*(10*lambda + 5*lambda^2 - 72*im*omega) - 
            4*z^3*((26 - 16*m^2)*lambda + 15*lambda^2 + lambda^3 - 240*im*omega) + 12*(-2 + 2*m^2 - lambda)*omega^2 + 
            z^2*(3*lambda^3 + lambda^2*(21 - m^2 + 16*im*omega) + lambda*(30 - 34*m^2 - 28*im*omega) + 12*im*(-35 + 5*m^2 + 12*im*omega)*omega) + 
            8*z*omega*(6*im - im*lambda^2 + 15*omega + lambda*(im + 2*im*m^2 + 3*omega))) + 
            a^4*z^2*(96*z^5 + 16*z^4*(-9 + 12*m^2 - 2*lambda - lambda^2 + 36*im*omega) + 60*(3 - 3*m^2 + lambda)*omega^2 + 
            2*z^3*(36 + 66*lambda + 35*lambda^2 + lambda^3 - 4*m^2*(27 + 16*lambda) - 1068*im*omega + 84*im*lambda*omega) - 
            8*z*omega*(39*im - 3*im*lambda^2 + im*m^2*(3 + 7*lambda) + 84*omega + 2*lambda*(-im + 6*omega)) - 
            z^2*(24*m^4 + 3*lambda^3 - 2*m^2*(30 + 50*lambda + lambda^2) + lambda^2*(38 + 32*im*omega) + 4*lambda*(16 + 33*im*omega) - 
            12*(-1 + 134*im*omega + 48*omega^2))) - a^6*z^4*(48*(-3 + 2*m^2)*z^4 + 192*z^5 + 
            8*z^3*(m^2*(30 + 8*lambda) - 3*(lambda^2 + lambda*(2 + 2*im*omega) - 54*im*omega)) + 12*(-42 + 34*m^2 - 9*lambda)*omega^2 + 
            8*z*omega*(87*im - 3*im*lambda^2 + im*m^2*(9 + 8*lambda) + 159*omega + lambda*(-11*im + 15*omega)) + 
            z^2*(48*m^4 + lambda^3 + lambda^2*(29 + 16*im*omega) - m^2*(156 + 98*lambda + lambda^2 + 180*im*omega) + lambda*(54 + 236*im*omega) - 
            12*(-1 + 179*im*omega + 48*omega^2)))))/((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*
            (2*lambda + lambda^2 - 12*im*omega + 24*a^3*m*z^2*((-im)*z + 2*omega) + 4*a*m*(6*im*z^2 + 2*im*z*lambda + 3*omega) + 
            12*a^4*z^2*(z^2 + 2*im*z*omega - 2*omega^2) - 4*a^2*(6*z^3 + z^2*(-3 + 6*m^2 + 6*im*omega) + 2*im*z*(-3 + lambda)*omega + 3*omega^2)))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_P_out_inf(s::Int, m::Int, a, omega, lambda, z)
    if s == 0
        return begin
            (2*im*(omega + a^2*z^4*(im + a^2*omega) + z^2*(-im + 2*a^2*omega)))/
            ((1 + a^2*z^2)*(1 - 2*z + a^2*z^2))
        end
    elseif s == 1
        return begin
            ((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - (2*z*(1 - 2*z + a^2*z^2))/
            (1 + a^2*z^2)^2 + (2*a*z^2*(1 - 2*z + a^2*z^2)*((-im)*m*(1 + 3*a^2*z^2) + 
            a*z*(3 + 3*z + 2*lambda + 2*a^2*z^2*(1 + lambda))))/((1 + a^2*z^2)*
            (2 - 2*im*a*m*z - 2*im*a^3*m*z^3 + lambda + a^4*z^4*(1 + lambda) + 
            a^2*z^2*(3 + 2*z + 2*lambda))) + 2*im*omega))/(1 - 2*z + a^2*z^2)
        end
    elseif s == -1
        return begin
            (1/(1 - 2*z + a^2*z^2))*((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - 
            (2*z*(1 - 2*z + a^2*z^2))/(1 + a^2*z^2)^2 + 
            (2*a*z^2*(1 - 2*z + a^2*z^2)*(im*m*(1 + 3*a^2*z^2) + 
            a*z*(-1 + 3*z + 2*a^2*z^2*(-1 + lambda) + 2*lambda)))/
            ((1 + a^2*z^2)*(2*im*a*m*z + 2*im*a^3*m*z^3 + a^4*z^4*(-1 + lambda) + lambda + 
            a^2*z^2*(-1 + 2*z + 2*lambda))) + 2*im*omega))
        end
    elseif s == 2
        return begin
            (1/(1 - 2*z + a^2*z^2))*((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - 
            (2*z*(1 - 2*z + a^2*z^2))/(1 + a^2*z^2)^2 + 2*im*omega + 
            (8*a*z^2*(1 - 2*z + a^2*z^2)*(-6*a*m^2*z + im*m*(-4 + 9*a^2*z^2 - lambda + 
            z*(-6 - 12*im*a^2*omega)) + a*(6*a^2*z^3 + im*(1 + lambda)*omega + z^2*(-9 - 9*im*a^2*omega) + 
            z*(3 + 6*im*omega - 6*a^2*omega^2))))/((1 + a^2*z^2)*(24 + 10*lambda + lambda^2 - 
            4*im*a*m*(6*z^2 + 2*z*(4 + lambda) + 3*im*omega) + 12*im*omega + 24*a^3*m*z^2*(im*z + 2*omega) + 
            12*a^4*z^2*(z^2 - 2*im*z*omega - 2*omega^2) - 4*a^2*(6*z^3 + z^2*(-3 + 6*m^2 - 6*im*omega) - 
            2*im*z*(1 + lambda)*omega + 3*omega^2)))))
        end
    elseif s == -2
        return begin
            (1/(1 - 2*z + a^2*z^2))*((1 + a^2*z^2)*(-((2*(-1 + z)*z)/(1 + a^2*z^2)) - 
            (2*z*(1 + z*(-2 + a^2*z)))/(1 + a^2*z^2)^2 + 2*im*omega + 
            (8*a*z^2*(1 + z*(-2 + a^2*z))*(-3*a*z*(-1 + 2*m^2 + 3*z) + im*m*(6*z + lambda) - 
            im*a*(-3 + 6*z + lambda)*omega + 3*a^2*m*z*(-3*im*z + 4*omega) + 
            3*a^3*z*(2*z^2 + 3*im*z*omega - 2*omega^2)))/((1 + a^2*z^2)*(lambda*(2 + lambda) - 12*im*omega + 
            12*a^4*z^2*(z - (1 - im)*omega)*(z + (1 + im)*omega) + 24*a^3*m*z^2*((-im)*z + 2*omega) + 
            4*a*m*(2*im*z*(3*z + lambda) + 3*omega) + 4*a^2*(-3*z^2*(-1 + 2*m^2 + 2*z) - 
            2*im*z*(-3 + 3*z + lambda)*omega - 3*omega^2)))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_Q_out_inf(s::Int, m::Int, a, omega, lambda, z)
    if s == 0
        return begin
            -((z^2*(-2*z*(-1 + lambda) + lambda - 4*a^2*z^3*(1 + lambda) - 2*a^4*z^5*(3 + lambda) + 2*a*m*omega + 
            a^6*z^6*(1 - m^2 + lambda + 2*a*m*omega) + a^2*z^4*(8 + a^2*(2 - 2*m^2 + 3*lambda) + 
            6*a^3*m*omega) + z^2*(-4 + a^2*(1 - m^2 + 3*lambda) + 6*a^3*m*omega)))/
            ((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2))
        end
    elseif s == 1
        return begin
            (z^2*(2*a^12*z^12*(1 + lambda) + 4*im*a^10*z^11*(1 + lambda)*(3*im + a^2*omega) - 
            (2 + lambda)*(2 + lambda + 2*a*m*omega) - a^6*z^8*(-4 - 8*lambda + 2*im*a*m*(9 + 2*lambda) + 
            a^2*(31 + 5*lambda^2 - m^2*(1 + 4*lambda) + 8*lambda*(3 + 4*im*omega) + 16*im*omega) + 
            10*a^3*m*(-2 + lambda)*omega) - a^4*z^6*(52 + 32*lambda + 2*im*a*m*(31 + 6*lambda) + 
            a^2*(71 + 10*lambda^2 - 3*m^2*(3 + 2*lambda) + 8*lambda*(7 + 6*im*omega) + 48*im*omega) + 
            20*a^3*m*(-1 + lambda)*omega) - a*z^2*(4*im*m*(2 + lambda) - a*m^2*(4 + lambda) + 
            a*(2 + lambda)*(14 + 5*lambda + 8*im*omega) + 10*a^2*m*(1 + lambda)*omega) + 
            2*z*(im*a*m*(2 + lambda) + (2 + lambda)^2 + 2*im*a^2*(2 + m^2 + lambda)*omega) + 
            2*a^6*z^9*(8 - 6*im*a*m + a^2*(-3 + 3*m^2 - 4*lambda + lambda^2 - 8*im*omega) - 
            im*a^3*m*(-4 + m^2 - lambda - 10*im*omega) + 2*im*a^4*(6 + m^2 + 5*lambda)*omega) + 
            2*a^4*z^7*(-4 + 8*im*a*m + a^3*m*(14*im - 3*im*m^2 + 4*im*lambda - 30*omega) + 
            4*a^2*(10 + m^2 + 7*lambda + lambda^2 - 4*im*omega) + 4*im*a^4*(7 + 2*m^2 + 5*lambda)*omega) + 
            2*im*a^3*z^5*(14*m + a^2*(-3*m^3 + 6*m*(3 + lambda + 5*im*omega)) + 
            im*a*(-69 + m^2 - 48*lambda - 6*lambda^2 + 8*im*omega) + 4*a^3*(8 + 3*m^2 + 5*lambda)*omega) + 
            2*a^2*z^3*(36 - 2*m^2 + 26*lambda + 4*lambda^2 - im*a*(m^3 - 2*m*(5 + 2*lambda + 5*im*omega)) + 
            2*im*a^2*(9 + 4*m^2 + 5*lambda)*omega) - a^8*z^10*(-6*im*a*m - 8*(1 + 2*lambda) + 
            2*a^3*m*(-3 + lambda)*omega + a^2*(1 - m^2*(-1 + lambda) + lambda^2 + 8*im*lambda*omega)) - 
            a^2*z^4*(24*(2 + lambda) + 2*im*a*m*(23 + 6*lambda) + 20*a^3*m*lambda*omega + 
            a^2*(67 + 54*lambda + 10*lambda^2 - m^2*(11 + 4*lambda) + 48*im*omega + 32*im*lambda*omega))))/
            ((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*(2 - 2*im*a*m*z - 2*im*a^3*m*z^3 + lambda + 
            a^4*z^4*(1 + lambda) + a^2*z^2*(3 + 2*z + 2*lambda)))
        end
    elseif s == -1
        return begin
            (z^2*((-1 + 2*z)*lambda^2 + 2*a^12*z^11*(-1 + lambda)*(z + 2*im*omega) - 2*a^11*m*z^10*(1 + lambda)*omega + 
            2*a^9*m*z^8*(-3*im*z^2 + im*z*(-2 + m^2 - lambda - 2*im*omega) - (6 + 5*lambda)*omega) + 
            2*im*a^7*m*z^6*(3*m^2*z + 6*z^3 + z^2*(5 + 2*lambda) - 2*z*(3 + 2*lambda + 5*im*omega) + 
            2*im*(7 + 5*lambda)*omega) + 2*im*a^5*m*z^4*(-8*z^3 + z^2*(19 + 6*lambda) + 
            3*z*(m^2 - 2*(1 + lambda + 3*im*omega)) + 2*im*(8 + 5*lambda)*omega) + 
            2*im*a^3*m*z^2*(-14*z^3 + z^2*(11 + 6*lambda) + z*(m^2 - 2*(1 + 2*lambda + 7*im*omega)) + 
            im*(9 + 5*lambda)*omega) + 2*a*m*(2*im*z^2*lambda - (2 + lambda)*omega + z*((-im)*lambda + 4*omega)) - 
            a^10*z^9*(12*z^2*(-1 + lambda) + 4*im*(5 + m^2 - 5*lambda)*omega + 
            z*(5 - m^2*(-3 + lambda) - 4*lambda + lambda^2 - 12*im*omega + 8*im*lambda*omega)) + 
            a^2*z*(-24*z^3*lambda + m^2*(-4*z^2 + z*(2 + lambda) - 4*im*omega) + 
            4*z^2*(5*lambda + 2*lambda^2 - 6*im*omega) + 4*im*(-1 + lambda)*omega - 
            z*(4*lambda + 5*lambda^2 - 20*im*omega + 8*im*lambda*omega)) + 
            a^8*z^7*(8*z^3*(-3 + 2*lambda) + 2*z^2*(9 + 3*m^2 - 8*lambda + lambda^2 - 4*im*omega) - 
            8*im*(5 + 2*m^2 - 5*lambda)*omega + z*(-3 - 4*lambda - 5*lambda^2 + m^2*(-7 + 4*lambda) + 56*im*omega - 
            32*im*lambda*omega)) + a^4*z^3*(-8*z^4 - 4*z^3*(-3 + 8*lambda) - 
            2*z^2*(3 + m^2 - 24*lambda - 6*lambda^2 + 28*im*omega) - 4*im*(5 + 4*m^2 - 5*lambda)*omega + 
            z*(1 - 14*lambda - 10*lambda^2 + m^2*(3 + 4*lambda) + 72*im*omega - 32*im*lambda*omega)) + 
            a^6*z^5*(16*z^4 + 4*z^3*(-3 + 2*lambda) + 8*z^2*(m^2 + 3*lambda + lambda^2 - 5*im*omega) - 
            8*im*(5 + 3*m^2 - 5*lambda)*omega + z*(1 - 16*lambda - 10*lambda^2 + m^2*(-3 + 6*lambda) + 96*im*omega - 
            48*im*lambda*omega))))/((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*
            (2*im*a*m*z + 2*im*a^3*m*z^3 + a^4*z^4*(-1 + lambda) + lambda + a^2*z^2*(-1 + 2*z + 2*lambda)))
        end
    elseif s == 2
        return begin
            (z^2*(24*a^12*z^12 - lambda^3 + 48*im*a^10*z^11*(3*im + a^2*omega) - 2*lambda^2*(8 + a*m*omega) + 
            8*im*a^6*z^9*(24*im + 18*a*m + 3*im*a^2*(9 + m^2) + 
            a^3*m*(19 - 3*m^2 + 4*lambda - 30*im*omega) + a^4*(23 + 3*m^2 - 4*lambda + 24*im*omega)*omega) + 
            4*lambda*(-21 - (3*im + 8*a*m)*omega + 3*a^2*omega^2) + 
            12*a^8*z^10*(24 - 6*im*a*m + a^2*(5 + 3*m^2 - 8*im*omega) - 12*a^3*m*omega + 8*a^4*omega^2) - 
            2*im*a^3*z^5*(-24*m*(9 + 4*lambda) + a*(-4*im*m^2*(91 + 16*lambda) + 
            im*(924 + 47*lambda^2 + lambda^3 + lambda*(394 - 84*im*omega) + 732*im*omega)) + 
            4*a^2*m*(-143 - 3*lambda^2 + m^2*(23 + 2*lambda) + lambda*(-44 + 31*im*omega) + 313*im*omega) + 
            4*a^3*(5 + 3*lambda^2 - m^2*(41 + 8*lambda) + 5*lambda*(7 - 3*im*omega) - 219*im*omega)*omega + 
            24*a^4*m*(1 + lambda)*omega^2) + 24*(-6 - (3*im + 5*a*m)*omega - 
            a*(im*m + a*(-3 + m^2))*omega^2 + a^3*m*omega^3) - 
            4*im*a*z^3*(-12*m*(3 + lambda) + im*a*(408 + 194*lambda + 27*lambda^2 + lambda^3 - 16*m^2*(4 + lambda) + 
            240*im*omega) + 12*a^4*m*(1 + lambda)*omega^2 + 2*a^3*omega*(17 + 26*lambda + 3*lambda^2 - 
            m^2*(31 + 7*lambda) - 132*im*omega - 12*im*lambda*omega) + 
            2*a^2*m*(-98 - 35*lambda - 3*lambda^2 + m^2*(10 + lambda) + 164*im*omega + 20*im*lambda*omega)) + 
            4*a^6*z^8*(36 - 24*m^2 + 2*a^3*m*(-100 + 12*m^2 - 7*lambda + 6*im*omega)*omega + 
            6*a^4*(21 - 5*m^2 + lambda)*omega^2 + 12*a^5*m*omega^3 + 
            2*a*m*(-44*im + 9*im*m^2 - 11*im*lambda + 12*omega) - 
            a^2*(39 + 6*m^4 + 20*lambda + 2*lambda^2 + m^2*(-65 - 8*lambda + 30*im*omega) + 170*im*omega - 
            22*im*lambda*omega)) + 2*z*((9 + lambda)*(24 + 10*lambda + lambda^2 + 12*im*omega) - 
            8*im*a^3*m*(1 + lambda)*omega^2 + 4*a*m*(24*im + 10*im*lambda + im*lambda^2 + 27*omega + 3*lambda*omega) - 
            4*a^2*omega*(6*im + 7*im*lambda + im*lambda^2 - 2*im*m^2*(4 + lambda) + 27*omega + 3*lambda*omega)) + 
            a^3*z^6*(384*im*m + 16*a*(-33 + 12*m^2 - 10*lambda - lambda^2 - 36*im*omega) + 
            8*im*a^2*m*(-177 + 18*m^2 - 53*lambda - 2*lambda^2 + 96*im*omega) + 
            2*a^4*m*(-810 + 126*m^2 - 100*lambda - lambda^2 + 12*im*omega)*omega + 
            12*a^5*(83 - 31*m^2 + 7*lambda)*omega^2 + 168*a^6*m*omega^3 - 
            a^3*(756 + 48*m^4 + lambda^3 + lambda^2*(41 - 16*im*omega) - 
            m^2*(564 + 106*lambda + lambda^2 - 180*im*omega) + lambda*(334 - 364*im*omega) + 948*im*omega - 
            576*omega^2)) + a*z^4*(96*im*m + 8*im*a^2*m*(-202 + 9*m^2 - 65*lambda - 4*lambda^2 + 
            72*im*omega) - 8*a*(120 + 50*lambda + 5*lambda^2 + 72*im*omega) + 
            6*a^4*m*(-256 + 36*m^2 - 44*lambda - lambda^2 - 20*im*omega)*omega - 12*a^5*(-78 + 34*m^2 - 9*lambda)*
            omega^2 + 216*a^6*m*omega^3 + a^3*(-1068 - 24*m^4 - 3*lambda^3 + 
            2*m^2*(246 + 58*lambda + lambda^2) + lambda^2*(-74 + 32*im*omega) + lambda*(-512 + 388*im*omega) - 
            568*im*omega + 576*omega^2)) - 8*im*a^4*z^7*(12*im - 6*a*m*(-5 + lambda) + 
            2*a^5*m*(1 + lambda)*omega^2 + a^4*omega*(-23 + 20*lambda + lambda^2 - 3*m^2*(7 + lambda) - 138*im*omega - 
            6*im*lambda*omega) + a^3*m*(-88 - 23*lambda - lambda^2 + m^2*(16 + lambda) + 206*im*omega + 
            14*im*lambda*omega) + a^2*(-2*im*m^2*(31 + 4*lambda) + 3*im*(24 + 10*lambda + lambda^2 + 46*im*omega - 
            2*im*lambda*omega))) + z^2*(-12*(24 + 10*lambda + lambda^2 + 12*im*omega) - 
            8*im*a*m*(60 + 23*lambda + 2*lambda^2 - 18*im*omega) + 2*a^3*m*(-346 + 30*m^2 - 76*lambda - 
            3*lambda^2 - 60*im*omega)*omega - 60*a^4*(-7 + 3*m^2 - lambda)*omega^2 + 120*a^5*m*omega^3 + 
            a^2*(-3*lambda^3 + lambda^2*(-57 + 16*im*omega) + m^2*(152 + 42*lambda + lambda^2 + 60*im*omega) + 
            lambda*(-342 + 100*im*omega) + 12*(-54 - 23*im*omega + 12*omega^2)))))/
            ((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*(24 + 10*lambda + lambda^2 - 
            4*im*a*m*(6*z^2 + 2*z*(4 + lambda) + 3*im*omega) + 12*im*omega + 24*a^3*m*z^2*(im*z + 2*omega) + 
            12*a^4*z^2*(z^2 - 2*im*z*omega - 2*omega^2) - 4*a^2*(6*z^3 + z^2*(-3 + 6*m^2 - 6*im*omega) - 
            2*im*z*(1 + lambda)*omega + 3*omega^2)))
        end
    elseif s == -2
        return begin
            (z^2*(-((-1 + 2*z)*(-2 + 6*z - lambda)*(2*lambda + lambda^2 - 12*im*omega)) + 
            48*a^11*m*z^8*omega^2*(4*im*z + omega) + 24*a^12*z^9*(z^3 + 2*im*z^2*omega - 2*z*omega^2 - 
            4*im*omega^3) + 8*im*a^9*m*z^6*(9*z^4 + z^3*(-3 + 3*m^2 - 4*lambda + 18*im*omega) + 
            3*z^2*(-4*im*m^2 + 3*im*lambda - 18*omega)*omega + 2*z*(45 + lambda)*omega^2 - 21*im*omega^3) + 
            a^4*z*(96*z^6 + 16*z^5*(-9 + 12*m^2 - 2*lambda - lambda^2 + 36*im*omega) - 
            z^3*(12 + 24*m^4 + 3*lambda^3 + lambda^2*(38 + 32*im*omega) + 4*lambda*(16 + 33*im*omega) - 
            2*m^2*(30 + 50*lambda + lambda^2 + 288*im*omega) - 552*im*omega) + 
            2*z^4*(36 + 35*lambda^2 + lambda^3 - 4*m^2*(27 + 16*lambda) + lambda*(66 + 84*im*omega) - 636*im*omega) - 
            8*im*z^2*(3 - 3*lambda^2 + m^2*(51 + 7*lambda) + lambda*(-2 - 24*im*omega) + 36*im*omega)*omega + 
            4*z*(-39 - 45*m^2 + 31*lambda + 48*im*omega)*omega^2 - 96*im*omega^3) - 
            2*im*a^7*m*z^4*(72*z^5 + 4*z^3*(12 + lambda^2 - m^2*(12 + lambda) + 3*lambda*(5 + 6*im*omega) - 
            6*im*omega) + 4*z^4*(9*m^2 - 11*lambda + 36*im*omega) + 
            im*z^2*(6 + 126*m^2 - 124*lambda - lambda^2 - 588*im*omega)*omega - 24*z*(21 + lambda)*omega^2 + 
            108*im*omega^3) - a^6*z^3*(48*(-3 + 2*m^2)*z^5 + 192*z^6 + 
            8*z^4*(m^2*(30 + 8*lambda) - 3*(lambda^2 + lambda*(2 + 2*im*omega) - 18*im*omega)) + 
            z^3*(12 + 48*m^4 + lambda^3 + lambda^2*(29 + 16*im*omega) + lambda*(54 + 236*im*omega) - 
            m^2*(156 + 98*lambda + lambda^2 + 756*im*omega) - 420*im*omega) + 
            8*im*z^2*(3 - 3*lambda^2 + m^2*(81 + 8*lambda) + lambda*(-11 - 27*im*omega) + 57*im*omega)*omega + 
            12*z*(30 + 34*m^2 - 17*lambda - 48*im*omega)*omega^2 + 384*im*omega^3) - 
            4*a^10*z^7*(36*z^4 - 3*z^3*(5 + 3*m^2 - 20*im*omega) + 
            2*z*(27 + 15*m^2 - 5*lambda - 24*im*omega)*omega^2 + 96*im*omega^3 - 
            2*z^2*omega*(15*im - 15*im*m^2 + 4*im*lambda + 24*omega)) - 
            2*im*a^3*m*(192*z^6 + 24*z^5*(-7 + 4*lambda) + 4*z^4*(-6 + 9*m^2 - 33*lambda - 4*lambda^2) + 
            3*im*z^2*(-6 + 10*m^2 - 28*lambda - lambda^2 - 44*im*omega)*omega - 8*z*(9 + lambda)*omega^2 + 
            12*im*omega^3 - 4*z^3*(-6 - 11*lambda - 3*lambda^2 + m^2*(6 + lambda) - 32*im*lambda*omega)) - 
            2*im*a^5*m*z^2*(24*z^5*(-9 + lambda) + 4*z^4*(3 + 18*m^2 - 37*lambda - 2*lambda^2 - 24*im*omega) + 
            3*z^2*(36*im*m^2 - im*(52*lambda + lambda^2 + 172*im*omega))*omega - 24*z*(13 + lambda)*omega^2 + 
            60*im*omega^3 + 4*z^3*(15 + 20*lambda + 3*lambda^2 - m^2*(15 + 2*lambda) + 9*im*omega + 43*im*lambda*omega)) - 
            a^2*(8*z^4*(10*lambda + 5*lambda^2 - 72*im*omega) - 4*z^3*((26 - 16*m^2)*lambda + 15*lambda^2 + lambda^3 - 
            168*im*omega) + 4*(6 + 6*m^2 - 7*lambda)*omega^2 + 8*im*z*omega*(lambda - lambda^2 + 2*m^2*(6 + lambda) + 
            9*im*omega - 7*im*lambda*omega) + z^2*(3*lambda^3 + lambda^2*(21 - m^2 + 16*im*omega) + 
            lambda*(30 - 34*m^2 - 28*im*omega) - 12*im*(15 + 11*m^2 + 4*im*omega)*omega)) + 
            4*a^8*z^5*(72*z^5 - 6*z^4*(9 + m^2 - 12*im*omega) + z*(-99 - 93*m^2 + 37*lambda + 144*im*omega)*
            omega^2 - 144*im*omega^3 - 2*im*z^2*omega*(-9 - 12*lambda - lambda^2 + 3*m^2*(19 + lambda) + 54*im*omega - 
            10*im*lambda*omega) - z^3*(-9 + 6*m^4 + 4*lambda + 2*lambda^2 - m^2*(33 + 8*lambda + 78*im*omega) + 
            54*im*omega + 22*im*lambda*omega + 48*omega^2)) + 2*a*m*(-48*im*z^4 - 24*im*z^3*(-1 + lambda) - 
            (12 + 16*lambda + lambda^2 - 12*im*omega)*omega + 4*z^2*(7*im*lambda + 2*im*lambda^2 + 6*omega) + 
            4*z*((-im)*lambda^2 + 3*omega + lambda*(-2*im + 7*omega)))))/
            ((1 + a^2*z^2)^2*(1 - 2*z + a^2*z^2)^2*(2*lambda + lambda^2 - 12*im*omega + 
            24*a^3*m*z^2*((-im)*z + 2*omega) + 4*a*m*(6*im*z^2 + 2*im*z*lambda + 3*omega) + 
            12*a^4*z^2*(z^2 + 2*im*z*omega - 2*omega^2) - 4*a^2*(6*z^3 + z^2*(-3 + 6*m^2 + 6*im*omega) + 
            2*im*z*(-3 + lambda)*omega + 3*omega^2)))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_P_out_hor(s::Int, m::Int, a, omega, lambda, x)
    if s == 0
        return begin
            ((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            2*im*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)))/(2*sqrt(1 - a^2) + x)
        end
    elseif s == 1
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            (2*a*x*(2*sqrt(1 - a^2) + x)*(-3*im*a^2*m*(1 + sqrt(1 - a^2) + x) - 
            im*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(-2*im*a^3*m*(1 + sqrt(1 - a^2) + x) - 
            2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda)))) + 2*im*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)))
        end
    elseif s == -1
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            (2*a*x*(2*sqrt(1 - a^2) + x)*(3*im*a^2*m*(1 + sqrt(1 - a^2) + x) + 
            im*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(-1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(2*im*a^3*m*(1 + sqrt(1 - a^2) + x) + 
            2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            lambda + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))) + 
            2*im*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)))
        end
    elseif s == 2
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            2*im*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega) + 
            (8*a*x*(2*sqrt(1 - a^2) + x)*(im*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*(4 + lambda)) - 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*im + 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*im*omega) - im*(1 + sqrt(1 - a^2) + x)^2*
            (1 + lambda)*omega) + a^3*(-6 + 9*im*(1 + sqrt(1 - a^2) + x)*omega + 
            6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*((-(1 + sqrt(1 - a^2) + x)^4)*
            (24 + 10*lambda + lambda^2 + 12*im*omega) - 24*a^3*m*(1 + sqrt(1 - a^2) + x)*
            (im + 2*(1 + sqrt(1 - a^2) + x)*omega) + 4*a*m*(1 + sqrt(1 - a^2) + x)^2*
            (6*im + 2*im*(1 + sqrt(1 - a^2) + x)*(4 + lambda) - 3*(1 + sqrt(1 - a^2) + x)^2*omega) + 
            12*a^4*(-1 + 2*im*(1 + sqrt(1 - a^2) + x)*omega + 2*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 
            4*a^2*(1 + sqrt(1 - a^2) + x)*(6 + (1 + sqrt(1 - a^2) + x)*
            (-3 + 6*m^2 - 6*im*omega) - 2*im*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2)))))
        end
    elseif s == -2
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            2*im*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega) - 
            (8*a*x*(2*sqrt(1 - a^2) + x)*((-im)*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*lambda) + 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*im - 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*im*omega) + im*(1 + sqrt(1 - a^2) + x)^2*
            (-3 + lambda)*omega) + a^3*(-6 - 9*im*(1 + sqrt(1 - a^2) + x)*omega + 
            6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*((1 + sqrt(1 - a^2) + x)^4*
            (2*lambda + lambda^2 - 12*im*omega) + 24*a^3*m*(1 + sqrt(1 - a^2) + x)*
            (-im + 2*(1 + sqrt(1 - a^2) + x)*omega) + 4*a*m*(1 + sqrt(1 - a^2) + x)^2*
            (6*im + 2*im*(1 + sqrt(1 - a^2) + x)*lambda + 3*(1 + sqrt(1 - a^2) + x)^2*omega) + 
            a^4*(12 + 24*im*(1 + sqrt(1 - a^2) + x)*omega - 24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 
            4*a^2*(1 + sqrt(1 - a^2) + x)*(6 + (1 + sqrt(1 - a^2) + x)*
            (-3 + 6*m^2 + 6*im*omega) + 2*im*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2)))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_Q_out_hor(s::Int, m::Int, a, omega, lambda, x)
    if s == 0
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            (x^2*(1 + sqrt(1 - a^2) + x)^2*(2*sqrt(1 - a^2) + x)^2 + 
            x*(2*sqrt(1 - a^2) + x)*(a^4 - 4*a^2*(1 + sqrt(1 - a^2) + x) - 
            (-3 + sqrt(1 - a^2) + x)*(1 + sqrt(1 - a^2) + x)^3) + 
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*(x*(2*sqrt(1 - a^2) + x)*lambda - 
            (a*m - (a^2 + (1 + sqrt(1 - a^2) + x)^2)*omega)^2))))
        end
    elseif s == 1
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            ((2*im*a*x*(2*sqrt(1 - a^2) + x)*(-3*im*a^2*m*(1 + sqrt(1 - a^2) + x) - 
            im*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda)))*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + 
            omega))/((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-2*im*a^3*m*(1 + sqrt(1 - a^2) + x) - 2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*(2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda)))) - 
            (-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            (2*x*(2*sqrt(1 - a^2) + x)*(a^4 - a^2*(1 + sqrt(1 - a^2) + x) - 
            (-2 + sqrt(1 - a^2) + x)*(1 + sqrt(1 - a^2) + x)^3) + 
            ((1 + sqrt(1 - a^2) + x)^2*(-1 + 2*sqrt(1 - a^2) + 2*x) + 
            a^2*(1 + 2*sqrt(1 - a^2) + 2*x))^2 + (2*a*x*(2*sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*((1 + sqrt(1 - a^2) + x)^2*
            (-1 + 2*sqrt(1 - a^2) + 2*x) + a^2*(1 + 2*sqrt(1 - a^2) + 2*x))*
            (-3*im*a^2*m*(1 + sqrt(1 - a^2) + x) - im*m*(1 + sqrt(1 - a^2) + x)^3 + 
            2*a^3*(1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*(3 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*(-2*im*a^3*m*(1 + sqrt(1 - a^2) + x) - 
            2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda)))) - ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (2*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*(-lambda - im*(1 + sqrt(1 - a^2) + x)*omega) + 
            a^8*(1 + lambda)*(2 + (1 + sqrt(1 - a^2) + x)^2*omega^2) - 
            2*im*a^5*m*(1 + sqrt(1 - a^2) + x)^2*(-6 + (1 + sqrt(1 - a^2) + x)*
            (4 + m^2 - lambda - 4*im*omega) - 3*im*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2) - 2*im*a*m*(1 + sqrt(1 - a^2) + x)^5*
            (3 + (1 + sqrt(1 - a^2) + x)*(-5 + 2*lambda) + (1 + sqrt(1 - a^2) + x)^2*
            (2 - lambda + 2*im*omega) - im*(1 + sqrt(1 - a^2) + x)^3*(3 + lambda)*omega + 
            (1 + sqrt(1 - a^2) + x)^4*omega^2) - 2*im*a^3*m*(1 + sqrt(1 - a^2) + x)^3*
            (13 + (1 + sqrt(1 - a^2) + x)*(-17 + 2*lambda) + (1 + sqrt(1 - a^2) + x)^2*
            (6 + m^2 - 2*lambda - 2*im*omega) - 3*im*(1 + sqrt(1 - a^2) + x)^3*(2 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^4*omega^2) + (1 + sqrt(1 - a^2) + x)^6*
            (-3*(2 + lambda) - (1 + sqrt(1 - a^2) + x)^2*lambda*(2 + lambda) + 
            2*(1 + sqrt(1 - a^2) + x)*(2 + 3*lambda + lambda^2) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda)*omega^2) + a^6*(1 + sqrt(1 - a^2) + x)*(-16*(1 + lambda) + 
            (1 + sqrt(1 - a^2) + x)*(5 + m^2*(-1 + lambda) + 6*lambda - lambda^2 + 2*im*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(5 + 4*lambda)*omega^2 + 2*(1 + sqrt(1 - a^2) + x)^2*omega*
            (im + 2*im*m^2 + omega)) + a^4*(1 + sqrt(1 - a^2) + x)^2*
            (11 + 25*lambda + 2*(1 + sqrt(1 - a^2) + x)*(-5 + 3*m^2 - 13*lambda + lambda^2 - 2*im*omega) + 
            (1 + sqrt(1 - a^2) + x)^2*(5 + 4*lambda - 3*lambda^2 + m^2*(3 + 2*lambda) - 8*im*omega) + 
            3*(1 + sqrt(1 - a^2) + x)^4*(3 + 2*lambda)*omega^2 + 4*(1 + sqrt(1 - a^2) + x)^3*omega*
            (im + 2*im*m^2 + omega)) + a^2*(1 + sqrt(1 - a^2) + x)^3*
            (30 + 5*(1 + sqrt(1 - a^2) + x)*(-7 + 2*lambda) - 4*(1 + sqrt(1 - a^2) + x)^2*
            (-2 + m^2 + lambda - lambda^2 - 3*im*omega) + (1 + sqrt(1 - a^2) + x)^3*
            (2 - 2*lambda - 3*lambda^2 + m^2*(4 + lambda) - 10*im*omega) + (1 + sqrt(1 - a^2) + x)^5*
            (7 + 4*lambda)*omega^2 + 2*(1 + sqrt(1 - a^2) + x)^4*omega*(im + 2*im*m^2 + omega))))/
            ((1 + sqrt(1 - a^2) + x)^2*(-2*im*a^3*m*(1 + sqrt(1 - a^2) + x) - 
            2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda)))))))
        end
    elseif s == -1
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            ((2*im*a*x*(2*sqrt(1 - a^2) + x)*(3*im*a^2*m*(1 + sqrt(1 - a^2) + x) + 
            im*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(-1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + 
            omega))/((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (2*im*a^3*m*(1 + sqrt(1 - a^2) + x) + 2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*lambda + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))) - 
            (-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            ((a^2 - (1 + sqrt(1 - a^2) + x)^2)^2 + 2*x*(1 + sqrt(1 - a^2) + x)*
            (2*sqrt(1 - a^2) + x)*(-3*a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            (2*a*x*(2*sqrt(1 - a^2) + x)*(a^2 - (1 + sqrt(1 - a^2) + x)^2)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(3*im*a^2*m*(1 + sqrt(1 - a^2) + x) + 
            im*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(-1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*
            (2*im*a^3*m*(1 + sqrt(1 - a^2) + x) + 2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*lambda + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))) - 
            ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*(2*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*
            (2 - lambda + im*(1 + sqrt(1 - a^2) + x)*omega) + a^8*(-1 + lambda)*
            (2 + (1 + sqrt(1 - a^2) + x)^2*omega^2) + 2*im*a^5*m*(1 + sqrt(1 - a^2) + x)^3*
            (-2 + m^2 - lambda + 4*im*omega - 3*im*(1 + sqrt(1 - a^2) + x)*omega + 
            3*im*(1 + sqrt(1 - a^2) + x)*lambda*omega + 3*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 
            (1 + sqrt(1 - a^2) + x)^6*lambda*(-3 - (1 + sqrt(1 - a^2) + x)^2*lambda + 
            2*(1 + sqrt(1 - a^2) + x)*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            2*im*a*m*(1 + sqrt(1 - a^2) + x)^5*(sqrt(1 - a^2) + x + 
            2*(1 + sqrt(1 - a^2) + x)*lambda - (1 + sqrt(1 - a^2) + x)^2*(lambda + 2*im*omega) + 
            im*(1 + sqrt(1 - a^2) + x)^3*(1 + lambda)*omega + (1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            2*im*a^3*m*(1 + sqrt(1 - a^2) + x)^3*(1 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda) + 
            (1 + sqrt(1 - a^2) + x)^2*(-2 + m^2 - 2*lambda + 2*im*omega) + 
            3*im*(1 + sqrt(1 - a^2) + x)^3*lambda*omega + 3*(1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            a^6*(1 + sqrt(1 - a^2) + x)*(8 - 8*lambda + (1 + sqrt(1 - a^2) + x)*
            (-1 + m^2*(-3 + lambda) - lambda^2 - 2*im*omega) + (1 + sqrt(1 - a^2) + x)^3*(-3 + 4*lambda)*
            omega^2 + 2*(1 + sqrt(1 - a^2) + x)^2*omega*(-im - 2*im*m^2 + omega)) + 
            a^4*(1 + sqrt(1 - a^2) + x)^2*(-11 + 9*lambda + 2*(1 + sqrt(1 - a^2) + x)*
            (1 + 3*m^2 + lambda + lambda^2 + 2*im*omega) + (1 + sqrt(1 - a^2) + x)^2*
            (1 - 6*lambda - 3*lambda^2 + m^2*(-1 + 2*lambda) + 8*im*omega) + 3*(1 + sqrt(1 - a^2) + x)^4*
            (-1 + 2*lambda)*omega^2 + 4*(1 + sqrt(1 - a^2) + x)^3*omega*(-im - 2*im*m^2 + omega)) + 
            a^2*(1 + sqrt(1 - a^2) + x)^3*(6 - 3*(1 + sqrt(1 - a^2) + x)*(1 + 2*lambda) - 
            4*(1 + sqrt(1 - a^2) + x)^2*(m^2 - 3*lambda - lambda^2 + 3*im*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(-4*lambda - 3*lambda^2 + m^2*(2 + lambda) + 10*im*omega) + 
            (1 + sqrt(1 - a^2) + x)^5*(-1 + 4*lambda)*omega^2 + 2*(1 + sqrt(1 - a^2) + x)^4*omega*
            (-im - 2*im*m^2 + omega))))/((1 + sqrt(1 - a^2) + x)^2*
            (2*im*a^3*m*(1 + sqrt(1 - a^2) + x) + 2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*lambda + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))))))
        end
    elseif s == 2
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - 
            (8*im*a*x*(2*sqrt(1 - a^2) + x)*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)*
            (im*m*(1 + sqrt(1 - a^2) + x)^2*(6 + (1 + sqrt(1 - a^2) + x)*(4 + lambda)) - 
            3*a^2*m*(1 + sqrt(1 - a^2) + x)*(3*im + 4*(1 + sqrt(1 - a^2) + x)*omega) + 
            a*(1 + sqrt(1 - a^2) + x)*(9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*im*omega) - 
            im*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega) + 
            a^3*(-6 + 9*im*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            ((1 + sqrt(1 - a^2) + x)^4*(24 + 10*lambda + lambda^2 + 12*im*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(im + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(-6*im - 2*im*(1 + sqrt(1 - a^2) + x)*(4 + lambda) + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) - 12*a^4*(-1 + 2*im*(1 + sqrt(1 - a^2) + x)*omega + 
            2*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*im*omega) - 
            2*im*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))) - (((1 + sqrt(1 - a^2) + x)^2*(-1 + 3*sqrt(1 - a^2) + 3*x) + 
            a^2*(1 + 3*sqrt(1 - a^2) + 3*x))^2 + x*(2*sqrt(1 - a^2) + x)*
            (3*a^4 + (1 + sqrt(1 - a^2) + x)^3*(8 - 3*(1 + sqrt(1 - a^2) + x))) + 
            (8*a*x*(2*sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            ((1 + sqrt(1 - a^2) + x)^2*(-1 + 3*sqrt(1 - a^2) + 3*x) + 
            a^2*(1 + 3*sqrt(1 - a^2) + 3*x))*(im*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*(4 + lambda)) - 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*im + 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*im*omega) - 
            im*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega) + 
            a^3*(-6 + 9*im*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*((-(1 + sqrt(1 - a^2) + x)^4)*(24 + 10*lambda + lambda^2 + 
            12*im*omega) - 24*a^3*m*(1 + sqrt(1 - a^2) + x)*(im + 2*(1 + sqrt(1 - a^2) + x)*
            omega) + 4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*im + 2*im*(1 + sqrt(1 - a^2) + x)*
            (4 + lambda) - 3*(1 + sqrt(1 - a^2) + x)^2*omega) + 
            12*a^4*(-1 + 2*im*(1 + sqrt(1 - a^2) + x)*omega + 2*(1 + sqrt(1 - a^2) + x)^2*
            omega^2) + 4*a^2*(1 + sqrt(1 - a^2) + x)*(6 + (1 + sqrt(1 - a^2) + x)*
            (-3 + 6*m^2 - 6*im*omega) - 2*im*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2))) - ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (24*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*(3 + 3*im*(1 + sqrt(1 - a^2) + x)*omega - 
            4*(1 + sqrt(1 - a^2) + x)^2*omega^2) - (1 + sqrt(1 - a^2) + x)^6*
            (24 + 10*lambda + lambda^2 + 12*im*omega)*(-12 - (1 + sqrt(1 - a^2) + x)^2*lambda + 
            2*(1 + sqrt(1 - a^2) + x)*(4 + lambda) + (1 + sqrt(1 - a^2) + x)^4*omega^2) - 
            4*a^5*m*(1 + sqrt(1 - a^2) + x)^2*(-54*im - 2*im*(1 + sqrt(1 - a^2) + x)*
            (-55 + 3*m^2 - 4*lambda - 18*im*omega) + 2*(1 + sqrt(1 - a^2) + x)^2*
            (4 + 12*m^2 - 8*lambda + 27*im*omega)*omega - 6*im*(1 + sqrt(1 - a^2) + x)^3*(8 + lambda)*
            omega^2 + 45*(1 + sqrt(1 - a^2) + x)^4*omega^3) + 
            12*a^8*(-2 - 3*(1 + sqrt(1 - a^2) + x)^2*omega^2 - 2*im*(1 + sqrt(1 - a^2) + x)^3*
            omega^3 + 2*(1 + sqrt(1 - a^2) + x)^4*omega^4) + 4*a^6*(1 + sqrt(1 - a^2) + x)*
            (60 - 3*(1 + sqrt(1 - a^2) + x)*(16 + 3*m^2 + 18*im*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(1 + 36*m^2 - 8*lambda + 18*im*omega)*omega^2 - 
            2*im*(1 + sqrt(1 - a^2) + x)^4*(7 + lambda)*omega^3 + 15*(1 + sqrt(1 - a^2) + x)^5*
            omega^4 - 2*(1 + sqrt(1 - a^2) + x)^2*omega*(-40*im + 9*im*m^2 - 4*im*lambda + 9*omega)) + 
            2*a*m*(1 + sqrt(1 - a^2) + x)^4*(144*im + 8*im*(1 + sqrt(1 - a^2) + x)*
            (-38 + lambda) + (1 + sqrt(1 - a^2) + x)^4*(40 + 20*lambda + lambda^2 + 24*im*omega)*omega + 
            4*im*(1 + sqrt(1 - a^2) + x)^5*(4 + lambda)*omega^2 - 6*(1 + sqrt(1 - a^2) + x)^6*
            omega^3 + 4*(1 + sqrt(1 - a^2) + x)^2*(16*im + 3*im*lambda + 2*im*lambda^2 + 6*omega) - 
            4*im*(1 + sqrt(1 - a^2) + x)^3*(-12 + lambda + lambda^2 - 14*im*omega - 5*im*lambda*omega)) - 
            2*a^3*m*(1 + sqrt(1 - a^2) + x)^3*(312*im + (1 + sqrt(1 - a^2) + x)^3*
            (84 + 30*m^2 - 52*lambda - lambda^2 + 36*im*omega)*omega - 4*im*(1 + sqrt(1 - a^2) + x)^4*
            (19 + 4*lambda)*omega^2 + 48*(1 + sqrt(1 - a^2) + x)^5*omega^3 + 
            12*(1 + sqrt(1 - a^2) + x)*(-63*im + 3*im*m^2 - 3*im*lambda + 32*omega) + 
            4*(1 + sqrt(1 - a^2) + x)^2*(71*im + im*lambda^2 - im*m^2*(10 + lambda) - 104*omega + 
            lambda*(9*im + 16*omega))) + a^2*(1 + sqrt(1 - a^2) + x)^3*
            (576 + 96*(1 + sqrt(1 - a^2) + x)*(-12 + 4*m^2 - 3*im*omega) + 
            8*(1 + sqrt(1 - a^2) + x)^2*(18 - 3*lambda^2 + m^2*(-58 + 8*lambda) + 
            lambda*(-30 - 2*im*omega) + 46*im*omega) - 8*im*(1 + sqrt(1 - a^2) + x)^4*
            (-lambda^2 + 2*m^2*(7 + lambda) + lambda*(2 + 5*im*omega) + 11*im*omega)*omega + 
            2*(1 + sqrt(1 - a^2) + x)^5*(-34 + 24*m^2 - 20*lambda - lambda^2 - 24*im*omega)*omega^2 - 
            8*im*(1 + sqrt(1 - a^2) + x)^6*(1 + lambda)*omega^3 + 12*(1 + sqrt(1 - a^2) + x)^7*
            omega^4 + (1 + sqrt(1 - a^2) + x)^3*(lambda^3 + 36*lambda*(4 + im*omega) + 
            2*lambda^2*(11 - 8*im*omega) - m^2*(-136 + 42*lambda + lambda^2 - 36*im*omega) - 
            8*(-18 + 19*im*omega + 6*omega^2))) + a^4*(1 + sqrt(1 - a^2) + x)^2*
            (-672 + (1 + sqrt(1 - a^2) + x)*(960 - 72*m^2 + 624*im*omega) + 
            (1 + sqrt(1 - a^2) + x)^4*(44 + 180*m^2 - 62*lambda - lambda^2 + 36*im*omega)*omega^2 - 
            8*im*(1 + sqrt(1 - a^2) + x)^5*(5 + 2*lambda)*omega^3 + 48*(1 + sqrt(1 - a^2) + x)^6*
            omega^4 + 8*(1 + sqrt(1 - a^2) + x)^3*omega*(56*im + im*lambda^2 - 3*im*m^2*(9 + lambda) - 
            46*omega + lambda*(6*im + 8*omega)) + 4*(1 + sqrt(1 - a^2) + x)^2*
            (6*m^4 + m^2*(7 - 8*lambda + 54*im*omega) + 2*(-15 + 10*lambda + lambda^2 - 144*im*omega - 9*im*
            lambda*omega + 48*omega^2)))))/((1 + sqrt(1 - a^2) + x)^2*
            ((-(1 + sqrt(1 - a^2) + x)^4)*(24 + 10*lambda + lambda^2 + 12*im*omega) - 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(im + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*im + 2*im*(1 + sqrt(1 - a^2) + x)*(4 + lambda) - 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + 12*a^4*(-1 + 2*im*(1 + sqrt(1 - a^2) + x)*
            omega + 2*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*im*omega) - 
            2*im*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))))/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4))
        end
    elseif s == -2
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - 
            (8*im*a*x*(2*sqrt(1 - a^2) + x)*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)*
            ((-im)*m*(1 + sqrt(1 - a^2) + x)^2*(6 + (1 + sqrt(1 - a^2) + x)*lambda) + 
            3*a^2*m*(1 + sqrt(1 - a^2) + x)*(3*im - 4*(1 + sqrt(1 - a^2) + x)*omega) + 
            a*(1 + sqrt(1 - a^2) + x)*(9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*im*omega) + 
            im*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega) + 
            a^3*(-6 - 9*im*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            ((1 + sqrt(1 - a^2) + x)^4*(2*lambda + lambda^2 - 12*im*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(-im + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*im + 2*im*(1 + sqrt(1 - a^2) + x)*lambda + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + a^4*(12 + 24*im*(1 + sqrt(1 - a^2) + x)*omega - 
            24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*im*omega) + 
            2*im*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))) - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            ((a^2*(-1 + sqrt(1 - a^2) + x) + (1 + sqrt(1 - a^2) + x)^3)^2 + 
            x*(2*sqrt(1 - a^2) + x)*(-a^4 - 8*a^2*(1 + sqrt(1 - a^2) + x) + 
            (1 + sqrt(1 - a^2) + x)^4) + (8*a*x*(2*sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(a^2*(-1 + sqrt(1 - a^2) + x) + 
            (1 + sqrt(1 - a^2) + x)^3)*((-im)*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*lambda) + 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*im - 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*im*omega) + 
            im*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega) + 
            a^3*(-6 - 9*im*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*((1 + sqrt(1 - a^2) + x)^4*(2*lambda + lambda^2 - 12*im*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(-im + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*im + 2*im*(1 + sqrt(1 - a^2) + x)*lambda + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + a^4*(12 + 24*im*(1 + sqrt(1 - a^2) + x)*omega - 
            24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*im*omega) + 
            2*im*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))) - ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-24*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*(-3 + 3*im*(1 + sqrt(1 - a^2) + x)*omega + 
            4*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 2*a^3*m*(1 + sqrt(1 - a^2) + x)^3*
            (72*im - 4*im*(1 + sqrt(1 - a^2) + x)*(15 + 9*m^2 - 13*lambda) + 
            4*(1 + sqrt(1 - a^2) + x)^2*(im*m^2*(6 + lambda) - im*(3 + lambda^2 + lambda*
            (9 + 16*im*omega) + 24*im*omega)) + (1 + sqrt(1 - a^2) + x)^3*
            (-36 + 30*m^2 - 44*lambda - lambda^2 - 36*im*omega)*omega + 4*im*(1 + sqrt(1 - a^2) + x)^4*
            (3 + 4*lambda)*omega^2 + 48*(1 + sqrt(1 - a^2) + x)^5*omega^3) + 
            12*a^8*(-2 - 3*(1 + sqrt(1 - a^2) + x)^2*omega^2 + 2*im*(1 + sqrt(1 - a^2) + x)^3*
            omega^3 + 2*(1 + sqrt(1 - a^2) + x)^4*omega^4) + 4*a^6*(1 + sqrt(1 - a^2) + x)*
            (12 - 3*(1 + sqrt(1 - a^2) + x)*(-4 + 3*m^2 + 6*im*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(-39 + 36*m^2 - 8*lambda - 18*im*omega)*omega^2 + 
            2*im*(1 + sqrt(1 - a^2) + x)^4*(3 + lambda)*omega^3 + 15*(1 + sqrt(1 - a^2) + x)^5*
            omega^4 + 2*(1 + sqrt(1 - a^2) + x)^2*omega*(24*im + 9*im*m^2 - 4*im*lambda + 15*omega)) - 
            4*a^5*m*(1 + sqrt(1 - a^2) + x)^2*(-18*im + 2*(1 + sqrt(1 - a^2) + x)^2*
            (-36 + 12*m^2 - 8*lambda - 27*im*omega)*omega + 6*im*(1 + sqrt(1 - a^2) + x)^3*(4 + lambda)*
            omega^2 + 45*(1 + sqrt(1 - a^2) + x)^4*omega^3 + 2*(1 + sqrt(1 - a^2) + x)*
            (9*im + 3*im*m^2 - 4*im*lambda + 30*omega)) - 2*a*m*(1 + sqrt(1 - a^2) + x)^4*
            (-48*im - 24*im*(1 + sqrt(1 - a^2) + x)*(-2 + lambda) + 
            4*im*(1 + sqrt(1 - a^2) + x)^2*(7*lambda + 2*lambda^2 + 6*im*omega) - 
            (1 + sqrt(1 - a^2) + x)^4*(12*lambda + lambda^2 - 24*im*omega)*omega + 
            4*im*(1 + sqrt(1 - a^2) + x)^5*lambda*omega^2 + 6*(1 + sqrt(1 - a^2) + x)^6*omega^3 - 
            4*im*(1 + sqrt(1 - a^2) + x)^3*(lambda + lambda^2 + 6*im*omega + 5*im*lambda*omega)) + 
            a^4*(1 + sqrt(1 - a^2) + x)^3*(24*m^4*(1 + sqrt(1 - a^2) + x) + 
            48*(-4 + 3*im*omega) + 8*(1 + sqrt(1 - a^2) + x)*(9 + lambda^2 + lambda*(2 + 13*im*omega) - 
            72*im*omega) - (1 + sqrt(1 - a^2) + x)^3*(60 + 54*lambda + lambda^2 + 36*im*omega)*omega^2 + 
            8*im*(1 + sqrt(1 - a^2) + x)^4*(-3 + 2*lambda)*omega^3 + 48*(1 + sqrt(1 - a^2) + x)^5*
            omega^4 + 8*(1 + sqrt(1 - a^2) + x)^2*omega*(24*im - 6*im*lambda - im*lambda^2 + 18*omega + 
            8*lambda*omega) + 4*m^2*(30 - (1 + sqrt(1 - a^2) + x)*(33 + 8*lambda + 54*im*omega) + 
            6*im*(1 + sqrt(1 - a^2) + x)^2*(5 + lambda)*omega + 45*(1 + sqrt(1 - a^2) + x)^3*
            omega^2)) + (1 + sqrt(1 - a^2) + x)^6*((-1 + sqrt(1 - a^2) + x)*
            (1 + sqrt(1 - a^2) + x)*lambda^3 + lambda^2*(12 - 12*(1 + sqrt(1 - a^2) + x) + 
            2*(1 + sqrt(1 - a^2) + x)^2 - (1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            12*im*omega*(-12 + 8*(1 + sqrt(1 - a^2) + x) + (1 + sqrt(1 - a^2) + x)^4*omega^2) - 
            2*lambda*(-12 + 4*(1 + sqrt(1 - a^2) + x)*(2 - 3*im*omega) + 
            6*im*(1 + sqrt(1 - a^2) + x)^2*omega + (1 + sqrt(1 - a^2) + x)^4*omega^2)) + 
            a^2*(1 + sqrt(1 - a^2) + x)^4*(96 + 8*(1 + sqrt(1 - a^2) + x)*
            (m^2*(6 + 8*lambda) - 3*(2 + lambda^2 + lambda*(2 + 2*im*omega) - 22*im*omega)) - 96*im*omega + 
            2*(1 + sqrt(1 - a^2) + x)^4*(6 + 24*m^2 - 12*lambda - lambda^2 + 24*im*omega)*omega^2 + 
            8*im*(1 + sqrt(1 - a^2) + x)^5*(-3 + lambda)*omega^3 + 12*(1 + sqrt(1 - a^2) + x)^6*
            omega^4 + (1 + sqrt(1 - a^2) + x)^2*(lambda^3 + lambda*(24 - 34*m^2 - 4*im*omega) + 
            lambda^2*(14 - m^2 + 16*im*omega) - 12*im*(22 + 3*m^2 - 4*im*omega)*omega) + 
            8*(1 + sqrt(1 - a^2) + x)^3*omega*((-im)*lambda^2 + 2*im*m^2*(3 + lambda) + 3*omega + 
            lambda*(2*im + 5*omega)))))/((1 + sqrt(1 - a^2) + x)^2*
            ((-(1 + sqrt(1 - a^2) + x)^4)*(2*lambda + lambda^2 - 12*im*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(im - 2*(1 + sqrt(1 - a^2) + x)*omega) - 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*im + 2*im*(1 + sqrt(1 - a^2) + x)*lambda + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + 12*a^4*(-1 - 2*im*(1 + sqrt(1 - a^2) + x)*
            omega + 2*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*im*omega) + 
            2*im*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_P_in_hor(s::Int, m::Int, a, omega, lambda, x)
    if s == 0
        return begin
            ((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            im*((a*m)/(1 + sqrt(1 - a^2)) - 2*omega)))/(2*sqrt(1 - a^2) + x)
        end
    elseif s == 1
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            (2*a*x*(2*sqrt(1 - a^2) + x)*(-3*im*a^2*m*(1 + sqrt(1 - a^2) + x) - 
            im*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(-2*im*a^3*m*(1 + sqrt(1 - a^2) + x) - 
            2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda)))) + im*((a*m)/(1 + sqrt(1 - a^2)) - 2*omega)))
        end
    elseif s == -1
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            (2*a*x*(2*sqrt(1 - a^2) + x)*(3*im*a^2*m*(1 + sqrt(1 - a^2) + x) + 
            im*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(-1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(2*im*a^3*m*(1 + sqrt(1 - a^2) + x) + 
            2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            lambda + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))) + 
            im*((a*m)/(1 + sqrt(1 - a^2)) - 2*omega)))
        end
    elseif s == 2
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            im*((a*m)/(1 + sqrt(1 - a^2)) - 2*omega) + 
            (8*a*x*(2*sqrt(1 - a^2) + x)*(im*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*(4 + lambda)) - 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*im + 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*im*omega) - im*(1 + sqrt(1 - a^2) + x)^2*
            (1 + lambda)*omega) + a^3*(-6 + 9*im*(1 + sqrt(1 - a^2) + x)*omega + 
            6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*((-(1 + sqrt(1 - a^2) + x)^4)*
            (24 + 10*lambda + lambda^2 + 12*im*omega) - 24*a^3*m*(1 + sqrt(1 - a^2) + x)*
            (im + 2*(1 + sqrt(1 - a^2) + x)*omega) + 4*a*m*(1 + sqrt(1 - a^2) + x)^2*
            (6*im + 2*im*(1 + sqrt(1 - a^2) + x)*(4 + lambda) - 3*(1 + sqrt(1 - a^2) + x)^2*omega) + 
            12*a^4*(-1 + 2*im*(1 + sqrt(1 - a^2) + x)*omega + 2*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 
            4*a^2*(1 + sqrt(1 - a^2) + x)*(6 + (1 + sqrt(1 - a^2) + x)*
            (-3 + 6*m^2 - 6*im*omega) - 2*im*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2)))))
        end
    elseif s == -2
        return begin
            (1/(2*sqrt(1 - a^2) + x))*((a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-((2*x*(1 + sqrt(1 - a^2) + x)*(2*sqrt(1 - a^2) + x))/
            (2*(1 + sqrt(1 - a^2)) + 2*(1 + sqrt(1 - a^2))*x + x^2)^2) + 
            (2*(sqrt(1 - a^2) + x))/(a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            im*((a*m)/(1 + sqrt(1 - a^2)) - 2*omega) - 
            (8*a*x*(2*sqrt(1 - a^2) + x)*((-im)*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*lambda) + 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*im - 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*im*omega) + im*(1 + sqrt(1 - a^2) + x)^2*
            (-3 + lambda)*omega) + a^3*(-6 - 9*im*(1 + sqrt(1 - a^2) + x)*omega + 
            6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/((1 + sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*((1 + sqrt(1 - a^2) + x)^4*
            (2*lambda + lambda^2 - 12*im*omega) + 24*a^3*m*(1 + sqrt(1 - a^2) + x)*
            (-im + 2*(1 + sqrt(1 - a^2) + x)*omega) + 4*a*m*(1 + sqrt(1 - a^2) + x)^2*
            (6*im + 2*im*(1 + sqrt(1 - a^2) + x)*lambda + 3*(1 + sqrt(1 - a^2) + x)^2*omega) + 
            a^4*(12 + 24*im*(1 + sqrt(1 - a^2) + x)*omega - 24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 
            4*a^2*(1 + sqrt(1 - a^2) + x)*(6 + (1 + sqrt(1 - a^2) + x)*
            (-3 + 6*m^2 + 6*im*omega) + 2*im*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2)))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end


function _gsn_Q_in_hor(s::Int, m::Int, a, omega, lambda, x)
    if s == 0
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            (x^2*(1 + sqrt(1 - a^2) + x)^2*(2*sqrt(1 - a^2) + x)^2 + 
            x*(2*sqrt(1 - a^2) + x)*(a^4 - 4*a^2*(1 + sqrt(1 - a^2) + x) - 
            (-3 + sqrt(1 - a^2) + x)*(1 + sqrt(1 - a^2) + x)^3) + 
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*(x*(2*sqrt(1 - a^2) + x)*lambda - 
            (a*m - (a^2 + (1 + sqrt(1 - a^2) + x)^2)*omega)^2))))
        end
    elseif s == 1
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-((2*im*a*x*(2*sqrt(1 - a^2) + x)*(-3*im*a^2*m*(1 + sqrt(1 - a^2) + x) - 
            im*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda)))*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + 
            omega))/((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (-2*im*a^3*m*(1 + sqrt(1 - a^2) + x) - 2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*(2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda))))) - 
            (-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            (2*x*(2*sqrt(1 - a^2) + x)*(a^4 - a^2*(1 + sqrt(1 - a^2) + x) - 
            (-2 + sqrt(1 - a^2) + x)*(1 + sqrt(1 - a^2) + x)^3) + 
            ((1 + sqrt(1 - a^2) + x)^2*(-1 + 2*sqrt(1 - a^2) + 2*x) + 
            a^2*(1 + 2*sqrt(1 - a^2) + 2*x))^2 + (2*a*x*(2*sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*((1 + sqrt(1 - a^2) + x)^2*
            (-1 + 2*sqrt(1 - a^2) + 2*x) + a^2*(1 + 2*sqrt(1 - a^2) + 2*x))*
            (-3*im*a^2*m*(1 + sqrt(1 - a^2) + x) - im*m*(1 + sqrt(1 - a^2) + x)^3 + 
            2*a^3*(1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*(3 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*(-2*im*a^3*m*(1 + sqrt(1 - a^2) + x) - 
            2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda)))) - ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (2*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*(-lambda - im*(1 + sqrt(1 - a^2) + x)*omega) + 
            a^8*(1 + lambda)*(2 + (1 + sqrt(1 - a^2) + x)^2*omega^2) - 
            2*im*a^5*m*(1 + sqrt(1 - a^2) + x)^2*(-6 + (1 + sqrt(1 - a^2) + x)*
            (4 + m^2 - lambda - 4*im*omega) - 3*im*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2) - 2*im*a*m*(1 + sqrt(1 - a^2) + x)^5*
            (3 + (1 + sqrt(1 - a^2) + x)*(-5 + 2*lambda) + (1 + sqrt(1 - a^2) + x)^2*
            (2 - lambda + 2*im*omega) - im*(1 + sqrt(1 - a^2) + x)^3*(3 + lambda)*omega + 
            (1 + sqrt(1 - a^2) + x)^4*omega^2) - 2*im*a^3*m*(1 + sqrt(1 - a^2) + x)^3*
            (13 + (1 + sqrt(1 - a^2) + x)*(-17 + 2*lambda) + (1 + sqrt(1 - a^2) + x)^2*
            (6 + m^2 - 2*lambda - 2*im*omega) - 3*im*(1 + sqrt(1 - a^2) + x)^3*(2 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^4*omega^2) + (1 + sqrt(1 - a^2) + x)^6*
            (-3*(2 + lambda) - (1 + sqrt(1 - a^2) + x)^2*lambda*(2 + lambda) + 
            2*(1 + sqrt(1 - a^2) + x)*(2 + 3*lambda + lambda^2) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda)*omega^2) + a^6*(1 + sqrt(1 - a^2) + x)*(-16*(1 + lambda) + 
            (1 + sqrt(1 - a^2) + x)*(5 + m^2*(-1 + lambda) + 6*lambda - lambda^2 + 2*im*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(5 + 4*lambda)*omega^2 + 2*(1 + sqrt(1 - a^2) + x)^2*omega*
            (im + 2*im*m^2 + omega)) + a^4*(1 + sqrt(1 - a^2) + x)^2*
            (11 + 25*lambda + 2*(1 + sqrt(1 - a^2) + x)*(-5 + 3*m^2 - 13*lambda + lambda^2 - 2*im*omega) + 
            (1 + sqrt(1 - a^2) + x)^2*(5 + 4*lambda - 3*lambda^2 + m^2*(3 + 2*lambda) - 8*im*omega) + 
            3*(1 + sqrt(1 - a^2) + x)^4*(3 + 2*lambda)*omega^2 + 4*(1 + sqrt(1 - a^2) + x)^3*omega*
            (im + 2*im*m^2 + omega)) + a^2*(1 + sqrt(1 - a^2) + x)^3*
            (30 + 5*(1 + sqrt(1 - a^2) + x)*(-7 + 2*lambda) - 4*(1 + sqrt(1 - a^2) + x)^2*
            (-2 + m^2 + lambda - lambda^2 - 3*im*omega) + (1 + sqrt(1 - a^2) + x)^3*
            (2 - 2*lambda - 3*lambda^2 + m^2*(4 + lambda) - 10*im*omega) + (1 + sqrt(1 - a^2) + x)^5*
            (7 + 4*lambda)*omega^2 + 2*(1 + sqrt(1 - a^2) + x)^4*omega*(im + 2*im*m^2 + omega))))/
            ((1 + sqrt(1 - a^2) + x)^2*(-2*im*a^3*m*(1 + sqrt(1 - a^2) + x) - 
            2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + a^4*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*
            (2 + lambda) + a^2*(1 + sqrt(1 - a^2) + x)*(2 + (1 + sqrt(1 - a^2) + x)*
            (3 + 2*lambda)))))))
        end
    elseif s == -1
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-((2*im*a*x*(2*sqrt(1 - a^2) + x)*(3*im*a^2*m*(1 + sqrt(1 - a^2) + x) + 
            im*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(-1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + 
            omega))/((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            (2*im*a^3*m*(1 + sqrt(1 - a^2) + x) + 2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*lambda + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda))))) - 
            (-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            ((a^2 - (1 + sqrt(1 - a^2) + x)^2)^2 + 2*x*(1 + sqrt(1 - a^2) + x)*
            (2*sqrt(1 - a^2) + x)*(-3*a^2 + (1 + sqrt(1 - a^2) + x)^2) + 
            (2*a*x*(2*sqrt(1 - a^2) + x)*(a^2 - (1 + sqrt(1 - a^2) + x)^2)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(3*im*a^2*m*(1 + sqrt(1 - a^2) + x) + 
            im*m*(1 + sqrt(1 - a^2) + x)^3 + 2*a^3*(-1 + lambda) + a*(1 + sqrt(1 - a^2) + x)*
            (3 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda))))/((1 + sqrt(1 - a^2) + x)*
            (2*im*a^3*m*(1 + sqrt(1 - a^2) + x) + 2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*lambda + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))) - 
            ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*(2*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*
            (2 - lambda + im*(1 + sqrt(1 - a^2) + x)*omega) + a^8*(-1 + lambda)*
            (2 + (1 + sqrt(1 - a^2) + x)^2*omega^2) + 2*im*a^5*m*(1 + sqrt(1 - a^2) + x)^3*
            (-2 + m^2 - lambda + 4*im*omega - 3*im*(1 + sqrt(1 - a^2) + x)*omega + 
            3*im*(1 + sqrt(1 - a^2) + x)*lambda*omega + 3*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 
            (1 + sqrt(1 - a^2) + x)^6*lambda*(-3 - (1 + sqrt(1 - a^2) + x)^2*lambda + 
            2*(1 + sqrt(1 - a^2) + x)*(1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            2*im*a*m*(1 + sqrt(1 - a^2) + x)^5*(sqrt(1 - a^2) + x + 
            2*(1 + sqrt(1 - a^2) + x)*lambda - (1 + sqrt(1 - a^2) + x)^2*(lambda + 2*im*omega) + 
            im*(1 + sqrt(1 - a^2) + x)^3*(1 + lambda)*omega + (1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            2*im*a^3*m*(1 + sqrt(1 - a^2) + x)^3*(1 + (1 + sqrt(1 - a^2) + x)*(3 + 2*lambda) + 
            (1 + sqrt(1 - a^2) + x)^2*(-2 + m^2 - 2*lambda + 2*im*omega) + 
            3*im*(1 + sqrt(1 - a^2) + x)^3*lambda*omega + 3*(1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            a^6*(1 + sqrt(1 - a^2) + x)*(8 - 8*lambda + (1 + sqrt(1 - a^2) + x)*
            (-1 + m^2*(-3 + lambda) - lambda^2 - 2*im*omega) + (1 + sqrt(1 - a^2) + x)^3*(-3 + 4*lambda)*
            omega^2 + 2*(1 + sqrt(1 - a^2) + x)^2*omega*(-im - 2*im*m^2 + omega)) + 
            a^4*(1 + sqrt(1 - a^2) + x)^2*(-11 + 9*lambda + 2*(1 + sqrt(1 - a^2) + x)*
            (1 + 3*m^2 + lambda + lambda^2 + 2*im*omega) + (1 + sqrt(1 - a^2) + x)^2*
            (1 - 6*lambda - 3*lambda^2 + m^2*(-1 + 2*lambda) + 8*im*omega) + 3*(1 + sqrt(1 - a^2) + x)^4*
            (-1 + 2*lambda)*omega^2 + 4*(1 + sqrt(1 - a^2) + x)^3*omega*(-im - 2*im*m^2 + omega)) + 
            a^2*(1 + sqrt(1 - a^2) + x)^3*(6 - 3*(1 + sqrt(1 - a^2) + x)*(1 + 2*lambda) - 
            4*(1 + sqrt(1 - a^2) + x)^2*(m^2 - 3*lambda - lambda^2 + 3*im*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(-4*lambda - 3*lambda^2 + m^2*(2 + lambda) + 10*im*omega) + 
            (1 + sqrt(1 - a^2) + x)^5*(-1 + 4*lambda)*omega^2 + 2*(1 + sqrt(1 - a^2) + x)^4*omega*
            (-im - 2*im*m^2 + omega))))/((1 + sqrt(1 - a^2) + x)^2*
            (2*im*a^3*m*(1 + sqrt(1 - a^2) + x) + 2*im*a*m*(1 + sqrt(1 - a^2) + x)^3 + 
            a^4*(-1 + lambda) + (1 + sqrt(1 - a^2) + x)^4*lambda + a^2*(1 + sqrt(1 - a^2) + x)*
            (2 + (1 + sqrt(1 - a^2) + x)*(-1 + 2*lambda)))))))
        end
    elseif s == 2
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 + 
            (8*im*a*x*(2*sqrt(1 - a^2) + x)*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)*
            (im*m*(1 + sqrt(1 - a^2) + x)^2*(6 + (1 + sqrt(1 - a^2) + x)*(4 + lambda)) - 
            3*a^2*m*(1 + sqrt(1 - a^2) + x)*(3*im + 4*(1 + sqrt(1 - a^2) + x)*omega) + 
            a*(1 + sqrt(1 - a^2) + x)*(9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*im*omega) - 
            im*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega) + 
            a^3*(-6 + 9*im*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            ((1 + sqrt(1 - a^2) + x)^4*(24 + 10*lambda + lambda^2 + 12*im*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(im + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(-6*im - 2*im*(1 + sqrt(1 - a^2) + x)*(4 + lambda) + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) - 12*a^4*(-1 + 2*im*(1 + sqrt(1 - a^2) + x)*omega + 
            2*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*im*omega) - 
            2*im*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))) - (((1 + sqrt(1 - a^2) + x)^2*(-1 + 3*sqrt(1 - a^2) + 3*x) + 
            a^2*(1 + 3*sqrt(1 - a^2) + 3*x))^2 + x*(2*sqrt(1 - a^2) + x)*
            (3*a^4 + (1 + sqrt(1 - a^2) + x)^3*(8 - 3*(1 + sqrt(1 - a^2) + x))) + 
            (8*a*x*(2*sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            ((1 + sqrt(1 - a^2) + x)^2*(-1 + 3*sqrt(1 - a^2) + 3*x) + 
            a^2*(1 + 3*sqrt(1 - a^2) + 3*x))*(im*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*(4 + lambda)) - 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*im + 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*im*omega) - 
            im*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega) + 
            a^3*(-6 + 9*im*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*((-(1 + sqrt(1 - a^2) + x)^4)*(24 + 10*lambda + lambda^2 + 
            12*im*omega) - 24*a^3*m*(1 + sqrt(1 - a^2) + x)*(im + 2*(1 + sqrt(1 - a^2) + x)*
            omega) + 4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*im + 2*im*(1 + sqrt(1 - a^2) + x)*
            (4 + lambda) - 3*(1 + sqrt(1 - a^2) + x)^2*omega) + 
            12*a^4*(-1 + 2*im*(1 + sqrt(1 - a^2) + x)*omega + 2*(1 + sqrt(1 - a^2) + x)^2*
            omega^2) + 4*a^2*(1 + sqrt(1 - a^2) + x)*(6 + (1 + sqrt(1 - a^2) + x)*
            (-3 + 6*m^2 - 6*im*omega) - 2*im*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 
            3*(1 + sqrt(1 - a^2) + x)^3*omega^2))) - ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (24*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*(3 + 3*im*(1 + sqrt(1 - a^2) + x)*omega - 
            4*(1 + sqrt(1 - a^2) + x)^2*omega^2) - (1 + sqrt(1 - a^2) + x)^6*
            (24 + 10*lambda + lambda^2 + 12*im*omega)*(-12 - (1 + sqrt(1 - a^2) + x)^2*lambda + 
            2*(1 + sqrt(1 - a^2) + x)*(4 + lambda) + (1 + sqrt(1 - a^2) + x)^4*omega^2) - 
            4*a^5*m*(1 + sqrt(1 - a^2) + x)^2*(-54*im - 2*im*(1 + sqrt(1 - a^2) + x)*
            (-55 + 3*m^2 - 4*lambda - 18*im*omega) + 2*(1 + sqrt(1 - a^2) + x)^2*
            (4 + 12*m^2 - 8*lambda + 27*im*omega)*omega - 6*im*(1 + sqrt(1 - a^2) + x)^3*(8 + lambda)*
            omega^2 + 45*(1 + sqrt(1 - a^2) + x)^4*omega^3) + 
            12*a^8*(-2 - 3*(1 + sqrt(1 - a^2) + x)^2*omega^2 - 2*im*(1 + sqrt(1 - a^2) + x)^3*
            omega^3 + 2*(1 + sqrt(1 - a^2) + x)^4*omega^4) + 4*a^6*(1 + sqrt(1 - a^2) + x)*
            (60 - 3*(1 + sqrt(1 - a^2) + x)*(16 + 3*m^2 + 18*im*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(1 + 36*m^2 - 8*lambda + 18*im*omega)*omega^2 - 
            2*im*(1 + sqrt(1 - a^2) + x)^4*(7 + lambda)*omega^3 + 15*(1 + sqrt(1 - a^2) + x)^5*
            omega^4 - 2*(1 + sqrt(1 - a^2) + x)^2*omega*(-40*im + 9*im*m^2 - 4*im*lambda + 9*omega)) + 
            2*a*m*(1 + sqrt(1 - a^2) + x)^4*(144*im + 8*im*(1 + sqrt(1 - a^2) + x)*
            (-38 + lambda) + (1 + sqrt(1 - a^2) + x)^4*(40 + 20*lambda + lambda^2 + 24*im*omega)*omega + 
            4*im*(1 + sqrt(1 - a^2) + x)^5*(4 + lambda)*omega^2 - 6*(1 + sqrt(1 - a^2) + x)^6*
            omega^3 + 4*(1 + sqrt(1 - a^2) + x)^2*(16*im + 3*im*lambda + 2*im*lambda^2 + 6*omega) - 
            4*im*(1 + sqrt(1 - a^2) + x)^3*(-12 + lambda + lambda^2 - 14*im*omega - 5*im*lambda*omega)) - 
            2*a^3*m*(1 + sqrt(1 - a^2) + x)^3*(312*im + (1 + sqrt(1 - a^2) + x)^3*
            (84 + 30*m^2 - 52*lambda - lambda^2 + 36*im*omega)*omega - 4*im*(1 + sqrt(1 - a^2) + x)^4*
            (19 + 4*lambda)*omega^2 + 48*(1 + sqrt(1 - a^2) + x)^5*omega^3 + 
            12*(1 + sqrt(1 - a^2) + x)*(-63*im + 3*im*m^2 - 3*im*lambda + 32*omega) + 
            4*(1 + sqrt(1 - a^2) + x)^2*(71*im + im*lambda^2 - im*m^2*(10 + lambda) - 104*omega + 
            lambda*(9*im + 16*omega))) + a^2*(1 + sqrt(1 - a^2) + x)^3*
            (576 + 96*(1 + sqrt(1 - a^2) + x)*(-12 + 4*m^2 - 3*im*omega) + 
            8*(1 + sqrt(1 - a^2) + x)^2*(18 - 3*lambda^2 + m^2*(-58 + 8*lambda) + 
            lambda*(-30 - 2*im*omega) + 46*im*omega) - 8*im*(1 + sqrt(1 - a^2) + x)^4*
            (-lambda^2 + 2*m^2*(7 + lambda) + lambda*(2 + 5*im*omega) + 11*im*omega)*omega + 
            2*(1 + sqrt(1 - a^2) + x)^5*(-34 + 24*m^2 - 20*lambda - lambda^2 - 24*im*omega)*omega^2 - 
            8*im*(1 + sqrt(1 - a^2) + x)^6*(1 + lambda)*omega^3 + 12*(1 + sqrt(1 - a^2) + x)^7*
            omega^4 + (1 + sqrt(1 - a^2) + x)^3*(lambda^3 + 36*lambda*(4 + im*omega) + 
            2*lambda^2*(11 - 8*im*omega) - m^2*(-136 + 42*lambda + lambda^2 - 36*im*omega) - 
            8*(-18 + 19*im*omega + 6*omega^2))) + a^4*(1 + sqrt(1 - a^2) + x)^2*
            (-672 + (1 + sqrt(1 - a^2) + x)*(960 - 72*m^2 + 624*im*omega) + 
            (1 + sqrt(1 - a^2) + x)^4*(44 + 180*m^2 - 62*lambda - lambda^2 + 36*im*omega)*omega^2 - 
            8*im*(1 + sqrt(1 - a^2) + x)^5*(5 + 2*lambda)*omega^3 + 48*(1 + sqrt(1 - a^2) + x)^6*
            omega^4 + 8*(1 + sqrt(1 - a^2) + x)^3*omega*(56*im + im*lambda^2 - 3*im*m^2*(9 + lambda) - 
            46*omega + lambda*(6*im + 8*omega)) + 4*(1 + sqrt(1 - a^2) + x)^2*
            (6*m^4 + m^2*(7 - 8*lambda + 54*im*omega) + 2*(-15 + 10*lambda + lambda^2 - 144*im*omega - 9*im*
            lambda*omega + 48*omega^2)))))/((1 + sqrt(1 - a^2) + x)^2*
            ((-(1 + sqrt(1 - a^2) + x)^4)*(24 + 10*lambda + lambda^2 + 12*im*omega) - 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(im + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*im + 2*im*(1 + sqrt(1 - a^2) + x)*(4 + lambda) - 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + 12*a^4*(-1 + 2*im*(1 + sqrt(1 - a^2) + x)*
            omega + 2*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 - 6*im*omega) - 
            2*im*(1 + sqrt(1 - a^2) + x)^2*(1 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))))/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4))
        end
    elseif s == -2
        return begin
            (1/(2*sqrt(1 - a^2) + x)^2)*((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)^2 + 
            (8*im*a*x*(2*sqrt(1 - a^2) + x)*(-((a*m)/(2*(1 + sqrt(1 - a^2)))) + omega)*
            ((-im)*m*(1 + sqrt(1 - a^2) + x)^2*(6 + (1 + sqrt(1 - a^2) + x)*lambda) + 
            3*a^2*m*(1 + sqrt(1 - a^2) + x)*(3*im - 4*(1 + sqrt(1 - a^2) + x)*omega) + 
            a*(1 + sqrt(1 - a^2) + x)*(9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*im*omega) + 
            im*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega) + 
            a^3*(-6 - 9*im*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*(a^2 + (1 + sqrt(1 - a^2) + x)^2)*
            ((1 + sqrt(1 - a^2) + x)^4*(2*lambda + lambda^2 - 12*im*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(-im + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*im + 2*im*(1 + sqrt(1 - a^2) + x)*lambda + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + a^4*(12 + 24*im*(1 + sqrt(1 - a^2) + x)*omega - 
            24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*im*omega) + 
            2*im*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))) - (1/(a^2 + (1 + sqrt(1 - a^2) + x)^2)^4)*
            ((a^2*(-1 + sqrt(1 - a^2) + x) + (1 + sqrt(1 - a^2) + x)^3)^2 + 
            x*(2*sqrt(1 - a^2) + x)*(-a^4 - 8*a^2*(1 + sqrt(1 - a^2) + x) + 
            (1 + sqrt(1 - a^2) + x)^4) + (8*a*x*(2*sqrt(1 - a^2) + x)*
            (a^2 + (1 + sqrt(1 - a^2) + x)^2)*(a^2*(-1 + sqrt(1 - a^2) + x) + 
            (1 + sqrt(1 - a^2) + x)^3)*((-im)*m*(1 + sqrt(1 - a^2) + x)^2*
            (6 + (1 + sqrt(1 - a^2) + x)*lambda) + 3*a^2*m*(1 + sqrt(1 - a^2) + x)*
            (3*im - 4*(1 + sqrt(1 - a^2) + x)*omega) + a*(1 + sqrt(1 - a^2) + x)*
            (9 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*im*omega) + 
            im*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega) + 
            a^3*(-6 - 9*im*(1 + sqrt(1 - a^2) + x)*omega + 6*(1 + sqrt(1 - a^2) + x)^2*omega^2)))/
            ((1 + sqrt(1 - a^2) + x)*((1 + sqrt(1 - a^2) + x)^4*(2*lambda + lambda^2 - 12*im*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(-im + 2*(1 + sqrt(1 - a^2) + x)*omega) + 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*im + 2*im*(1 + sqrt(1 - a^2) + x)*lambda + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + a^4*(12 + 24*im*(1 + sqrt(1 - a^2) + x)*omega - 
            24*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*im*omega) + 
            2*im*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))) - ((a^2 + (1 + sqrt(1 - a^2) + x)^2)^2*
            (-24*a^7*m*(1 + sqrt(1 - a^2) + x)^2*omega*(-3 + 3*im*(1 + sqrt(1 - a^2) + x)*omega + 
            4*(1 + sqrt(1 - a^2) + x)^2*omega^2) - 2*a^3*m*(1 + sqrt(1 - a^2) + x)^3*
            (72*im - 4*im*(1 + sqrt(1 - a^2) + x)*(15 + 9*m^2 - 13*lambda) + 
            4*(1 + sqrt(1 - a^2) + x)^2*(im*m^2*(6 + lambda) - im*(3 + lambda^2 + lambda*
            (9 + 16*im*omega) + 24*im*omega)) + (1 + sqrt(1 - a^2) + x)^3*
            (-36 + 30*m^2 - 44*lambda - lambda^2 - 36*im*omega)*omega + 4*im*(1 + sqrt(1 - a^2) + x)^4*
            (3 + 4*lambda)*omega^2 + 48*(1 + sqrt(1 - a^2) + x)^5*omega^3) + 
            12*a^8*(-2 - 3*(1 + sqrt(1 - a^2) + x)^2*omega^2 + 2*im*(1 + sqrt(1 - a^2) + x)^3*
            omega^3 + 2*(1 + sqrt(1 - a^2) + x)^4*omega^4) + 4*a^6*(1 + sqrt(1 - a^2) + x)*
            (12 - 3*(1 + sqrt(1 - a^2) + x)*(-4 + 3*m^2 + 6*im*omega) + 
            (1 + sqrt(1 - a^2) + x)^3*(-39 + 36*m^2 - 8*lambda - 18*im*omega)*omega^2 + 
            2*im*(1 + sqrt(1 - a^2) + x)^4*(3 + lambda)*omega^3 + 15*(1 + sqrt(1 - a^2) + x)^5*
            omega^4 + 2*(1 + sqrt(1 - a^2) + x)^2*omega*(24*im + 9*im*m^2 - 4*im*lambda + 15*omega)) - 
            4*a^5*m*(1 + sqrt(1 - a^2) + x)^2*(-18*im + 2*(1 + sqrt(1 - a^2) + x)^2*
            (-36 + 12*m^2 - 8*lambda - 27*im*omega)*omega + 6*im*(1 + sqrt(1 - a^2) + x)^3*(4 + lambda)*
            omega^2 + 45*(1 + sqrt(1 - a^2) + x)^4*omega^3 + 2*(1 + sqrt(1 - a^2) + x)*
            (9*im + 3*im*m^2 - 4*im*lambda + 30*omega)) - 2*a*m*(1 + sqrt(1 - a^2) + x)^4*
            (-48*im - 24*im*(1 + sqrt(1 - a^2) + x)*(-2 + lambda) + 
            4*im*(1 + sqrt(1 - a^2) + x)^2*(7*lambda + 2*lambda^2 + 6*im*omega) - 
            (1 + sqrt(1 - a^2) + x)^4*(12*lambda + lambda^2 - 24*im*omega)*omega + 
            4*im*(1 + sqrt(1 - a^2) + x)^5*lambda*omega^2 + 6*(1 + sqrt(1 - a^2) + x)^6*omega^3 - 
            4*im*(1 + sqrt(1 - a^2) + x)^3*(lambda + lambda^2 + 6*im*omega + 5*im*lambda*omega)) + 
            a^4*(1 + sqrt(1 - a^2) + x)^3*(24*m^4*(1 + sqrt(1 - a^2) + x) + 
            48*(-4 + 3*im*omega) + 8*(1 + sqrt(1 - a^2) + x)*(9 + lambda^2 + lambda*(2 + 13*im*omega) - 
            72*im*omega) - (1 + sqrt(1 - a^2) + x)^3*(60 + 54*lambda + lambda^2 + 36*im*omega)*omega^2 + 
            8*im*(1 + sqrt(1 - a^2) + x)^4*(-3 + 2*lambda)*omega^3 + 48*(1 + sqrt(1 - a^2) + x)^5*
            omega^4 + 8*(1 + sqrt(1 - a^2) + x)^2*omega*(24*im - 6*im*lambda - im*lambda^2 + 18*omega + 
            8*lambda*omega) + 4*m^2*(30 - (1 + sqrt(1 - a^2) + x)*(33 + 8*lambda + 54*im*omega) + 
            6*im*(1 + sqrt(1 - a^2) + x)^2*(5 + lambda)*omega + 45*(1 + sqrt(1 - a^2) + x)^3*
            omega^2)) + (1 + sqrt(1 - a^2) + x)^6*((-1 + sqrt(1 - a^2) + x)*
            (1 + sqrt(1 - a^2) + x)*lambda^3 + lambda^2*(12 - 12*(1 + sqrt(1 - a^2) + x) + 
            2*(1 + sqrt(1 - a^2) + x)^2 - (1 + sqrt(1 - a^2) + x)^4*omega^2) + 
            12*im*omega*(-12 + 8*(1 + sqrt(1 - a^2) + x) + (1 + sqrt(1 - a^2) + x)^4*omega^2) - 
            2*lambda*(-12 + 4*(1 + sqrt(1 - a^2) + x)*(2 - 3*im*omega) + 
            6*im*(1 + sqrt(1 - a^2) + x)^2*omega + (1 + sqrt(1 - a^2) + x)^4*omega^2)) + 
            a^2*(1 + sqrt(1 - a^2) + x)^4*(96 + 8*(1 + sqrt(1 - a^2) + x)*
            (m^2*(6 + 8*lambda) - 3*(2 + lambda^2 + lambda*(2 + 2*im*omega) - 22*im*omega)) - 96*im*omega + 
            2*(1 + sqrt(1 - a^2) + x)^4*(6 + 24*m^2 - 12*lambda - lambda^2 + 24*im*omega)*omega^2 + 
            8*im*(1 + sqrt(1 - a^2) + x)^5*(-3 + lambda)*omega^3 + 12*(1 + sqrt(1 - a^2) + x)^6*
            omega^4 + (1 + sqrt(1 - a^2) + x)^2*(lambda^3 + lambda*(24 - 34*m^2 - 4*im*omega) + 
            lambda^2*(14 - m^2 + 16*im*omega) - 12*im*(22 + 3*m^2 - 4*im*omega)*omega) + 
            8*(1 + sqrt(1 - a^2) + x)^3*omega*((-im)*lambda^2 + 2*im*m^2*(3 + lambda) + 3*omega + 
            lambda*(2*im + 5*omega)))))/((1 + sqrt(1 - a^2) + x)^2*
            ((-(1 + sqrt(1 - a^2) + x)^4)*(2*lambda + lambda^2 - 12*im*omega) + 
            24*a^3*m*(1 + sqrt(1 - a^2) + x)*(im - 2*(1 + sqrt(1 - a^2) + x)*omega) - 
            4*a*m*(1 + sqrt(1 - a^2) + x)^2*(6*im + 2*im*(1 + sqrt(1 - a^2) + x)*lambda + 
            3*(1 + sqrt(1 - a^2) + x)^2*omega) + 12*a^4*(-1 - 2*im*(1 + sqrt(1 - a^2) + x)*
            omega + 2*(1 + sqrt(1 - a^2) + x)^2*omega^2) + 4*a^2*(1 + sqrt(1 - a^2) + x)*
            (6 + (1 + sqrt(1 - a^2) + x)*(-3 + 6*m^2 + 6*im*omega) + 
            2*im*(1 + sqrt(1 - a^2) + x)^2*(-3 + lambda)*omega + 3*(1 + sqrt(1 - a^2) + x)^3*
            omega^2))))))
        end
    else
        throw(DomainError(s, "Currently only spin weight s of 0, +/-1, +/-2 are supported"))
    end
end

