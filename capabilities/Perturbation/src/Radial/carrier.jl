# Carrier for a normalized homogeneous radial Teukolsky solution (one boundary condition). Shared
# across the radial backends; this round only the static (ω=0) backend populates it. Callable
# `R(r)` dispatches on `R.method`; the static branch routes to `_static_radial`. The float type T
# is inferred from the inputs (Float64 by default; BigFloat preserved end-to-end so the radial
# value is evaluated at full precision — no ComplexF64 downcast).

"""
    TeukolskyRadialFunction

Carrier for one normalized homogeneous radial Teukolsky solution (one boundary condition `bc`).
Callable as `R(r)`; `R.amplitudes` holds the asymptotic amplitudes; `R.method` selects the backend
(`:Static` this round). Type parameter `T<:Real` carries the working precision.
"""
mutable struct TeukolskyRadialFunction{T<:Real}
    const s::Int
    const l::Int
    const m::Int
    const a::T
    const omega::T
    const nu::Complex{T}         # renorm angular momentum (= λ on the static path)
    const lambda::Complex{T}     # SWSH eigenvalue (= l(l+1)-s(s+1) statically)
    const bc::Symbol             # :In | :Up
    const norm::Complex{T}       # raw transmission amplitude used to normalize
    const amplitudes::NamedTuple # static ℋ/ℐ amplitudes (or the unit-transmission set)
    const method::Symbol         # :Static | :MST | :HeunC | :NumericalIntegration | :SasakiNakamura | :GSN
    cache::Any                   # lazily-built dense ODE solution (numerical backends; later rounds)
    function TeukolskyRadialFunction(s::Integer, l::Integer, m::Integer, a::Real, omega::Real,
            nu::Number, lambda::Number, bc::Symbol, norm::Number, amplitudes, method::Symbol)
        T = promote_type(typeof(float(a)), typeof(float(omega)),
                         typeof(float(real(nu))), typeof(float(real(lambda))), typeof(float(real(norm))))
        new{T}(s, l, m, convert(T, a), convert(T, omega), convert(Complex{T}, nu),
               convert(Complex{T}, lambda), bc, convert(Complex{T}, norm), amplitudes, method, nothing)
    end
end

"Normalized radial value R(r). Dispatches on `R.method`: static (ω=0) → closed-form 2F1; :MST →
series; :HeunC → confluent-Heun; :SasakiNakamura → SN X→R back-transform (s=-2 Up)."
function (R::TeukolskyRadialFunction)(r::Real)
    R.method === :Static && return _static_radial(R.s, R.l, R.m, R.a, R.bc, r; norm=R.norm)
    R.method === :MST && return R.bc === :In ?
        _mst_radial_in(R.s, R.m, R.a, R.omega, R.nu, R.lambda, r; norm=R.norm) :
        _mst_radial_up(R.s, R.m, R.a, R.omega, R.nu, R.lambda, r; norm=R.norm)
    R.method === :HeunC && return _heunc_radial_value(R, r)
    R.method === :SasakiNakamura && return _sn_value(R, r, 0)
    R.method === :GSN && return _gsn_value(R, r, 0)
    error("TeukolskyRadialFunction: method :$(R.method) is not built in this round (:Static, :MST, :HeunC, :SasakiNakamura, :GSN).")
end
