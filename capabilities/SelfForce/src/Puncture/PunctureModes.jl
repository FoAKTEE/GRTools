# First-order puncture-mode amplitudes hP[i][l,m,r0,M,Δr] and their one-sided radial derivatives
# dhP. Port of Punctures/puncture_modes.jl (h1P-exact-order-eps2.wl, Wardell et al.
# arXiv:2508.00087). Each mode is `prefactor(i,l,m,r0,M,Δr) × ∫_{-π}^{π} integrand_i(β) dβ`, where the
# integrand is a machine-generated rational function of (r0,M,Δr,sinβ) contracted against the
# AlphaIntegration α-moments (`alphaint`/`alphaintv*`/`alphaintt*`) at the substituted radial argument
# δ(β). The large auto-generated `_hP_prefactor`/`_hP_integrand` tables live in
# `puncture_modes_generated.jl` (vendored verbatim, emitted by oracle/gen_puncture_modes.wl); this
# leaf supplies the Gauss-Kronrod β-quadrature (consuming `QuadGK`) and the public entry points.
#
# Type-generic {T}: r0/M/Δr promote to `float(promote_type(...))`, so BigFloat precision is inherited
# (the generated rationals stay exact via `rat`). dhP differentiates the (pref·integrand) w.r.t. Δr
# EXACTLY by forward-mode AD (`ForwardDiff`) inside the β-integral — the same generated integrand as
# hP — for n ∈ {1,2}.

import QuadGK
import ForwardDiff

"""
    rat(a, b) -> Rational

Exact rational literal `a // b`, used pervasively by the machine-generated puncture-mode
prefactor/integrand tables (`puncture_modes_generated.jl`). Keeping the generated coefficients and
fractional powers as exact `Rational`s (instead of pre-evaluated floats) is what lets a BigFloat
evaluation of `hP`/`dhP` inherit full working precision rather than collapsing to Float64.
"""
rat(a::Integer, b::Integer) = a // b

# auto-generated _hP_prefactor / _hP_integrand (the inert AINT/AV*/AT* tags bind to the
# AlphaIntegration alphaint* calls already in scope in this Puncture module).
include("puncture_modes_generated.jl")

"""
    hP(i, l, m, r0, M, Δr; order=8, rtol=…)

First-order puncture mode `hP[i][l,m,r0,M,Δr]`. `i ∈ 1..10` selects the tensor/parity component; the
β-integral is evaluated by Gauss-Kronrod quadrature. Real for `i ∈ 1..7`, complex for `i ∈ 8,9,10`.
Inherits the float precision of `r0`,`M`,`Δr`.
"""
function hP(i::Integer, l::Integer, m::Integer, r0::Real, M::Real, Δr::Real;
            order::Integer = 8, rtol::Union{Nothing,Real} = nothing)
    T = float(promote_type(typeof(r0), typeof(M), typeof(Δr)))
    r0T, MT, ΔrT = T(r0), T(M), T(Δr)
    tol = rtol === nothing ? sqrt(eps(T)) : T(rtol)
    pref = _hP_prefactor(i, l, m, r0T, MT, ΔrT)
    val, _ = QuadGK.quadgk(β -> _hP_integrand(i, l, m, r0T, MT, ΔrT, β),
                           -T(π), T(π); order = order, rtol = tol)
    return pref * val
end

"""
    dhP(i, l, m, r0, M, Δr, n=1; order=8, rtol=…)

`n`-th one-sided radial derivative `∂_{Δr}^n hP[i]`, `n ∈ {1,2}`. The β-prefactor is Δr-only and pulls
inside the β-integral, so `dhP(n) = ∫ ∂_{Δr}^n(pref·integrand) dβ`; the inner Δr-derivative is taken
exactly by forward-mode AD over the same generated integrand/prefactor as `hP`. Real for `i ∈ 1..7`,
complex for `i ∈ 8,9,10`.
"""
function dhP(i::Integer, l::Integer, m::Integer, r0::Real, M::Real, Δr::Real,
             n::Integer = 1; order::Integer = 8, rtol::Union{Nothing,Real} = nothing)
    n in (1, 2) || throw(ArgumentError("dhP: only n ∈ {1,2} supported (got $n)"))
    T = float(promote_type(typeof(r0), typeof(M), typeof(Δr)))
    r0T, MT, Δr0 = T(r0), T(M), T(Δr)
    tol = rtol === nothing ? sqrt(eps(T)) : T(rtol)
    f(d, β) = _hP_prefactor(i, l, m, r0T, MT, d) * _hP_integrand(i, l, m, r0T, MT, d, β)
    dintegrand = if n == 1
        β -> ForwardDiff.derivative(d -> f(d, β), Δr0)
    else
        β -> ForwardDiff.derivative(x -> ForwardDiff.derivative(d -> f(d, β), x), Δr0)
    end
    val, _ = QuadGK.quadgk(dintegrand, -T(π), T(π); order = order, rtol = tol)
    return val
end
