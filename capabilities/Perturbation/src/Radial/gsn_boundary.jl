# :GSN backend — boundary values (X, dX/dr*) of the GSN homogeneous solutions from the asymptotic
# ansatzes (GSN paper arXiv:2306.16469 Eqs. (33)-(41) region):
#   Up at ∞     : X ~ e^{+iωr*}·f(r),  f = Σᵢ aᵢ^out/(ωr)^i     (outgoing-∞ coefficients)
#   In at horizon: X ~ e^{−ipr*}·g(r),  g = Σᵢ aᵢ^in·(r−r₊)^i,  p = ω − mΩ_H
# Port of src/Homogeneous/InitialConditions.jl (v0.8.5 commit 982b087) with the audited typed-
# zero accumulator fix (`ans = 0.0` → zero(Complex{T}), incl. the empty-loop order=0 derivative
# return — GSN_PLAN §3c-3). The dr* → dr conversion factor is dr/dr* = Δ/(r²+a²).

"Asymptotic ∞-series template Σᵢ₌₀ func(i)/(ωr)^i. [GSN InitialConditions.jl fansatz]"
function _gsn_fansatz(func, omega, r; order::Int=3)
    ans = zero(func(0)*one(omega*r))
    for i in 0:order
        ans += func(i)/((omega*r)^i)
    end
    return ans
end

"d/dr of the ∞-series template. [GSN InitialConditions.jl dfansatz_dr]"
function _gsn_dfansatz_dr(func, omega, r; order::Int=3)
    ans = zero(func(0)*one(omega*r))
    for i in 1:order
        ans += -i*func(i)/(omega^i * r^(i + 1))
    end
    return ans
end

"Horizon-series template Σᵢ₌₀ func(i)·(r−r₊)^i. [GSN InitialConditions.jl gansatz]"
function _gsn_gansatz(func, a, r; order::Int=1)
    ans = zero(func(0)*one(r))
    for i in 0:order
        ans += func(i)*(r - _gsn_rp(a))^i
    end
    return ans
end

"d/dr of the horizon-series template. [GSN InitialConditions.jl dgansatz_dr]"
function _gsn_dgansatz_dr(func, a, r; order::Int=1)
    ans = zero(func(0)*one(r))
    for i in 1:order
        ans += i*func(i)*(r - _gsn_rp(a))^(i - 1)
    end
    return ans
end

"""
    _gsn_Xup_bc(s, m, a, omega, lambda, rsout, order, ::Type{T}) -> (X, dX/dr*)

Purely-outgoing boundary values at r* = rsout: X = e^{iωr*}f(r), dX/dr* = e^{iωr*}(iωf +
(Δ/(r²+a²))f′). [GSN InitialConditions.jl Xup_initialconditions; upstream default order 3]
"""
function _gsn_Xup_bc(s::Int, m::Int, a, omega, lambda, rsout, order::Int, ::Type{T}) where {T<:Real}
    func(i) = _gsn_coeff_out_inf(s, m, a, omega, lambda, i, T)
    rout = _gsn_r_from_rstar(a, rsout)
    f = _gsn_fansatz(func, omega, rout; order=order)
    dfdr = _gsn_dfansatz_dr(func, omega, rout; order=order)
    phase = exp(im*T(omega)*T(rsout))
    return phase*f, phase*(im*omega*f + (_gsn_Delta(a, rout)/(rout^2 + a^2))*dfdr)
end

"""
    _gsn_Xin_bc(s, m, a, omega, lambda, rsin, order, ::Type{T}) -> (X, dX/dr*)

Purely-ingoing boundary values at r* = rsin: X = e^{−ipr*}g(r), p = ω − mΩ_H, dX/dr* =
e^{−ipr*}(−ipg + (Δ/(r²+a²))g′). [GSN InitialConditions.jl Xin_initialconditions; upstream
default order 0]
"""
function _gsn_Xin_bc(s::Int, m::Int, a, omega, lambda, rsin, order::Int, ::Type{T}) where {T<:Real}
    func(i) = _gsn_coeff_in_hor(s, m, a, omega, lambda, i, T)
    rin = _gsn_r_from_rstar(a, rsin)
    g = _gsn_gansatz(func, a, rin; order=order)
    dgdr = _gsn_dgansatz_dr(func, a, rin; order=order)
    p = omega - m*_gsn_Omega_h(a)
    phase = exp(-im*T(p)*T(rsin))
    return phase*g, phase*(-im*p*g + (_gsn_Delta(a, rin)/(rin^2 + a^2))*dgdr)
end
