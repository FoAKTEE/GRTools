# First-law binding-energy closed forms over the redshift interpolant `z` (1PAT1.m L14;
# 1PAT1R.m L84-89). `z` is the order-8 redshift [`Interp1`](@ref); all expressions take `M = 1`
# as upstream. `x = M/r0` is the dimensionless inverse radius. Pure algebra leaf consumed by the
# 1PA forcing.

# Algebraic (z-independent) tail of EFLx and its derivative.
_eflx_tail(x) = -1 + sqrt(1 - 3x) + (x / 6) * (7 - 24x) / (1 - 3x)^(3 / 2)
function _eflx_tail_deriv(x)
    s = 1 - 3x
    return -3 / (2 * sqrt(s)) + ((7 - 48x) / 6) * s^(-3 / 2) +
           ((7x - 24x^2) / 6) * (9 / 2) * s^(-5 / 2)
end

"""
    E_FLx(z, x)

First-law binding-energy flux function `EFLx[x] = ½ z(x) - (x/3) z'(x) - 1 + √(1-3x)
+ (x/6)(7-24x)/(1-3x)^{3/2}` (1PAT1.m L14), where `z` is the redshift [`Interp1`](@ref).
"""
E_FLx(z::Interp1, x::Real) = (1 // 2) * z(x) - (x / 3) * derivative(z, x) + _eflx_tail(x)

"""
    E_FLx_deriv(z, x)

`EFLx'[x]` (used by the inspiral forcing and by [`dEdΩ`](@ref)): `(1/6) z'(x) - (x/3) z''(x)
+ tail'(x)`, where `z''` is the second derivative of the interpolant.
"""
E_FLx_deriv(z::Interp1, x::Real) =
    (1 // 6) * derivative(z, x) - (x / 3) * derivative2(z, x) + _eflx_tail_deriv(x)

"""
    E_dm(y)

Secondary-mass binding-energy correction `EΔm[y] = -(y/3)(1-6y)/(1-3y)^{3/2}` (1PAT1R.m). Pure
algebra (no self-force data).
"""
E_dm(y::Real) = -(y / 3) * (1 - 6y) / (1 - 3y)^(3 / 2)

"""
    E_dm_deriv(y)

Analytic derivative `EΔm'[y]`.
"""
function E_dm_deriv(y::Real)
    s = 1 - 3y
    u = -(y - 6y^2) / 3          # EΔm = u · s^{-3/2}
    du = -(1 - 12y) / 3
    return du * s^(-3 / 2) + u * (9 / 2) * s^(-5 / 2)
end

"""
    dEdΩ(z, r0, ν, at1, at2, δm1)

Energy flux per unit orbital frequency `dEdΩ` (1PAT1R.m L84-87): the geodesic term plus
secondary-spin (`at1`, `at2`) and self-force/`δm1` corrections, with `x = 1/r0`,
`dx/dΩ = (2/3)√x` (M = 1). Reuses [`E_FLx_deriv`](@ref) and [`E_dm_deriv`](@ref).
"""
function dEdΩ(z::Interp1, r0::Real, ν::Real, at1::Real, at2::Real, δm1::Real)
    M = one(float(r0))
    x = M / r0
    dxdΩ = (2M) / (3 * sqrt(x))
    geo = (M^2 * (-1 + 6x)) / (3 * (1 - 3x)^(3 / 2) * sqrt(x))
    spin = (at2 * M^2 * x * (-5 + 12x)) / (3 * (1 - 3x)^(3 / 2)) -
           (2 * at1 * M^2 * x * (10 - 33x + 36x^2)) / (9 * (1 - 3x)^(5 / 2))
    sf = ν * (E_FLx_deriv(z, x) + δm1 * E_dm_deriv(x)) * dxdΩ
    return geo + spin + sf
end
