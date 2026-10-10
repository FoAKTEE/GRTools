# ∂(r1, r2, z1²)/∂(En, Lz, K) — implicit differentiation of the turning-point conditions
# R(r1)=R(r2)=0 and Θ(z1)=0 (SpinningTrajectory.wl `dpexdCFun`, first block).
"""
    turning_point_jacobian(a, p, e, x)

Return `∂(r_apo,r_peri,z_max²)/∂(E,Lz,K)` at fixed spin, where
`K = Q + (Lz-aE)²`. Requires distinct radial turning points (e > 0),
a bound orbit away from the separatrix, and a nondegenerate polar root.
Equatorial negative-a constants are normalized by the constants owner; the
derivatives here are with respect to the original signed-a coordinates.
"""
function turning_point_jacobian(a, p, e, x)
    En, Lz, Q = constants_of_motion(a, p, e, x)
    Δ(r) = r^2 - 2r + a^2
    dΔ(r) = 2r - 2
    dRdr(r) = 4 * En * r * (-a * Lz + En * (a^2 + r^2)) - 2 * r * Δ(r) - ((-a * En + Lz)^2 + Q + r^2) * dΔ(r)
    dRdEn(r) = 2 * (a^2 + r^2) * (a^2 * En - a * Lz + En * r^2)
    dRdLz(r) = 2 * a^2 * Lz - 2 * a * En * (a^2 + r^2)
    dRdK(r) = -Δ(r)
    dZdzsq(z) = -(Lz^2 + Q + a^2 * (-1 + En^2) * (-1 + 2 * z^2))
    dZdEn(z) = 2 * a * (1 - z^2) * (Lz - a * En * (1 - z^2))
    dZdLz(z) = -2 * (Lz - a * En * (1 - z^2))
    dZdK(z) = 1 - z^2
    r1 = p / (1 - e); r2 = p / (1 + e); z1 = sqrt(1 - x^2)
    return -vcat([dRdEn(r1) dRdLz(r1) dRdK(r1)] / dRdr(r1),
                 [dRdEn(r2) dRdLz(r2) dRdK(r2)] / dRdr(r2),
                 [dZdEn(z1) dZdLz(z1) dZdK(z1)] / dZdzsq(z1))
end

"""
    orbital_params_jacobian(a, p, e, x) -> ∂(p,e,x)/∂(En,Lz,K)

Analytic turning-point Jacobian (3×3), SpinningTrajectory.wl `dpexdCFun`: the chain of
`∂(p,e,x)/∂(r1,r2,z1²)` with the implicit `∂(r1,r2,z1²)/∂C`. Replaces `SpinningOrbit._dpexdC`.
"""
function orbital_params_jacobian(a, p, e, x)
    z = zero(float(p))
    dpexdr12z1sq = [0.5*(-1 + e)^2      0.5*(1 + e)^2      z;
                    ((-1 + e)^2*(1 + e))/(2p)  ((-1 + e)*(1 + e)^2)/(2p)  z;
                    z                   z                  -(1/(2x))]
    return dpexdr12z1sq * turning_point_jacobian(a, p, e, x)
end

