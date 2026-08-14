# App-H jump conditions — spin-general source jumps feeding the linear-system closure (eq:jumps).
#
# Consolidated from KerrLorenzCirc jumps.jl into ONE spin-general jumps_pm(bg,l,sgn) (the a=0 forms
# are the literal a==0 limit) + hj1(bg,l) trace jump + jumps_eigenvalue.
#   jumps_pm(bg,l,sgn): sgn=−1 → J0/J1(P_{-2}); sgn=+1 → J0/J1(P_{+2}). S0=SWSH(2sgn,l,m,aω).
#   hj1(bg,l):          J1(h)=−16π/(ut·Δ0)·{}_0S_{lm}(π/2);  J0(h)=0.
#
# PRECISION NOTE: the spheroidal harmonic S(π/2) and its θ-derivative are Float64-capped (the SWSH
# SphericalExpansion LAPACK path). The background scalars (calA,calB,ut,Δ0 in `bg`) carry the input
# precision; the angular S-factors are evaluated once per (l,m) at the equator and well-conditioned.

import Harmonics: SpinWeightedSpheroidalHarmonicS, SpinWeightedSpheroidalEigenvalue

"Spheroidal s=−2 eigenvalue Λ at this orbit's aω (Float64-capped angular value)."
jumps_eigenvalue(bg::OrbitBackground, l) = real(SpinWeightedSpheroidalEigenvalue(-2, l, bg.m, Float64(bg.aω)))

# Equatorial spheroidal harmonic value + θ-derivative (Float64 angular cap, accepted).
function _S_and_deriv(s, l, m, aω)
    aωf = Float64(aω)
    S = SpinWeightedSpheroidalHarmonicS(s, l, m, aωf)
    S0  = real(S(π/2, 0.0))
    S0p = (real(S(π/2+1e-6, 0.0)) - real(S(π/2-1e-6, 0.0)))/2e-6
    (S0, S0p)
end

"""
    jumps_pm(bg::OrbitBackground, l, sgn) -> (J0, J1)

Spin-general s=∓2 Teukolsky SOURCE jumps at r0 (eq:jumps). `sgn=-1` → (J0,J1)(P_{-2}); `sgn=+1` →
(J0,J1)(P_{+2}). Returns Complex{T} jumps where T is the background element type (BigFloat-capable for
the radial coefficients; the entering S-factors are Float64). Validated vs App-H upper tables (~1.8e-8).
"""
function jumps_pm(bg::OrbitBackground{T}, l::Integer, sgn::Integer) where {T}
    a = bg.a; r0 = bg.r0; m = bg.m; Ω = bg.Ω; ut = bg.ut
    Δ = bg.Δ0; Δp = 2r0 - 2; calA = bg.calA; calB = bg.calB
    cD = 4*T(π)/(r0^2*ut)
    S0f, S0pf = _S_and_deriv(2sgn, l, m, bg.aω)
    S0 = T(S0f); S0p = T(S0pf)
    Λ = T(jumps_eigenvalue(bg, l))
    meff = sgn*m
    W0  = meff*(Ω*(r0^2+a^2)-a)/Δ
    W0p = meff*(Ω*2r0*Δ - (Ω*(r0^2+a^2)-a)*Δp)/Δ^2
    Q0  = meff*(1 - a*Ω)
    c2 = calB^2*S0
    c1 = 2calB*(im*calA*(S0p+Q0*S0) + calB*(-im*W0 + 1/r0)*S0)
    c0 = 2calA*(im*calB*(2/r0 - im*W0) + calA*(-Q0 + im*a/r0))*(S0p+Q0*S0) +
         (calA^2*Λ + calB^2*(im*W0p - 2im*W0/r0 - W0^2))*S0
    V0m4 = Δ*(W0^2 + im*W0p) - im*Δp*W0 + 6im*meff*Ω*r0 - Λ
    (cD*Δ*(c1 - Δp/Δ*c2), cD*Δ*(c0 - V0m4/Δ*c2))
end

"""
    hj1(bg::OrbitBackground, l) -> J1

Trace-derivative jump J1(h) = −16π/(u^t·Δ0)·{}_0S_{lm}(π/2). J0(h)=0. Vanishes for l+m odd.
Validated vs App-H J1(h) column, both spins (~1e-8).
"""
function hj1(bg::OrbitBackground{T}, l::Integer) where {T}
    S0 = T(real(SpinWeightedSpheroidalHarmonicS(0, l, bg.m, Float64(bg.aω))(π/2, 0.0)))
    -16*T(π)/(bg.ut*bg.Δ0)*S0
end
