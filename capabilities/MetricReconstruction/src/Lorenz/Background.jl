# Orbit-kinematics leaf — circular-equatorial Kerr geodesic scalars bundled in OrbitBackground{T}.
#
# Consolidated circular-equatorial orbit scalars at r0 (these were re-derived in ~6 legacy scripts).
# Returns E, L, Ω, ut, calA, calB, Δ0, ω=mΩ, aω=amΩ. Type-generic (BigFloat r0,a → BigFloat scalars).
# Uses the INLINED equatorial inverse-metric (kerr_metric_inv_eq from Operators) for u^t — NO DiscModes.
#
# Forms (M=1, prograde):  v=1/√r0; den=√(1−3v²+2a v³);  E=(1−2v²+a v³)/den;
#   L=r0 v(1−2a v³+a²v⁴)/den;  Ω=v³/(1+a v³)=1/(√(r0³)+a);  u^t=−E g^{tt}+L g^{tφ};  Δ0=r0²−2r0+a²;
#   calA=(E(r0²+a²)−aL)/Δ0;  calB=L−aE.

"""
    OrbitBackground{T}

Circular-equatorial Kerr orbit-kinematics bundle at radius `r0`, spin `a`, azimuthal `m`: the
constants of motion `E,L`, orbital frequency `Ω`, time-component `ut`, the horizon function `Δ0`, the
Hertz-amplitude prefactors `calA,calB`, and the mode frequency `ω=mΩ` with `aω`. Consumed by
Jumps/Closure/KerrRadial. Type-generic `{T}` (BigFloat-inherited).
"""
struct OrbitBackground{T<:Real}
    a::T; r0::T; m::Int
    E::T; L::T; Ω::T; ut::T
    Δ0::T; calA::T; calB::T
    ω::T; aω::T
end

"""
    orbit_background(a, r0, m) -> OrbitBackground

Circular-equatorial orbit scalars at radius `r0` for spin `a`, azimuthal `m`. Element type inherits
from `promote(a,r0)` (BigFloat in → BigFloat out). At a=0 the standard Schwarzschild forms result
(v=1/√r0, Ω=v³, u^t=1/√(1−3/r0), …).
"""
function orbit_background(a::Real, r0::Real, m::Integer)
    T = typeof(float(promote(a, r0)[1]))
    aT = T(a); r0T = T(r0)
    v = one(T)/sqrt(r0T)
    den = sqrt(1 - 3v^2 + 2aT*v^3)
    E = (1 - 2v^2 + aT*v^3)/den
    L = r0T*v*(1 - 2aT*v^3 + aT^2*v^4)/den
    Ω = one(T)/(sqrt(r0T^3) + aT)
    gi = kerr_metric_inv_eq(aT, r0T)
    ut = -E*gi[1,1] + L*gi[1,4]
    Δ0 = r0T^2 - 2r0T + aT^2
    calA = (E*(r0T^2 + aT^2) - aT*L)/Δ0
    calB = L - aT*E
    ω = m*Ω
    aω = aT == 0 ? zero(T) : aT*ω
    OrbitBackground{T}(aT, r0T, Int(m), E, L, Ω, ut, Δ0, calA, calB, ω, aω)
end
