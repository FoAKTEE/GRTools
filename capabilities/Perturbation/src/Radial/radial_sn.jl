# SasakiNakamura radial backend (:SasakiNakamura, s=-2 Up only). Spec: BHPToolkit
# Kernel/SasakiNakamura.m. Transforms the Teukolsky equation to the regular SN variable X,
# integrates the SN ODE inward from rout=1000 (asymptotic OuterBCs recurrence) to the target,
# and back-transforms X→R^up via the χ relation, then scale-matches the raw solution to the
# already-green :MST unit-transmission R^up at a deep-field anchor. Only s=-2 Up is implemented
# upstream (SasakiNakamura.m l.121 prints "Only s=-2 is currently implemented"); any other s/BC
# throws `NotYetPorted`. M=1.
#
# FIREWALL: the inward integration is driven by the VENDORED, in-firewall `NumericsCore.dp5_step!`
# (never OrdinaryDiffEq). The complex SN state (X, X') is carried as a real 4-vector
# [Re X, Im X, Re X', Im X'] so the real-typed DP5 primitive applies; the dense mesh + Hermite
# evaluation is preserved. This is a Float64 numerical backend (M=1, s=-2).

"Faithful-restriction marker: the SN backend implements only s=-2 Up upstream; other s / the In BC
are not ported. Matches the legacy `Teukolsky.NotYetPorted` boundary."
struct NotYetPorted <: Exception
    what::Symbol
end
Base.showerror(io::IO, e::NotYetPorted) =
    print(io, "NotYetPorted: SasakiNakamura backend supports only s=-2 Up (got $(e.what)).")

# ---------------------------------------------------------------------------------------------
# SN equation coefficients and their analytic r-derivatives (SasakiNakamura.m l.74-100).
# The MMA code builds α,β,η,V,U1,G,F,U with symbolic D[·,r]; here every derivative is in closed
# form. Constants c0..c4 are r-independent.
# ---------------------------------------------------------------------------------------------

"""
    _sn_coeffs(lambda, m, a, omega, r) -> NamedTuple

All SN ODE building blocks at radius r (M=1, s=-2), with analytic r-derivatives:
returns η, η', Δ, K, K', V, β, β', α, α', U1, G, G', F, U, f, f'. Mirrors
SasakiNakamuraEquation (SasakiNakamura.m l.74-100). R3.
"""
function _sn_coeffs(lambda, m, a, omega, r)
    am = a * omega - m
    # c0..c4 (constants in r)
    c0 = -12im * omega + lambda * (lambda + 2) - 12 * a * omega * am
    c1 = 8im * a * (3 * a * omega - lambda * am)
    c2 = -24im * a * (am) + 12 * a^2 * (1 - 2 * am^2)
    c3 = 24im * a^3 * am - 24 * a^2
    c4 = 12 * a^4

    Δ  = r^2 - 2r + a^2
    dΔ = 2r - 2
    r2a2 = r^2 + a^2
    f  = Δ / r2a2
    # f' = (dΔ (r²+a²) - Δ·2r)/(r²+a²)²
    df = (dΔ * r2a2 - Δ * 2r) / r2a2^2

    η  = c0 + c1/r + c2/r^2 + c3/r^3 + c4/r^4
    dη = -c1/r^2 - 2c2/r^3 - 3c3/r^4 - 4c4/r^5

    K  = r2a2 * omega - m * a
    dK = 2r * omega                          # D[K,r]

    # V = -((K²+4i(r-M)K)/Δ) + 8 i ω r + λ
    V = -((K^2 + 4im * (r - 1) * K) / Δ) + 8im * omega * r + lambda

    # β = 2Δ(-iK + r - M - 2Δ/r)
    g = -im * K + r - 1 - 2Δ/r               # inner factor
    dg = -im * dK + 1 - 2 * (dΔ * r - Δ) / r^2   # d/dr[-iK + r -1 -2Δ/r]
    β  = 2Δ * g
    dβ = 2 * (dΔ * g + Δ * dg)               # β' = 2(Δ' g + Δ g')

    # α = -iKβ/Δ² + 3i K' + λ + 6Δ/r²
    α = -im * K * β / Δ^2 + 3im * dK + lambda + 6Δ/r^2
    # α' : term by term.
    #  d/dr[-iK β/Δ²] = -i[ (K'β + K β')/Δ² - K β·2Δ'/Δ³ ]
    t1 = -im * ((dK * β + K * dβ) / Δ^2 - K * β * 2 * dΔ / Δ^3)
    #  d/dr[3i K'] = 3i K''  with K''=2ω
    t2 = 3im * (2 * omega)
    #  d/dr[6Δ/r²] = 6(Δ' r² - Δ·2r)/r⁴
    t3 = 6 * (dΔ * r^2 - Δ * 2r) / r^4
    dα = t1 + t2 + t3

    # Let P ≡ 2α + β'/Δ.  U1 = V + Δ²/β ( P' - (η'/η)(α + β'/Δ) ).
    # need β'' for P' :  β = 2Δ g  ⇒ β'' = 2(Δ'' g + 2Δ' g' + Δ g''), Δ''=2.
    #  g'' = d/dr[-i dK + 1 - 2Δ'/r + 2Δ/r²]  (rewrote dg = -i dK +1 -2(Δ' r-Δ)/r² = -i dK+1-2Δ'/r+2Δ/r²)
    #   d/dr[-i dK] = -i K'' = -i·2ω
    #   d/dr[1] = 0
    #   d/dr[-2Δ'/r] = -2(Δ'' r - Δ')/r² = -2(2r - dΔ)/r²
    #   d/dr[2Δ/r²] = 2(Δ' r² - Δ·2r)/r⁴ = 2(dΔ r - 2Δ)/r³
    d2g = -im * (2 * omega) - 2 * (2r - dΔ) / r^2 + 2 * (dΔ * r - 2Δ) / r^3
    d2β = 2 * (2 * g + 2 * dΔ * dg + Δ * d2g)
    # P = 2α + β'/Δ ; P' = 2α' + (β''Δ - β'Δ')/Δ²
    dP = 2 * dα + (d2β * Δ - dβ * dΔ) / Δ^2
    AB = α + dβ / Δ                          # (α + β'/Δ)
    U1 = V + (Δ^2 / β) * (dP - (dη / η) * AB)

    # G = -(2(r-M)/(r²+a²)) + rΔ/(r²+a²)²
    G = -(2 * (r - 1) / r2a2) + r * Δ / r2a2^2
    # G' :  d/dr[-2(r-1)/(r²+a²)] = -2[(r²+a²) - (r-1)2r]/(r²+a²)²
    g1 = -2 * (r2a2 - (r - 1) * 2r) / r2a2^2
    #       d/dr[rΔ/(r²+a²)²] = [(Δ + rΔ')(r²+a²)² - rΔ·2(r²+a²)2r]/(r²+a²)⁴
    #                        = [(Δ + rΔ')(r²+a²) - 4r²Δ]/(r²+a²)³
    g2 = ((Δ + r * dΔ) * r2a2 - 4 * r^2 * Δ) / r2a2^3
    dG = g1 + g2

    # F = (η'/η) Δ/(r²+a²)  = (η'/η) f
    F = (dη / η) * f

    # U = Δ U1/(r²+a²)² + G² + Δ G'/(r²+a²) - F G
    U = Δ * U1 / r2a2^2 + G^2 + Δ * dG / r2a2 - F * G

    return (η=η, dη=dη, Δ=Δ, dΔ=dΔ, K=K, dK=dK, β=β, dβ=dβ, d2β=d2β, α=α, dα=dα,
            f=f, df=df, F=F, U=U, r2a2=r2a2)
end

"""
    _sn_rhs(lambda, m, a, omega, r, X, dX) -> d2X

Second r-derivative X''(r) of the SN equation (SasakiNakamura.m l.98):
f² X'' + f(f' - F) X' - U X = 0  ⇒  X'' = (U X - f(f'-F) X') / f². M=1, s=-2. R3.
"""
@inline function _sn_rhs(lambda, m, a, omega, r, X, dX)
    co = _sn_coeffs(lambda, m, a, omega, r)
    return (co.U * X - co.f * (co.df - co.F) * dX) / co.f^2
end

# ---------------------------------------------------------------------------------------------
# Outer boundary conditions: asymptotic series X(r) = e^{iω r*} Σ_{j=0}^{20} c_j (ω r)^{-j}.
# The 13-term recurrence for c_k is OuterBCs (SasakiNakamura.m l.35-43), Gralla-Porfyriadis-
# Warburton (arXiv:1506.08496). c_j=0 for j<0, c_0=1; c_k from c_{k-13..k-1}. M=1.
# ---------------------------------------------------------------------------------------------

"""
    _sn_outer_c(lambda, m, a, omega, kmax=20) -> Vector{ComplexF64}

Asymptotic-series coefficients c_0..c_kmax for the SN Up outer BC (SasakiNakamura.m OuterBCs
recurrence, M=1). c[1]=c_0=1. R3.
"""
function _sn_outer_c(lambda, m, a, omega; kmax::Int=20)
    M = 1.0
    c = zeros(ComplexF64, kmax + 14)         # index shift: cidx(j) = j + 14 (so j=-13 → 1)
    cidx(j) = j + 14
    for j in -13:-1
        c[cidx(j)] = 0
    end
    c[cidx(0)] = 1
    den = 2*lambda + lambda^2 - 12*omega*(-a*m + im*M + a^2*omega)
    for k in 1:kmax
        cm(j) = c[cidx(k + j)]                # c[-n+k] in MMA = c at index k+(-n)
        term =
        -((6im*a^12*(210 - 29k + k^2)*omega^12*cm(-13))/(k*den)) +
        (12*a^10*(-14 + k)*omega^11*(-a*(-12 + k)*m + 3im*(-13 + k)*M + a^2*(-13 + k)*omega)*cm(-12))/(k*den) +
        (6im*a^8*omega^10*(-2im*a*(729 - 121k + 5k^2)*m*M - 12*(156 - 25k + k^2)*M^2 - 4*a^3*(143 - 24k + k^2)*m*omega + 2*a^4*(156 - 25k + k^2)*omega^2 + a^2*(k^2*(-5 + 2*m^2 + 10im*M*omega) + k*(123 - 46*m^2 - 250im*M*omega) + 3*(-251 + 87*m^2 + 520im*M*omega)))*cm(-11))/(k*den) +
        (4*a^6*omega^9*(-6*a*(493 - 89k + 4k^2)*m*M^2 + 12im*(132 - 23k + k^2)*M^3 - 12*a^5*(-12 + k)*m*omega^2 + 6*a^6*(-12 + k)*omega^3 - 3im*a^2*M*(4*k^2*(-3 + m^2 + 2im*M*omega) - 2k*(-135 + 43*m^2 + 91im*M*omega) + 3*(-503 + 153*m^2 + 344im*M*omega)) + a^4*omega*(-k^2*(-15 + lambda + 12im*M*omega) + k*(-339 + 6*m^2 + 20*lambda + 276im*M*omega) - 3*(-635 + 25*m^2 + 32*lambda + 528im*M*omega)) + a^3*m*(k^2*(-12 + lambda + 24im*M*omega) + k*(264 - 20*lambda - 534im*M*omega) + 3*(-481 + m^2 + 32*lambda + 986im*M*omega)))*cm(-10))/(k*den) -
        (im*a^5*omega^8*(96im*(-11 + k)*(-9 + k)*m*M^3 + 48*a^6*m*omega^3 + 48*a*M^2*(3 - 2*m^2 + (-9 + k)^2*(9 - 2*m^2 - 2im*M*omega) + (-9 + k)*(-21 + 4*m^2 + 4im*M*omega)) + 4*a^5*omega^2*(-2970 - 27*k^2 - 30*m^2 + 46*lambda + 264im*M*omega + k*(567 - 4*lambda - 24im*M*omega)) - 8im*a^2*m*M*(9*m^2 - 11*lambda + 36im*M*omega + 4*(-9 + k)^2*(-12 + lambda + 6im*M*omega) - 2*(-9 + k)*(-45 + lambda + 30im*M*omega)) + 4*a^4*m*omega*(51*k^2 + k*(-1011 + 4*lambda + 48im*M*omega) + 6*(4*m^2 - 9*(-92 + lambda + 10im*M*omega))) + a^3*((-9 + k)^2*(120 - 96*m^2 + 2*lambda + lambda^2 - 492im*M*omega + 32im*M*lambda*omega - 96*M^2*omega^2) - 4*(-9 + 6*m^4 + 2*lambda^2 + 54im*M*omega + 48*M^2*omega^2 + lambda*(4 + 22im*M*omega) - m^2*(33 + 8*lambda + 78im*M*omega)) + (-9 + k)*(lambda^2 + 96*m^2*(1 - im*M*omega) + 2*lambda*(1 - 8im*M*omega) + 12*(-22 + 95im*M*omega + 24*M^2*omega^2))))*cm(-9))/(2k*den) +
        (a^4*omega^7*(48im*(159 - 36k + 2k^2)*M^3 + 8*a*m*M^2*(k*(525 - 34*lambda) + 2*k^2*(-15 + lambda) + 3*(-751 + 47*lambda)) - 4*a^5*m*(51k - 2*(249 + lambda))*omega^2 + 12*a^6*(-88 + 9k)*omega^3 + 4*a^3*m*(-12 - 15*lambda - lambda^2 + m^2*(12 + lambda) + 6im*M*omega - 18im*M*lambda*omega + 2*(-8 + k)^2*(-9 + 2*lambda + 39im*M*omega) + (-8 + k)*(36 + im*M*(-141 + 4*lambda)*omega)) + a^4*omega*(-8*(-8 + k)^2*(-15 + 2*lambda + 21im*M*omega) + 4*(9 + lambda^2 - 3*m^2*(19 + lambda) - 54im*M*omega + 2*lambda*(6 + 5im*M*omega)) + (-8 + k)*(-264 + 96*m^2 - lambda^2 + 456im*M*omega + lambda*(-2 - 16im*M*omega))) - im*a^2*M*(2*(-8 + k)^2*(-108 + 72*m^2 - lambda^2 + 156im*M*omega + lambda*(-2 - 8im*M*omega)) - (-8 + k)*(-432 + 168*m^2 + lambda^2 + 540im*M*omega + 2*lambda*(1 - 8im*M*omega)) - 4*(m^2*(30 + 8*lambda) - 3*(lambda^2 - 18im*M*omega + lambda*(2 + 2im*M*omega)))))*cm(-8))/(k*den) +
        (im*a^3*omega^6*(-192im*(54 - 15k + k^2)*m*M^3 - 168*a^6*m*omega^3 + 4*a^5*omega^2*(3459 + 48*k^2 + 93*m^2 - 149*lambda - 732im*M*omega + 4k*(-204 + 4*lambda + 21im*M*omega)) - 2*a^4*m*omega*(10086 + 168*k^2 + 126*m^2 - 348*lambda - lambda^2 - 2772im*M*omega + 8k*(-327 + 4*lambda + 39im*M*omega)) + 8im*a^2*m*M*(3 + 18*m^2 - 37*lambda - 2*lambda^2 - 24im*M*omega + 6*(-7 + k)^2*(-9 + 2*lambda + 9im*M*omega) - 2*(-7 + k)*(-42 + lambda + 42im*M*omega)) + 4*a*M^2*(-12*(-7 + k)*(-15 + 4*m^2 + 4im*M*omega) + 4*(9 - 12*m^2 + 2*lambda + lambda^2 - 36im*M*omega) + (-7 + k)^2*(-108 + 48*m^2 - 2*lambda - lambda^2 + 60im*M*omega)) + a^3*(12 + 48*m^4 + 54*lambda + 29*lambda^2 + lambda^3 - 420im*M*omega + 236im*M*lambda*omega + 16im*M*lambda^2*omega - m^2*(156 + 98*lambda + lambda^2 + 756im*M*omega) + 4*(-7 + k)^2*(-30 + 36*m^2 - 2*lambda - lambda^2 + 192im*M*omega - 24im*M*lambda*omega + 60*M^2*omega^2) + 4*(-7 + k)*(lambda^2*(-1 - im*M*omega) + lambda*(-2 + 2im*M*omega) + m^2*(-36 + 72im*M*omega) - 6*(-9 + 52im*M*omega + 22*M^2*omega^2))))*cm(-7))/(2k*den) +
        (a^2*omega^5*(48im*(-7 + k)*(-6 + k)*M^3 + 8*a*m*M^2*(21 + 21*(-6 + k) + 4*(-6 + k)^2*(-6 + lambda) - 12*lambda) + 24*a^5*m*(105 - 14k + lambda)*omega^2 + 96*a^6*(-15 + 2k)*omega^3 + 4*a^3*m*(-15 - 20*lambda - 3*lambda^2 + m^2*(15 + 2*lambda) - 9im*M*omega - 43im*M*lambda*omega + 6*(-6 + k)^2*(-2 + lambda + 15im*M*omega) + 3*(-6 + k)*(8 + im*M*(-39 + 4*lambda)*omega)) + 4*a^4*omega*(-3 + 11*lambda + 3*lambda^2 - m^2*(81 + 8*lambda) - 57im*M*omega + 27im*M*lambda*omega - 6*(-6 + k)^2*(-5 + lambda + 9im*M*omega) + (-6 + k)*(-54 + 36*m^2 - lambda^2 + 135im*M*omega + lambda*(-2 - 12im*M*omega))) - im*a^2*M*(36 + 66*lambda + 35*lambda^2 + lambda^3 - 4*m^2*(27 + 16*lambda) - 636im*M*omega + 84im*M*lambda*omega - (-6 + k)*(120*m^2 + 10*lambda + 5*lambda^2 + 36*(-6 + 7im*M*omega)) + 2*(-6 + k)^2*(-72 + 72*m^2 - 3*lambda^2 + 180im*M*omega + lambda*(-6 - 16im*M*omega))))*cm(-6))/(k*den) -
        (im*a*omega^4*(96im*(-1 + (-5 + k)^2)*m*M^3 + 216*a^6*m*omega^3 - 8im*a^2*m*M*(-6 + 9*m^2 - 33*lambda - 4*lambda^2 + 12*(-5 + k)^2*(-2 + lambda + 3im*M*omega) + 2*(-5 + k)*(15 + lambda - 12im*M*omega)) + 12*a^5*omega^2*(-590 - 14*k^2 - 34*m^2 + 57*lambda + 228im*M*omega - 2k*(-91 + 4*lambda + 18im*M*omega)) + 8*a*M^2*(-10*lambda - 5*lambda^2 + 72im*M*omega + (-5 + k)*(-18 + 2*lambda + lambda^2 - 12im*M*omega) + (-5 + k)^2*(18 - 12*m^2 + 2*lambda + lambda^2 - 24im*M*omega)) + 6*a^4*m*omega*(1360 + 44*k^2 + 36*m^2 - 132*lambda - lambda^2 - 772im*M*omega + 4k*(-123 + 4*lambda + 30im*M*omega)) + a^3*(-12 - 24*m^4 - 64*lambda - 38*lambda^2 - 3*lambda^3 + 552im*M*omega - 132im*M*lambda*omega - 32im*M*lambda^2*omega + 2*m^2*(30 + 50*lambda + lambda^2 + 288im*M*omega) - 6*(-5 + k)^2*(-10 + 16*m^2 - lambda^2 + 92im*M*omega + 32*M^2*omega^2 + lambda*(-2 - 16im*M*omega)) + 2*(-5 + k)*(48*m^2*(1 - 3im*M*omega) + lambda^2*(3 + 6im*M*omega) + lambda*(6 + 20im*M*omega) + 42*(-1 + 6im*M*omega + 4*M^2*omega^2))))*cm(-5))/(2k*den) +
        (a*omega^3*(8*m*M^2*(k*(51 - 14*lambda) + 21*(-5 + lambda) + 2*k^2*(-3 + lambda)) + 24*a^4*m*(57 - 11k + lambda)*omega^2 + 24*a^5*(-36 + 7k)*omega^3 + 2*a^3*omega*(k^2*(30 - 8*lambda - 60im*M*omega) + k*(-282 + 48*m^2 + 58*lambda - 3*lambda^2 + 630im*M*omega - 24im*M*lambda*omega) - 2*(-321 + 50*lambda - 9*lambda^2 + 7*m^2*(21 + lambda) + 816im*M*omega - 72im*M*lambda*omega)) + 4*a^2*m*(-6 - 11*lambda - 3*lambda^2 + m^2*(6 + lambda) - 32im*M*lambda*omega + (-4 + k)^2*(-3 + 4*lambda + 42im*M*omega) + 3*(-4 + k)*(2 + im*M*(-9 + 4*lambda)*omega)) - im*a*M*(2*((26 - 16*m^2)*lambda + 15*lambda^2 + lambda^3 - 168im*M*omega) + 2*(-4 + k)^2*(-18 + 24*m^2 - 3*lambda^2 + 84im*M*omega + lambda*(-6 - 8im*M*omega)) - (-4 + k)*(-36 + 24*m^2 + 7*lambda^2 - 60im*M*omega + 2*lambda*(7 + 8im*M*omega))))*cm(-4))/(k*den) -
        (im*omega^2*(120*a^5*m*omega^3 + 4*(-4 + k)*k*M^2*(2*lambda + lambda^2 - 12im*M*omega) - 8im*a*m*M*(-7*lambda - 2*lambda^2 + (-3 + k)*(3 + 2*lambda) + 6im*M*omega + (-3 + k)^2*(-3 + 4*lambda + 6im*M*omega)) + 4*a^4*omega^2*(-363 - 18*k^2 - 45*m^2 + 79*lambda + 228im*M*omega - 2k*(-81 + 8*lambda + 30im*M*omega)) + 2*a^3*m*omega*(48*k^2 + 30*m^2 + 8k*(-39 + 4*lambda + 21im*M*omega) - 3*(-162 + 60*lambda + lambda^2 + 212im*M*omega)) + a^2*(-3*lambda^3 + 12im*M*omega*(15 + 11*m^2 + 4im*M*omega) + lambda^2*(-21 + m^2 - 16im*M*omega) + lambda*(-30 + 34*m^2 + 28im*M*omega) + 4*(-3 + k)*(-3 + 6im*M*omega + 36*M^2*omega^2 + lambda^2*(1 + 3im*M*omega) + 6*m^2*(1 - 4im*M*omega) + 2*lambda*(1 + 5im*M*omega)) - 4*(-3 + k)^2*(-3 + 6*m^2 - lambda^2 + 42im*M*omega + 12*M^2*omega^2 + lambda*(-2 - 8im*M*omega))))*cm(-3))/(2k*den) -
        (omega*(8*a^3*m*(-33 + 12k - lambda)*omega^2 - 24*a^4*(-8 + 3k)*omega^3 + M*(3 + 5k - 2k^2 + lambda)*(2im*lambda + im*lambda^2 + 12*M*omega) + 4*a^2*omega*(-18 + lambda - 3*lambda^2 + 2*m^2*(12 + lambda) + 75im*M*omega - 15im*M*lambda*omega + k^2*(-3 + lambda + 6im*M*omega) + k*(15 - 6*m^2 - 2*lambda + lambda^2 - 45im*M*omega + 4im*M*lambda*omega)) + 4*a*m*(lambda^2 - 3im*(5 - 7k + 2k^2)*M*omega - lambda*(2 + k^2 - 15im*M*omega + k*(-4 + 4im*M*omega))))*cm(-2))/(k*den) -
        (im*(-lambda^3 + lambda^2*(k^2 - 2*(2 + a*m*omega + 2im*M*omega) + k*(-1 + 4im*M*omega)) + 12*omega*(2*a^3*m*omega^2 + a^2*omega*(-6 + 5k - k^2 - 2*m^2 + 4im*M*omega - 4im*k*M*omega) + a*m*(-2 - k + k^2 - 2im*M*omega + 4im*k*M*omega) + M*(2im + im*k - im*k^2 - 4*M*omega + 4k*M*omega)) + 2*lambda*(-2 + k^2 - 24*a*m*omega + 2im*M*omega + 22*a^2*omega^2 + k*(-1 + 8*a*m*omega + 4im*M*omega - 8*a^2*omega^2)))*cm(-1))/(2k*den)
        c[cidx(k)] = term
    end
    return c[cidx(0):cidx(kmax)]             # length kmax+1, index 1 → c_0
end

"""
    _sn_outer_X(c, omega, a, r) -> (X, dX, d2X)

Asymptotic SN Up solution X(r) and its first two r-derivatives at r (SasakiNakamura.m l.48-51):
X = e^{iω r*} Σ c_j (ω r)^{-j}, r* = r + 2rp/(rp-rm) log((r-rp)/2) - 2rm/(rp-rm) log((r-rm)/2)
(M=1). R3.
"""
function _sn_outer_X(c::Vector{<:Complex}, omega, a, r)
    M = 1.0
    rp = M + sqrt(M^2 - a^2)
    rm = M - sqrt(M^2 - a^2)
    # tortoise r* and its derivatives
    rstar = r + (2M*rp)/(rp - rm) * log((r - rp)/(2M)) - (2M*rm)/(rp - rm) * log((r - rm)/(2M))
    drs   = 1 + (2M*rp)/(rp - rm) / (r - rp) - (2M*rm)/(rp - rm) / (r - rm)
    d2rs  = -(2M*rp)/(rp - rm) / (r - rp)^2 + (2M*rm)/(rp - rm) / (r - rm)^2
    # S(r) = Σ_{j} c_j (ω r)^{-j}, S', S''
    S = zero(ComplexF64); dS = zero(ComplexF64); d2S = zero(ComplexF64)
    for j in 0:(length(c) - 1)
        cj = c[j + 1]
        wr = omega * r
        S   += cj * wr^(-j)
        if j >= 1
            dS  += cj * (-j) * omega^(-j) * r^(-j - 1)
            d2S += cj * (-j) * (-j - 1) * omega^(-j) * r^(-j - 2)
        end
    end
    E = exp(im * omega * rstar)
    # X = E S ;  X' = E(iω r*' S + S') ;  X'' = E[(iω r*'' + (iω r*')²) S + 2 iω r*' S' + S'']
    X   = E * S
    dX  = E * (im * omega * drs * S + dS)
    d2X = E * ((im * omega * d2rs + (im * omega * drs)^2) * S + 2 * im * omega * drs * dS + d2S)
    return X, dX, d2X
end

# ---------------------------------------------------------------------------------------------
# Integrate the SN equation inward and back-transform to R^up.
# ---------------------------------------------------------------------------------------------

const _SN_ROUT = 1000.0                       # outer BC radius (SasakiNakamura.m l.109)

"""
    _SNSolution

Immutable dense-output carrier for the s=-2 SN Up solution: stores the adaptive mesh of X and X'
(ascending r). Evaluation Hermite-interpolates; back-transform to R^up uses the χ relation. R3.
"""
struct _SNSolution
    rs::Vector{Float64}
    X::Vector{ComplexF64}
    dX::Vector{ComplexF64}
    d2X::Vector{ComplexF64}
end

"Adaptive inward integration of the SN equation from rout=1000 to rnear via the vendored,
 in-firewall `NumericsCore.dp5_step!`. The complex state (X,X') is carried as a real 4-vector
 [Re X, Im X, Re X', Im X']; the adaptive mesh (X, X', X'') is stored ascending for Hermite eval."
function _sn_integrate(lambda, m, a, omega, rnear; tol::Real=1e-11)
    c = _sn_outer_c(lambda, m, a, omega)
    X0, dX0, _ = _sn_outer_X(c, omega, a, _SN_ROUT)
    # real 4-vector RHS: y = [ReX, ImX, ReX', ImX'];  g = [X', X''] split into re/im
    f!(g, r, y) = begin
        X = complex(y[1], y[2]); dX = complex(y[3], y[4])
        d2X = _sn_rhs(lambda, m, a, omega, r, X, dX)
        g[1] = y[3]; g[2] = y[4]; g[3] = real(d2X); g[4] = imag(d2X)
        return g
    end

    r0 = _SN_ROUT; r1 = float(rnear)
    span = r0 - r1
    h = -min(span, max(span / 200, 1e-2))         # inward (h<0)
    hmin = span * 1e-12 + 1e-13
    tolT = float(tol)

    rcur = r0
    y = Float64[real(X0), imag(X0), real(dX0), imag(dX0)]
    g = similar(y); f!(g, rcur, y); k1 = copy(g)
    rmesh = Float64[rcur]; Xmesh = ComplexF64[X0]; dXmesh = ComplexF64[dX0]
    d2mesh = ComplexF64[complex(g[3], g[4])]

    maxsteps = 5_000_000; nsteps = 0
    while rcur > r1 && nsteps < maxsteps
        nsteps += 1
        rcur + h < r1 && (h = r1 - rcur)          # don't overshoot the inner endpoint
        err, y5, k7 = dp5_step!(f!, rcur, y, k1, h, tolT)
        if err <= 1.0 || abs(h) <= hmin
            rcur += h; y = y5; k1 = k7
            push!(rmesh, rcur)
            push!(Xmesh, complex(y5[1], y5[2])); push!(dXmesh, complex(y5[3], y5[4]))
            push!(d2mesh, complex(k7[3], k7[4]))  # X'' at the new node (FSAL final stage)
            h *= clamp(0.9 * err^(-1 / 5), 0.2, 5.0)
        else
            h *= max(0.2, 0.9 * err^(-1 / 5))
        end
    end
    nsteps >= maxsteps && error("SasakiNakamura: step budget exhausted (rnear=$rnear).")
    reverse!(rmesh); reverse!(Xmesh); reverse!(dXmesh); reverse!(d2mesh)
    return _SNSolution(rmesh, Xmesh, dXmesh, d2mesh)
end

"Hermite eval of X (deriv=0) or X' (deriv=1) on the SN mesh — domain check + the mesh
`NumericsCore.cubic_hermite` (values X with slopes X', or X' with slopes X'')."
function _sn_hermite(sol::_SNSolution, r::Real, deriv::Int)
    rs = sol.rs
    (r >= rs[1] - 1e-9 && r <= rs[end] + 1e-9) ||
        error("SasakiNakamura: r=$r outside integrated domain [$(rs[1]), $(rs[end])].")
    return deriv == 0 ? cubic_hermite(rs, sol.X, sol.dX, r) : cubic_hermite(rs, sol.dX, sol.d2X, r)
end

"""
    _sn_to_R(lambda, m, a, omega, r, X, dX, d2X) -> (R, dR)

Back-transform SN X→Teukolsky R^up at r (s=-2, M=1), TeukolskyRadialUp (SasakiNakamura.m
l.123-150). χ = X Δ/√(r²+a²);  R = (1/η)[(α + β'/Δ)χ - (β/Δ)χ'];  derivative analytic. R3.
"""
function _sn_to_R(lambda, m, a, omega, r, X, dX, d2X)
    co = _sn_coeffs(lambda, m, a, omega, r)
    Δ = co.Δ; dΔ = co.dΔ; η = co.η; dη = co.dη; β = co.β; dβ = co.dβ; α = co.α; dα = co.dα
    r2a2 = co.r2a2
    sq = sqrt(r2a2)
    dsq = r / sq                              # d/dr √(r²+a²)
    d2sq = (sq - r * dsq) / r2a2             # d/dr (r/√) = (√ - r·r/√)/(r²+a²)
    # χ = X Δ / √(r²+a²)
    χ   = X * Δ / sq
    # χ' = (X'Δ + X Δ')/√ - XΔ √'/(r²+a²)  = derivative of XΔ/√
    χp  = (dX * Δ + X * dΔ) / sq - X * Δ * dsq / r2a2
    # χ'' for dR/dr
    #  numerator N = XΔ ;  χ = N/√ ;  χ' = N'/√ - N√'/(r²+a²) ... use full product rule
    N   = X * Δ
    Np  = dX * Δ + X * dΔ
    Npp = d2X * Δ + 2 * dX * dΔ + X * 2      # Δ''=2
    # χ = N sq^{-1}; χ' = N' sq^{-1} - N sq^{-2} sq'; χ'' = N'' sq^{-1} - 2N' sq^{-2} sq' + N(2 sq^{-3} sq'^2 - sq^{-2} sq'')
    χpp = Npp/sq - 2*Np*dsq/r2a2 + N*(2*dsq^2/sq^3 - d2sq/sq^2)

    AB = α + dβ / Δ                          # (α + β'/Δ)
    # R = (1/η)[ AB χ - (β/Δ) χ' ]
    R = (1/η) * (AB * χ - (β/Δ) * χp)
    # dR/dr = (-η'/η²)[...] + (1/η) d/dr[ AB χ - (β/Δ)χ' ]
    dAB = dα + co.d2β / Δ - dβ * dΔ / Δ^2    # AB' = α' + (β''/Δ - β'Δ'/Δ²)
    bracket = AB * χ - (β/Δ) * χp
    dbracket = dAB * χ + AB * χp -
               ((dβ * Δ - β * dΔ)/Δ^2) * χp - (β/Δ) * χpp
    dR = (-dη/η^2) * bracket + (1/η) * dbracket
    return R, dR
end

# ---------------------------------------------------------------------------------------------
# Carrier integration + caching for the :SasakiNakamura Up backend (s=-2 only).
# ---------------------------------------------------------------------------------------------

"""
    _SNCache

Lazy cache of the integrated SN solution plus the raw-to-unit-transmission `scale`, held in a
carrier's `cache` field. `scale` matches the raw SN R^up (X→R back-transform, asymptotically
e^{iω r*}) to the carrier's unit-transmission normalization by comparing against the ported MST
R^up at a deep-field anchor. R3.
"""
mutable struct _SNCache
    sol::_SNSolution
    scale::ComplexF64
    have_scale::Bool
end

const _SN_LOCK = ReentrantLock()

# Anchor radius where the raw SN R^up is matched to the MST unit-transmission R^up to fix the
# overall complex normalization. Both are well-conditioned there (Up decays outward but is finite).
const _SN_ANCHOR = 50.0

"""
    _sn_value(R, r, deriv) -> Complex

s=-2 Up SasakiNakamura radial value (`deriv=0`) / first derivative (`deriv=1`) at `r`, normalized
to unit transmission (matched to the MST Up backend). Lazily integrates the SN equation inward,
caches the dense mesh + scale in the carrier's `cache` field, Hermite-interpolates X, and
back-transforms. Thread-safe via a process lock guarding the cache field. R3.
"""
function _sn_value(R::TeukolskyRadialFunction, r::Real, deriv::Int)
    R.s == -2 || throw(NotYetPorted(:SasakiNakamura_s))      # only s=-2 implemented upstream
    R.bc === :Up || throw(NotYetPorted(:SasakiNakamura_In))  # only Up implemented upstream
    s, m, a, omega, nu, lambda, norm = R.s, R.m, R.a, R.omega, R.nu, R.lambda, R.norm
    need = float(r)
    rnear = min(need, _SN_ANCHOR)
    lock(_SN_LOCK)
    try
        cache = R.cache
        if !(cache isa _SNCache) || need < cache.sol.rs[1] - 1e-9 || need > cache.sol.rs[end] + 1e-9
            sol = _sn_integrate(lambda, m, a, omega, rnear)
            cache = _SNCache(sol, one(ComplexF64), false)
            R.cache = cache
        end
        if !cache.have_scale
            # raw SN R^up at the anchor, then match to MST unit-transmission R^up there.
            Xa  = _sn_hermite(cache.sol, _SN_ANCHOR, 0)
            dXa = _sn_hermite(cache.sol, _SN_ANCHOR, 1)
            d2Xa = _sn_rhs(lambda, m, a, omega, _SN_ANCHOR, Xa, dXa)
            Rraw, _ = _sn_to_R(lambda, m, a, omega, _SN_ANCHOR, Xa, dXa, d2Xa)
            Rmst = _mst_radial_up(s, m, a, omega, nu, lambda, _SN_ANCHOR; norm=norm)
            cache.scale = Rmst / Rraw
            cache.have_scale = true
        end
        X   = _sn_hermite(cache.sol, need, 0)
        dX  = _sn_hermite(cache.sol, need, 1)
        d2X = _sn_rhs(lambda, m, a, omega, need, X, dX)
        Rr, dRr = _sn_to_R(lambda, m, a, omega, need, X, dX, d2X)
        return cache.scale * (deriv == 0 ? Rr : dRr)
    finally
        unlock(_SN_LOCK)
    end
end
