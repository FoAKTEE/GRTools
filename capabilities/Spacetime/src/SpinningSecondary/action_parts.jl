# Linear-in-spin ACTION parts + the J(action) / pex(turning-point) gauge transforms.
# Faithful 1:1 of SpinningTrajectory.wl `J1CFun` (L75-87), `dΩJdCFun` (L94-187: the analytic
# ∂J/∂C and ∂Ω/∂C Jacobian rows), `dpexdCFun` (L189-210) and the `J1pex` / `Ω1J` transforms
# (L236) — the machinery behind arXiv:2506.20726 fig:deltaJr and the J-gauge curve of
# arXiv:2603.13482 fig:delta_Omega.
#
# `C = (En, Lz, K)` with `K = Q + (Lz − aEn)²` the squared Carter–Killing constant; `J = (Jr, Jθ, Jφ)`
# the action variables; `pex = (p, e, x)` the turning-point parametrisation. M = 1, bound orbits.
# Elliptic layer = NumericsCore's complex Carlson `ellipK/ellipE/ellipPi/RD` in the Mathematica
# parameter convention `m = k²` (identical to the convention the legacy code used), coerced with
# `real` exactly where the legacy `_EK/_EE/_EP` wrappers did.

@inline _EK(m)    = real(ellipK(m))
@inline _EE(m)    = real(ellipE(m))
@inline _EP(n, m) = real(ellipPi(n, m))

# real radial roots of a bound orbit (r3, r4 come back as a Complex carrier)
@inline function _real_radial_roots(a, p, e, x, En, Q)
    r1, r2, r3, r4 = radial_roots(a, p, e, x, En, Q)
    return real(r1), real(r2), real(r3), real(r4)
end

"""
    spin_linear_actions(a, p, e, x) -> (δJr, δJθ, δJφ)

Linear-in-spin correction (coefficient of `s∥`) to the action variables at FIXED constants of
motion `(En, Lz, K)` — SpinningTrajectory.wl `J1CFun`. `δJφ = 0` identically (Jφ = Lz is held
fixed with the constants). Replaces `SpinningOrbit._linear_actions_CTilde` (`KerrSpinActionParts.J1C`).
"""
function spin_linear_actions(a, p, e, x)
    En, Lz, Q, K = _ctilde_K(a, p, e, x)
    r1, r2, r3, r4 = _real_radial_roots(a, p, e, x, En, Q)
    krsq = ((r1 - r2) * (r3 - r4)) / ((r1 - r3) * (r2 - r4))
    sK = sqrt(K)
    char = ((r1 - r2) * (sK + im * r3)) / ((sK + im * r2) * (r1 - r3))
    imterm = imag(1 / ((r2 - im * sK) * (r3 - im * sK)) * ellipPi(char, krsq))
    Jr1 = 2 / pi * sqrt(K / ((1 - En^2) * (r1 - r3) * (r2 - r4))) *
          ((a^2 * En - a * Lz + En * r3^2) / (K + r3^2) * _EK(krsq) +
           ((r2 - r3) * (-a^2 * En + a * Lz + En * K)) / sK * imterm)
    z1 = sqrt(1 - x^2)
    az2 = sqrt(a^2 * (1 - z1^2) + (Lz^2 + Q) / (1 - En^2))
    mz = (a * z1 / az2)^2
    Jz1 = -2 / pi * 1 / (az2 * sK * sqrt(1 - En^2)) *
          (En * K * _EK(mz) + (En * (a^2 - K) - a * Lz) * _EP((a^2 * z1^2) / K, mz))
    return (Jr1, Jz1, zero(Jr1))
end

"""
    spin_action_jacobian(a, p, e, x) -> ∂(Jr,Jθ,Jφ)/∂(En,Lz,K)

Analytic action Jacobian (3×3 matrix), SpinningTrajectory.wl `dΩJdCFun` ∂J/∂C block. The third
row is `∂Jφ/∂C = (0, 1, 0)` (`Jφ ≡ Lz`). Replaces `SpinningOrbit._dJdC`.
"""
function spin_action_jacobian(a, p, e, x)
    En, Lz, Q, K = _ctilde_K(a, p, e, x)
    r1, r2, r3, r4 = _real_radial_roots(a, p, e, x, En, Q)
    krsq = ((r1 - r2) * (r3 - r4)) / ((r1 - r3) * (r2 - r4))
    Kkr = _EK(krsq); Ekr = _EE(krsq); Piarkr = _EP((r1 - r2) / (r1 - r3), krsq)
    rminus = 1 - sqrt(1 - a^2); rplus = 1 + sqrt(1 - a^2)
    Pirpkr = _EP((r1 - r2) / (r1 - r3) * (r3 - rplus) / (r2 - rplus), krsq)
    Pirmkr = _EP((r1 - r2) / (r1 - r3) * (r3 - rminus) / (r2 - rminus), krsq)
    C = 2 / sqrt((1 - En^2) * (r1 - r3) * (r2 - r4))
    dJrdK = -((C * Kkr) / (2pi))
    dJrdEn = C / pi * (Kkr * ((4 + a^2) * En - a * Lz) + En * (0.5 * (Kkr * (-r1 * r2 + r3 * (4 + r1 + r2 + r3)) + Ekr * (r1 - r3) * (r2 - r4) + Piarkr * (r2 - r3) * (4 + r1 + r2 + r3 + r4)) + 2 / (-rminus + rplus) * (-(((Kkr - (Pirmkr * (r2 - r3)) / (r2 - rminus)) * (-2 * a^2 + (4 - (a * Lz) / En) * rminus)) / (r3 - rminus)) + ((Kkr - (Pirpkr * (r2 - r3)) / (r2 - rplus)) * (-2 * a^2 + (4 - (a * Lz) / En) * rplus)) / (r3 - rplus))))
    dJrdLz = -(C / pi) * (a * En * Kkr + a / (-rminus + rplus) * (-(((Kkr - (Pirmkr * (r2 - r3)) / (r2 - rminus)) * (-a * Lz + 2 * En * rminus)) / (r3 - rminus)) + ((Kkr - (Pirpkr * (r2 - r3)) / (r2 - rplus)) * (-a * Lz + 2 * En * rplus)) / (r3 - rplus)))
    z1 = sqrt(1 - x^2)
    az2 = sqrt(a^2 * (1 - z1^2) + (Lz^2 + Q) / (1 - En^2))
    mz = (a * z1 / az2)^2
    Kkz = _EK(mz); Pikz = _EP(z1^2, mz)
    dJzdK = Kkz / (sqrt(1 - En^2) * pi * az2)
    dJzdEn = 1 / (sqrt(1 - En^2) * pi * az2) * 2 * (-En * az2^2 * _EE(mz) + (a * Lz - En * (a^2 - az2^2)) * Kkz)
    dJzdLz = -((2 * (-a * En * Kkz + Lz * Pikz)) / (sqrt(1 - En^2) * pi * az2))
    z = zero(dJrdK); o = one(dJrdK)
    return [dJrdEn dJrdLz dJrdK; dJzdEn dJzdLz dJzdK; z o z]
end

# ∂(r1, r2, z1²)/∂(En, Lz, K) — implicit differentiation of the turning-point conditions
# R(r1)=R(r2)=0 and Θ(z1)=0 (SpinningTrajectory.wl `dpexdCFun`, first block).
function _dr12z1sq_dC(a, p, e, x)
    En, Lz, Q, _ = _ctilde_K(a, p, e, x)
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
    spin_orbital_params_jacobian(a, p, e, x) -> ∂(p,e,x)/∂(En,Lz,K)

Analytic turning-point Jacobian (3×3), SpinningTrajectory.wl `dpexdCFun`: the chain of
`∂(p,e,x)/∂(r1,r2,z1²)` with the implicit `∂(r1,r2,z1²)/∂C`. Replaces `SpinningOrbit._dpexdC`.
"""
function spin_orbital_params_jacobian(a, p, e, x)
    z = zero(float(p))
    dpexdr12z1sq = [0.5*(-1 + e)^2      0.5*(1 + e)^2      z;
                    ((-1 + e)^2*(1 + e))/(2p)  ((-1 + e)*(1 + e)^2)/(2p)  z;
                    z                   z                  -(1/(2x))]
    return dpexdr12z1sq * _dr12z1sq_dC(a, p, e, x)
end

"""
    spin_frequency_jacobian(a, p, e, x) -> ∂(Ωr,Ωθ,Ωφ)/∂(En,Lz,K)

Analytic Boyer–Lindquist frequency Jacobian (3×3), SpinningTrajectory.wl `dΩJdCFun` ∂Ω/∂C block:
built from the first and SECOND derivatives of the radial and polar actions with respect to
`(En, Lz, K)` (the Mino-frequency factors are `Υr = −1/(2 ∂Jr/∂K)`, `Υθ = 1/(2 ∂Jθ/∂K)`, …).
Replaces `SpinningOrbit._dOmegadC`.
"""
function spin_frequency_jacobian(a, p, e, x)
    M = one(float(p))
    En, Lz, Q, K = _ctilde_K(a, p, e, x)
    r1, r2, r3, r4 = _real_radial_roots(a, p, e, x, En, Q)
    krsq = ((r1 - r2) * (r3 - r4)) / ((r1 - r3) * (r2 - r4))
    Kkr = _EK(krsq); Ekr = _EE(krsq); Piarkr = _EP((r1 - r2) / (r1 - r3), krsq)
    rminus = 1 - sqrt(1 - a^2); rplus = 1 + sqrt(1 - a^2)
    Pirpkr = _EP((r1 - r2) / (r1 - r3) * (r3 - rplus) / (r2 - rplus), krsq)
    Pirmkr = _EP((r1 - r2) / (r1 - r3) * (r3 - rminus) / (r2 - rminus), krsq)
    C = 2 / sqrt((1 - En^2) * (r1 - r3) * (r2 - r4))
    dJrdK = -((C * Kkr) / (2pi))
    dJrdEn = C / pi * (Kkr * ((4 + a^2) * En - a * Lz) + En * (0.5 * (Kkr * (-r1 * r2 + r3 * (4 + r1 + r2 + r3)) + Ekr * (r1 - r3) * (r2 - r4) + Piarkr * (r2 - r3) * (4 + r1 + r2 + r3 + r4)) + 2 / (-rminus + rplus) * (-(((Kkr - (Pirmkr * (r2 - r3)) / (r2 - rminus)) * (-2 * a^2 + (4 - (a * Lz) / En) * rminus)) / (r3 - rminus)) + ((Kkr - (Pirpkr * (r2 - r3)) / (r2 - rplus)) * (-2 * a^2 + (4 - (a * Lz) / En) * rplus)) / (r3 - rplus))))
    dJrdLz = -(C / pi) * (a * En * Kkr + a / (-rminus + rplus) * (-(((Kkr - (Pirmkr * (r2 - r3)) / (r2 - rminus)) * (-a * Lz + 2 * En * rminus)) / (r3 - rminus)) + ((Kkr - (Pirpkr * (r2 - r3)) / (r2 - rplus)) * (-a * Lz + 2 * En * rplus)) / (r3 - rplus)))
    int1 = ((r2 - r4) * Ekr + (-r1 + r4) * Kkr) / ((r1 - r2) * (r1 - r4))
    int2 = ((-r1 + r3) * Ekr + (r2 - r3) * Kkr) / ((r1 - r2) * (r2 - r3))
    int3 = ((r2 - r4) * Ekr + (-r2 + r3) * Kkr) / ((r2 - r3) * (r3 - r4))
    int4 = ((r1 - r3) * Ekr + (-r1 + r4) * Kkr) / ((r3 - r4) * (-r1 + r4))
    d2JrdK2 = C / (4pi * (1 - En^2)) * ((a^2 - 2 * M * r1 + r1^2) / ((r1 - r2) * (r1 - r3) * (r1 - r4)) * int1 - (a^2 - 2 * M * r2 + r2^2) / ((r1 - r2) * (r2 - r3) * (r2 - r4)) * int2 - (a^2 - 2 * M * r3 + r3^2) / ((r1 - r3) * (-r2 + r3) * (r3 - r4)) * int3 - (a^2 - 2 * M * r4 + r4^2) / ((r1 - r4) * (-r2 + r4) * (-r3 + r4)) * int4)
    d2JrdKdEn = -C / (2pi * (1 - En^2)) * (En * Kkr + ((a^2 + r1^2) * (a^2 * En - a * Lz + En * r1^2)) / ((r1 - r2) * (r1 - r3) * (r1 - r4)) * int1 - ((a^2 + r2^2) * (a^2 * En - a * Lz + En * r2^2)) / ((r1 - r2) * (r2 - r3) * (r2 - r4)) * int2 - ((a^2 + r3^2) * (a^2 * En - a * Lz + En * r3^2)) / ((r1 - r3) * (-r2 + r3) * (r3 - r4)) * int3 - ((a^2 + r4^2) * (a^2 * En - a * Lz + En * r4^2)) / ((r1 - r4) * (-r2 + r4) * (-r3 + r4)) * int4)
    d2JrdKdLz = a * C / (2pi * (1 - En^2)) * ((a^2 * En - a * Lz + En * r1^2) / ((r1 - r2) * (r1 - r3) * (r1 - r4)) * int1 - (a^2 * En - a * Lz + En * r2^2) / ((r1 - r2) * (r2 - r3) * (r2 - r4)) * int2 - (a^2 * En - a * Lz + En * r3^2) / ((r1 - r3) * (-r2 + r3) * (r3 - r4)) * int3 - (a^2 * En - a * Lz + En * r4^2) / ((r1 - r4) * (-r2 + r4) * (-r3 + r4)) * int4)
    intr = r3 * Kkr + (r2 - r3) * Piarkr
    intr2 = 0.5 * (r1 * (-r2 + r3) + r3 * (r2 + r3)) * Kkr + 0.5 * (r1 - r3) * (r2 - r4) * Ekr + 0.5 * (r2 - r3) * (r1 + r2 + r3 + r4) * Piarkr
    d2JrdEn2 = C / (pi * (1 - En^2)) * ((4 + (-1 + En^2) * (-a^2 + En^2 * K + 2 * a * En * Lz)) / (-1 + En^2)^2 * Kkr + 2 / (1 - En^2) * intr + intr2 + ((a^2 + r1^2)^2 * (K + r1^2)) / ((r1 - r2) * (r1 - r3) * (r1 - r4)) * int1 - ((a^2 + r2^2)^2 * (K + r2^2)) / ((r1 - r2) * (r2 - r3) * (r2 - r4)) * int2 - ((a^2 + r3^2)^2 * (K + r3^2)) / ((r1 - r3) * (-r2 + r3) * (r3 - r4)) * int3 - ((a^2 + r4^2)^2 * (K + r4^2)) / ((r1 - r4) * (-r2 + r4) * (-r3 + r4)) * int4)
    d2JrdEndLz = C * a / (pi * (1 - En^2)) * (-Kkr - (int1 * (a^2 + r1^2) * (K + r1^2)) / ((r1 - r2) * (r1 - r3) * (r1 - r4)) + (int2 * (a^2 + r2^2) * (K + r2^2)) / ((r1 - r2) * (r2 - r3) * (r2 - r4)) + (int3 * (a^2 + r3^2) * (K + r3^2)) / ((r1 - r3) * (-r2 + r3) * (r3 - r4)) + (int4 * (a^2 + r4^2) * (K + r4^2)) / ((r1 - r4) * (-r2 + r4) * (-r3 + r4)))
    d2JrdLz2 = C * a^2 / (pi * (1 - En^2)) * ((K + r1^2) / ((r1 - r2) * (r1 - r3) * (r1 - r4)) * int1 + (-K - r2^2) / ((r1 - r2) * (r2 - r3) * (r2 - r4)) * int2 + (-K - r3^2) / ((r1 - r3) * (-r2 + r3) * (r3 - r4)) * int3 + (-K - r4^2) / ((r1 - r4) * (-r2 + r4) * (-r3 + r4)) * int4)
    Upr = -1 / (2 * dJrdK)
    dUprdK = d2JrdK2 / (2 * dJrdK^2); dUprdEn = d2JrdKdEn / (2 * dJrdK^2); dUprdLz = d2JrdKdLz / (2 * dJrdK^2)
    Uptr = -dJrdEn / (2 * dJrdK); Upphr = dJrdLz / (2 * dJrdK)
    dUptrdK = -d2JrdKdEn / (2 * dJrdK) + dJrdEn / (2 * dJrdK^2) * d2JrdK2
    dUptrdEn = -d2JrdEn2 / (2 * dJrdK) + dJrdEn / (2 * dJrdK^2) * d2JrdKdEn
    dUptrdLz = -d2JrdEndLz / (2 * dJrdK) + dJrdEn / (2 * dJrdK^2) * d2JrdKdLz
    dUpphrdK = d2JrdKdLz / (2 * dJrdK) - dJrdLz / (2 * dJrdK^2) * d2JrdK2
    dUpphrdEn = d2JrdEndLz / (2 * dJrdK) - dJrdLz / (2 * dJrdK^2) * d2JrdKdEn
    dUpphrdLz = d2JrdLz2 / (2 * dJrdK) - dJrdLz / (2 * dJrdK^2) * d2JrdKdLz
    z1 = sqrt(1 - x^2)
    az2 = sqrt(a^2 * (1 - z1^2) + (Lz^2 + Q) / (1 - En^2))
    kzsq = (a * z1 / az2)^2
    Kkz = _EK(kzsq); Ekz = _EE(kzsq)
    Dkz = 1 / 3 * real(RD(complex(zero(kzsq)), complex(1 - kzsq), complex(one(kzsq))))
    Pikz = _EP(z1^2, kzsq)
    dJzdK = Kkz / (sqrt(1 - En^2) * pi * az2)
    Iz2 = (-(Ekz / (1 - kzsq))); Iz1 = ((-Dkz + Kkz) / (-1 + kzsq))
    d2JzdK2 = -1 / (2pi * az2^3 * (1 - En^2)^1.5 * (az2^2 - a^2 * z1^2)) * ((a^2 - az2^2) * Iz2 + (1 - z1^2) * a^2 * Iz1)
    d2JzdKdEn = 1 / (pi * az2 * (1 - En^2)^1.5) * (En * Kkz - 1 / ((az2^2 - a^2 * z1^2) * az2^2) * ((a^2 - az2^2) * (a * Lz + En * (-a^2 + az2^2)) * Iz2 + a^3 * (1 - z1^2) * (Lz + a * En * (-1 + z1^2)) * Iz1))
    d2JzdKdLz = 1 / (pi * az2^3 * (1 - En^2)^1.5 * (az2^2 - a^2 * z1^2)) * ((a^2 * Lz + a * En * (-a^2 + az2^2)) * Iz2 + (Lz + a * En * (-1 + z1^2)) * a^2 * Iz1)
    d2JzdEn2 = 2 / (pi * az2 * (1 - En^2)^1.5) * (-(K - (-2 * a^2 + az2^2 + a^2 * z1^2)) * Kkz + az2^2 * (Kkz - Ekz) - 1 / ((az2^2 - a^2 * z1^2) * az2^2) * ((a^2 - az2^2)^2 * (K - az2^2) * Iz2 + a^4 * (1 - z1^2)^2 * (K - a^2 * z1^2) * Iz1))
    d2JzdEndLz = 2 * a / (pi * az2 * (1 - En^2)^1.5) * (Kkz + 1 / ((az2^2 - a^2 * z1^2) * az2^2) * ((a^2 - az2^2) * (K - az2^2) * Iz2 + a^2 * (1 - z1^2) * (K - a^2 * z1^2) * Iz1))
    d2JzdLz2 = -2 * a^2 / (pi * az2^3 * (1 - En^2)^1.5 * (az2^2 - a^2 * z1^2)) * ((K - az2^2) * Iz2 + (K - a^2 * z1^2) * Iz1)
    dJzdEn = 1 / (sqrt(1 - En^2) * pi * az2) * 2 * (-En * az2^2 * Ekz + (a * Lz - En * (a^2 - az2^2)) * Kkz)
    dJzdLz = -((2 * (-a * En * Kkz + Lz * Pikz)) / (sqrt(1 - En^2) * pi * az2))
    Upz = 1 / (2 * dJzdK)
    dUpzdK = -d2JzdK2 / (2 * dJzdK^2); dUpzdEn = -d2JzdKdEn / (2 * dJzdK^2); dUpzdLz = -d2JzdKdLz / (2 * dJzdK^2)
    Uptz = dJzdEn / (2 * dJzdK); Upphz = -dJzdLz / (2 * dJzdK)
    dUptzdK = d2JzdKdEn / (2 * dJzdK) - dJzdEn / (2 * dJzdK^2) * d2JzdK2
    dUptzdEn = d2JzdEn2 / (2 * dJzdK) - dJzdEn / (2 * dJzdK^2) * d2JzdKdEn
    dUptzdLz = d2JzdEndLz / (2 * dJzdK) - dJzdEn / (2 * dJzdK^2) * d2JzdKdLz
    dUpphzdK = -d2JzdKdLz / (2 * dJzdK) + dJzdLz / (2 * dJzdK^2) * d2JzdK2
    dUpphzdEn = -d2JzdEndLz / (2 * dJzdK) + dJzdLz / (2 * dJzdK^2) * d2JzdKdEn
    dUpphzdLz = -d2JzdLz2 / (2 * dJzdK) + dJzdLz / (2 * dJzdK^2) * d2JzdKdLz
    Upt = Uptr + Uptz; Upph = Upphr + Upphz
    dUptdK = dUptrdK + dUptzdK; dUptdEn = dUptrdEn + dUptzdEn; dUptdLz = dUptrdLz + dUptzdLz
    dUpphdK = dUpphrdK + dUpphzdK; dUpphdEn = dUpphrdEn + dUpphzdEn; dUpphdLz = dUpphrdLz + dUpphzdLz
    dOrdK = dUprdK / Upt - dUptdK * Upr / Upt^2
    dOrdEn = dUprdEn / Upt - dUptdEn * Upr / Upt^2
    dOrdLz = dUprdLz / Upt - dUptdLz * Upr / Upt^2
    dOzdK = dUpzdK / Upt - dUptdK * Upz / Upt^2
    dOzdEn = dUpzdEn / Upt - dUptdEn * Upz / Upt^2
    dOzdLz = dUpzdLz / Upt - dUptdLz * Upz / Upt^2
    dOphdK = dUpphdK / Upt - dUptdK * Upph / Upt^2
    dOphdEn = dUpphdEn / Upt - dUptdEn * Upph / Upt^2
    dOphdLz = dUpphdLz / Upt - dUptdLz * Upph / Upt^2
    return [dOrdEn dOrdLz dOrdK; dOzdEn dOzdLz dOzdK; dOphdEn dOphdLz dOphdK]
end

"""
    spin_actions_pex(a, p, e, x) -> [δJr, δJθ, δJφ]

Linear-in-spin action correction at FIXED orbital parameters `(p, e, x)` — the published
`δJ_r` of arXiv:2506.20726 (fig:deltaJr, O(e²) near circularity):
`J1pex = J1C + (∂J/∂C)·(C1CTilde − (∂pex/∂C)⁻¹·pex1CTilde)`.
Replaces `SpinningOrbit.linear_actions_pex` (`KerrSpinActionParts.J1pex`).
"""
function spin_actions_pex(a, p, e, x)
    J1C = collect(spin_linear_actions(a, p, e, x))
    C1pex = collect(spin_linear_constants(a, p, e, x)) -
            spin_orbital_params_jacobian(a, p, e, x) \ collect(spin_linear_orbital_params(a, p, e, x))
    return J1C + spin_action_jacobian(a, p, e, x) * C1pex
end

"""
    spin_frequencies_action_gauge(a, p, e, x) -> [δΩr, δΩθ, δΩφ]

Linear-in-spin Boyer–Lindquist frequency correction in the ACTION (J) gauge — the fourth curve of
arXiv:2603.13482 fig:delta_Omega: `Ω1J = Ω1CTilde − (∂Ω/∂C)·(C1CTilde + (∂J/∂C)⁻¹·J1C)`. Diverges
near the separatrix like the constants (C) gauge, while the CTilde gauge stays finite.
Replaces `SpinningOrbit.linear_frequencies_J` (`KerrSpinActionParts.Omega1J`).
"""
function spin_frequencies_action_gauge(a, p, e, x)
    J1C = collect(spin_linear_actions(a, p, e, x))
    return collect(spin_linear_frequencies(a, p, e, x)) -
           spin_frequency_jacobian(a, p, e, x) *
           (collect(spin_linear_constants(a, p, e, x)) + spin_action_jacobian(a, p, e, x) \ J1C)
end

"""
    spin_frequencies_constants_gauge(a, p, e, x) -> [δΩr, δΩθ, δΩφ]

Linear-in-spin Boyer–Lindquist frequency correction in the CONSTANTS (C) gauge — the reference
geodesic keeps `(En, Lz, K)` fixed: `Ω1C = Ω1CTilde − (∂Ω/∂C)·C1CTilde`. Diverges near the
separatrix (Υr → 0 degenerates ∂Ω/∂C), unlike the CTilde gauge. arXiv:2603.13482 fig:delta_Omega.
"""
function spin_frequencies_constants_gauge(a, p, e, x)
    return collect(spin_linear_frequencies(a, p, e, x)) -
           spin_frequency_jacobian(a, p, e, x) * collect(spin_linear_constants(a, p, e, x))
end

"""
    spin_frequencies_turning_point_gauge(a, p, e, x) -> [δΩr, δΩθ, δΩφ]

Linear-in-spin Boyer–Lindquist frequency correction in the TURNING-POINT (p, "pex") gauge — the
reference geodesic keeps `(p, e, x)` fixed: `Ω1p = Ω1CTilde − (∂Ω/∂(p,e,x))·pex1CTilde`, with
`∂Ω/∂(p,e,x) = (∂Ω/∂C)·(∂(p,e,x)/∂C)⁻¹` from the analytic Jacobians (the legacy figure script used a
finite-difference Jacobian). Diverges only mildly (∼(p−p_s)^(−1/2)). arXiv:2603.13482 fig:delta_Omega.
"""
function spin_frequencies_turning_point_gauge(a, p, e, x)
    dOmega_dpex = spin_frequency_jacobian(a, p, e, x) / spin_orbital_params_jacobian(a, p, e, x)
    return collect(spin_linear_frequencies(a, p, e, x)) -
           dOmega_dpex * collect(spin_linear_orbital_params(a, p, e, x))
end
