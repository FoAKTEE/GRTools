# MST asymptotic amplitudes (Sasaki–Tagoshi / Casals–Ottewill). The renormalized-angular-momentum
# ν + the MST minimal-solution f_n (from the Radial submodule) → the In-side transmission/incidence/
# reflection (Btrans/Binc/Bref) and the Up-side (Ctrans/Cinc/Cref via the direct D1/D2 forms).
# Pure-function tier consumed by the radial-normalization, the fluxes and the point-particle
# convolution. Pochhammer delegates to BigSpecial; Γ/loggamma via SpecialFunctions (BigFloat
# real-negative args routed through the reflection formula so the inner Γ sees a positive base).

import SpecialFunctions

# complex Γ branches (args may be negative real, e.g. K_ν2 at ν→−1−ν); exp(_lg(z)) == _g(z).
# BigFloat note: `SpecialFunctions.{gamma,loggamma}(::Complex{BigFloat})` with a *purely real*
# non-positive argument throws (falls through to the REAL branch), whereas Complex{Float64}
# analytically continues. For the extended-precision (high-ℓ small-ω) amplitude chain we route such
# real-negative Complex{BigFloat} args through the reflection formula so the inner Γ calls see a
# positive-real argument (1−z has Re≥1); this matches the Float64 complex result to machine
# precision modulo the 2πi branch (the residual phase cancels in the loggamma sums of _knu/_d1).
_needs_refl(z) = z isa Complex{BigFloat} && imag(z) == 0 && real(z) <= 0
function _lg(z0)
    z = complex(z0)
    _needs_refl(z) || return SpecialFunctions.loggamma(z)
    return log(complex(big(pi))) - log(sin(big(pi) * z)) - SpecialFunctions.loggamma(1 - z)
end
function _g(z0)
    z = complex(z0)
    _needs_refl(z) || return SpecialFunctions.gamma(complex(z))
    return complex(big(pi)) / (sin(big(pi) * z) * SpecialFunctions.gamma(1 - z))
end

"Sum a coefficient function over n (up from 0, down from −1) until the partial sum is unchanged."
function _coeff_sum(cf, Nmax::Integer)
    up = zero(cf(0)); n = 0           # zero of the coeff type (ComplexF64 or Complex{BigFloat})
    while n <= Nmax; t = cf(n); prev = up; up += t; up == prev && break; n += 1 end
    dn = zero(up); n = -1
    while n >= -Nmax; t = cf(n); prev = dn; dn += t; dn == prev && break; n -= 1 end
    return up + dn
end

"""
    _knu(s, m, a, omega, nuv, lambda, Nmax) -> Complex

Sasaki–Tagoshi `K_ν` (Eq. 165): prefactor × (upward Γ-series / downward Pochhammer-series), with
`f_n` recomputed at the passed `nuv` (K_ν2 calls this at nuv = −1−ν). `loggamma` keeps the up-series
Γ's from overflowing at large n. epsp = (τ+ε)/2.
"""
function _knu(s::Integer, m::Integer, a::Real, omega::Real, nuv, lambda, Nmax::Integer)
    q = a; eps = 2 * omega; kappa = sqrt(1 - q^2); tau = (eps - m * q) / kappa; epsp = (tau + eps) / 2
    farr, _, _, _ = _mst_fn_array(s, m, a, omega, nuv, lambda, Nmax)
    upterm = n -> (iseven(n) ? 1 : -1) * exp(
        _lg(1 + n + s + im * eps + nuv) + _lg(1 + n + 2nuv) + _lg(1 + n + nuv + im * tau)
        - _lg(n + 1) - _lg(1 + n - s - im * eps + nuv) - _lg(1 + n + nuv - im * tau)) * farr[n]
    downterm = n -> (iseven(n) ? 1 : -1) * pochhammer(1 + s - im * eps + nuv, n) /
        (_g(-n + 1) * pochhammer(1 - s + im * eps + nuv, n) * pochhammer(2 + 2nuv, n)) * farr[n]
    sumUp = zero(ComplexF64); n = 0
    while n <= Nmax; t = upterm(n); p = sumUp; sumUp += t; sumUp == p && break; n += 1 end
    sumDown = zero(ComplexF64); n = 0
    while n >= -Nmax; t = downterm(n); p = sumDown; sumDown += t; sumDown == p && break; n -= 1 end
    pref = 2.0^(-nuv) * exp(im * eps * kappa) * complex(eps * kappa)^(s - nuv) *
           _g(1 - s - 2im * epsp) * _g(2 + 2nuv) /
           (_g(1 - s + im * eps + nuv) * _g(1 + s + im * eps + nuv) * _g(1 + nuv + im * tau))
    return pref * sumUp / sumDown
end

"D2 factor for renorm angmom `nuv` (D22 = D2 at nuv=−1−ν). Uses Σ f_n(nuv)."
function _d2(s::Integer, m::Integer, a::Real, omega::Real, nuv, lambda, Nmax::Integer)
    q = a; eps = 2 * omega; kappa = sqrt(1 - q^2); tau = (eps - m * q) / kappa
    farr, _, _, _ = _mst_fn_array(s, m, a, omega, nuv, lambda, Nmax)
    fSum = _coeff_sum(n -> farr[n], Nmax)
    return -exp(im * kappa * (eps + tau) * (1 + kappa + 2 * log(kappa)) / (2 * (1 + kappa))) *
           (2 * kappa)^(2s) *
           (sin(pi * (nuv - im * eps)) * sin(pi * (nuv - im * tau))) /
           (sin(2 * pi * nuv) * sin(pi * im * (eps + tau))) * fSum
end

"D1 factor for renorm angmom `nuv` (D12 = D1 at nuv=−1−ν). Uses Σ (Γ-ratio) f_n(nuv)."
function _d1(s::Integer, m::Integer, a::Real, omega::Real, nuv, lambda, Nmax::Integer)
    q = a; eps = 2 * omega; kappa = sqrt(1 - q^2); tau = (eps - m * q) / kappa
    farr, _, _, _ = _mst_fn_array(s, m, a, omega, nuv, lambda, Nmax)
    d1term = n -> exp(_lg(1 + nuv + n + s + im * eps) + _lg(1 + nuv + n + im * tau)
                      - _lg(1 + nuv + n - s - im * eps) - _lg(1 + nuv + n - im * tau)) * farr[n]
    fSumD1 = _coeff_sum(d1term, Nmax)
    return exp(-im * kappa * (eps + tau) * (1 + kappa + 2 * log(kappa)) / (2 * (1 + kappa))) *
           (sin(pi * (nuv + im * eps)) * sin(pi * (nuv + im * tau)) * _g(1 - s - im * (eps + tau))) /
           (sin(2 * pi * nuv) * sin(pi * im * (eps + tau)) * _g(1 + s + im * eps + im * tau)) * fSumD1
end

"""
    _mst_amplitudes(s, m, a, omega, nu, lambda; Nmax=80) -> NamedTuple

All asymptotic amplitudes: In (Btrans, Binc, Bref) + Up (Ctrans, Cinc, Cref) + Aplus/Aminus. The
Up incidence/reflection use the DIRECT Teukolsky D1/D2 forms (not the flux-balance shortcut).
"""
function _mst_amplitudes(s::Integer, m::Integer, a::Real, omega::Real, nu, lambda; Nmax::Integer=80)
    q = a; eps = 2 * omega; kappa = sqrt(1 - q^2); tau = (eps - m * q) / kappa
    farr, _, _, _ = _mst_fn_array(s, m, a, omega, nu, lambda, Nmax)
    a1 = nu + 1 + s - im * eps; a2 = nu + 1 - s + im * eps
    fSum = _coeff_sum(n -> farr[n], Nmax)
    fSumAminus = _coeff_sum(n -> (iseven(n) ? 1 : -1) * pochhammer(a1, n) / pochhammer(a2, n) * farr[n], Nmax)

    prefacInTrans = 4.0^s * kappa^(2s) * exp(im * (eps + tau) * kappa * (0.5 + log(kappa) / (1 + kappa)))
    prefacUpTrans = (eps / 2)^(-1 - 2s) * exp(im * eps * (log(complex(eps)) - (1 - kappa) / 2))
    InTrans = prefacInTrans * fSum                                                  # Btrans
    Aminus = 2.0^(-s - 1 + im * eps) * exp(-pi * eps / 2 - im * pi * (nu + 1 + s) / 2) * fSumAminus
    UpTrans = prefacUpTrans * Aminus                                                # Ctrans

    Knu1 = _knu(s, m, a, omega, nu, lambda, Nmax)
    Knu2 = _knu(s, m, a, omega, -1 - nu, lambda, Nmax)
    prefacAplus = 2.0^(-1 + s - im * eps) * exp(-pi * eps / 2 + im * pi * (1 - s + nu) / 2) *
                  _g(1 - s + im * eps + nu) / _g(1 + s - im * eps + nu)
    Aplus = prefacAplus * fSum
    prefacInInc = (2 / eps) *
        (Knu1 - im * exp(-im * pi * nu) * sin(pi * (nu - s + im * eps)) / sin(pi * (nu + s - im * eps)) * Knu2) *
        exp(-im * eps * (log(complex(eps)) - (1 - kappa) / 2))
    InInc = prefacInInc * Aplus                                                     # Binc
    InRef = UpTrans * (Knu1 + im * exp(im * pi * nu) * Knu2)                         # Bref

    # Up incidence/reflection — direct Teukolsky D1/D2 forms
    D1 = _d1(s, m, a, omega, nu, lambda, Nmax); D12 = _d1(s, m, a, omega, -1 - nu, lambda, Nmax)
    D2 = _d2(s, m, a, omega, nu, lambda, Nmax); D22 = _d2(s, m, a, omega, -1 - nu, lambda, Nmax)
    outer = exp(-pi * eps - im * pi * s) / sin(2 * pi * nu)
    c1 = exp(-im * pi * nu) * sin(pi * (nu - s + im * eps)) / Knu1
    c2 = im * sin(pi * (nu + s - im * eps)) / Knu2
    UpInc = outer * (c1 * D1 - c2 * D12)                                            # Cinc
    UpRef = outer * (c1 * D2 - c2 * D22)                                            # Cref

    return (InTrans=InTrans, InInc=InInc, InRef=InRef,
            UpTrans=UpTrans, UpInc=UpInc, UpRef=UpRef, Aplus=Aplus, Aminus=Aminus)
end
