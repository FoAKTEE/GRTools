# Energy / angular-momentum fluxes from the asymptotic amplitudes Z^∞ (ℐ) / Z^H (ℋ), for all
# spin weights s∈{−2,−1,0,+1,+2} (Weyl/Maxwell/scalar/ψ0). A thin closed-form layer over the
# amplitudes: takes Z as input (decoupled from the source convolution). M=1. The s-specific FluxInf
# and FluxHor normalizations follow the BHPToolkit convention (arXiv:2008.12703 spin-1,
# arXiv:1003.1860 scalar).

"Invariant Wronskian, including the BHPToolkit static-amplitude normalization."
function _teuk_wronskian(Rin::TeukolskyRadialFunction)
    if Rin.method === :Static
        incidence = iszero(Rin.a * Rin.m) ? Rin.amplitudes.I : Rin.amplitudes.Ip
        return -(2 * Rin.l + 1) * incidence
    end
    return 2 * im * Rin.omega * Rin.amplitudes.Incidence
end

"""
    energy_flux(s, l, m, a, omega, lambda, ZI, ZH) -> (I, H)

Energy flux to infinity (I) and through the horizon (H) from Z^∞=ZI, Z^H=ZH. s∈{−2,−1,0,+1,+2}.
The ℐ part carries `ω^(2(1−|s|))`; the ℋ part the s-specific horizon normalization. Returns (0,0)
for ω=0.
"""
function energy_flux(s::Integer, l::Integer, m::Integer, a::Real, omega::Real, lambda::Real, ZI, ZH)
    iszero(omega) && return (I=0.0, H=0.0)
    M = 1.0
    rh = M + sqrt(M^2 - a^2); Omh = a / (2 * M * rh)
    kf = omega - m * Omh; epsf = sqrt(M^2 - a^2) / (4 * M * rh)
    wfac = omega^(2 * (1 - abs(s)))                  # s=−2:ω^−2, s=±1:ω^0, s=0:ω^2
    # --- FluxInf (switch on s) ---
    FI = if s == -2 || s == 0
        abs2(ZI) * wfac / (4 * pi)
    elseif s == -1
        2 * abs2(ZI) * wfac / (4 * pi)
    elseif s == 1
        AbsCSq = (lambda + 2)^2 + 4 * a * omega * (m - a * omega)
        2 * (4 * omega^4 * abs2(ZI)) / AbsCSq * wfac / (4 * pi)
    else  # s == 2
        AbsCSq = (4 + lambda)^2 * (6 + lambda)^2 + 144 * M^2 * omega^2 +
                 8 * a * (4 + lambda) * (-4 + 5 * (6 + lambda)) * omega * (m - a * omega) +
                 48 * a^2 * omega^2 * (2 * (4 + lambda) + 3 * (m - a * omega)^2)
        (16 * omega^8 * abs2(ZI)) / AbsCSq * wfac / (4 * pi)
    end
    # --- FluxHor (switch on s) ---
    twoMrh = 2 * M * rh
    FH = if s == -2
        AbsCSq = ((lambda + 2)^2 + 4 * a * m * omega - 4 * a^2 * omega^2) *
                 (lambda^2 + 36 * m * a * omega - 36 * a^2 * omega^2) +
                 (2 * lambda + 3) * (96 * a^2 * omega^2 - 48 * m * a * omega) +
                 144 * omega^2 * (M^2 - a^2)
        alphah = (256 * twoMrh^5 * kf * (kf^2 + 4 * epsf^2) * (kf^2 + 16 * epsf^2) * omega^3) / AbsCSq
        alphah * abs2(ZH) / (4 * pi * omega^2)
    elseif s == -1
        ph = lambda^2 + 4 * a * omega * (m - a * omega)
        2 * omega * abs2(ZH) * (twoMrh * kf) * 4 * ((twoMrh * kf)^2 + (M^2 - a^2)) / (ph * pi)
    elseif s == 0
        # rh^2 factor (vs arXiv:1003.1860 Eq.55) because ψ = r R
        1 / (2 * pi * rh) * omega * (omega - m * Omh) * abs2(ZH) * rh^2
    elseif s == 1
        2 * (omega * abs2(ZH)) / (32 * pi * kf * rh)
    else  # s == 2
        (omega * abs2(ZH)) / (512 * pi * rh^3 * kf * (kf^2 + 4 * epsf^2))
    end
    return (I=FI, H=FH)
end

"""
    angular_momentum_flux(s, l, m, a, omega, lambda, ZI, ZH) -> (I, H)

Angular-momentum flux = energy flux × m/ω. Returns (0,0) for ω=0.
"""
function angular_momentum_flux(s::Integer, l::Integer, m::Integer, a::Real, omega::Real, lambda::Real, ZI, ZH)
    iszero(omega) && return (I=0.0, H=0.0)
    ef = energy_flux(s, l, m, a, omega, lambda, ZI, ZH)
    return (I=ef.I * m / omega, H=ef.H * m / omega)
end
