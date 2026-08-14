# Point-particle GHP source convolution → Z^∞/Z^H + the public TeukolskyPointParticleMode*
# entrypoints. The source is the GHP-projected stress-energy of a Kerr geodesic, convolved with the
# normalized radial R^in/R^up (+ Teukolsky-ODE d²R) and the SWSH angular operators. Orbit scalars
# (E, Lz, Υt, the BL/Mino frequencies, the libration trajectories) come from `Spacetime`; the SWSH
# from `Harmonics`. The per-spin GHP projections (s∈{−2,−1,0,+1,+2}) are kept as separate leaves and
# routed by a single dispatch on `s` (the projection FORMULA genuinely differs per spin weight).
# Cross-assignment: Z^∞ ← In-radial, Z^H ← Up-radial. M=1.

"d²R from the radial Teukolsky ODE given R, R'."
function _d2R_ode(R, dR, s::Integer, m::Integer, a::Real, r0::Real, omega::Real, lambda::Real)
    Δ = a^2 - 2r0 + r0^2
    K = -a * m + (a^2 + r0^2) * omega
    return (-(-lambda + 4 * im * r0 * s * omega +
              (-2 * im * (-1 + r0) * s * K + K^2) / Δ) * R - (-2 + 2r0) * (1 + s) * dR) / Δ
end

"composite trapezoid of f over [0,π] with N panels (endpoints included)."
function _trapz_halfperiod(f, N::Integer)
    h = pi / N
    acc = (f(0.0) + f(pi)) / 2
    for i in 1:(N-1)
        acc += f(i * h)
    end
    return acc * h
end

"Richardson-extrapolated trapezoid `2·I(2N)−I(N)` over [0,π] (removes the c/N polar endpoint term)."
_trapz_romberg(f, N::Integer) = 2 * _trapz_halfperiod(f, 2N) - _trapz_halfperiod(f, N)

"""
    _trapz2d_romberg_sep(rblock, θblock, combine, Nr, Nθ)

Richardson-extrapolated tensor trapezoid over [0,π]² for an integrand that SEPARATES its expensive
node work: `rblock(qr)` evaluated once per radial node (the SWSH-free radial data), `θblock(qθ)`
once per polar node (the SWSH data), and the cheap per-point GHP algebra in `combine`. The blocks
are cached on the 2N grids so the four trapezoid rules of the per-dimension Richardson step reuse
the same blocks — SWSH is evaluated 2Nθ+1 times, not ~4 Nr Nθ.
"""
function _trapz2d_romberg_sep(rblock, θblock, combine, Nr::Integer, Nθ::Integer)
    rnodes = [rblock(i * pi / (2Nr)) for i in 0:(2Nr)]
    θnodes = [θblock(j * pi / (2Nθ)) for j in 0:(2Nθ)]
    w(idx, n, s) = (idx == 0 || idx == n * s) ? 0.5 : 1.0
    function trap(nr, nθ)
        sr = (2Nr) ÷ nr; sθ = (2Nθ) ÷ nθ
        acc = zero(ComplexF64)
        for i in 0:nr
            wr = w(i * sr, nr, sr); rb = rnodes[i*sr+1]
            row = zero(ComplexF64)
            for j in 0:nθ
                row += w(j * sθ, nθ, sθ) * combine(rb, θnodes[j*sθ+1])
            end
            acc += wr * row
        end
        return acc * (pi / nr) * (pi / nθ)
    end
    Iθ(nr) = 2 * trap(nr, 2Nθ) - trap(nr, Nθ)
    return 2 * Iθ(2Nr) - Iθ(Nr)
end

# ============================================================================
# CIRCULAR-EQUATORIAL source convolutions (e=0, |x|=1, n=k=0).
# ============================================================================

"Z^∞/Z^H for a circular-equatorial s=−2 point particle. Z^∞=−8π αIn/(W Υt), Z^H=−8π αUp/(W Υt)."
function _convolve_source_circular_s2(l::Integer, m::Integer, a::Real, r0::Real,
        En::Real, Lz::Real, Upst::Real, omega::Real, lambda::Real; Nmax::Integer=80)
    s = -2; theta0 = pi / 2; sin0 = sin(theta0); cos0 = cos(theta0)
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    S0 = Sc(theta0, 0.0)
    dS0 = swsh_derivative(Sc, 1, 0, theta0, 0.0)
    d2S0 = swsh_derivative(Sc, 2, 0, theta0, 0.0)
    L1 = -m / sin0 + a * omega * sin0 + cos0 / sin0
    L2 = -m / sin0 + a * omega * sin0 + 2 * cos0 / sin0
    L2S = dS0 + L2 * S0
    L2p = m * cos0 / sin0^2 + a * omega * cos0 - 2 / sin0^2
    L1Sp = d2S0 + L1 * dS0
    L1L2S = L1Sp + L2p * S0 + L2 * dS0 + L1 * L2 * S0
    rho = -1 / (r0 - im * a * cos0); rhobar = -1 / (r0 + im * a * cos0)
    Sigma = 1 / (rho * rhobar); Delta = r0^2 - 2r0 + a^2; Kt = (r0^2 + a^2) * omega - m * a
    sq2 = sqrt(2)
    Ann0 = -rho^(-2) * rhobar^(-1) * (sq2 * Delta)^(-2) *
           (rho^(-1) * L1L2S + 3 * im * a * sin0 * L1 * S0 + 3 * im * a * cos0 * S0 +
            2 * im * a * sin0 * dS0 - im * a * sin0 * L2 * S0)
    Anmbar0 = rho^(-3) * (sq2 * Delta)^(-1) *
              ((rho + rhobar - im * Kt / Delta) * L2S + (rho - rhobar) * a * sin0 * Kt / Delta * S0)
    Anmbar1 = -rho^(-3) * (sq2 * Delta)^(-1) * (L2S + im * (rho - rhobar) * a * sin0 * S0)
    Ambarmbar0 = (Kt^2 * S0 * rhobar) / (4 * Delta^2 * rho^3) +
                 (im * Kt * S0 * (1 - r0 + Delta * rho) * rhobar) / (2 * Delta^2 * rho^3) +
                 (im * r0 * S0 * rhobar * omega) / (2 * Delta * rho^3)
    Ambarmbar1 = -rho^(-3) * rhobar * S0 / 2 * (im * Kt / Delta - rho)
    Ambarmbar2 = -rho^(-3) * rhobar * S0 / 4
    rcomp = (En * (r0^2 + a^2) - a * Lz) / (2 * Sigma)
    thetacomp = rho * (im * sin0 * (a * En - Lz / sin0^2)) / sq2
    Cnn = rcomp^2; Cnmbar = rcomp * thetacomp; Cmbarmbar = thetacomp^2
    alphafn(R0, dR0, d2R0) =
        (Ann0 * Cnn + Anmbar0 * Cnmbar + Ambarmbar0 * Cmbarmbar) * R0 -
        (Anmbar1 * Cnmbar + Ambarmbar1 * Cmbarmbar) * dR0 + Ambarmbar2 * Cmbarmbar * d2R0
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    Rin = R.In(r0); dRin = teukolsky_radial_deriv(R.In, r0); d2Rin = _d2R_ode(Rin, dRin, s, m, a, r0, omega, lambda)
    Rup = R.Up(r0); dRup = teukolsky_radial_deriv(R.Up, r0); d2Rup = _d2R_ode(Rup, dRup, s, m, a, r0, omega, lambda)
    alphaIn = alphafn(Rin, dRin, d2Rin); alphaUp = alphafn(Rup, dRup, d2Rup)
    W = _teuk_wronskian(R.In)
    Zinf = -8 * pi * alphaIn / W / Upst
    ZH = -8 * pi * alphaUp / W / Upst
    return (Zinf=Zinf, ZH=ZH)
end

"Z^∞/Z^H for a circular-equatorial s=+2 (ψ0) point particle (the Δ²-rescaled projection)."
function _convolve_source_circular_s2p(l::Integer, m::Integer, a::Real, r0::Real,
        En::Real, Lz::Real, Upst::Real, omega::Real, lambda::Real; Nmax::Integer=80)
    s = 2; theta0 = pi / 2; sin0 = sin(theta0); cos0 = cos(theta0)
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    S0 = Sc(theta0, 0.0)
    dS0 = swsh_derivative(Sc, 1, 0, theta0, 0.0)
    d2S0 = swsh_derivative(Sc, 2, 0, theta0, 0.0)
    rho = -1 / (r0 - im * a * cos0); rhobar = -1 / (r0 + im * a * cos0)
    Delta = r0^2 - 2r0 + a^2; Kt = (r0^2 + a^2) * omega - m * a
    sq2 = sqrt(2)
    deltaLd1 = m / sin0 - a * omega * sin0 + cos0 / sin0
    deltaLd2 = m / sin0 - a * omega * sin0 + 2 * cos0 / sin0
    ddeltaLd2 = -m * cos0 / sin0^2 - a * omega * cos0 - 2 / sin0^2
    drhoOverrho = im * a * rho * sin0
    d2rhoOverrho = im * a * rho * (cos0 + 2 * sin0 * drhoOverrho)
    All0 = -(1 / 2) * rho^(-1) * rhobar *
           (d2S0 + (deltaLd1 + deltaLd2 + 2 * drhoOverrho) * dS0 +
            (ddeltaLd2 + deltaLd1 * deltaLd2 - 6 * drhoOverrho^2 + 3 * d2rhoOverrho +
             (3 * deltaLd1 - deltaLd2) * drhoOverrho) * S0)
    Alm0 = (2 / sq2) * rho^(-1) *
           (-(rho + rhobar + im * Kt / Delta) * (dS0 + deltaLd2 * S0) +
            (rho - rhobar) * a * sin0 * Kt / Delta * S0)
    Alm1 = (2 / sq2) * rho^(-1) * ((dS0 + deltaLd2 * S0) + im * (rho - rhobar) * a * sin0 * S0)
    Amm0 = (Kt^2 * S0) / (Delta^2 * rho * rhobar) +
           (2 * im * Kt * S0 * (-1 + r0 - Delta * rho)) / (Delta^2 * rho * rhobar) -
           (2 * im * r0 * S0 * omega) / (Delta * rho * rhobar)
    Amm1 = 2 * rho^(-1) * rhobar^(-1) * S0 * (im * Kt / Delta + rho)
    Amm2 = -rho^(-1) * rhobar^(-1) * S0
    rcomp = (En * (r0^2 + a^2) - a * Lz) / Delta
    thetacomp = -rhobar * (im * sin0 * (a * En - Lz / sin0^2)) / sq2
    Cll = rcomp^2; Clm = rcomp * thetacomp; Cmm = thetacomp^2
    dDelta = 2 * (r0 - 1); d2Delta = 2.0
    alphafn(R0, dR0, d2R0) = begin
        D2R = Delta^2 * R0
        dD2R = Delta^2 * dR0 + 2 * Delta * dDelta * R0
        d2D2R = Delta^2 * d2R0 + 4 * Delta * dDelta * dR0 + (2 * dDelta^2 + 2 * Delta * d2Delta) * R0
        (All0 * Cll + Alm0 * Clm + Amm0 * Cmm) * D2R -
        (Alm1 * Clm + Amm1 * Cmm) * dD2R + Amm2 * Cmm * d2D2R
    end
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    Rin = R.In(r0); dRin = teukolsky_radial_deriv(R.In, r0); d2Rin = _d2R_ode(Rin, dRin, s, m, a, r0, omega, lambda)
    Rup = R.Up(r0); dRup = teukolsky_radial_deriv(R.Up, r0); d2Rup = _d2R_ode(Rup, dRup, s, m, a, r0, omega, lambda)
    alphaIn = alphafn(Rin, dRin, d2Rin); alphaUp = alphafn(Rup, dRup, d2Rup)
    W = _teuk_wronskian(R.In)
    Zinf = -8 * pi * alphaIn / W / Upst
    ZH = -8 * pi * alphaUp / W / Upst
    return (Zinf=Zinf, ZH=ZH)
end

"Z^∞/Z^H for a circular-equatorial s=0 SCALAR point particle: α=−4π r0² S (no 8π, no phase)."
function _convolve_source_circular_s0(l::Integer, m::Integer, a::Real, r0::Real,
        Upst::Real, omega::Real, lambda::Real; Nmax::Integer=80)
    s = 0; theta0 = pi / 2
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    S0 = Sc(theta0, 0.0)
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    Rin = R.In(r0); Rup = R.Up(r0)
    W = _teuk_wronskian(R.In)
    alpha = -4 * pi * r0^2 * S0
    Zinf = alpha * Rin / W / Upst
    ZH = alpha * Rup / W / Upst
    return (Zinf=Zinf, ZH=ZH)
end

"Z^∞/Z^H for a circular-equatorial s=∓1 MAXWELL point particle (closed form; s=−1: P=R, s=+1: P=ΔR)."
function _convolve_source_circular_maxwell(s::Integer, l::Integer, m::Integer, a::Real, r0::Real,
        Upst::Real, omega::Real, lambda::Real; Nmax::Integer=80)
    theta0 = pi / 2
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    S0 = Sc(theta0, 0.0)
    dS0 = swsh_derivative(Sc, 1, 0, theta0, 0.0)
    Delta = r0^2 - 2r0 + a^2; Deltap = 2r0 - 2
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    Rin = R.In(r0); dRin = teukolsky_radial_deriv(R.In, r0)
    Rout = R.Up(r0); dRout = teukolsky_radial_deriv(R.Up, r0)
    if s == -1
        Pin, Pout, dPin, dPout = Rin, Rout, dRin, dRout
    else  # s == +1
        Pin = Delta * Rin; Pout = Delta * Rout
        dPin = Delta * dRin + Deltap * Rin; dPout = Delta * dRout + Deltap * Rout
    end
    W = _teuk_wronskian(R.In)
    Om = 1 / (sqrt(r0^3) + a)
    Scal = (4 * pi) / (sqrt(2) * r0) * (s == -1 ? 1 / 2 : 1)
    Bcal = Delta * ((r0^2 + a^2) * Om - a)
    Ar = r0 * (r0 * ((r0^2 + a^2) * Om^2 - 1) + 2 * (1 - a * Om)^2)
    Ati = r0 * Delta * Om
    c = -Delta * (1 - a * Om)
    A = Scal * ((m * Ar + s * im * Ati) * S0 + s * c * dS0)
    B = s * im * Scal * Bcal * S0
    ZIn = (Pout * A - dPout * B) / (Delta * W)
    ZOut = (Pin * A - dPin * B) / (Delta * W)
    return (Zinf=ZOut, ZH=ZIn)
end

"Dispatch the circular-equatorial source convolution on spin weight `s`."
function _convolve_source_circular(s::Integer, l::Integer, m::Integer, a::Real, p::Real,
        En::Real, Lz::Real, Upst::Real, omega::Real, lambda::Real; Nmax::Integer=80)
    s == -2 && return _convolve_source_circular_s2(l, m, a, p, En, Lz, Upst, omega, lambda; Nmax=Nmax)
    s == 2 && return _convolve_source_circular_s2p(l, m, a, p, En, Lz, Upst, omega, lambda; Nmax=Nmax)
    s == 0 && return _convolve_source_circular_s0(l, m, a, p, Upst, omega, lambda; Nmax=Nmax)
    return _convolve_source_circular_maxwell(s, l, m, a, p, Upst, omega, lambda; Nmax=Nmax)
end

# ============================================================================
# ECCENTRIC-EQUATORIAL source convolutions (e>0, |x|=1, k=0): radial-libration Mino integral.
# ============================================================================

"Z^∞/Z^H for an ECCENTRIC-EQUATORIAL s=−2 point particle (radial-libration Mino integral)."
function _convolve_source_eccentric_s2(l::Integer, m::Integer, n::Integer, a::Real,
        p::Real, e::Real, En::Real, Lz::Real, Upst::Real, Upsr::Real,
        omega::Real, lambda::Real, rad; Nmax::Integer=80, N::Integer=64)
    s = -2; theta0 = pi / 2; sin0 = sin(theta0); cos0 = cos(theta0)
    sq2 = sqrt(2)
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    S0 = Sc(theta0, 0.0)
    dS0 = swsh_derivative(Sc, 1, 0, theta0, 0.0)
    d2S0 = swsh_derivative(Sc, 2, 0, theta0, 0.0)
    L1 = -m / sin0 + a * omega * sin0 + cos0 / sin0
    L2 = -m / sin0 + a * omega * sin0 + 2 * cos0 / sin0
    L2S = dS0 + L2 * S0
    L2p = m * cos0 / sin0^2 + a * omega * cos0 - 2 / sin0^2
    L1Sp = d2S0 + L1 * dS0
    L1L2S = L1Sp + L2p * S0 + L2 * dS0 + L1 * L2 * S0
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    W = _teuk_wronskian(R.In)
    function integrand(qr, RF)
        d = rad(qr)
        r0 = d.r
        R0 = RF(r0); dR0 = teukolsky_radial_deriv(RF, r0)
        d2R0 = _d2R_ode(R0, dR0, s, m, a, r0, omega, lambda)
        Δ = r0^2 + a^2 - 2 * r0
        Kt = (r0^2 + a^2) * omega - m * a
        rho = -1 / (r0 - im * a * cos0); rhobar = -1 / (r0 + im * a * cos0)
        Sigma = 1 / (rho * rhobar)
        Ann0 = -rho^(-2) * rhobar^(-1) * (sq2 * Δ)^(-2) *
               (rho^(-1) * L1L2S + 3 * im * a * sin0 * L1 * S0 + 3 * im * a * cos0 * S0 +
                2 * im * a * sin0 * dS0 - im * a * sin0 * L2 * S0)
        Anmbar0 = rho^(-3) * (sq2 * Δ)^(-1) *
                  ((rho + rhobar - im * Kt / Δ) * L2S + (rho - rhobar) * a * sin0 * Kt / Δ * S0)
        Anmbar1 = -rho^(-3) * (sq2 * Δ)^(-1) * (L2S + im * (rho - rhobar) * a * sin0 * S0)
        Ambarmbar0 = (Kt^2 * S0 * rhobar) / (4 * Δ^2 * rho^3) +
                     (im * Kt * S0 * (1 - r0 + Δ * rho) * rhobar) / (2 * Δ^2 * rho^3) +
                     (im * r0 * S0 * rhobar * omega) / (2 * Δ * rho^3)
        Ambarmbar1 = -rho^(-3) * rhobar * S0 / 2 * (im * Kt / Δ - rho)
        Ambarmbar2 = -rho^(-3) * rhobar * S0 / 4
        thetacomp = rho * (im * sin0 * (a * En - Lz / sin0^2)) / sq2
        urp = d.Ur
        urm = rad(2 * pi - qr).Ur
        rcompP = (En * (r0^2 + a^2) - a * Lz + urp) / (2 * Sigma)
        rcompM = (En * (r0^2 + a^2) - a * Lz + urm) / (2 * Sigma)
        CnnP = rcompP^2; CnmbarP = rcompP * thetacomp; CmbarmbarP = thetacomp^2
        CnnM = rcompM^2; CnmbarM = rcompM * thetacomp; CmbarmbarM = thetacomp^2
        rphase = omega * d.Δtr - m * d.Δφr + n * qr
        ph = cis(rphase)
        resP = (Ann0 * CnnP + Anmbar0 * CnmbarP + Ambarmbar0 * CmbarmbarP) * R0 -
               (Anmbar1 * CnmbarP + Ambarmbar1 * CmbarmbarP) * dR0 + Ambarmbar2 * CmbarmbarP * d2R0
        resM = (Ann0 * CnnM + Anmbar0 * CnmbarM + Ambarmbar0 * CmbarmbarM) * R0 -
               (Anmbar1 * CnmbarM + Ambarmbar1 * CmbarmbarM) * dR0 + Ambarmbar2 * CmbarmbarM * d2R0
        return resP * ph + resM * conj(ph)
    end
    αIn = _trapz_halfperiod(qr -> integrand(qr, R.In), N) / (2 * pi)
    αUp = _trapz_halfperiod(qr -> integrand(qr, R.Up), N) / (2 * pi)
    d0 = rad(0.0)
    xi = m * d0.Δφr - omega * d0.Δtr
    eixi = cis(xi)
    Zinf = -8 * pi * αIn / W / Upst * eixi
    ZH = -8 * pi * αUp / W / Upst * eixi
    return (Zinf=Zinf, ZH=ZH, N=N)
end

"Z^∞/Z^H for an ECCENTRIC-EQUATORIAL s=+2 (ψ0) point particle (Δ²-rescaled, radial libration)."
function _convolve_source_eccentric_s2p(l::Integer, m::Integer, n::Integer, a::Real,
        p::Real, e::Real, En::Real, Lz::Real, Upst::Real, Upsr::Real,
        omega::Real, lambda::Real, rad; Nmax::Integer=80, N::Integer=64)
    s = 2; theta0 = pi / 2; sin0 = sin(theta0); cos0 = cos(theta0)
    sq2 = sqrt(2)
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    S0 = Sc(theta0, 0.0)
    dS0 = swsh_derivative(Sc, 1, 0, theta0, 0.0)
    d2S0 = swsh_derivative(Sc, 2, 0, theta0, 0.0)
    deltaLd1 = m / sin0 - a * omega * sin0 + cos0 / sin0
    deltaLd2 = m / sin0 - a * omega * sin0 + 2 * cos0 / sin0
    ddeltaLd2 = -m * cos0 / sin0^2 - a * omega * cos0 - 2 / sin0^2
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    W = _teuk_wronskian(R.In)
    function integrand(qr, RF)
        d = rad(qr)
        r0 = d.r
        R0 = RF(r0); dR0 = teukolsky_radial_deriv(RF, r0)
        d2R0 = _d2R_ode(R0, dR0, s, m, a, r0, omega, lambda)
        Δ = r0^2 + a^2 - 2 * r0; dΔ = 2 * (r0 - 1); d2Δ = 2.0
        Kt = (r0^2 + a^2) * omega - m * a
        rho = -1 / (r0 - im * a * cos0); rhobar = -1 / (r0 + im * a * cos0)
        drhoOverrho = im * a * rho * sin0
        d2rhoOverrho = im * a * rho * (cos0 + 2 * sin0 * drhoOverrho)
        All0 = -(1 / 2) * rho^(-1) * rhobar *
               (d2S0 + (deltaLd1 + deltaLd2 + 2 * drhoOverrho) * dS0 +
                (ddeltaLd2 + deltaLd1 * deltaLd2 - 6 * drhoOverrho^2 + 3 * d2rhoOverrho +
                 (3 * deltaLd1 - deltaLd2) * drhoOverrho) * S0)
        Alm0 = (2 / sq2) * rho^(-1) *
               (-(rho + rhobar + im * Kt / Δ) * (dS0 + deltaLd2 * S0) +
                (rho - rhobar) * a * sin0 * Kt / Δ * S0)
        Alm1 = (2 / sq2) * rho^(-1) * ((dS0 + deltaLd2 * S0) + im * (rho - rhobar) * a * sin0 * S0)
        Amm0 = (Kt^2 * S0) / (Δ^2 * rho * rhobar) +
               (2 * im * Kt * S0 * (-1 + r0 - Δ * rho)) / (Δ^2 * rho * rhobar) -
               (2 * im * r0 * S0 * omega) / (Δ * rho * rhobar)
        Amm1 = 2 * rho^(-1) * rhobar^(-1) * S0 * (im * Kt / Δ + rho)
        Amm2 = -rho^(-1) * rhobar^(-1) * S0
        thetacomp = -rhobar * (im * sin0 * (a * En - Lz / sin0^2)) / sq2
        urp = d.Ur
        urm = rad(2 * pi - qr).Ur
        rcompP = (En * (r0^2 + a^2) - a * Lz - urp) / Δ
        rcompM = (En * (r0^2 + a^2) - a * Lz - urm) / Δ
        CllP = rcompP^2; ClmP = rcompP * thetacomp; CmmP = thetacomp^2
        CllM = rcompM^2; ClmM = rcompM * thetacomp; CmmM = thetacomp^2
        D2R = Δ^2 * R0
        dD2R = Δ^2 * dR0 + 2 * Δ * dΔ * R0
        d2D2R = Δ^2 * d2R0 + 4 * Δ * dΔ * dR0 + (2 * dΔ^2 + 2 * Δ * d2Δ) * R0
        rphase = omega * d.Δtr - m * d.Δφr + n * qr
        ph = cis(rphase)
        resP = (All0 * CllP + Alm0 * ClmP + Amm0 * CmmP) * D2R -
               (Alm1 * ClmP + Amm1 * CmmP) * dD2R + Amm2 * CmmP * d2D2R
        resM = (All0 * CllM + Alm0 * ClmM + Amm0 * CmmM) * D2R -
               (Alm1 * ClmM + Amm1 * CmmM) * dD2R + Amm2 * CmmM * d2D2R
        return resP * ph + resM * conj(ph)
    end
    αIn = _trapz_halfperiod(qr -> integrand(qr, R.In), N) / (2 * pi)
    αUp = _trapz_halfperiod(qr -> integrand(qr, R.Up), N) / (2 * pi)
    d0 = rad(0.0)
    xi = m * d0.Δφr - omega * d0.Δtr
    eixi = cis(xi)
    Zinf = -8 * pi * αIn / W / Upst * eixi
    ZH = -8 * pi * αUp / W / Upst * eixi
    return (Zinf=Zinf, ZH=ZH, N=N)
end

"Z^∞/Z^H for an ECCENTRIC-EQUATORIAL s=∓1 MAXWELL point particle (radial libration)."
function _convolve_source_eccentric_maxwell(s::Integer, l::Integer, m::Integer, n::Integer,
        a::Real, p::Real, e::Real, En::Real, Lz::Real, Upst::Real, Upsr::Real,
        omega::Real, lambda::Real, rad; Nmax::Integer=80, N::Integer=64)
    theta0 = pi / 2; sin0 = sin(theta0); cos0 = cos(theta0)
    sq2 = sqrt(2); csc0 = 1 / sin0
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    S0 = Sc(theta0, 0.0)
    dS0 = swsh_derivative(Sc, 1, 0, theta0, 0.0)
    L1 = (s == -1 ? -1 : 1) * m / sin0 + (s == -1 ? 1 : -1) * a * omega * sin0 + cos0 / sin0
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    W = _teuk_wronskian(R.In)
    function integrand(qr, RF)
        d = rad(qr)
        r0 = d.r
        R0 = RF(r0); dR0 = teukolsky_radial_deriv(RF, r0)
        Δ = r0^2 + a^2 - 2 * r0
        Kt = (r0^2 + a^2) * omega - m * a
        rho = -1 / (r0 - im * a * cos0); rhobar = -1 / (r0 + im * a * cos0)
        Sigma = 1 / (rho * rhobar)
        urp = d.Ur
        urm = rad(2 * pi - qr).Ur
        if s == -1
            Ar0 = -((dS0 + L1 * S0 + im * a * S0 * rho * sin0) / (2 * sq2 * Δ * rho^2 * rhobar))
            Am0 = (S0 * (-(im * Kt / Δ) + rho)) / (4 * rho^2)
            Am1 = -(S0 / (4 * rho^2))
            thetacomp = -(rho * (im * Lz * csc0 - im * a * En * sin0)) / sq2
            rcompP = (urp + (a^2 + r0^2) * En - a * Lz) / (2 * Sigma)
            rcompM = (urm + (a^2 + r0^2) * En - a * Lz) / (2 * Sigma)
        else  # s == +1
            Δp = 2 * (r0 - 1)
            Ar0 = -((Δ * (dS0 + L1 * S0 + im * a * S0 * rho * sin0)) / (2 * sq2 * rho))
            Am0 = -((S0 * (im * Kt + Δp + Δ * rho)) / (2 * rho * rhobar))
            Am1 = (S0 * Δ) / (2 * rho * rhobar)
            thetacomp = (rhobar * (-im * Lz * csc0 + im * a * En * sin0)) / sq2
            rcompP = (urp - (a^2 + r0^2) * En + a * Lz) / Δ
            rcompM = (urm - (a^2 + r0^2) * En + a * Lz) / Δ
        end
        resP = (Ar0 * rcompP + Am0 * thetacomp) * R0 - (Am1 * thetacomp) * dR0
        resM = (Ar0 * rcompM + Am0 * thetacomp) * R0 - (Am1 * thetacomp) * dR0
        rphase = omega * d.Δtr - m * d.Δφr + n * qr
        ph = cis(rphase)
        return resP * ph + resM * conj(ph)
    end
    αIn = _trapz_halfperiod(qr -> integrand(qr, R.In), N) / (2 * pi)
    αUp = _trapz_halfperiod(qr -> integrand(qr, R.Up), N) / (2 * pi)
    d0 = rad(0.0)
    xi = m * d0.Δφr - omega * d0.Δtr
    eixi = cis(xi)
    Zinf = 8 * pi * αIn / W / Upst * eixi
    ZH = 8 * pi * αUp / W / Upst * eixi
    return (Zinf=Zinf, ZH=ZH, N=N)
end

"Z^∞/Z^H for an ECCENTRIC-EQUATORIAL s=0 SCALAR point particle (cos-fold radial libration)."
function _convolve_source_eccentric_s0(l::Integer, m::Integer, n::Integer, a::Real,
        p::Real, e::Real, Upst::Real, Upsr::Real,
        omega::Real, lambda::Real, rad; Nmax::Integer=80, N::Integer=64)
    s = 0; theta0 = pi / 2; cos0 = cos(theta0)
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    S0 = Sc(theta0, 0.0)
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    W = _teuk_wronskian(R.In)
    function integrand(qr, RF)
        d = rad(qr)
        r0 = d.r
        rphase = omega * d.Δtr - m * d.Δφr + n * qr
        return -8 * pi * (r0^2 + a^2 * cos0^2) * RF(r0) * cos(rphase)
    end
    αIn = S0 * _trapz_halfperiod(qr -> integrand(qr, R.In), N) / (2 * pi)
    αUp = S0 * _trapz_halfperiod(qr -> integrand(qr, R.Up), N) / (2 * pi)
    d0 = rad(0.0)
    xi = m * d0.Δφr - omega * d0.Δtr
    eixi = cis(xi)
    Zinf = αIn / W / Upst * eixi
    ZH = αUp / W / Upst * eixi
    return (Zinf=Zinf, ZH=ZH, N=N)
end

"Dispatch the eccentric-equatorial source convolution on spin weight `s`."
function _convolve_source_eccentric(s::Integer, l::Integer, m::Integer, n::Integer, a::Real,
        p::Real, e::Real, En::Real, Lz::Real, Upst::Real, Upsr::Real,
        omega::Real, lambda::Real, rad; Nmax::Integer=80, N::Integer=64)
    s == -2 && return _convolve_source_eccentric_s2(l, m, n, a, p, e, En, Lz, Upst, Upsr, omega, lambda, rad; Nmax=Nmax, N=N)
    s == 2 && return _convolve_source_eccentric_s2p(l, m, n, a, p, e, En, Lz, Upst, Upsr, omega, lambda, rad; Nmax=Nmax, N=N)
    s == 0 && return _convolve_source_eccentric_s0(l, m, n, a, p, e, Upst, Upsr, omega, lambda, rad; Nmax=Nmax, N=N)
    return _convolve_source_eccentric_maxwell(s, l, m, n, a, p, e, En, Lz, Upst, Upsr, omega, lambda, rad; Nmax=Nmax, N=N)
end

# ============================================================================
# SPHERICAL (e=0, inclined |x|<1, n=0) source convolutions: polar-libration Mino integral at r=p.
# ============================================================================

"Z^∞/Z^H for a SPHERICAL s=−2 point particle (polar-libration Mino integral)."
function _convolve_source_spherical_s2(l::Integer, m::Integer, k::Integer, a::Real,
        p::Real, En::Real, Lz::Real, Upst::Real,
        omega::Real, lambda::Real, pol; Nmax::Integer=80, N::Integer=64)
    s = -2; sq2 = sqrt(2); r0 = float(p)
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    Δ = r0^2 - 2r0 + a^2
    Kt = (r0^2 + a^2) * omega - m * a
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    W = _teuk_wronskian(R.In)
    RinV = R.In(r0); dRinV = teukolsky_radial_deriv(R.In, r0)
    d2RinV = _d2R_ode(RinV, dRinV, s, m, a, r0, omega, lambda)
    RupV = R.Up(r0); dRupV = teukolsky_radial_deriv(R.Up, r0)
    d2RupV = _d2R_ode(RupV, dRupV, s, m, a, r0, omega, lambda)
    function integrand(qθ, R0, dR0, d2R0)
        d = pol(qθ)
        theta0 = d.θ; sin0 = sin(theta0); cos0 = cos(theta0)
        S0 = Sc(theta0, 0.0)
        dS0 = swsh_derivative(Sc, 1, 0, theta0, 0.0)
        d2S0 = swsh_derivative(Sc, 2, 0, theta0, 0.0)
        L1 = -m / sin0 + a * omega * sin0 + cos0 / sin0
        L2 = -m / sin0 + a * omega * sin0 + 2 * cos0 / sin0
        L2S = dS0 + L2 * S0
        L2p = m * cos0 / sin0^2 + a * omega * cos0 - 2 / sin0^2
        L1Sp = d2S0 + L1 * dS0
        L1L2S = L1Sp + L2p * S0 + L2 * dS0 + L1 * L2 * S0
        rho = -1 / (r0 - im * a * cos0); rhobar = -1 / (r0 + im * a * cos0)
        Sigma = 1 / (rho * rhobar)
        Ann0 = -rho^(-2) * rhobar^(-1) * (sq2 * Δ)^(-2) *
               (rho^(-1) * L1L2S + 3 * im * a * sin0 * L1 * S0 + 3 * im * a * cos0 * S0 +
                2 * im * a * sin0 * dS0 - im * a * sin0 * L2 * S0)
        Anmbar0 = rho^(-3) * (sq2 * Δ)^(-1) *
                  ((rho + rhobar - im * Kt / Δ) * L2S + (rho - rhobar) * a * sin0 * Kt / Δ * S0)
        Anmbar1 = -rho^(-3) * (sq2 * Δ)^(-1) * (L2S + im * (rho - rhobar) * a * sin0 * S0)
        Ambarmbar0 = (Kt^2 * S0 * rhobar) / (4 * Δ^2 * rho^3) +
                     (im * Kt * S0 * (1 - r0 + Δ * rho) * rhobar) / (2 * Δ^2 * rho^3) +
                     (im * r0 * S0 * rhobar * omega) / (2 * Δ * rho^3)
        Ambarmbar1 = -rho^(-3) * rhobar * S0 / 2 * (im * Kt / Δ - rho)
        Ambarmbar2 = -rho^(-3) * rhobar * S0 / 4
        uθP = -d.Uz / sin0
        uθM = -pol(2 * pi - qθ).Uz / sin0
        rcomp = (En * (r0^2 + a^2) - a * Lz) / (2 * Sigma)
        θcompP = rho * (im * sin0 * (a * En - Lz / sin0^2) + uθP) / sq2
        θcompM = rho * (im * sin0 * (a * En - Lz / sin0^2) + uθM) / sq2
        CnnP = rcomp^2; CnmbarP = rcomp * θcompP; CmbarmbarP = θcompP^2
        CnnM = rcomp^2; CnmbarM = rcomp * θcompM; CmbarmbarM = θcompM^2
        θphase = omega * d.Δtz - m * d.Δφz + k * qθ
        ph = cis(θphase)
        resP = (Ann0 * CnnP + Anmbar0 * CnmbarP + Ambarmbar0 * CmbarmbarP) * R0 -
               (Anmbar1 * CnmbarP + Ambarmbar1 * CmbarmbarP) * dR0 + Ambarmbar2 * CmbarmbarP * d2R0
        resM = (Ann0 * CnnM + Anmbar0 * CnmbarM + Ambarmbar0 * CmbarmbarM) * R0 -
               (Anmbar1 * CnmbarM + Ambarmbar1 * CmbarmbarM) * dR0 + Ambarmbar2 * CmbarmbarM * d2R0
        return resP * ph + resM * conj(ph)
    end
    αIn = _trapz_romberg(qθ -> integrand(qθ, RinV, dRinV, d2RinV), N) / (2 * pi)
    αUp = _trapz_romberg(qθ -> integrand(qθ, RupV, dRupV, d2RupV), N) / (2 * pi)
    d0 = pol(0.0)
    xi = m * d0.Δφz - omega * d0.Δtz
    eixi = cis(xi)
    Zinf = -8 * pi * αIn / W / Upst * eixi
    ZH = -8 * pi * αUp / W / Upst * eixi
    return (Zinf=Zinf, ZH=ZH, N=N)
end

"Z^∞/Z^H for a SPHERICAL s=+2 (ψ0) point particle (Δ²-rescaled, polar libration)."
function _convolve_source_spherical_s2p(l::Integer, m::Integer, k::Integer, a::Real,
        p::Real, En::Real, Lz::Real, Upst::Real,
        omega::Real, lambda::Real, pol; Nmax::Integer=80, N::Integer=64)
    s = 2; sq2 = sqrt(2); r0 = float(p)
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    Δ = r0^2 + a^2 - 2r0; dΔ = 2 * (r0 - 1); d2Δ = 2.0
    Kt = (r0^2 + a^2) * omega - m * a
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    W = _teuk_wronskian(R.In)
    d2rfn(R0, dR0) = _d2R_ode(R0, dR0, s, m, a, r0, omega, lambda)
    function D2pack(RF)
        R0 = RF(r0); dR0 = teukolsky_radial_deriv(RF, r0); d2R0 = d2rfn(R0, dR0)
        D2R = Δ^2 * R0
        dD2R = Δ^2 * dR0 + 2 * Δ * dΔ * R0
        d2D2R = Δ^2 * d2R0 + 4 * Δ * dΔ * dR0 + (2 * dΔ^2 + 2 * Δ * d2Δ) * R0
        return (D2R, dD2R, d2D2R)
    end
    D2In = D2pack(R.In); D2Up = D2pack(R.Up)
    function integrand(qθ, D2R, dD2R, d2D2R)
        d = pol(qθ)
        theta0 = d.θ; sin0 = sin(theta0); cos0 = cos(theta0)
        S0 = Sc(theta0, 0.0)
        dS0 = swsh_derivative(Sc, 1, 0, theta0, 0.0)
        d2S0 = swsh_derivative(Sc, 2, 0, theta0, 0.0)
        deltaLd1 = m / sin0 - a * omega * sin0 + cos0 / sin0
        deltaLd2 = m / sin0 - a * omega * sin0 + 2 * cos0 / sin0
        ddeltaLd2 = -m * cos0 / sin0^2 - a * omega * cos0 - 2 / sin0^2
        rho = -1 / (r0 - im * a * cos0); rhobar = -1 / (r0 + im * a * cos0)
        drhoOverrho = im * a * rho * sin0
        d2rhoOverrho = im * a * rho * (cos0 + 2 * sin0 * drhoOverrho)
        All0 = -(1 / 2) * rho^(-1) * rhobar *
               (d2S0 + (deltaLd1 + deltaLd2 + 2 * drhoOverrho) * dS0 +
                (ddeltaLd2 + deltaLd1 * deltaLd2 - 6 * drhoOverrho^2 + 3 * d2rhoOverrho +
                 (3 * deltaLd1 - deltaLd2) * drhoOverrho) * S0)
        Alm0 = (2 / sq2) * rho^(-1) *
               (-(rho + rhobar + im * Kt / Δ) * (dS0 + deltaLd2 * S0) +
                (rho - rhobar) * a * sin0 * Kt / Δ * S0)
        Alm1 = (2 / sq2) * rho^(-1) * ((dS0 + deltaLd2 * S0) + im * (rho - rhobar) * a * sin0 * S0)
        Amm0 = (Kt^2 * S0) / (Δ^2 * rho * rhobar) +
               (2 * im * Kt * S0 * (-1 + r0 - Δ * rho)) / (Δ^2 * rho * rhobar) -
               (2 * im * r0 * S0 * omega) / (Δ * rho * rhobar)
        Amm1 = 2 * rho^(-1) * rhobar^(-1) * S0 * (im * Kt / Δ + rho)
        Amm2 = -rho^(-1) * rhobar^(-1) * S0
        uθP = -d.Uz / sin0
        uθM = -pol(2 * pi - qθ).Uz / sin0
        rcomp = (En * (r0^2 + a^2) - a * Lz) / Δ
        θcompP = -rhobar * (im * sin0 * (a * En - Lz / sin0^2) - uθP) / sq2
        θcompM = -rhobar * (im * sin0 * (a * En - Lz / sin0^2) - uθM) / sq2
        CllP = rcomp^2; ClmP = rcomp * θcompP; CmmP = θcompP^2
        CllM = rcomp^2; ClmM = rcomp * θcompM; CmmM = θcompM^2
        θphase = omega * d.Δtz - m * d.Δφz + k * qθ
        ph = cis(θphase)
        resP = (All0 * CllP + Alm0 * ClmP + Amm0 * CmmP) * D2R -
               (Alm1 * ClmP + Amm1 * CmmP) * dD2R + Amm2 * CmmP * d2D2R
        resM = (All0 * CllM + Alm0 * ClmM + Amm0 * CmmM) * D2R -
               (Alm1 * ClmM + Amm1 * CmmM) * dD2R + Amm2 * CmmM * d2D2R
        return resP * ph + resM * conj(ph)
    end
    αIn = _trapz_romberg(qθ -> integrand(qθ, D2In...), N) / (2 * pi)
    αUp = _trapz_romberg(qθ -> integrand(qθ, D2Up...), N) / (2 * pi)
    d0 = pol(0.0)
    xi = m * d0.Δφz - omega * d0.Δtz
    eixi = cis(xi)
    Zinf = -8 * pi * αIn / W / Upst * eixi
    ZH = -8 * pi * αUp / W / Upst * eixi
    return (Zinf=Zinf, ZH=ZH, N=N)
end

"Z^∞/Z^H for a SPHERICAL s=∓1 MAXWELL point particle (polar libration)."
function _convolve_source_spherical_maxwell(s::Integer, l::Integer, m::Integer, k::Integer,
        a::Real, p::Real, En::Real, Lz::Real, Upst::Real,
        omega::Real, lambda::Real, pol; Nmax::Integer=80, N::Integer=64)
    sq2 = sqrt(2); r0 = float(p)
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    Δ = r0^2 - 2r0 + a^2; Δp = 2 * (r0 - 1)
    Kt = (r0^2 + a^2) * omega - m * a
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    W = _teuk_wronskian(R.In)
    RinV = R.In(r0); dRinV = teukolsky_radial_deriv(R.In, r0)
    RupV = R.Up(r0); dRupV = teukolsky_radial_deriv(R.Up, r0)
    function integrand(qθ, R0, dR0)
        d = pol(qθ)
        theta0 = d.θ; sin0 = sin(theta0); cos0 = cos(theta0)
        S0 = Sc(theta0, 0.0)
        dS0 = swsh_derivative(Sc, 1, 0, theta0, 0.0)
        L1 = (s == -1 ? -1 : 1) * m / sin0 + (s == -1 ? 1 : -1) * a * omega * sin0 + cos0 / sin0
        rho = -1 / (r0 - im * a * cos0); rhobar = -1 / (r0 + im * a * cos0)
        Sigma = 1 / (rho * rhobar)
        uθP = -d.Uz / sin0
        uθM = -pol(2 * pi - qθ).Uz / sin0
        if s == -1
            Ar0 = -((dS0 + L1 * S0 + im * a * S0 * rho * sin0) / (2 * sq2 * Δ * rho^2 * rhobar))
            Am0 = (S0 * (-(im * Kt / Δ) + rho)) / (4 * rho^2)
            Am1 = -(S0 / (4 * rho^2))
            rcomp = ((a^2 + r0^2) * En - a * Lz) / (2 * Sigma)
            θcompP = rho * (im * sin0 * (a * En - Lz / sin0^2) + uθP) / sq2
            θcompM = rho * (im * sin0 * (a * En - Lz / sin0^2) + uθM) / sq2
        else  # s == +1
            Ar0 = -((Δ * (dS0 + L1 * S0 + im * a * S0 * rho * sin0)) / (2 * sq2 * rho))
            Am0 = -((S0 * (im * Kt + Δp + Δ * rho)) / (2 * rho * rhobar))
            Am1 = (S0 * Δ) / (2 * rho * rhobar)
            rcomp = -((a^2 + r0^2) * En - a * Lz) / Δ
            θcompP = rhobar * (im * sin0 * (a * En - Lz / sin0^2) - uθP) / sq2
            θcompM = rhobar * (im * sin0 * (a * En - Lz / sin0^2) - uθM) / sq2
        end
        θphase = omega * d.Δtz - m * d.Δφz + k * qθ
        ph = cis(θphase)
        resP = (Ar0 * rcomp + Am0 * θcompP) * R0 - (Am1 * θcompP) * dR0
        resM = (Ar0 * rcomp + Am0 * θcompM) * R0 - (Am1 * θcompM) * dR0
        return resP * ph + resM * conj(ph)
    end
    αIn = _trapz_romberg(qθ -> integrand(qθ, RinV, dRinV), N) / (2 * pi)
    αUp = _trapz_romberg(qθ -> integrand(qθ, RupV, dRupV), N) / (2 * pi)
    d0 = pol(0.0)
    xi = m * d0.Δφz - omega * d0.Δtz
    eixi = cis(xi)
    Zinf = 8 * pi * αIn / W / Upst * eixi
    ZH = 8 * pi * αUp / W / Upst * eixi
    return (Zinf=Zinf, ZH=ZH, N=N)
end

"Z^∞/Z^H for a SPHERICAL s=0 SCALAR point particle (separable r0²/a²cos²θ polar integrals)."
function _convolve_source_spherical_s0(l::Integer, m::Integer, k::Integer, a::Real,
        p::Real, Upst::Real, omega::Real, lambda::Real, pol; Nmax::Integer=80, N::Integer=64)
    s = 0; r0 = float(p)
    if a != 0
        isodd(l + m + k) && return (Zinf=zero(ComplexF64), ZH=zero(ComplexF64), N=N)
    else
        (isodd(l + m + k) || abs(m + k) > l) && return (Zinf=zero(ComplexF64), ZH=zero(ComplexF64), N=N)
    end
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    W = _teuk_wronskian(R.In)
    RinV = R.In(r0); RupV = R.Up(r0)
    II2(qθ) = begin
        d = pol(qθ)
        -8 * pi * Sc(d.θ, 0.0) * cos(omega * d.Δtz - m * d.Δφz + k * qθ)
    end
    α2 = _trapz_halfperiod(II2, N) / (2 * pi)
    α1In = r0^2 * RinV; α1Up = r0^2 * RupV
    if a == 0
        α3In = zero(ComplexF64); α3Up = zero(ComplexF64); α4 = zero(ComplexF64)
    else
        II4(qθ) = begin
            d = pol(qθ)
            -8 * pi * a^2 * cos(d.θ)^2 * Sc(d.θ, 0.0) * cos(omega * d.Δtz - m * d.Δφz + k * qθ)
        end
        α4 = _trapz_halfperiod(II4, N) / (2 * pi)
        α3In = RinV; α3Up = RupV
    end
    d0 = pol(0.0)
    xi = m * d0.Δφz - omega * d0.Δtz
    eixi = cis(xi)
    Zinf = (α1In * α2 + α3In * α4) / W / Upst * eixi
    ZH = (α1Up * α2 + α3Up * α4) / W / Upst * eixi
    return (Zinf=Zinf, ZH=ZH, N=N)
end

"Dispatch the spherical (inclined e=0) source convolution on spin weight `s`."
function _convolve_source_spherical(s::Integer, l::Integer, m::Integer, k::Integer, a::Real,
        p::Real, En::Real, Lz::Real, Upst::Real, omega::Real, lambda::Real, pol; Nmax::Integer=80, N::Integer=64)
    s == -2 && return _convolve_source_spherical_s2(l, m, k, a, p, En, Lz, Upst, omega, lambda, pol; Nmax=Nmax, N=N)
    s == 2 && return _convolve_source_spherical_s2p(l, m, k, a, p, En, Lz, Upst, omega, lambda, pol; Nmax=Nmax, N=N)
    s == 0 && return _convolve_source_spherical_s0(l, m, k, a, p, Upst, omega, lambda, pol; Nmax=Nmax, N=N)
    return _convolve_source_spherical_maxwell(s, l, m, k, a, p, En, Lz, Upst, omega, lambda, pol; Nmax=Nmax, N=N)
end

# ============================================================================
# GENERIC (inclined eccentric: e>0 AND |x|<1) source convolutions: (qr,qθ) double Mino integral.
# ============================================================================

"Z^∞/Z^H for a GENERIC s=−2 point particle (double Mino integral, four-term velocity fold)."
function _convolve_source_generic_s2(l::Integer, m::Integer, n::Integer, k::Integer, a::Real,
        p::Real, e::Real, x::Real, En::Real, Lz::Real, Upst::Real,
        omega::Real, lambda::Real, rad, pol; Nmax::Integer=80, Nr::Integer=64, Nθ::Integer=64)
    s = -2; sq2 = sqrt(2)
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    W = _teuk_wronskian(R.In)
    function θblock(qθ)
        dθ_ = pol(qθ); theta0 = dθ_.θ; sin0 = sin(theta0); cos0 = cos(theta0)
        S0 = Sc(theta0, 0.0)
        dS0 = swsh_derivative(Sc, 1, 0, theta0, 0.0)
        d2S0 = swsh_derivative(Sc, 2, 0, theta0, 0.0)
        uθp = -dθ_.Uz / sin0; uθm = -pol(2 * pi - qθ).Uz / sin0
        cθ = cis(omega * dθ_.Δtz - m * dθ_.Δφz + k * qθ)
        return (sin0=sin0, cos0=cos0, S0=S0, dS0=dS0, d2S0=d2S0, uθp=uθp, uθm=uθm, cθ=cθ)
    end
    function rblock(qr, RF)
        dr_ = rad(qr); r0 = dr_.r
        R0 = RF(r0); dR0 = teukolsky_radial_deriv(RF, r0); d2R0 = _d2R_ode(R0, dR0, s, m, a, r0, omega, lambda)
        Δ = r0^2 + a^2 - 2 * r0; Kt = (r0^2 + a^2) * omega - m * a
        urp = dr_.Ur; urm = rad(2 * pi - qr).Ur
        cr = cis(omega * dr_.Δtr - m * dr_.Δφr + n * qr)
        return (r0=r0, R0=R0, dR0=dR0, d2R0=d2R0, Δ=Δ, Kt=Kt, urp=urp, urm=urm, cr=cr)
    end
    function combine(rb, θb)
        r0 = rb.r0; Δ = rb.Δ; Kt = rb.Kt; sin0 = θb.sin0; cos0 = θb.cos0
        S0 = θb.S0; dS0 = θb.dS0; d2S0 = θb.d2S0
        L1 = -m / sin0 + a * omega * sin0 + cos0 / sin0
        L2 = -m / sin0 + a * omega * sin0 + 2 * cos0 / sin0
        L2S = dS0 + L2 * S0
        L2p = m * cos0 / sin0^2 + a * omega * cos0 - 2 / sin0^2
        L1Sp = d2S0 + L1 * dS0
        L1L2S = L1Sp + L2p * S0 + L2 * dS0 + L1 * L2 * S0
        rho = -1 / (r0 - im * a * cos0); rhobar = -1 / (r0 + im * a * cos0)
        Sigma = 1 / (rho * rhobar)
        Ann0 = -rho^(-2) * rhobar^(-1) * (sq2 * Δ)^(-2) *
               (rho^(-1) * L1L2S + 3 * im * a * sin0 * L1 * S0 + 3 * im * a * cos0 * S0 +
                2 * im * a * sin0 * dS0 - im * a * sin0 * L2 * S0)
        Anmbar0 = rho^(-3) * (sq2 * Δ)^(-1) *
                  ((rho + rhobar - im * Kt / Δ) * L2S + (rho - rhobar) * a * sin0 * Kt / Δ * S0)
        Anmbar1 = -rho^(-3) * (sq2 * Δ)^(-1) * (L2S + im * (rho - rhobar) * a * sin0 * S0)
        Ambarmbar0 = (Kt^2 * S0 * rhobar) / (4 * Δ^2 * rho^3) +
                     (im * Kt * S0 * (1 - r0 + Δ * rho) * rhobar) / (2 * Δ^2 * rho^3) +
                     (im * r0 * S0 * rhobar * omega) / (2 * Δ * rho^3)
        Ambarmbar1 = -rho^(-3) * rhobar * S0 / 2 * (im * Kt / Δ - rho)
        Ambarmbar2 = -rho^(-3) * rhobar * S0 / 4
        rcomp(ur0) = (En * (r0^2 + a^2) - a * Lz + ur0) / (2 * Sigma)
        θcomp(uθ0) = rho * (im * sin0 * (a * En - Lz / sin0^2) + uθ0) / sq2
        proj(rc, θc) = (Ann0 * rc^2 + Anmbar0 * rc * θc + Ambarmbar0 * θc^2) * rb.R0 -
                       (Anmbar1 * rc * θc + Ambarmbar1 * θc^2) * rb.dR0 + Ambarmbar2 * θc^2 * rb.d2R0
        cr = rb.cr; cθ = θb.cθ
        return proj(rcomp(rb.urp), θcomp(θb.uθp)) * cr * cθ +
               proj(rcomp(rb.urp), θcomp(θb.uθm)) * cr * conj(cθ) +
               proj(rcomp(rb.urm), θcomp(θb.uθp)) * conj(cr) * cθ +
               proj(rcomp(rb.urm), θcomp(θb.uθm)) * conj(cr) * conj(cθ)
    end
    αIn = _trapz2d_romberg_sep(qr -> rblock(qr, R.In), θblock, combine, Nr, Nθ) / (2 * pi)^2
    αUp = _trapz2d_romberg_sep(qr -> rblock(qr, R.Up), θblock, combine, Nr, Nθ) / (2 * pi)^2
    Zinf = -8 * pi * αIn / W / Upst
    ZH = -8 * pi * αUp / W / Upst
    return (Zinf=Zinf, ZH=ZH, Nr=Nr, Nθ=Nθ)
end

"Z^∞/Z^H for a GENERIC s=+2 (ψ0) point particle (Δ²-rescaled, double Mino integral)."
function _convolve_source_generic_s2p(l::Integer, m::Integer, n::Integer, k::Integer, a::Real,
        p::Real, e::Real, x::Real, En::Real, Lz::Real, Upst::Real,
        omega::Real, lambda::Real, rad, pol; Nmax::Integer=80, Nr::Integer=64, Nθ::Integer=64)
    s = 2; sq2 = sqrt(2)
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    W = _teuk_wronskian(R.In)
    function θblock(qθ)
        dθ_ = pol(qθ); theta0 = dθ_.θ; sin0 = sin(theta0); cos0 = cos(theta0)
        S0 = Sc(theta0, 0.0)
        dS0 = swsh_derivative(Sc, 1, 0, theta0, 0.0)
        d2S0 = swsh_derivative(Sc, 2, 0, theta0, 0.0)
        uθp = -dθ_.Uz / sin0; uθm = -pol(2 * pi - qθ).Uz / sin0
        cθ = cis(omega * dθ_.Δtz - m * dθ_.Δφz + k * qθ)
        return (sin0=sin0, cos0=cos0, S0=S0, dS0=dS0, d2S0=d2S0, uθp=uθp, uθm=uθm, cθ=cθ)
    end
    function rblock(qr, RF)
        dr_ = rad(qr); r0 = dr_.r
        R0 = RF(r0); dR0 = teukolsky_radial_deriv(RF, r0); d2R0 = _d2R_ode(R0, dR0, s, m, a, r0, omega, lambda)
        Δ = r0^2 + a^2 - 2 * r0; dΔ = 2 * (r0 - 1); d2Δ = 2.0; Kt = (r0^2 + a^2) * omega - m * a
        urp = dr_.Ur; urm = rad(2 * pi - qr).Ur
        cr = cis(omega * dr_.Δtr - m * dr_.Δφr + n * qr)
        D2R = Δ^2 * R0
        dD2R = Δ^2 * dR0 + 2 * Δ * dΔ * R0
        d2D2R = Δ^2 * d2R0 + 4 * Δ * dΔ * dR0 + (2 * dΔ^2 + 2 * Δ * d2Δ) * R0
        return (r0=r0, Δ=Δ, Kt=Kt, urp=urp, urm=urm, cr=cr, D2R=D2R, dD2R=dD2R, d2D2R=d2D2R)
    end
    function combine(rb, θb)
        r0 = rb.r0; Δ = rb.Δ; Kt = rb.Kt; sin0 = θb.sin0; cos0 = θb.cos0
        S0 = θb.S0; dS0 = θb.dS0; d2S0 = θb.d2S0
        deltaLd1 = m / sin0 - a * omega * sin0 + cos0 / sin0
        deltaLd2 = m / sin0 - a * omega * sin0 + 2 * cos0 / sin0
        ddeltaLd2 = -m * cos0 / sin0^2 - a * omega * cos0 - 2 / sin0^2
        rho = -1 / (r0 - im * a * cos0); rhobar = -1 / (r0 + im * a * cos0)
        drhoOverrho = im * a * rho * sin0
        d2rhoOverrho = im * a * rho * (cos0 + 2 * sin0 * drhoOverrho)
        All0 = -(1 / 2) * rho^(-1) * rhobar *
               (d2S0 + (deltaLd1 + deltaLd2 + 2 * drhoOverrho) * dS0 +
                (ddeltaLd2 + deltaLd1 * deltaLd2 - 6 * drhoOverrho^2 + 3 * d2rhoOverrho +
                 (3 * deltaLd1 - deltaLd2) * drhoOverrho) * S0)
        Alm0 = (2 / sq2) * rho^(-1) *
               (-(rho + rhobar + im * Kt / Δ) * (dS0 + deltaLd2 * S0) +
                (rho - rhobar) * a * sin0 * Kt / Δ * S0)
        Alm1 = (2 / sq2) * rho^(-1) * ((dS0 + deltaLd2 * S0) + im * (rho - rhobar) * a * sin0 * S0)
        Amm0 = (Kt^2 * S0) / (Δ^2 * rho * rhobar) +
               (2 * im * Kt * S0 * (-1 + r0 - Δ * rho)) / (Δ^2 * rho * rhobar) -
               (2 * im * r0 * S0 * omega) / (Δ * rho * rhobar)
        Amm1 = 2 * rho^(-1) * rhobar^(-1) * S0 * (im * Kt / Δ + rho)
        Amm2 = -rho^(-1) * rhobar^(-1) * S0
        rcomp(ur0) = (En * (r0^2 + a^2) - a * Lz - ur0) / Δ
        θcomp(uθ0) = -rhobar * (im * sin0 * (a * En - Lz / sin0^2) - uθ0) / sq2
        proj(rc, θc) = (All0 * rc^2 + Alm0 * rc * θc + Amm0 * θc^2) * rb.D2R -
                       (Alm1 * rc * θc + Amm1 * θc^2) * rb.dD2R + Amm2 * θc^2 * rb.d2D2R
        cr = rb.cr; cθ = θb.cθ
        return proj(rcomp(rb.urp), θcomp(θb.uθp)) * cr * cθ +
               proj(rcomp(rb.urp), θcomp(θb.uθm)) * cr * conj(cθ) +
               proj(rcomp(rb.urm), θcomp(θb.uθp)) * conj(cr) * cθ +
               proj(rcomp(rb.urm), θcomp(θb.uθm)) * conj(cr) * conj(cθ)
    end
    αIn = _trapz2d_romberg_sep(qr -> rblock(qr, R.In), θblock, combine, Nr, Nθ) / (2 * pi)^2
    αUp = _trapz2d_romberg_sep(qr -> rblock(qr, R.Up), θblock, combine, Nr, Nθ) / (2 * pi)^2
    Zinf = -8 * pi * αIn / W / Upst
    ZH = -8 * pi * αUp / W / Upst
    return (Zinf=Zinf, ZH=ZH, Nr=Nr, Nθ=Nθ)
end

"Z^∞/Z^H for a GENERIC s=∓1 MAXWELL point particle (double Mino integral)."
function _convolve_source_generic_maxwell(s::Integer, l::Integer, m::Integer, n::Integer,
        k::Integer, a::Real, p::Real, e::Real, x::Real, En::Real, Lz::Real, Upst::Real,
        omega::Real, lambda::Real, rad, pol; Nmax::Integer=80, Nr::Integer=64, Nθ::Integer=64)
    sq2 = sqrt(2)
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    W = _teuk_wronskian(R.In)
    function θblock(qθ)
        dθ_ = pol(qθ); theta0 = dθ_.θ; sin0 = sin(theta0); cos0 = cos(theta0)
        S0 = Sc(theta0, 0.0)
        dS0 = swsh_derivative(Sc, 1, 0, theta0, 0.0)
        uθp = -dθ_.Uz / sin0; uθm = -pol(2 * pi - qθ).Uz / sin0
        cθ = cis(omega * dθ_.Δtz - m * dθ_.Δφz + k * qθ)
        return (sin0=sin0, cos0=cos0, S0=S0, dS0=dS0, uθp=uθp, uθm=uθm, cθ=cθ)
    end
    function rblock(qr, RF)
        dr_ = rad(qr); r0 = dr_.r
        R0 = RF(r0); dR0 = teukolsky_radial_deriv(RF, r0)
        Δ = r0^2 + a^2 - 2 * r0; Kt = (r0^2 + a^2) * omega - m * a
        urp = dr_.Ur; urm = rad(2 * pi - qr).Ur
        cr = cis(omega * dr_.Δtr - m * dr_.Δφr + n * qr)
        return (r0=r0, R0=R0, dR0=dR0, Δ=Δ, Kt=Kt, urp=urp, urm=urm, cr=cr)
    end
    function combine(rb, θb)
        r0 = rb.r0; Δ = rb.Δ; Kt = rb.Kt; sin0 = θb.sin0; cos0 = θb.cos0; csc0 = 1 / sin0
        S0 = θb.S0; dS0 = θb.dS0
        L1 = (s == -1 ? -1 : 1) * m / sin0 + (s == -1 ? 1 : -1) * a * omega * sin0 + cos0 / sin0
        rho = -1 / (r0 - im * a * cos0); rhobar = -1 / (r0 + im * a * cos0)
        Sigma = 1 / (rho * rhobar)
        local Ar0, Am0, Am1, rcomp, θcomp
        if s == -1
            Ar0 = -((dS0 + L1 * S0 + im * a * S0 * rho * sin0) / (2 * sq2 * Δ * rho^2 * rhobar))
            Am0 = (S0 * (-(im * Kt / Δ) + rho)) / (4 * rho^2)
            Am1 = -(S0 / (4 * rho^2))
            rcomp = ur0 -> (ur0 + (a^2 + r0^2) * En - a * Lz) / (2 * Sigma)
            θcomp = uθ0 -> -(rho * (-uθ0 + im * Lz * csc0 - im * a * En * sin0)) / sq2
        else  # s == +1
            Δp = 2 * (r0 - 1)
            Ar0 = -((Δ * (dS0 + L1 * S0 + im * a * S0 * rho * sin0)) / (2 * sq2 * rho))
            Am0 = -((S0 * (im * Kt + Δp + Δ * rho)) / (2 * rho * rhobar))
            Am1 = (S0 * Δ) / (2 * rho * rhobar)
            rcomp = ur0 -> (ur0 - (a^2 + r0^2) * En + a * Lz) / Δ
            θcomp = uθ0 -> (rhobar * (-uθ0 - im * Lz * csc0 + im * a * En * sin0)) / sq2
        end
        proj = (rc, θc) -> (Ar0 * rc + Am0 * θc) * rb.R0 - (Am1 * θc) * rb.dR0
        cr = rb.cr; cθ = θb.cθ
        return proj(rcomp(rb.urp), θcomp(θb.uθp)) * cr * cθ +
               proj(rcomp(rb.urp), θcomp(θb.uθm)) * cr * conj(cθ) +
               proj(rcomp(rb.urm), θcomp(θb.uθp)) * conj(cr) * cθ +
               proj(rcomp(rb.urm), θcomp(θb.uθm)) * conj(cr) * conj(cθ)
    end
    αIn = _trapz2d_romberg_sep(qr -> rblock(qr, R.In), θblock, combine, Nr, Nθ) / (2 * pi)^2
    αUp = _trapz2d_romberg_sep(qr -> rblock(qr, R.Up), θblock, combine, Nr, Nθ) / (2 * pi)^2
    Zinf = 8 * pi * αIn / W / Upst
    ZH = 8 * pi * αUp / W / Upst
    return (Zinf=Zinf, ZH=ZH, Nr=Nr, Nθ=Nθ)
end

"Z^∞/Z^H for a GENERIC s=0 SCALAR point particle (separable double Mino integral)."
function _convolve_source_generic_s0(l::Integer, m::Integer, n::Integer, k::Integer, a::Real,
        p::Real, e::Real, x::Real, Upst::Real, omega::Real, lambda::Real, rad, pol;
        Nmax::Integer=80, Nr::Integer=64, Nθ::Integer=64)
    s = 0
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, a * omega; numterms=Nmax > 60 ? :auto : Nmax)
    nu = RenormalizedAngularMomentum(s, l, m, a, omega, lambda)
    R = TeukolskyRadial(s, l, m, a, omega; lambda=lambda, nu=nu, Nmax=Nmax)
    W = _teuk_wronskian(R.In)
    II1(qr, RF) = begin
        d = rad(qr); r0 = d.r
        2 * r0^2 * RF(r0) * cos(omega * d.Δtr - m * d.Δφr + n * qr)
    end
    II3(qr, RF) = begin
        d = rad(qr); r0 = d.r
        2 * RF(r0) * cos(omega * d.Δtr - m * d.Δφr + n * qr)
    end
    α1In = _trapz_romberg(qr -> II1(qr, R.In), Nr) / (2 * pi)
    α1Up = _trapz_romberg(qr -> II1(qr, R.Up), Nr) / (2 * pi)
    II2(qθ) = begin
        d = pol(qθ)
        -8 * pi * Sc(d.θ, 0.0) * cos(omega * d.Δtz - m * d.Δφz + k * qθ)
    end
    α2 = _trapz_romberg(II2, Nθ) / (2 * pi)
    if a == 0
        α3In = zero(ComplexF64); α3Up = zero(ComplexF64); α4 = zero(ComplexF64)
    else
        α3In = _trapz_romberg(qr -> II3(qr, R.In), Nr) / (2 * pi)
        α3Up = _trapz_romberg(qr -> II3(qr, R.Up), Nr) / (2 * pi)
        II4(qθ) = begin
            d = pol(qθ)
            -8 * pi * a^2 * cos(d.θ)^2 * Sc(d.θ, 0.0) * cos(omega * d.Δtz - m * d.Δφz + k * qθ)
        end
        α4 = _trapz_romberg(II4, Nθ) / (2 * pi)
    end
    Zinf = (α1In * α2 + α3In * α4) / W / Upst
    ZH = (α1Up * α2 + α3Up * α4) / W / Upst
    return (Zinf=Zinf, ZH=ZH, Nr=Nr, Nθ=Nθ)
end

"Dispatch the generic (inclined eccentric) source convolution on spin weight `s`."
function _convolve_source_generic(s::Integer, l::Integer, m::Integer, n::Integer, k::Integer,
        a::Real, p::Real, e::Real, x::Real, En::Real, Lz::Real, Upst::Real,
        omega::Real, lambda::Real, rad, pol; Nmax::Integer=80, Nr::Integer=64, Nθ::Integer=64)
    s == -2 && return _convolve_source_generic_s2(l, m, n, k, a, p, e, x, En, Lz, Upst, omega, lambda, rad, pol; Nmax=Nmax, Nr=Nr, Nθ=Nθ)
    s == 2 && return _convolve_source_generic_s2p(l, m, n, k, a, p, e, x, En, Lz, Upst, omega, lambda, rad, pol; Nmax=Nmax, Nr=Nr, Nθ=Nθ)
    s == 0 && return _convolve_source_generic_s0(l, m, n, k, a, p, e, x, Upst, omega, lambda, rad, pol; Nmax=Nmax, Nr=Nr, Nθ=Nθ)
    return _convolve_source_generic_maxwell(s, l, m, n, k, a, p, e, x, En, Lz, Upst, omega, lambda, rad, pol; Nmax=Nmax, Nr=Nr, Nθ=Nθ)
end

# ============================================================================
# Public point-particle mode entrypoints (orbit scalars from Spacetime, λ from Harmonics).
# ============================================================================

_check_spin(s) = s in (-2, -1, 0, 1, 2) ||
    throw(ArgumentError("TeukolskyPointParticleMode: only s∈{−2,−1,0,+1,+2} are built (got $s)."))

"""
    TeukolskyPointParticleModeCircular(s, l, m, a, p; Nmax=80) -> NamedTuple

Full (l,m) Teukolsky mode of a circular, equatorial point particle on a Kerr geodesic (M=1), built
end-to-end from `(a, p)`. Orbit scalars (E, Lz, Υt, Ω_φ) from `Spacetime`, ω=mΩ_φ, λ from the
spin-weighted spheroidal eigenvalue, source convolution → Z^∞/Z^H, then the energy/L fluxes.
s∈{−2,−1,0,+1,+2}.
"""
function TeukolskyPointParticleModeCircular(s::Integer, l::Integer, m::Integer,
        a::Real, p::Real; Nmax::Integer=80)
    _check_spin(s)
    En = geodesic_energy(a, p, 0.0, 1.0)
    Lz = geodesic_angular_momentum(a, p, 0.0, 1.0)
    Omega = kerr_frequencies(a, p, 0.0, 1.0).Omegaph
    Upst = kerr_frequencies(a, p, 0.0, 1.0; Time="Mino").Ut
    omega = m * Omega
    lambda = SpinWeightedSpheroidalEigenvalue(s, l, m, a * omega)
    Z = _convolve_source_circular(s, l, m, a, p, En, Lz, Upst, omega, lambda; Nmax=Nmax)
    ef = energy_flux(s, l, m, a, omega, lambda, Z.Zinf, Z.ZH)
    jf = angular_momentum_flux(s, l, m, a, omega, lambda, Z.Zinf, Z.ZH)
    return (s=s, l=l, m=m, a=float(a), p=float(p),
            En=En, Lz=Lz, Upst=Upst, Omega=Omega, omega=omega, lambda=lambda,
            Zinf=Z.Zinf, ZH=Z.ZH, EnergyFlux=ef, AngularMomentumFlux=jf)
end

"""
    TeukolskyPointParticleModeEccentric(l, m, n, a, p, e; s=-2, Nmax=80, N=64) -> NamedTuple

Full (l,m,n) Teukolsky mode of an eccentric, equatorial point particle (M=1, e>0, |x|=1, k=0),
end-to-end from `(a, p, e)`. ω=mΩ_φ+nΩ_r; the source is convolved over the radial libration.
s∈{−2,−1,0,+1,+2}.
"""
function TeukolskyPointParticleModeEccentric(l::Integer, m::Integer, n::Integer,
        a::Real, p::Real, e::Real; s::Integer=-2, Nmax::Integer=80, N::Integer=64)
    _check_spin(s)
    e > 0 || throw(ArgumentError("TeukolskyPointParticleModeEccentric requires e>0 (use Circular for e=0)"))
    En = geodesic_energy(a, p, e, 1.0)
    Lz = geodesic_angular_momentum(a, p, e, 1.0)
    bl = kerr_frequencies(a, p, e, 1.0)
    Omegar = bl.Omegar; Omegaph = bl.Omegaph
    mino = kerr_frequencies(a, p, e, 1.0; Time="Mino")
    Upsr = mino.Ur; Upst = mino.Ut
    omega = m * Omegaph + n * Omegar
    lambda = SpinWeightedSpheroidalEigenvalue(s, l, m, a * omega)
    rad = geodesic_trajectory_radial(a, p, e, 1.0)
    Z = _convolve_source_eccentric(s, l, m, n, a, p, e, En, Lz, Upst, Upsr, omega, lambda, rad; Nmax=Nmax, N=N)
    ef = energy_flux(s, l, m, a, omega, lambda, Z.Zinf, Z.ZH)
    jf = angular_momentum_flux(s, l, m, a, omega, lambda, Z.Zinf, Z.ZH)
    return (s=s, l=l, m=m, n=n, a=float(a), p=float(p), e=float(e),
            En=En, Lz=Lz, Upst=Upst, Upsr=Upsr, Omegar=Omegar, Omegaph=Omegaph,
            omega=omega, lambda=lambda, Zinf=Z.Zinf, ZH=Z.ZH,
            EnergyFlux=ef, AngularMomentumFlux=jf, N=Z.N)
end

"""
    TeukolskyPointParticleModeSpherical(l, m, k, a, p, x; s=-2, Nmax=80, N=64) -> NamedTuple

Full (l,m,k) Teukolsky mode of a spherical (e=0, inclined |x|<1) point particle (M=1, n=0),
end-to-end from `(a, p, x)`. ω=mΩ_φ+kΩ_θ; the source is convolved over the polar libration.
s∈{−2,−1,0,+1,+2}.
"""
function TeukolskyPointParticleModeSpherical(l::Integer, m::Integer, k::Integer,
        a::Real, p::Real, x::Real; s::Integer=-2, Nmax::Integer=80, N::Integer=64)
    _check_spin(s)
    abs(x) < 1 || throw(ArgumentError("TeukolskyPointParticleModeSpherical requires |x|<1 (use Circular for |x|=1)"))
    En = geodesic_energy(a, p, 0.0, x)
    Lz = geodesic_angular_momentum(a, p, 0.0, x)
    bl = kerr_frequencies(a, p, 0.0, x)
    Omegath = bl.Omegath; Omegaph = bl.Omegaph
    Upst = kerr_frequencies(a, p, 0.0, x; Time="Mino").Ut
    omega = m * Omegaph + k * Omegath
    lambda = SpinWeightedSpheroidalEigenvalue(s, l, m, a * omega)
    pol = geodesic_trajectory_polar(a, p, 0.0, x)
    Z = _convolve_source_spherical(s, l, m, k, a, p, En, Lz, Upst, omega, lambda, pol; Nmax=Nmax, N=N)
    ef = energy_flux(s, l, m, a, omega, lambda, Z.Zinf, Z.ZH)
    jf = angular_momentum_flux(s, l, m, a, omega, lambda, Z.Zinf, Z.ZH)
    return (s=s, l=l, m=m, k=k, a=float(a), p=float(p), x=float(x),
            En=En, Lz=Lz, Upst=Upst, Omegath=Omegath, Omegaph=Omegaph,
            omega=omega, lambda=lambda, Zinf=Z.Zinf, ZH=Z.ZH,
            EnergyFlux=ef, AngularMomentumFlux=jf, N=Z.N)
end

"""
    TeukolskyPointParticleModeGeneric(l, m, n, k, a, p, e, x; s=-2, Nmax=80, Nr=64, Nθ=64) -> NamedTuple

Full (l,m,n,k) Teukolsky mode of a generic (inclined eccentric: e>0 AND |x|<1) point particle
(M=1), end-to-end from `(a, p, e, x)`. ω=mΩ_φ+nΩ_r+kΩ_θ; the source is convolved over BOTH the
radial and polar librations. s∈{−2,−1,0,+1,+2}.
"""
function TeukolskyPointParticleModeGeneric(l::Integer, m::Integer, n::Integer, k::Integer,
        a::Real, p::Real, e::Real, x::Real; s::Integer=-2,
        Nmax::Integer=80, Nr::Integer=64, Nθ::Integer=64)
    _check_spin(s)
    e > 0 || throw(ArgumentError("TeukolskyPointParticleModeGeneric requires e>0 (use Spherical for e=0)"))
    abs(x) < 1 || throw(ArgumentError("TeukolskyPointParticleModeGeneric requires |x|<1 (use Eccentric for |x|=1)"))
    En = geodesic_energy(a, p, e, x)
    Lz = geodesic_angular_momentum(a, p, e, x)
    bl = kerr_frequencies(a, p, e, x)
    Omegar = bl.Omegar; Omegath = bl.Omegath; Omegaph = bl.Omegaph
    Upst = kerr_frequencies(a, p, e, x; Time="Mino").Ut
    omega = m * Omegaph + n * Omegar + k * Omegath
    lambda = SpinWeightedSpheroidalEigenvalue(s, l, m, a * omega)
    rad = geodesic_trajectory_radial(a, p, e, x)
    pol = geodesic_trajectory_polar(a, p, e, x)
    Z = _convolve_source_generic(s, l, m, n, k, a, p, e, x, En, Lz, Upst, omega, lambda, rad, pol; Nmax=Nmax, Nr=Nr, Nθ=Nθ)
    ef = energy_flux(s, l, m, a, omega, lambda, Z.Zinf, Z.ZH)
    jf = angular_momentum_flux(s, l, m, a, omega, lambda, Z.Zinf, Z.ZH)
    return (s=s, l=l, m=m, n=n, k=k, a=float(a), p=float(p), e=float(e), x=float(x),
            En=En, Lz=Lz, Upst=Upst, Omegar=Omegar, Omegath=Omegath, Omegaph=Omegaph,
            omega=omega, lambda=lambda, Zinf=Z.Zinf, ZH=Z.ZH,
            EnergyFlux=ef, AngularMomentumFlux=jf, Nr=Z.Nr, Nθ=Z.Nθ)
end
