# DG-P5: SOURCED generic-spin circular-equatorial point-particle sector at a = 0
# (s ∈ {0, −2, +2}) for the RT-chart DG evolution — the sourced companion of the DG-P4
# vacuum machinery (spinweight.jl / rtchart.jl / spinmodes.jl). The distributional
# δ/δ′/δ″ point-particle source enters the DG scheme PURELY through the analytically
# known time-dependent field jumps at the particle interface (the jump-conditioned
# numerical-flux treatment of δ-derivative sources, [arXiv:2211.14390, §2.1–§3.2]:
# the field ansatz ψ = ψ⁺Θ + ψ⁻Θ̄ + ψ_δ·δ makes the smooth pieces homogeneous solutions
# stitched by closed-form jumps; the δ-layer ψ_δ is pinned at r₀ and never enters
# element interiors). At s = 0 this reduces EXACTLY to the locked Flux1D δ-injection
# (gate S0 below, closed-form identity).
#
# PROVENANCE CHAIN (every step gated in test_timedomain_spinsource.jl; measured
# agreements quoted there):
#
#  (1) Radial FD source. The circular-equatorial point-particle source of the radial
#      Teukolsky equation, in the form used by the in-package FD oracle
#      (AmplitudesFluxes/pointparticle.jl `_convolve_source_circular_{s2,s2p,s0}` —
#      the GHP projection of Hughes [arXiv:gr-qc/9910091, Eqs. (4.39)–(4.40)]),
#      specialized to a = 0, θ₀ = π/2:  the self-adjoint-form source is
#         𝔗 = (−8π/Υt)·Δ^{w_s}·(a₀ δ + a₁ δ′ + a₂ δ″),  δ = δ(r − r₀),
#      with ∫R_in·𝔗-pairings equal to the oracle's α_In by construction. At a = 0 the
#      spheroidal harmonics are EXACTLY the spin-weighted spherical ₛYₗₘ (no basis
#      re-projection; the Kerr a ≠ 0 source in the ₛYₗₘ basis needs the
#      spheroidal→spherical re-expansion not in print — dossier DG_SPIN_FORMULATION.md
#      §8.3 — and is an explicit [FUTURE]; constructors below enforce a = 0).
#
#  (2) Jump algebra (DERIVED-J). For X″ + p(r)X′ + q(r)X = S₀δ + S₁δ′ + S₂δ″ the ansatz
#      X = X⁺Θ(r−r₀) + X⁻Θ(r₀−r) + X_δ·δ(r−r₀) and the distributional identities
#      F(r)δ′ = F₀δ′ − F₀′δ, F(r)δ″ = F₀δ″ − 2F₀′δ′ + F₀″δ [arXiv:1702.04204, App. A,
#      display after Eq. (A13)] give, by matching δ″/δ′/δ coefficients:
#         X_δ  = S₂
#         [X]  = S₁ − p₀S₂
#         [X′] = S₀ − p₀[X] + (p₀′ − q₀)S₂
#      (p, q from the radial Teukolsky ODE, ≡ AmplitudesFluxes `_d2R_ode`). GATED: with
#      the s = +2 source of (1) this reproduces the INDEPENDENT Barack–Giudice closed
#      forms ψ_δ, [ψ], [ψ′] [arXiv:1702.04204, Eqs. (A15)–(A17)] to ≤ 5e-16 (test gate
#      J2), and the FD Green's-function jumps Z^∞R̃_up − Z^H R̃_in to ≤ 1e-12 (gate J3).
#
#  (3) Variable map. ψ̂ = rΔˢX (the Barack–Giudice / RT mode variable; at a = 0
#      ψ_RT,ℓm = φ_BG,sℓm — rtchart.jl header):  [ψ̂] = (rΔˢ)₀[X],
#      [ψ̂′_r] = (rΔˢ)′₀[X] + (rΔˢ)₀[X′].
#
#  (4) Chart push-forward (DERIVED-C, chain rule through the RT chart maps
#      [arXiv:1301.1591, Eqs. (3),(5),(6)] at a = 0: t = T + F(R), r = r(R) with
#      F = rt_F_of_R, r = rt_r_of_R). The steady state is monochromatic,
#      ψ(t,r) = ψ̂(r)e^{−iωt}, ω = mΩ, so at fixed RT time T:
#         J_ψ(T) = [ψ̂]·e^{−iω(T+F(R₀))}
#         J_π(T) = −∂_T J_ψ = iω·J_ψ(T)
#         J_χ(T) = ([ψ̂′_r]·(dr/dR)₀ − iωF′(R₀)[ψ̂])·e^{−iω(T+F(R₀))}
#      with dr/dR = 2(1+R²)/(1−R²)² and F′(R) = 8MR/(1−R²) + 4R/(1−R²)² + (1−1/f)dr/dR
#      (analytic derivatives of the chart maps; gated against central differences).
#
#  (5) DG interface enforcement. At the particle interface k_p (R₀ a forced mesh edge)
#      the LLF numerical flux and the ψ-jump penalty are evaluated on the SHIFTED right
#      state Ũ_R := U_R − J(T): for a linear flux F = A·U the left-face and right-face
#      modified fluxes F*_L = F*(U_L, Ũ_R), F*_R = F*_L + A·J differ by exactly A·J, so
#      BOTH face penalty terms equal the vacuum formulas applied to (U_L, Ũ_R) — the
#      A·J contributions cancel in F(U_R) − F*_R. On the exact solution (U_R − U_L = J)
#      all penalties vanish; consistency is preserved. This generalizes the s = 0
#      characteristic δ-flux injection [arXiv:2307.01349, Eqs. (42a,b)]: at s = 0 the
#      only nonzero jump is [χ] = G(t) (gate S0: [ψ̂′_r]·f₀ ≡ Flux1D's Gpre, closed form).
#      The jumps are ramped from zero with the shared smooth ramp `_ramp` (common.jl),
#      exactly like the Flux1D source.
#
#  (6) Flux extraction at ℐ⁺/H⁺. With unit-transmission In/Up normalization (Radial),
#      X → Z^∞ r^{−(2s+1)}e^{iωr*} at r → ∞ and X → Z^H Δ^{−s}e^{−iωr*} at r → r₊, so
#      the RT variable obeys |ψ̂| → |Z^∞| at ℐ⁺ and |ψ̂| → r₊|Z^H| at H⁺ (ANY s; the
#      rΔˢ weight cancels the asymptotics exactly — this is what the chart buys).
#      Steady state: |π| = ω|ψ̂|, so |Z^∞| = |π(R=1)|/ω, |Z^H| = |π(R₊)|/(ω r₊).
#      Ė^∞/Ė^H then follow from the in-package amplitude→flux layer
#      (AmplitudesFluxes/fluxes.jl `energy_flux`, BHPToolkit normalization —
#      equivalently Ė^∞ = lim r²/(4π)∮|∫ψ₄dt|²dΩ [arXiv:gr-qc/0703028, Eq. (2.3)]
#      for s = −2).
#
# The vacuum DG-P4 paths are untouched: this file only ADDS a problem type + RHS methods
# (`SpinModesProblem` gained the behavior-preserving `forced_edges` hook).

# ---- (1) GHP source coefficients at a = 0 --------------------------------------------------------

"ₛYₗₘ(π/2), ∂θ, ∂θθ values via the SWSH surface at aω = 0 (EXACT spherical limit; same
convention as the FD oracle and as Barack–Giudice's ₂Yₗₘ — gated vs closed forms)."
function _swsh_equator(s::Integer, l::Integer, m::Integer)
    Sc = SpinWeightedSpheroidalHarmonicS(s, l, m, 0.0)
    return (Sc(pi / 2, 0.0), swsh_derivative(Sc, 1, 0, pi / 2, 0.0),
            swsh_derivative(Sc, 2, 0, pi / 2, 0.0))
end

"""
    _ghp_source_a0(s, l, m, r0, En, Lz, omega) -> (a0, a1, a2, w)

Circular-equatorial GHP source triplet at a = 0 (M = 1), in the oracle's pairing form
`α_In = ∫ R̃_in·(a₀δ + a₁δ′ + a₂δ″) dr` with R̃ = R (s = −2, 0) or R̃ = Δ²R (s = +2, the
oracle's Chandrasekhar-rescaled projection). Equivalently, the source of the NORMALIZED
radial ODE X″ + pX′ + qX (i.e. 𝔗_selfadjoint/Δ^{s+1}) is
`(−8π/Υt)·Δ(r)^w·(a₀δ + a₁δ′ + a₂δ″)` with w = −(s+1) (s ≤ 0) or w = 2 − (s+1)
(s = +2), expanded distributionally by [`_ode_source_Sn`](@ref). Transcribed from the
in-package FD oracle `_convolve_source_circular_{s2,s2p,s0}` (AmplitudesFluxes/
pointparticle.jl; GHP projection of [arXiv:gr-qc/9910091, Eqs. (4.39)–(4.40)]),
specialized to a = 0, θ₀ = π/2 (cosθ₀ = 0 terms dropped; ρ = ρ̄ = −1/r₀). s = 0 uses the
scalar α = −4π r₀² S normalization rewritten as a₀ = r₀²S/2 under the common −8π
prefactor. GATED against the oracle's Z^∞/Z^H through the Green's-function jump identity
(gate J3) and against the independent Barack–Giudice closed forms (gate J2).
"""
function _ghp_source_a0(s::Integer, l::Integer, m::Integer, r0, En, Lz, omega)
    S0, dS0, d2S0 = _swsh_equator(s, l, m)
    rho = -1 / r0
    Δ0 = r0^2 - 2r0
    Kt = r0^2 * omega
    sq2 = sqrt(2.0)
    if s == -2
        # oracle `_convolve_source_circular_s2` body at a = 0 (L1† = L2† = −m at θ₀=π/2)
        L2S = dS0 - m * S0
        L1L2S = d2S0 - 2m * dS0 + (m^2 - 2) * S0
        Ann0 = -rho^(-2) * rho^(-1) * (sq2 * Δ0)^(-2) * (rho^(-1) * L1L2S)
        Anmbar0 = rho^(-3) * (sq2 * Δ0)^(-1) * ((2rho - im * Kt / Δ0) * L2S)
        Anmbar1 = -rho^(-3) * (sq2 * Δ0)^(-1) * L2S
        Ambarmbar0 = (Kt^2 * S0 * rho) / (4 * Δ0^2 * rho^3) +
                     (im * Kt * S0 * (1 - r0 + Δ0 * rho) * rho) / (2 * Δ0^2 * rho^3) +
                     (im * r0 * S0 * rho * omega) / (2 * Δ0 * rho^3)
        Ambarmbar1 = -rho^(-3) * rho * S0 / 2 * (im * Kt / Δ0 - rho)
        Ambarmbar2 = -rho^(-3) * rho * S0 / 4
        rcomp = En / 2                       # (En(r₀²+a²) − aLz)/(2Σ), Σ = r₀²
        thetacomp = rho * (-im * Lz) / sq2
        Cnn = rcomp^2; Cnm = rcomp * thetacomp; Cmm = thetacomp^2
        return (Ann0 * Cnn + Anmbar0 * Cnm + Ambarmbar0 * Cmm,
                Anmbar1 * Cnm + Ambarmbar1 * Cmm,
                Ambarmbar2 * Cmm, -(s + 1))
    elseif s == 2
        # oracle `_convolve_source_circular_s2p` body at a = 0 (δL₁ = δL₂ = m at θ₀=π/2)
        L2S = dS0 + m * S0
        All0 = -(1 / 2) * rho^(-1) * rho * (d2S0 + 2m * dS0 + (m^2 - 2) * S0)
        Alm0 = (2 / sq2) * rho^(-1) * (-(2rho + im * Kt / Δ0) * L2S)
        Alm1 = (2 / sq2) * rho^(-1) * L2S
        Amm0 = (Kt^2 * S0) / (Δ0^2 * rho^2) +
               (2 * im * Kt * S0 * (-1 + r0 - Δ0 * rho)) / (Δ0^2 * rho^2) -
               (2 * im * r0 * S0 * omega) / (Δ0 * rho^2)
        Amm1 = 2 * rho^(-2) * S0 * (im * Kt / Δ0 + rho)
        Amm2 = -rho^(-2) * S0
        rcomp = En * r0^2 / Δ0
        thetacomp = -rho * (-im * Lz) / sq2
        Cll = rcomp^2; Clm = rcomp * thetacomp; Cmm = thetacomp^2
        return (All0 * Cll + Alm0 * Clm + Amm0 * Cmm,
                Alm1 * Clm + Amm1 * Cmm,
                Amm2 * Cmm, 2 - (s + 1))
    elseif s == 0
        # oracle `_convolve_source_circular_s0`: α = −4π r₀² S ⇒ under the −8π prefactor
        # and w = −1 wrap: g₀a₀ = (−8π/Υt)Δ₀⁻¹·(r₀²S/2) = −4πr₀²S/(Υt·Δ₀). δ-only source.
        return (complex(r0^2 * S0 / 2), 0.0im, 0.0im, -1)
    else
        throw(ArgumentError("spin source: s ∈ {0, −2, +2} built (|s| = 1 is [FUTURE]; got s=$s)"))
    end
end

# ---- (2) jump algebra ----------------------------------------------------------------------------

"""
    _ode_source_Sn(a0, a1, a2, w, r0, Upst) -> (S0, S1, S2)

Distributional division: rewrite 𝔗/Δ^{s+1} = (−8π/Υt)·g(r)·(a₀δ + a₁δ′ + a₂δ″),
g(r) = Δ(r)^w, as S₀δ + S₁δ′ + S₂δ″ using g·δ′ = g₀δ′ − g₀′δ, g·δ″ = g₀δ″ − 2g₀′δ′ + g₀″δ
[arXiv:1702.04204, App. A distributional identities]:
S₂ = g₀a₂, S₁ = g₀a₁ − 2g₀′a₂, S₀ = g₀a₀ − g₀′a₁ + g₀″a₂ (all × −8π/Υt). M = 1.
"""
function _ode_source_Sn(a0, a1, a2, w::Integer, r0, Upst)
    pref = -8pi / Upst
    Δ0 = r0^2 - 2r0
    dΔ0 = 2r0 - 2
    g0 = pref * Δ0^w
    g1 = pref * w * Δ0^(w - 1) * dΔ0
    g2 = pref * (w * (w - 1) * Δ0^(w - 2) * dΔ0^2 + 2w * Δ0^(w - 1))
    return (g0 * a0 - g1 * a1 + g2 * a2, g0 * a1 - 2 * g1 * a2, g0 * a2)
end

"""
    _ode_jumps(S0, S1, S2, p0, dp0, q0) -> (Xδ, JX, JXp)

DERIVED-J (header (2)): jumps of X″ + pX′ + qX = S₀δ + S₁δ′ + S₂δ″ across r₀:
X_δ = S₂, [X] = S₁ − p₀S₂, [X′] = S₀ − p₀[X] + (p₀′ − q₀)S₂. Gated to ≤ 5e-16 against
the independent Barack–Giudice closed forms [arXiv:1702.04204, Eqs. (A15)–(A17)] (s=+2).
"""
_ode_jumps(S0, S1, S2, p0, dp0, q0) =
    (Xδ = S2, JX = S1 - p0 * S2, JXp = S0 - p0 * (S1 - p0 * S2) + (dp0 - q0) * S2)

"radial-Teukolsky ODE coefficients at r₀ (a = 0, M = 1): p = (s+1)Δ′/Δ, p′, and
q = (−λ + 4isrω + (K² − 2is(r−1)K)/Δ)/Δ, K = r²ω (≡ AmplitudesFluxes `_d2R_ode`;
λ = ℓ(ℓ+1) − s(s+1), the SWSH eigenvalue at aω = 0)."
function _ode_pq_a0(s::Integer, l::Integer, m::Integer, r0, omega)
    Δ0 = r0^2 - 2r0
    dΔ0 = 2r0 - 2
    λ = l * (l + 1) - s * (s + 1)
    K = r0^2 * omega
    p0 = (s + 1) * dΔ0 / Δ0
    dp0 = (s + 1) * (2 / Δ0 - dΔ0^2 / Δ0^2)
    q0 = (-λ + 4im * s * r0 * omega + (K^2 - 2im * s * (r0 - 1) * K) / Δ0) / Δ0
    return (p0, dp0, q0)
end

# ---- (3)+(4) field jumps and the RT-chart jump state ---------------------------------------------

"""
    spin_source_jumps(s, l, m, r0) -> NamedTuple

Complete t = 0 jump data of the sourced mode (s,ℓ,m) at circular orbit radius `r0`
(a = 0, M = 1, μ_particle = 1; time dependence e^{−iωt}, ω = mΩ): the radial-function
jumps `(Xδ, JX, JXp)` (header (1)+(2)) and the BG/RT-variable ψ̂ = rΔˢX jumps
`(psid, Jpsi, Jpsip)` (header (3)). `m = 0` is excluded (static modes carry no flux and
the monochromatic jump ansatz needs ω ≠ 0).
"""
function spin_source_jumps(s::Integer, l::Integer, m::Integer, r0::Real)
    m == 0 && throw(ArgumentError("spin_source_jumps: m ≠ 0 required (static m = 0 is fluxless)"))
    r0 > 3 || throw(ArgumentError("spin_source_jumps: r0 > 3M (timelike circular geodesic) required"))
    En = geodesic_energy(0.0, r0, 0.0, 1.0)
    Lz = geodesic_angular_momentum(0.0, r0, 0.0, 1.0)
    f = kerr_frequencies(0.0, r0, 0.0, 1.0; Time="Mino")
    Upst = f.Ut
    Omega = f.Uph / f.Ut
    omega = m * Omega
    a0, a1, a2, w = _ghp_source_a0(s, l, m, r0, En, Lz, omega)
    S0, S1, S2 = _ode_source_Sn(a0, a1, a2, w, r0, Upst)
    p0, dp0, q0 = _ode_pq_a0(s, l, m, r0, omega)
    J = _ode_jumps(S0, S1, S2, p0, dp0, q0)
    Δ0 = r0^2 - 2r0
    dΔ0 = 2r0 - 2
    W0 = r0 * float(Δ0)^s                       # (rΔˢ)₀
    W1 = float(Δ0)^s + s * r0 * float(Δ0)^(s - 1) * dΔ0   # (rΔˢ)′₀
    return (s=s, l=l, m=m, r0=float(r0), omega=omega, Omega=Omega, Upst=Upst,
            Xδ=J.Xδ, JX=J.JX, JXp=J.JXp,
            psid=W0 * J.Xδ, Jpsi=W0 * J.JX, Jpsip=W1 * J.JX + W0 * J.JXp)
end

"""
    bg_jumps_s2(l, m, r0) -> (psid, Jpsi, Jpsip)

INDEPENDENT s = +2 closed-form jump data at t = 0 (a = 0, M = 1, μ = 1):
Barack–Giudice ψ_δ = −4f₀⁻²s₂ [arXiv:1702.04204, Eq. (A15)], [ψ] [Eq. (A16)],
[ψ′] [Eq. (A17)] with s₂ = μγ₀π·Mf₀³r₀²𝒴 [Eqs. (A6),(A9)], 𝒴 = ₂Ȳₗₘ(π/2,0),
m̃ = m(1−3M/r₀), λ = (ℓ+2)(ℓ−1). Validation oracle for the generic jump algebra
(gate J2) — NOT used by the evolution path.
"""
function bg_jumps_s2(l::Integer, m::Integer, r0::Real)
    f0 = 1 - 2 / r0
    y = 1 / r0
    γ0 = 1 / sqrt(1 - 3 / r0)
    mt = m * (1 - 3 / r0)
    Ω = sqrt(1 / r0^3)
    λ = (l + 2) * (l - 1)
    Y, Yt, Ytt = _swsh_equator(2, l, m)      # 𝒴(0), 𝒴_θ(0), 𝒴_θθ(0) — real at (π/2, 0)
    s2 = pi * γ0 * (f0^3 * r0^2 * Y)
    psid = -4 * f0^(-2) * s2
    Jpsi = 8pi * γ0 * r0^2 * ((y^2 - im * mt * Ω * r0) * Y - im * f0 * r0 * Ω * Yt)
    Jpsip = 4pi * γ0 * r0 * ((2y^2 - 2 - y * (λ - 4) + mt^2 / f0 * (1 + y^2 * γ0^4)
                              - 2im * mt * r0 * Ω / f0 * (3 - 7y)) * Y
                             + 2 * (-3im * f0 * r0 * Ω + mt) * Yt + f0 * Ytt)
    return (psid=psid, Jpsi=Jpsi, Jpsip=Jpsip)
end

# rt_R_of_r (R = (√(1+r²)−1)/r, the inverse of rt_r_of_R) is a PURE chart map — RE-HOMED to
# Spacetime.RTChart (F5 Part B); re-imported+re-exported by the TimeDomain module header (same object,
# ===). The `rt_R_of_r` call sites below (and TimeDomain's `export rt_R_of_r`) resolve to it.

"analytic F′(R) (a = 0): d(t−T)/dR through [arXiv:1301.1591, Eqs. (3),(5),(6)] —
F′ = 8MR/(1−R²) + 4R/(1−R²)² + (1 − 1/f)·dr/dR, dr/dR = 2(1+R²)/(1−R²)²."
function _rt_dF_dR(R; M=1.0)
    drdR = 2 * (1 + R^2) / (1 - R^2)^2
    f = 1 - 2M / rt_r_of_R(R)
    return 8M * R / (1 - R^2) + 4R / (1 - R^2)^2 + (1 - 1 / f) * drdR
end

"""
    _rt_jump_state(J) -> (R0, Jψ, Jπ, Jχ)

DERIVED-C (header (4)): the (ψ, π, χ) jump amplitudes at the particle interface in the
RT chart, such that the enforced jumps are Ĵ·e^{−iωT}·ramp(T). From the monochromatic
steady state ψ(T,R) = ψ̂(r(R))e^{−iω(T+F(R))}: Ĵψ = [ψ̂]e^{−iωF(R₀)}, Ĵπ = iωĴψ,
Ĵχ = ([ψ̂′_r]·(dr/dR)₀ − iωF′(R₀)[ψ̂])e^{−iωF(R₀)}.
"""
function _rt_jump_state(J)
    R0 = rt_R_of_r(J.r0)
    phase = cis(-J.omega * rt_F_of_R(R0; M=1.0))
    drdR = 2 * (1 + R0^2) / (1 - R0^2)^2
    Jψ = J.Jpsi * phase
    Jπ = im * J.omega * Jψ
    Jχ = (J.Jpsip * drdR - im * J.omega * _rt_dF_dR(R0) * J.Jpsi) * phase
    return (R0=R0, Jψ=Jψ, Jπ=Jπ, Jχ=Jχ)
end

# ---- (5) sourced problem + RHS -------------------------------------------------------------------

"""
    SpinSourcedProblem{T}

Sourced RT-chart circular-orbit problem (DG-P5): a vacuum [`SpinModesProblem`](@ref)
(with the particle radius R₀ forced onto the mesh) plus the per-mode jump-state
amplitudes at the particle interface and the source ramp. a = 0, M = 1, s ∈ {0,−2,+2}.
Build with [`SpinSourcedProblem(s, ℓs, m, r0; ...)`](@ref); evolve with
[`solve_spin_sourced`](@ref); fluxes with [`spin_sourced_flux`](@ref).
"""
struct SpinSourcedProblem{T<:AbstractFloat}
    base::SpinModesProblem{T}
    r0::T
    R0::T
    kp::Int                       # element left of the particle interface
    omega::T
    lambda::Vector{T}             # SWSH eigenvalue per mode (aω = 0)
    Jψ::Vector{Complex{T}}        # per-mode jump amplitudes (× e^{−iωT}·ramp)
    Jπ::Vector{Complex{T}}
    Jχ::Vector{Complex{T}}
    ramp_τ::T
    ramp_δ::T
end

"""
    SpinSourcedProblem(s, ℓs, m, r0; N=8, K=48, excise=0.05, filter_alpha=36.0,
                       filter_power=16, constraint_kappa=0.3, psi_jump_toff=250.0,
                       ramp_τ=400.0, ramp_δ=0.00025, extract_R=Float64[], T=Float64)

Build the a = 0 (M = 1) sourced problem for spin weight `s ∈ {0, −2, +2}`, azimuthal
number `m ≠ 0`, coupled-mode list `ℓs` and circular orbit radius `r0` (each ℓ in `ℓs`
carries its own jump data; at a = 0 the ℓ-coupling matrices are diagonal so the modes
evolve independently but share the run). Kerr a ≠ 0 is a documented [FUTURE] (basis
re-projection, header (1)). Numerical-dissipation kwargs pass through to the vacuum
constructor; `ramp_τ`/`ramp_δ` are the Flux1D-standard smooth source ramp.
"""
function SpinSourcedProblem(s::Integer, ℓs::Vector{Int}, m::Integer, r0::Real;
                            N=8, K=48, excise=0.05, filter_alpha=36.0, filter_power=16,
                            constraint_kappa=0.3, psi_jump_toff=250.0,
                            ramp_τ=400.0, ramp_δ=0.00025, extract_R=Float64[],
                            T::Type{<:AbstractFloat}=Float64)
    s in (0, -2, 2) || throw(ArgumentError("SpinSourcedProblem: s ∈ {0, −2, +2} (got $s)"))
    R0 = T(rt_R_of_r(r0))
    base = SpinModesProblem(s, ℓs, m; N=N, K=K, M=1.0, a=0.0, excise=excise,
                            filter_alpha=filter_alpha, filter_power=filter_power,
                            constraint_kappa=constraint_kappa, psi_jump_toff=psi_jump_toff,
                            extract_R=extract_R, forced_edges=[Float64(R0)], T=T)
    kp = findmin(abs.(base.mesh.edges .- R0))[2] - 1
    (kp >= 1 && kp <= base.mesh.K - 1) || error("particle interface not interior to the mesh")
    n = length(ℓs)
    Jψ = Vector{Complex{T}}(undef, n)
    Jπ = Vector{Complex{T}}(undef, n)
    Jχ = Vector{Complex{T}}(undef, n)
    λv = Vector{T}(undef, n)
    local ω
    for (q, l) in enumerate(ℓs)
        J = spin_source_jumps(s, l, m, r0)
        st = _rt_jump_state(J)
        Jψ[q] = st.Jψ; Jπ[q] = st.Jπ; Jχ[q] = st.Jχ
        λv[q] = l * (l + 1) - s * (s + 1)
        ω = J.omega
    end
    return SpinSourcedProblem{T}(base, T(r0), R0, kp, T(ω), λv, Jψ, Jπ, Jχ, T(ramp_τ), T(ramp_δ))
end

function _rhs!(du, u, t::Real, p::SpinSourcedProblem{T}, sc) where {T}
    b = p.base
    el = b.el
    Np = el.Np
    K = b.mesh.K
    n = b.n
    U = reshape(u, Np, K, 3, n)
    dU = reshape(du, Np, K, 3, n)
    Dπ, Dχ, Dψ, tmpπ, tmpχ, Fπ, Fχ, FπR, FχR = sc
    D = el.D
    A4 = b.Aflux
    Gpi4 = b.Gpi
    Gch4 = b.Gchi
    Gps4 = b.Gpsi
    κc = b.κc
    # volume terms — identical to the vacuum RHS (spinmodes.jl)
    @inbounds for k in 1:K
        invJ = b.twoOverh[k]
        for q in 1:n, i in 1:Np
            aπ = zero(Complex{T}); aχ = zero(Complex{T}); aψ = zero(Complex{T})
            for j in 1:Np
                aπ += D[i, j] * U[j, k, 2, q]
                aχ += D[i, j] * U[j, k, 3, q]
                aψ += D[i, j] * U[j, k, 1, q]
            end
            Dπ[i, q] = invJ * aπ; Dχ[i, q] = invJ * aχ; Dψ[i, q] = invJ * aψ
        end
        for i in 1:Np
            for q in 1:n
                s1 = zero(Complex{T})
                for c in 1:n
                    s1 += A4[q, c, i, k] * Dπ[i, c] + A4[q, n+c, i, k] * Dχ[i, c]
                end
                w = zero(Complex{T})
                for c in 1:n
                    w += -Gpi4[q, c, i, k] * U[i, k, 2, c] + Gch4[q, c, i, k] * U[i, k, 3, c] +
                         Gps4[q, c, i, k] * U[i, k, 1, c]
                end
                dU[i, k, 1, q] = -U[i, k, 2, q]
                dU[i, k, 2, q] = -s1 + w
                dU[i, k, 3, q] = -Dπ[i, q] - κc * (U[i, k, 3, q] - Dψ[i, q])
            end
        end
    end
    # interface fluxes: vacuum LLF everywhere; at k = kp the right state is SHIFTED by the
    # ramped jump J(T) (header (5)) and the ψ-jump penalty is always on with target J_ψ(T).
    rmp = _ramp(t, p.ramp_τ, p.ramp_δ)
    ph = cis(-p.omega * t) * rmp
    @inbounds for k in 1:(K-1)
        A = view(b.Aint, :, :, k)
        α = b.αint[k]
        src = (k == p.kp)
        for q in 1:n
            tmpπ[q] = U[Np, k, 2, q]
            tmpχ[q] = U[Np, k, 3, q]
        end
        _AUvec!(Fπ, Fχ, A, tmpπ, tmpχ, n)          # F(U^L)
        for q in 1:n
            tmpπ[q] = U[1, k+1, 2, q] - (src ? p.Jπ[q] * ph : zero(Complex{T}))
            tmpχ[q] = U[1, k+1, 3, q] - (src ? p.Jχ[q] * ph : zero(Complex{T}))
        end
        _AUvec!(FπR, FχR, A, tmpπ, tmpχ, n)        # F(Ũ^R)
        ck = b.twoOverh[k]
        ck1 = b.twoOverh[k+1]
        for q in 1:n
            uLπ = U[Np, k, 2, q]; uLχ = U[Np, k, 3, q]
            uRπ = tmpπ[q]; uRχ = tmpχ[q]           # shifted right state at kp
            Fsπ = (Fπ[q] + FπR[q]) / 2 - (α / 2) * (uRπ - uLπ)
            Fsχ = (Fχ[q] + FχR[q]) / 2 - (α / 2) * (uRχ - uLχ)
            dRπ = Fπ[q] - Fsπ; dRχ = Fχ[q] - Fsχ
            dLπ = FπR[q] - Fsπ; dLχ = FχR[q] - Fsχ
            γψ = src ? one(T) :
                 t < b.ψjump_toff ? one(T) :
                 t < b.ψjump_toff + 30 ? one(T) - (T(t) - b.ψjump_toff) / 30 : zero(T)
            dψ = γψ * (α / 2) * (U[1, k+1, 1, q] - U[Np, k, 1, q] -
                                 (src ? p.Jψ[q] * ph : zero(Complex{T})))
            for i in 1:Np
                dU[i, k, 1, q] += ck * b.invMcN[i] * dψ
                dU[i, k+1, 1, q] -= ck1 * b.invMc1[i] * dψ
                dU[i, k, 2, q] += ck * b.invMcN[i] * dRπ
                dU[i, k, 3, q] += ck * b.invMcN[i] * dRχ
                dU[i, k+1, 2, q] -= ck1 * b.invMc1[i] * dLπ
                dU[i, k+1, 3, q] -= ck1 * b.invMc1[i] * dLχ
            end
        end
    end
    # physical boundaries: pure outflow (vacuum treatment, no condition imposed)
    return du
end

"""
    solve_spin_sourced(p; tfinal=1200.0, dt=nothing, record_every=100) -> (τs, Ψ, Π, maxu)

Trivial initial data; evolve the sourced system and record complex probe series of ψ AND
π at `p.base.probe_R` (defaults: horizon R₊, mid-domain, scri) — `Ψ[probe, mode, sample]`,
`Π[probe, mode, sample]` — plus the boundedness monitor `maxu`. The steady-state π
amplitudes at the scri/horizon probes are the flux amplitudes (header (6)).
"""
function solve_spin_sourced(p::SpinSourcedProblem{T}; tfinal=1200.0, dt=nothing,
                            record_every=100) where {T}
    b = p.base
    Np, K, n = b.el.Np, b.mesh.K, b.n
    h = dt === nothing ? spin_cfl_dt(b) : dt
    u = zeros(Complex{T}, Np * K * 3 * n)
    sc = (zeros(Complex{T}, Np, n), zeros(Complex{T}, Np, n), zeros(Complex{T}, Np, n),
          zeros(Complex{T}, n), zeros(Complex{T}, n), zeros(Complex{T}, n),
          zeros(Complex{T}, n), zeros(Complex{T}, n), zeros(Complex{T}, n))
    nprb = length(b.probe_R)
    τs = T[]
    colsψ = Array{Complex{T},2}[]
    colsπ = Array{Complex{T},2}[]
    maxu = T[]
    fbuf = zeros(Complex{T}, Np)
    rec(uu, t, step) = begin
        step > 0 && _apply_filter!(uu, b, fbuf)
        if step % record_every == 0
            UU = reshape(uu, Np, K, 3, n)
            push!(τs, t)
            vψ = Matrix{Complex{T}}(undef, nprb, n)
            vπ = Matrix{Complex{T}}(undef, nprb, n)
            for ip in 1:nprb, q in 1:n
                vψ[ip, q] = dot(view(b.probe_L, ip, :), view(UU, :, b.probe_k[ip], 1, q))
                vπ[ip, q] = dot(view(b.probe_L, ip, :), view(UU, :, b.probe_k[ip], 2, q))
            end
            push!(colsψ, vψ)
            push!(colsπ, vπ)
            push!(maxu, maximum(abs, uu))
        end
    end
    rhs!(du, uu, t) = _rhs!(du, uu, t, p, sc)
    rk4!(rhs!, u, (0.0, tfinal), h; callback=rec)
    Ψ = Array{Complex{T},3}(undef, nprb, n, length(τs))
    Π = Array{Complex{T},3}(undef, nprb, n, length(τs))
    for j in eachindex(colsψ), ip in 1:nprb, q in 1:n
        Ψ[ip, q, j] = colsψ[j][ip, q]
        Π[ip, q, j] = colsπ[j][ip, q]
    end
    return τs, Ψ, Π, maxu
end

# ---- (6) flux extraction -------------------------------------------------------------------------

"""
    spin_sourced_flux(p, τs, Π; twin=(0.75, 1.0)) -> Vector{NamedTuple}

Per-mode steady-state fluxes from the recorded π probe series (header (6)):
|Z^∞| = mean|π(R=1)|/ω and |Z^H| = mean|π(R₊)|/(ω r₊) over the late-time window
`twin·τ_end`, converted through the in-package amplitude→flux layer `energy_flux`
(AmplitudesFluxes/fluxes.jl). Returns, per mode ℓ: `(l, Zinf, ZH, EI, EH, devI, devH)`
where `devI/devH` are the max relative deviations of |π| over the window (steady-state
quality monitors — measured, they bound the flux systematic ≲ 2·dev).
"""
function spin_sourced_flux(p::SpinSourcedProblem{T}, τs, Π; twin=(0.75, 1.0)) where {T}
    b = p.base
    ω = p.omega
    rplus = 2.0                     # a = 0, M = 1
    tlo = twin[1] * τs[end]
    thi = twin[2] * τs[end]
    sel = [j for j in eachindex(τs) if tlo <= τs[j] <= thi]
    isempty(sel) && error("spin_sourced_flux: empty steady-state window")
    out = NamedTuple[]
    for (q, l) in enumerate(b.ℓs)
        aS = [abs(Π[3, q, j]) for j in sel]      # scri probe (default probe 3 = R = 1)
        aH = [abs(Π[1, q, j]) for j in sel]      # horizon probe (default probe 1 = R₊)
        Zi = sum(aS) / length(aS) / ω
        Zh = sum(aH) / length(aH) / (ω * rplus)
        devI = maximum(abs.(aS ./ (ω * Zi) .- 1))
        devH = maximum(abs.(aH ./ (ω * rplus * Zh) .- 1))
        ef = energy_flux(p.base.s, l, b.m, 0.0, ω, p.lambda[q], Zi, Zh)
        push!(out, (l=l, Zinf=Zi, ZH=Zh, EI=ef.I, EH=ef.H, devI=devI, devH=devH))
    end
    return out
end
