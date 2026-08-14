# Sourced |s|=2 / s=0 circular-equatorial point-particle sector at Kerr a ≠ 0.
#
# ROUTE: the FD source per mode (ℓ̂, m, ω = mΩ) lives on the spheroidal ₛS_{ℓ̂mω}(θ)
# [gr-qc/9910091 Eqs. (4.39)–(4.40)]; the TD evolution runs in the spherical basis. The
# re-projection ₛS_{ℓ̂mω} = Σ_ℓ b_ℓ^{(ℓ̂)} ₛYₗₘ uses the Harmonics SphericalExpansion
# coefficients (real for real aω); the source triplet is the FD oracle's own body
# (`_convolve_source_circular_*`) and the jump algebra is reused unchanged.
#
# THE CHAIN (formulas implemented below; gates in test_timedomain_spinsource_kerr.jl):
#
#  (1) Per-SPHEROIDAL-mode radial jumps at a ≠ 0. For each ℓ̂ the radial Teukolsky ODE is
#      still a scalar 2nd-order ODE, X″ + pX′ + qX = S₀δ + S₁δ′ + S₂δ″, with
#      p = (s+1)Δ′/Δ, q = (−λ_{ℓ̂}(aω) + 4isrω + (K² − 2is(r−1)K)/Δ)/Δ,
#      K = (r²+a²)ω − am, Δ = r² − 2r + a² (≡ AmplitudesFluxes `_d2R_ode`, M = 1), and the
#      GHP source triplet (a₀, a₁, a₂) is the general-a oracle body specialized to
#      θ₀ = π/2 (ρ = ρ̄ = −1/r₀; s = +2 in the oracle's Δ²-rescaled pairing). The
#      distributional division and jump algebra (`_ode_jumps`; [1702.04204 App. A]) are
#      unchanged — only the Δ(a), p, q, λ inputs generalize.
#
#  (2) Spheroidal→spherical projection. The
#      monochromatic physical solution is ψₘ(t, r, θ) = Σ_ℓ̂ X_ℓ̂(r) ₛS_{ℓ̂mω}(θ) e^{−iωt};
#      inserting ₛS_{ℓ̂mω} = Σ_ℓ b_ℓ^{(ℓ̂)} ₛYₗₘ (b real for real aω) and reading off the
#      ₛYₗₘ component gives ψ_ℓ(t, r) = Σ_ℓ̂ b_ℓ^{(ℓ̂)} X_ℓ̂(r) e^{−iωt}. Jumps are linear:
#         [ψ̂_ℓ] = Σ_ℓ̂ b_ℓ^{(ℓ̂)} [ψ̂_ℓ̂],   [ψ̂′_ℓ] = Σ_ℓ̂ b_ℓ^{(ℓ̂)} [ψ̂′_ℓ̂],
#      with ψ̂ = rΔˢX (mode variable, ℓ-independent weight — spinsource.jl header (3);
#      ψ_RT = φ_BG holds at ANY a: both are the Δ⁻ˢr⁻¹ rescaling [arXiv:1301.1591,
#      Eq. (12)]). The ℓ̂ sum runs over the evolved band PLUS `jext` guard modes; the b's
#      decay as O((aω)^{|ℓ−ℓ̂|}). At a = 0, b_ℓ^{(ℓ̂)} = δ_{ℓℓ̂} exactly.
#
#  (3) Kerr chart push-forward (generalizes spinsource.jl header (4)).
#      Chart maps [arXiv:1301.1591, Eqs. (3)–(6)]: τ = t − r + r*_a (ingoing Kerr;
#      r*_a = Kerr tortoise, Spacetime.tortoise_kerr), φ̃ = φ + a∫dr/Δ,
#      τ = T − 4M log(1−R²) + (1+R²)/(1−R²), r = 2R/(1−R²). Hence at fixed RT time T the
#      e^{imφ̃}-mode of the monochromatic solution carries the R-local phase
#         ψ_ℓ(T, R) = e^{−iωT} · Φ(R) · Σ_ℓ̂ b_ℓ^{(ℓ̂)} ψ̂_ℓ̂(r(R)),
#         Φ(R) = exp(−i[ω·H_a(R) + m·aG(r(R))]),
#         H_a(R) = −4M log(1−R²) + (1+R²)/(1−R²) + r(R) − r*_a(r(R)),
#         aG(r) = a∫dr/Δ = (a/(r₊−r₋))·log((r−r₊)/(r−r₋))   (φ̃ shift; e^{imφ} = e^{imφ̃}e^{−imaG}).
#      Jump state at the particle interface (jumps enforced as Ĵ·e^{−iωT}·ramp):
#         Ĵψ_ℓ = [ψ̂_ℓ]·Φ(R₀),   Ĵπ_ℓ = iω·Ĵψ_ℓ,
#         Ĵχ_ℓ = ([ψ̂′_ℓ]·(dr/dR)₀ − i(ω·H_a′(R₀) + m·a·(dr/dR)₀/Δ₀)·[ψ̂_ℓ])·Φ(R₀),
#         H_a′(R) = 8MR/(1−R²) + 4R/(1−R²)² + (1 − (r²+a²)/Δ)·dr/dR.
#      Additive constants in r*_a / aG are a common phase (drop out of |Z|). At a = 0:
#      H_a ≡ rt_F_of_R and aG ≡ 0.
#
#  (4) Evolution: the jump-conditioned interface fluxes, driver and ramp are reused
#      verbatim (`SpinSourcedProblem` + `solve_spin_sourced` + `_rhs!`); the A·J
#      cancellation argument of spinsource.jl header (5) is linear-algebraic and holds
#      for the coupled matrix flux unchanged.
#
#  (5) Flux extraction at a ≠ 0 (spheroidal back-projection). Unit-transmission asymptotics
#      (Radial): X_ℓ̂ → Z^∞_ℓ̂ r^{−(2s+1)}e^{iωr*} at ℐ⁺ and X_ℓ̂ → Z^H_ℓ̂ Δ^{−s}e^{−ikr*}
#      at H⁺ (k = ω − mΩ_H), so |ψ̂_ℓ̂| → |Z^∞_ℓ̂| at ℐ⁺ and → r₊|Z^H_ℓ̂| at H⁺ (the rΔˢ
#      weight cancels the asymptotics at ANY a; the ingoing-chart phase mΩ_H r* − maG is
#      horizon-regular since Ω_H·2Mr₊ = a). The T-frequency is ω EVERYWHERE (t − T and
#      φ̃ − φ are functions of R only), so |π_ℓ| = ω|ψ_ℓ| in steady state. The recorded
#      complex spherical amplitudes are re-projected per spheroidal mode with the SAME
#      (real, orthonormal) b's — the phase factors Φ(R_probe)e^{−iωT} are ℓ-independent, so
#         |Z^∞_ℓ̂| = |Σ_ℓ b_ℓ^{(ℓ̂)} π_ℓ(R = 1)| / ω,
#         |Z^H_ℓ̂| = |Σ_ℓ b_ℓ^{(ℓ̂)} π_ℓ(R₊)| / (ω r₊),   r₊ = 1 + √(1 − a²),
#      then Ė^∞/Ė^H per (ℓ̂, m) through `energy_flux` (a ≠ 0 horizon normalization incl.
#      superradiant k < 0).
#
# HONEST BOUNDARIES (in-code guards): |s| = 1 sourced remains [FUTURE] (as at a = 0);
# m = 0 excluded (fluxless static); eccentric/inclined Kerr sourced TD out of scope (the
# jump ansatz is monochromatic); |a| < 1 (subextremal).

# ---- (1) general-a GHP source triplet at θ₀ = π/2 ------------------------------------------------

"""
    _ghp_source_kerr(s, m, a, r0, En, Lz, omega, S0, dS0, d2S0) -> (a0, a1, a2, w)

Circular-equatorial GHP source triplet at general `a` (M = 1), θ₀ = π/2, in the FD oracle's
pairing form `α = ∫ R̃·(a₀δ + a₁δ′ + a₂δ″) dr`, R̃ = R (s = −2, 0) or Δ²R (s = +2).
Transcription of the in-package oracle bodies `_convolve_source_circular_{s2,s2p,s0}`
(AmplitudesFluxes/pointparticle.jl; GHP projection of [arXiv:gr-qc/9910091,
Eqs. (4.39)–(4.40)]) specialized to the equator: sinθ₀ = 1, cosθ₀ = 0, ρ = ρ̄ = −1/r₀,
Σ = r₀²; the SWSH values (S0, dS0, d2S0) are the SPHEROIDAL ₛS_{ℓ̂mω}(π/2)-values at
γ = aω. Reduces bit-tightly to `_ghp_source_a0` at a = 0 (gate K3 via the jump chain).
`w` is the Δ-power of the normalized-ODE source wrap (spinsource.jl `_ode_source_Sn` doc).
"""
function _ghp_source_kerr(s::Integer, m::Integer, a::Real, r0, En, Lz, omega, S0, dS0, d2S0)
    rho = -1 / r0
    Δ0 = r0^2 - 2r0 + a^2
    Kt = (r0^2 + a^2) * omega - m * a
    sq2 = sqrt(2.0)
    if s == -2
        # oracle `_convolve_source_circular_s2` at θ₀ = π/2 (L₁ = L₂ = −m + aω there)
        L1 = -m + a * omega
        L2 = L1
        L2S = dS0 + L2 * S0
        L1L2S = (d2S0 + L1 * dS0) + (-2) * S0 + L2 * dS0 + L1 * L2 * S0
        Ann0 = -rho^(-2) * rho^(-1) * (sq2 * Δ0)^(-2) *
               (rho^(-1) * L1L2S + 3im * a * L1 * S0 + 2im * a * dS0 - im * a * L2 * S0)
        Anmbar0 = rho^(-3) * (sq2 * Δ0)^(-1) * ((2rho - im * Kt / Δ0) * L2S)
        Anmbar1 = -rho^(-3) * (sq2 * Δ0)^(-1) * L2S
        Ambarmbar0 = (Kt^2 * S0 * rho) / (4 * Δ0^2 * rho^3) +
                     (im * Kt * S0 * (1 - r0 + Δ0 * rho) * rho) / (2 * Δ0^2 * rho^3) +
                     (im * r0 * S0 * rho * omega) / (2 * Δ0 * rho^3)
        Ambarmbar1 = -rho^(-3) * rho * S0 / 2 * (im * Kt / Δ0 - rho)
        Ambarmbar2 = -rho^(-3) * rho * S0 / 4
        rcomp = (En * (r0^2 + a^2) - a * Lz) / (2 * r0^2)
        thetacomp = rho * (im * (a * En - Lz)) / sq2
        Cnn = rcomp^2; Cnm = rcomp * thetacomp; Cmm = thetacomp^2
        return (Ann0 * Cnn + Anmbar0 * Cnm + Ambarmbar0 * Cmm,
                Anmbar1 * Cnm + Ambarmbar1 * Cmm,
                Ambarmbar2 * Cmm, -(s + 1))
    elseif s == 2
        # oracle `_convolve_source_circular_s2p` at θ₀ = π/2 (δL₁ = δL₂ = m − aω there;
        # ρ′/ρ = iaρ, ρ″/ρ = 2(iaρ)² ⇒ the −6(ρ′/ρ)² + 3(ρ″/ρ) group cancels exactly)
        dL = m - a * omega
        drr = im * a * rho
        All0 = -(1 / 2) * (d2S0 + (2 * dL + 2 * drr) * dS0 +
                           ((-2) + dL^2 + 2 * dL * drr) * S0)
        L2S = dS0 + dL * S0
        Alm0 = (2 / sq2) * rho^(-1) * (-(2rho + im * Kt / Δ0) * L2S)
        Alm1 = (2 / sq2) * rho^(-1) * L2S
        Amm0 = (Kt^2 * S0) / (Δ0^2 * rho^2) +
               (2 * im * Kt * S0 * (-1 + r0 - Δ0 * rho)) / (Δ0^2 * rho^2) -
               (2 * im * r0 * S0 * omega) / (Δ0 * rho^2)
        Amm1 = 2 * rho^(-2) * S0 * (im * Kt / Δ0 + rho)
        Amm2 = -rho^(-2) * S0
        rcomp = (En * (r0^2 + a^2) - a * Lz) / Δ0
        thetacomp = -rho * (im * (a * En - Lz)) / sq2
        Cll = rcomp^2; Clm = rcomp * thetacomp; Cmm = thetacomp^2
        return (All0 * Cll + Alm0 * Clm + Amm0 * Cmm,
                Alm1 * Clm + Amm1 * Cmm,
                Amm2 * Cmm, 2 - (s + 1))
    elseif s == 0
        # oracle `_convolve_source_circular_s0`: α = −4π r₀² S at ANY a (Σ = r₀² on the
        # equator) — same δ-only wrap as `_ghp_source_a0`.
        return (complex(r0^2 * S0 / 2), 0.0im, 0.0im, -1)
    else
        throw(ArgumentError("Kerr spin source: s ∈ {0, −2, +2} built (|s| = 1 is [FUTURE]; got s=$s)"))
    end
end

# ---- (1) Kerr ODE coefficients + distributional division -----------------------------------------

"Kerr `_ode_source_Sn` (spinsource.jl doc, Δ generalized): g(r) = (−8π/Υt)Δ(a,r)^w,
S₂ = g₀a₂, S₁ = g₀a₁ − 2g₀′a₂, S₀ = g₀a₀ − g₀′a₁ + g₀″a₂. M = 1. Bit-equal to the a = 0
function at a = 0 (Δ₀ + 0.0 is exact)."
function _ode_source_Sn_kerr(a0, a1, a2, w::Integer, a::Real, r0, Upst)
    pref = -8pi / Upst
    Δ0 = r0^2 - 2r0 + a^2
    dΔ0 = 2r0 - 2
    g0 = pref * Δ0^w
    g1 = pref * w * Δ0^(w - 1) * dΔ0
    g2 = pref * (w * (w - 1) * Δ0^(w - 2) * dΔ0^2 + 2w * Δ0^(w - 1))
    return (g0 * a0 - g1 * a1 + g2 * a2, g0 * a1 - 2 * g1 * a2, g0 * a2)
end

"radial-Teukolsky ODE coefficients at r₀ for general a (M = 1): p = (s+1)Δ′/Δ, p′, and
q = (−λ + 4isrω + (K² − 2is(r−1)K)/Δ)/Δ, K = (r²+a²)ω − am (≡ AmplitudesFluxes `_d2R_ode`;
λ = spheroidal eigenvalue λ_{ℓ̂}(aω), SpinWeightedSpheroidalEigenvalue convention)."
function _ode_pq_kerr(s::Integer, lambda::Real, m::Integer, a::Real, r0, omega)
    Δ0 = r0^2 - 2r0 + a^2
    dΔ0 = 2r0 - 2
    K = (r0^2 + a^2) * omega - a * m
    p0 = (s + 1) * dΔ0 / Δ0
    dp0 = (s + 1) * (2 / Δ0 - dΔ0^2 / Δ0^2)
    q0 = (-lambda + 4im * s * r0 * omega + (K^2 - 2im * s * (r0 - 1) * K) / Δ0) / Δ0
    return (p0, dp0, q0)
end

# ---- (1)+(2) per-spheroidal-mode jumps -------------------------------------------------------------

"""
    spin_source_jumps_kerr(s, lhat, m, a, r0; Sc=nothing) -> NamedTuple

Complete t = 0 jump data of the sourced SPHEROIDAL mode (s, ℓ̂, m, ω = mΩ) at circular
equatorial orbit radius `r0` around Kerr `a` (M = 1, μ_particle = 1): radial-function jumps
`(Xδ, JX, JXp)` and the mode-variable ψ̂ = rΔˢX jumps `(psid, Jpsi, Jpsip)` (provenance
header (1); the a = 0 committed `spin_source_jumps` is the tight reduction, gate K3).
Conventions match the FD oracle exactly: ω = m·Ω_φ (BL), Υt Mino, λ from the SWSHSFunction
eigen-solve. Pass `Sc` (the ₛS_{ℓ̂mω} SphericalExpansion carrier) to reuse it; `m = 0`
excluded (fluxless static; monochromatic ansatz needs ω ≠ 0).
"""
function spin_source_jumps_kerr(s::Integer, lhat::Integer, m::Integer, a::Real, r0::Real;
                                Sc=nothing)
    m == 0 && throw(ArgumentError("spin_source_jumps_kerr: m ≠ 0 required (static m = 0 is fluxless)"))
    abs(a) < 1 || throw(ArgumentError("spin_source_jumps_kerr: |a| < 1 (subextremal) required"))
    # timelike circular-geodesic existence (equatorial): 1 − 3/r₀ + 2a/r₀^{3/2} > 0
    1 - 3 / r0 + 2a / r0^(3 / 2) > 0 ||
        throw(ArgumentError("spin_source_jumps_kerr: no timelike circular geodesic at r0=$r0, a=$a"))
    En = geodesic_energy(a, r0, 0.0, 1.0)
    Lz = geodesic_angular_momentum(a, r0, 0.0, 1.0)
    Omega = kerr_frequencies(a, r0, 0.0, 1.0).Omegaph
    Upst = kerr_frequencies(a, r0, 0.0, 1.0; Time="Mino").Ut
    omega = m * Omega
    if Sc === nothing
        Sc = SpinWeightedSpheroidalHarmonicS(s, lhat, m, a * omega)
    end
    lambda = Sc.eigenvalue
    S0 = Sc(pi / 2, 0.0)
    dS0 = swsh_derivative(Sc, 1, 0, pi / 2, 0.0)
    d2S0 = swsh_derivative(Sc, 2, 0, pi / 2, 0.0)
    a0, a1, a2, w = _ghp_source_kerr(s, m, a, r0, En, Lz, omega, S0, dS0, d2S0)
    S0n, S1n, S2n = _ode_source_Sn_kerr(a0, a1, a2, w, a, r0, Upst)
    p0, dp0, q0 = _ode_pq_kerr(s, lambda, m, a, r0, omega)
    J = _ode_jumps(S0n, S1n, S2n, p0, dp0, q0)
    Δ0 = r0^2 - 2r0 + a^2
    dΔ0 = 2r0 - 2
    W0 = r0 * float(Δ0)^s                                  # (rΔˢ)₀
    W1 = float(Δ0)^s + s * r0 * float(Δ0)^(s - 1) * dΔ0    # (rΔˢ)′₀
    return (s=s, l=lhat, m=m, a=float(a), r0=float(r0), omega=omega, Omega=Omega,
            Upst=Upst, lambda=lambda,
            Xδ=J.Xδ, JX=J.JX, JXp=J.JXp,
            psid=W0 * J.Xδ, Jpsi=W0 * J.JX, Jpsip=W1 * J.JX + W0 * J.JXp)
end

# ---- (3) Kerr RT chart push-forward ---------------------------------------------------------------

"""
    rt_H_of_R_kerr(R; a=0.0) -> t − T   (Kerr, M = 1)

DERIVED-CK (header (3)): composition of the chart maps [arXiv:1301.1591, Eqs. (3),(5),(6)]
at general a: `H_a(R) = −4 log(1−R²) + (1+R²)/(1−R²) + r(R) − r*_a(r(R))` with the Kerr
tortoise `Spacetime.tortoise_kerr` (M = 1). At a = 0 this equals `rt_F_of_R` (identical
closed forms; `tortoise_kerr(r, 0) ≡ tortoise_schw(r)`). Additive tortoise constants are a
common phase of the monochromatic jump state (drop out of every |Z|)."""
function rt_H_of_R_kerr(R; a=0.0)
    r = rt_r_of_R(R)
    return -4 * log(1 - R^2) + (1 + R^2) / (1 - R^2) + r - tortoise_kerr(r, a)
end

"analytic H_a′(R) (DERIVED-CK, header (3)): `8R/(1−R²) + 4R/(1−R²)² + (1 − (r²+a²)/Δ)·dr/dR`,
dr/dR = 2(1+R²)/(1−R²)². Reduces to `_rt_dF_dR` at a = 0 (r²/Δ = 1/f). M = 1."
function _rt_dH_dR_kerr(R; a=0.0)
    drdR = 2 * (1 + R^2) / (1 - R^2)^2
    r = rt_r_of_R(R)
    Δ = r^2 - 2r + a^2
    return 8R / (1 - R^2) + 4R / (1 - R^2)^2 + (1 - (r^2 + a^2) / Δ) * drdR
end

"azimuthal chart shift aG(r) = a∫dr/Δ = (a/(r₊−r₋))·log((r−r₊)/(r−r₋)) (M = 1;
[arXiv:1301.1591, Eq. (4)]: φ̃ = φ + aG ⇒ the e^{imφ̃}-mode carries e^{−im·aG}).
Returns 0 at a = 0 (the shift enters only as m·aG)."
function _kerr_phi_shift(r, a)
    iszero(a) && return 0.0
    sq = sqrt(1 - a^2)
    rp = 1 + sq
    rm = 1 - sq
    return (a / (rp - rm)) * log((r - rp) / (r - rm))
end

# ---- (4) sourced Kerr problem ---------------------------------------------------------------------

"""
    SpinSourcedKerrProblem{T}

Sourced RT-chart circular-equatorial problem at Kerr a ≠ 0 (DG-S2): wraps a DG-P5
[`SpinSourcedProblem`](@ref) whose per-SPHERICAL-mode jump amplitudes are the b-projected
per-spheroidal jumps (provenance header (2)+(3)), plus the back-projection data for the
flux extraction (header (5)). Build with [`SpinSourcedKerrProblem(s, ℓs, m, a, r0; ...)`](@ref);
evolve with `solve_spin_sourced` (forwarded); fluxes with [`spin_sourced_flux_kerr`](@ref).
The committed a = 0 `SpinSourcedProblem` path is REUSED, not modified.
"""
struct SpinSourcedKerrProblem{T<:AbstractFloat}
    prob::SpinSourcedProblem{T}
    a::T
    omega::T
    rplus::T
    lamS::Vector{T}       # spheroidal eigenvalue λ_{ℓ̂}(aω) per band mode ℓ̂ = ℓs
    bmat::Matrix{T}       # bmat[q, j] = b_{ℓ_q}^{(ℓ̂_j)} over the evolved band (real)
end

"""
    SpinSourcedKerrProblem(s, ℓs, m, a, r0; jext=2, N=8, K=48, excise=0.05,
                           filter_alpha=36.0, filter_power=16, constraint_kappa=0.3,
                           psi_jump_toff=250.0, ramp_τ=400.0, ramp_δ=0.00025,
                           extract_R=Float64[], T=Float64)

Build the Kerr (M = 1, |a| < 1) sourced circular-equatorial problem for spin weight
`s ∈ {0, −2, +2}` (|s| = 1 remains [FUTURE], as at a = 0), azimuthal number `m ≠ 0`,
consecutive coupled-mode band `ℓs` and orbit radius `r0`. The per-spherical-ℓ jump data sum
the per-spheroidal jumps over ℓ̂ ∈ [max(|s|,|m|), max(ℓs)+jext] with Harmonics
b-coefficients (header (2)); `jext` guard modes bound the ℓ̂-truncation (b decay
O((aω)^{|ℓ−ℓ̂|}), default 2). Numerical kwargs pass through to the DG-P5/DG-P4
constructors. At a = 0 this reduces to the committed `SpinSourcedProblem` data (gate K3)."""
function SpinSourcedKerrProblem(s::Integer, ℓs::Vector{Int}, m::Integer, a::Real, r0::Real;
                                jext::Integer=2, N=8, K=48, excise=0.05, filter_alpha=36.0,
                                filter_power=16, constraint_kappa=0.3, psi_jump_toff=250.0,
                                ramp_τ=400.0, ramp_δ=0.00025, extract_R=Float64[],
                                T::Type{<:AbstractFloat}=Float64)
    s in (0, -2, 2) || throw(ArgumentError("SpinSourcedKerrProblem: s ∈ {0, −2, +2} (got $s)"))
    m == 0 && throw(ArgumentError("SpinSourcedKerrProblem: m ≠ 0 required"))
    abs(a) < 1 || throw(ArgumentError("SpinSourcedKerrProblem: |a| < 1 (subextremal) required"))
    jext >= 0 || throw(ArgumentError("SpinSourcedKerrProblem: jext ≥ 0 required"))
    R0 = T(rt_R_of_r(r0))
    base = SpinModesProblem(s, ℓs, m; N=N, K=K, M=1.0, a=float(a), excise=excise,
                            filter_alpha=filter_alpha, filter_power=filter_power,
                            constraint_kappa=constraint_kappa, psi_jump_toff=psi_jump_toff,
                            extract_R=extract_R, forced_edges=[Float64(R0)], T=T)
    kp = findmin(abs.(base.mesh.edges .- R0))[2] - 1
    (kp >= 1 && kp <= base.mesh.K - 1) || error("particle interface not interior to the mesh")
    n = length(ℓs)
    # ---- per-spheroidal jumps + b-coefficients over the guarded ℓ̂ set (headers (1)+(2))
    jmin = max(abs(s), abs(m))
    jset = collect(jmin:(maximum(ℓs) + jext))
    local ω
    JPsi = Dict{Int,Complex{Float64}}()
    JPsip = Dict{Int,Complex{Float64}}()
    bcof = Dict{Int,Dict{Int,Float64}}()
    lamhat = Dict{Int,Float64}()
    # ω needed for γ = aω before the loop: same convention as spin_source_jumps_kerr
    ω = m * kerr_frequencies(a, r0, 0.0, 1.0).Omegaph
    for ĵ in jset
        Sc = SpinWeightedSpheroidalHarmonicS(s, ĵ, m, a * ω)
        J = spin_source_jumps_kerr(s, ĵ, m, a, r0; Sc=Sc)
        ω = J.omega
        JPsi[ĵ] = J.Jpsi
        JPsip[ĵ] = J.Jpsip
        bcof[ĵ] = Dict(k => Float64(v) for (k, v) in expansion_coefficients(Sc))
        lamhat[ĵ] = J.lambda
    end
    # ---- chart push-forward factors at R₀ (header (3); common to all modes)
    Δ0 = r0^2 - 2r0 + a^2
    drdR = 2 * (1 + R0^2) / (1 - R0^2)^2
    phase = cis(-(ω * rt_H_of_R_kerr(R0; a=a) + m * _kerr_phi_shift(r0, a)))
    dchi = ω * _rt_dH_dR_kerr(R0; a=a) + m * a * drdR / Δ0
    Jψ = Vector{Complex{T}}(undef, n)
    Jπ = Vector{Complex{T}}(undef, n)
    Jχ = Vector{Complex{T}}(undef, n)
    λv = Vector{T}(undef, n)
    bmat = zeros(T, n, n)
    for (q, l) in enumerate(ℓs)
        jPsi = zero(Complex{Float64})
        jPsip = zero(Complex{Float64})
        for ĵ in jset
            b = get(bcof[ĵ], l, 0.0)
            jPsi += b * JPsi[ĵ]
            jPsip += b * JPsip[ĵ]
        end
        Jψ[q] = jPsi * phase
        Jπ[q] = im * ω * Jψ[q]
        Jχ[q] = (jPsip * drdR - im * dchi * jPsi) * phase
    end
    for (j, ĵ) in enumerate(ℓs)          # band back-projection data (header (5))
        λv[j] = lamhat[ĵ]
        for (q, l) in enumerate(ℓs)
            bmat[q, j] = get(bcof[ĵ], l, 0.0)
        end
    end
    prob = SpinSourcedProblem{T}(base, T(r0), R0, kp, T(ω), λv, Jψ, Jπ, Jχ, T(ramp_τ), T(ramp_δ))
    return SpinSourcedKerrProblem{T}(prob, T(a), T(ω), T(1 + sqrt(1 - a^2)), λv, bmat)
end

"forward the DG-P5 evolution driver (additive dispatch; the a = 0 path is untouched)."
solve_spin_sourced(p::SpinSourcedKerrProblem; kw...) = solve_spin_sourced(p.prob; kw...)

# ---- (5) flux extraction --------------------------------------------------------------------------

"""
    spin_sourced_flux_kerr(kp, τs, Π; twin=(0.75, 1.0)) -> Vector{NamedTuple}

Per-SPHEROIDAL-mode steady-state fluxes at a ≠ 0 from the recorded complex π probe series
(provenance header (5)): back-project the spherical band onto ₛS_{ℓ̂mω} with the (real,
orthonormal) b-matrix, then |Z^∞_ℓ̂| = mean|Σ_ℓ b π_ℓ(R=1)|/ω, |Z^H_ℓ̂| =
mean|Σ_ℓ b π_ℓ(R₊)|/(ω r₊) over the late-time window, converted through the in-package
`energy_flux` (a ≠ 0 normalization incl. superradiant modes). Returns per band mode ℓ̂:
`(l, Zinf, ZH, EI, EH, devI, devH)`; `devI/devH` are max relative |amplitude| deviations
over the window (steady-state monitors). Only the interior band modes ℓ̂ ≤ max(ℓs) − 1 are
gate-grade (the top band mode carries the ℓ-truncation edge)."""
function spin_sourced_flux_kerr(kp::SpinSourcedKerrProblem{T}, τs, Π;
                                twin=(0.75, 1.0)) where {T}
    p = kp.prob
    b = p.base
    ω = kp.omega
    tlo = twin[1] * τs[end]
    thi = twin[2] * τs[end]
    sel = [j for j in eachindex(τs) if tlo <= τs[j] <= thi]
    isempty(sel) && error("spin_sourced_flux_kerr: empty steady-state window")
    n = b.n
    out = NamedTuple[]
    for (j, ĵ) in enumerate(b.ℓs)
        aS = Vector{T}(undef, length(sel))
        aH = Vector{T}(undef, length(sel))
        for (i, jj) in enumerate(sel)
            zS = zero(Complex{T})
            zH = zero(Complex{T})
            for q in 1:n
                zS += kp.bmat[q, j] * Π[3, q, jj]      # scri probe (default probe 3 = R = 1)
                zH += kp.bmat[q, j] * Π[1, q, jj]      # horizon probe (default probe 1 = R₊)
            end
            aS[i] = abs(zS)
            aH[i] = abs(zH)
        end
        Zi = sum(aS) / length(aS) / ω
        Zh = sum(aH) / length(aH) / (ω * kp.rplus)
        devI = maximum(abs.(aS ./ (ω * Zi) .- 1))
        devH = maximum(abs.(aH ./ (ω * kp.rplus * Zh) .- 1))
        ef = energy_flux(b.s, ĵ, b.m, Float64(kp.a), Float64(ω), Float64(kp.lamS[j]), Zi, Zh)
        push!(out, (l=ĵ, Zinf=Zi, ZH=Zh, EI=ef.I, EH=ef.H, devI=devI, devH=devH))
    end
    return out
end
