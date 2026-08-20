# Newtonian planet/EMRI–disc background of Miranda & Rafikov (arXiv:1811.09628 "MR19" linear theory
# and its arXiv:2007.13787 "MR20" locally-cooled extension) — the model the MR19 figure set is built
# on. Ported from the legacy `DiscModes/newtonian_disc.jl` (itself the verified port of the authors'
# `planet-disk-modes-cooling/main.py`, commit 948d41c).
#
# SCOPE, and why it sits BESIDE the other three disc formulations of this package rather than inside
# one of them: `Modes.RelativisticDisc` is the Kerr Lin–Shu reduced-operator disc (arXiv:2601.19123),
# `SpiralWaves.DiscBackground` the Page–Thorne Lagrangian master system (arXiv:2608.17003) and its
# Newtonian LIMIT; this is the MR19/MR20 **locally-cooled, power-law, planet-forced** Newtonian disc:
# a different published model with its own cooling time t_c = β/Ω, its own Schwarzschild buoyancy
# N_r², and its own softened-planet forcing through Laplace coefficients. None of the three shares a
# formula with the others, so nothing is duplicated.
#
# Names here are DELIBERATELY unexported (`Modes.Omega`, `Modes.Sigma`, …): they are this model's
# profiles, and the exported disc generics of the package (`disc_OmegaK`, `csound`, `Omega` of
# `SpiralWaves`) belong to the other backgrounds. Only the solver API is exported.
#
# TYPE-GENERIC (Float64 default; BigFloat inherited).

"""
    NewtonianDisc(; h0, q, p, gamma, beta, soft, indirect, pole_disp, lap_nphi)

Newtonian locally-cooled accretion-disc background (Miranda & Rafikov 2019/2020), in units
r_p = Ω_p = M_* = 1 (the planet/secondary orbits at r = 1).

* `h0` — aspect ratio H/r at r = 1;
* `q`, `p` — temperature and surface-density power-law indices (c_s ∝ r^(−q/2), Σ ∝ r^(−p));
* `gamma` — adiabatic index Γ;
* `beta` — dimensionless cooling time β = Ω t_c. MR19 is the ADIABATIC (no-cooling) theory, i.e. the
  limit β → ∞ (β ≳ 10³ is converged); the MR20 extension's own `parameters.py` ships β = 0.1, which
  damps the outward-propagating waves — the two regimes give visibly different figures;
* `soft` — Plummer softening of the secondary's potential;
* `indirect` ∈ {0,1} — include the indirect (m = 1) term;
* `pole_disp` — imaginary displacement of the corotation singularity;
* `lap_nphi` — azimuthal quadrature points for the Laplace coefficients (4000 = legacy default; 800
  is already < 2e-11 relative at corotation and ~5× faster).
"""
Base.@kwdef struct NewtonianDisc{T<:Real}
    h0::T          = 0.1
    q::T           = 1.0
    p::T           = 1.0
    gamma::T       = 1.4
    beta::T        = 0.1
    soft::T        = 0.6 * 0.1
    indirect::Int  = 0
    pole_disp::T   = 1e-6
    lap_nphi::Int  = 4000
end

# --- background profiles (s = radius) -------------------------------------------------------------
Sigma(d::NewtonianDisc, s)     = 1 / s^d.p
dSigmadr(d::NewtonianDisc, s)  = -d.p / s^(d.p + 1)
OmegaK(d::NewtonianDisc, s)    = 1 / s^(3//2)
Omega(d::NewtonianDisc, s)     = OmegaK(d, s) * sqrt(1 - d.h0^2 * (d.p + d.q) * s^(1 - d.q))
function dOmegadr(d::NewtonianDisc, s)
    -3//2 * OmegaK(d, s) / s *
    (-(s^(-d.q)) * (d.h0^2 * (2 + d.q) * (d.p + d.q) * s - 3 * s^d.q) / 3 /
     sqrt(1 - d.h0^2 * (d.p + d.q) * s^(1 - d.q)))
end
kappa(d::NewtonianDisc, s) = OmegaK(d, s) * sqrt(1 + d.h0^2 * (d.q - 2) * (d.p + d.q) * s^(1 - d.q))
function dkappadr(d::NewtonianDisc, s)
    -3//2 * OmegaK(d, s) / s *
    (s^(-d.q) * (d.h0^2 * (d.p + d.q) * (d.q^2 - 4) * s + 3 * s^d.q) / 3 /
     sqrt(1 + d.h0^2 * (d.q - 2) * (d.p + d.q) * s^(1 - d.q)))
end
"Pattern frequency ω_p = mΩ(r=1) — corotation sits at r = 1 for every m."
omega_pat(d::NewtonianDisc, m) = m
"Resonance denominator D = κ² − (mΩ − ω_p)²; its Lindblad zeros bound the wave-propagating region."
Dfreq(d::NewtonianDisc, s, m)  = kappa(d, s)^2 - (m * Omega(d, s) - omega_pat(d, m))^2
dDdr(d::NewtonianDisc, s, m)   = 2 * kappa(d, s) * dkappadr(d, s) -
                                 2 * m * (m * Omega(d, s) - omega_pat(d, m)) * dOmegadr(d, s)
csiso(d::NewtonianDisc, s)     = d.h0 / s^(d.q / 2)
dcsisodr(d::NewtonianDisc, s)  = -(d.q / 2) * d.h0 / s^(d.q / 2 + 1)
csadi(d::NewtonianDisc, s)     = sqrt(d.gamma) * csiso(d, s)
dcsadidr(d::NewtonianDisc, s)  = sqrt(d.gamma) * dcsisodr(d, s)
"Radial Schwarzschild buoyancy frequency N_r² (negative for a convectively stable power-law disc)."
function Nr2(d::NewtonianDisc, s)
    zeta = (d.q + (1 - d.gamma) * d.p) * (d.q + d.p) / d.gamma^2
    -zeta * csadi(d, s)^2 / s^2
end
function dNr2dr(d::NewtonianDisc, s)
    zeta = (d.q + (1 - d.gamma) * d.p) * (d.q + d.p) / d.gamma^2
    zeta * (d.q + 2) * csadi(d, s)^2 / s^3
end
LSinv(d::NewtonianDisc, s)    = (eta = (d.q + (1 - d.gamma) * d.p) / d.gamma; -eta / s)
dLSinvdr(d::NewtonianDisc, s) = (eta = (d.q + (1 - d.gamma) * d.p) / d.gamma; eta / s^2)
LTinv(d::NewtonianDisc, s)    = -d.q / s
dLTinvdr(d::NewtonianDisc, s) = d.q / s^2
"Cooling time t_c = β/Ω."
tc(d::NewtonianDisc, s)       = d.beta / Omega(d, s)
dtcdr(d::NewtonianDisc, s)    = -d.beta / Omega(d, s)^2 * dOmegadr(d, s)

# --- WKB wavenumber + outgoing-wave coefficient (main.py k_func / wkb_func) ------------------------
"""
    kwave(d, r, m) -> k

WKB radial wavenumber of the m-th cooled density wave: `k = √(−D_eff)/c_eff` with the cooling-modified
sound speed and resonance denominator (MR20 §3). Complex inside the evanescent band.
"""
function kwave(d::NewtonianDisc, s, m)
    ot = omega_pat(d, m) - m * Omega(d, s)
    t = tc(d, s)
    cseff = ((1 + (d.gamma * ot * t)^2) / (1 + (ot * t)^2))^(one(s) / 4) * csiso(d, s)
    Deff = Dfreq(d, s, m) + (d.gamma * ot * t)^2 / (1 + (d.gamma * ot * t)^2) * Nr2(d, s)
    sqrt(complex(-Deff)) / cseff
end

"""
    wkbC(d, r, m) -> C

Outgoing-wave boundary coefficient: the WKB relation δh′ = C δh imposed at the domain edges.
Includes the cooling-induced imaginary part of `k` and the amplitude-transport terms.
"""
function wkbC(d::NewtonianDisc, s, m)
    ot = omega_pat(d, m) - m * Omega(d, s)
    t = tc(d, s)
    rek = kwave(d, s, m)
    imk = (d.gamma - 1) / 2 * ot * t / (1 + d.gamma * (ot * t)^2) * rek
    k = rek + im * imk
    w = im * k
    w += dDdr(d, s, m) / (4 * Dfreq(d, s, m)) - dSigmadr(d, s) / (2 * Sigma(d, s)) +
         dcsisodr(d, s) / (2 * csiso(d, s)) - 1 / (2 * s)
    w += (1 + im * d.gamma * ot * t) / 2 / (1 + (d.gamma * ot * t)^2) * LTinv(d, s)
    w
end

# --- softened Laplace coefficients b_s^(m)(α; ε) + derivatives ------------------------------------
"""
    laplace_b_soft(s, m, α, ε; nphi = 4000) -> b_s^{(m)}(α; ε)

SOFTENED Laplace coefficient
`b_s^{(m)}(α; ε) = (2/π)∫₀^π cos(mφ)(1 + α² − 2α cosφ + ε²)^{−s} dφ`, by the midpoint rule on the
full period (the integrand is smooth and periodic ⇒ spectral convergence). `ε` is the Plummer
softening of the secondary's potential, without which the m-sum of the planet forcing does not
converge at the secondary's own radius.

This is the ε ≠ 0 member of the family whose ε = 0 member is [`laplace_b`](@ref) (the analytic kernel
used by the relativistic `forcing_Fm`); the two agree to quadrature accuracy at ε = 0, which is a
gate in `test/Modes/test_newtonian_planet.jl`. They are kept as separate entry points because each
is pinned to the published numbers of a different code path (`laplace_b`: arXiv:2601.19123 reduced
ODE; here: the Miranda–Rafikov `main.py` quadrature, whose node count `nphi` is a disc field).
"""
function laplace_b_soft(svar, m, alpha, eps; nphi = 4000)
    T = typeof(float(alpha)); acc = zero(T); dphi = 2T(pi) / nphi
    for j in 0:nphi-1
        phi = (j + one(T) / 2) * dphi
        acc += cos(m * phi) * (1 + alpha^2 - 2 * alpha * cos(phi) + eps^2)^(-svar)
    end
    acc * dphi / T(pi)
end
"∂_α of [`laplace_b_soft`](@ref) (same quadrature, differentiated integrand)."
function laplace_b_soft_deriv(svar, m, alpha, eps; nphi = 4000)
    T = typeof(float(alpha)); acc = zero(T); dphi = 2T(pi) / nphi
    for j in 0:nphi-1
        phi = (j + one(T) / 2) * dphi
        acc += 2 * svar * cos(m * phi) * (cos(phi) - alpha) *
               (1 + alpha^2 - 2 * alpha * cos(phi) + eps^2)^(-(svar + 1))
    end
    acc * dphi / T(pi)
end
"∂²_α of [`laplace_b_soft`](@ref)."
function laplace_b_soft_2deriv(svar, m, alpha, eps; nphi = 4000)
    T = typeof(float(alpha)); acc = zero(T); dphi = 2T(pi) / nphi
    for j in 0:nphi-1
        phi = (j + one(T) / 2) * dphi
        acc += 2 * svar * cos(m * phi) *
               (svar + (1 + 2 * svar) * (alpha^2 - 2 * alpha * cos(phi)) +
                (1 + svar) * cos(2 * phi) - eps^2) *
               (1 + alpha^2 - 2 * alpha * cos(phi) + eps^2)^(-(svar + 2))
    end
    acc * dphi / T(pi)
end

"Azimuthal harmonic Φ_m(r) of the softened secondary potential (= −b_{1/2}^{(m)}), + indirect term."
function Phi_m(d::NewtonianDisc, s, m)
    phi = -laplace_b_soft(one(float(s)) / 2, m, s, d.soft; nphi = d.lap_nphi)
    (m == 1 && d.indirect == 1) && (phi += s)
    phi
end
function dPhidr_m(d::NewtonianDisc, s, m)
    dphi = -laplace_b_soft_deriv(one(float(s)) / 2, m, s, d.soft; nphi = d.lap_nphi)
    (m == 1 && d.indirect == 1) && (dphi += one(float(s)))
    dphi
end
d2Phidr2_m(d::NewtonianDisc, s, m) = -laplace_b_soft_2deriv(one(float(s)) / 2, m, s, d.soft; nphi = d.lap_nphi)
