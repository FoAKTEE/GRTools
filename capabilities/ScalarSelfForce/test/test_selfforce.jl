# test_selfforce.jl — ESSF node `mode-self-force-extraction` (GLOBAL_DAG row 44, wave 12,
# worker ESSF-W12-A). THIS IS THE NUMBER THE PAPER IS ABOUT: F_r^m at the particle, out of a
# converged m-mode solution.
#
# WHAT THE NODE DELIVERS (`Extraction/selfforce.jl`):
#   * `mode_self_force(orbit, m, psi, dpsi_r, dpsi_z; phi_p)` — the pure Osburn–Nishimura
#     Eq. (3.10) formula on raw ingredients (type-generic, BigFloat-clean);
#   * `mode_self_force(sol::ModeSolution)` — the same read through the SEALED row-42
#     `field_at` / `field_grad_at`, reduced by the MEAN over every element containing the
#     puncture (spectre parity); `element = e` gives one one-sided trace instead;
#   * `mode_self_force(hist::AmrHistory; level, reduction)` — per AMR level, `:trace` by
#     default (the sealed `AmrLevel` diagnostics) so a convergence sequence is homogeneous;
#   * `mode_self_force_history`, `mode_force_convergence` — the sequence and its
#     successive-difference table (rows 50 / 53's resolution floor);
#   * `regularized_mode_force(...; order = 4)` — `F_r^m − F^m_{r[order]}` through
#     `mmode_reg_term_Fr` (row 45's helper, which owns the m = 0 halving).
#
# THE ORACLES, IN ORDER OF STRENGTH (this ordering is the point — see gate (a0)):
#   1. INDEPENDENT, and the only one that is: a BigFloat central-difference of the FULL
#      FOLDED FIELD `Φ = 2Re[Ψ_m sin^{|m|}θ/r · e^{im(φ + χ(r) − Ωt)}]` (m = 0: the unfolded
#      `Re[Ψ_0/r]`) with respect to `(t, r, θ, φ)` at the particle. It never writes the
#      dephasing term, the ×2, the phase or the θ-chain down — they emerge from
#      differentiating `Φ`. It is what gives `F_θ` an oracle at all: spectre HARD-CODES
#      `F_θ = 0` (`ObserveSelfForce.cpp:49-51`) and publishes no value to compare against.
#   2. REGRESSION CHECKS, NOT INDEPENDENT ORACLES: `sf_spectre_force` (a line-by-line port of
#      `ObserveSelfForce.cpp:52-77`) and `sf_closed_form` (the same algebra transcribed from
#      `arXiv:2206.07031 main.tex:178-184`). They are transcriptions of the SAME algebra the
#      implementation transcribes, so agreement between them is a transcription check — it
#      catches a typo, never a shared misreading. Retained for exactly that.
#   3. THE EXTERNAL ANCHOR: spectre's own regression rows
#      `tests/InputFiles/SelfForce/ScalarCircularOrbit.yaml:19-23` at `a = 0.5`, `r₀ = 10`,
#      `m = 1`, column legend
#      `tests/Unit/Elliptic/Systems/SelfForce/Scalar/Events/Test_ObserveSelfForce.cpp:195-199`.
#      They exist at ONE configuration and for the RADIAL component only.
#
# PAPER ANCHORS: `arXiv:2206.07031 src/main.tex:178-184` (Eq. 3.10) and `:89-94` (`Δφ(r)`);
# CRP `arXiv:2606.04998 main.tex:428-430` (`F^self_α = q ∇_αΦ^ℛ`, the normalization),
# `:701-703` (the extraction follows Eq. 3.10), `:719-724` (the resolution floor);
# `src/MModeOperator/ansatz.jl` `mode_term` (the locked ansatz);
# `src/Extraction/regularization.jl` `mmode_reg_term_Fr` (the m = 0 halving rule).
#
# GATE CONSTANTS: named below, ALWAYS separately from the values observed — a threshold is a
# contract, an observation is a measurement of one mesh sequence. Every threshold was
# calibrated on throwaway prototypes in the node evidence
# (`proto{1,2,3}_w12a.jl` for the original landing, `proto4_w12a.jl` for this follow-up,
# re-measured on main ≥ 6df9f51 with the mean reduction) BEFORE this file was written or
# rewritten. ONE gate was deliberately WEAKENED after the original green run — gate (b)'s
# level-index statement became a contiguous-tail statement — because a golden AMR level index
# is not a contract; that change is recorded in the node's `process_incidents.md` §3. No
# threshold has ever been loosened to make a failing assertion pass.
#
# NO mode sum (row 46), NO figures (rows 48-53), NO parallel fan-out (row 47).

using Test
using Printf
using ScalarSelfForce
using EllipticCore
using ScalarSelfForce: kerr_varphi_offset, sinm_theta

# ------------------------------------ GATE CONSTANTS ---------------------------------------
const SF_FD_ORACLE_RTOL    = 1.0e-20   # (a0) vs the independent BigFloat FD of the full field
const SF_FORMULA_RTOL      = 1.0e-12   # (a1)/(a2) the DAG row-44 analytic-field tolerance
const SF_SPECTRE_ATOL      = 2.0e-7    # spectre's OWN tolerance, ScalarCircularOrbit.yaml:22
const SF_SPECTRE_BAND      = 5.0e-8    # (b) OUR regression band on the deviation, mean
                                       #     reduction, pinned config; 4x tighter than theirs
const SF_TRACE_SPREAD_ATOL = 1.0e-10   # (b) mean vs one-sided trace, tied to the DG spread
const SF_LINEARITY_RTOL    = 1.0e-14   # (b) mean-of-ingredients == mean-of-forces
const SF_TWO_STEP_MAX      = 0.5       # (c) two-level contraction factor
const SF_RATIO_BAND        = 2.0       # (c) broad non-divergence band on the 1-step ratios
const SF_RES_ERR_ATOL      = 1.0e-9    # (c) resolution floor, absolute
const SF_RES_ERR_RTOL      = 1.0e-5    # (c) resolution floor, relative
const SF_SLOPE_TARGET      = -6.0      # (e) the Heffernan-subtracted residual's decay
const SF_SLOPE_BAND        = 0.3       # (e) row 45's band
const SF_SLOPE_HEADROOM    = 3.0       # (e) |residual| / level-difference required of a point
const SF_SLOPE_SIGMA_MAX   = 0.1       # (e) propagated slope uncertainty required of a fit

# ---- configs (sf_-prefixed: sibling test files share the runtests.jl include scope) --------
sf_cfg(; spin = 0.5, r0 = 10.0, m = 1, L = 0, P = 4) =
    crp_paper_config(; spin = spin, orbital_radius = r0, m_mode = m,
                     inner_radius = CircularOrbit(spin, r0).rp, P = P, L = L)

# ---- ORACLE 1 (INDEPENDENT): the full folded field, differentiated numerically -------------
# An analytic Ψ_m with genuine r AND z dependence, so ∂_zΨ ≠ 0 and F_θ is a real number.
#
# TWO SCOPE LIMITS, stated so the oracle is not read as more than it is:
#   1. it CALLS `kerr_varphi_offset` and `sinm_theta`, the same sealed functions the
#      implementation calls. What it derives independently is the ASSEMBLY — that χ(r) belongs
#      inside the phase and therefore contributes `+im χ′` to ∂_r, that the fold is ×2, that
#      ∂_t → −imΩ. It does NOT re-derive the DEFINITION of χ(r) or of sin^{|m|}θ; a wrong
#      χ(r) would move `Φ` and `mode_self_force` together and pass. Those definitions are
#      gated where they are owned, in Spacetime and in `ansatz.jl`.
#   2. it differentiates AT θ = π/2, where sinθ = 1. So it confirms the magnitude and sign of
#      the polar term but cannot separate `∂_θ = −sinθ ∂_z` from `∂_θ = −∂_z`; only an
#      off-equatorial evaluation could, and the extraction never needs one.
const SFC = (Complex{BigFloat}(big"0.37", big"-0.21"), Complex{BigFloat}(big"-0.13", big"0.29"),
             Complex{BigFloat}(big"0.41", big"0.17"), Complex{BigFloat}(big"-0.07", big"0.23"),
             Complex{BigFloat}(big"0.19", big"-0.11"))
sf_psi(r, z, r0)    = SFC[1] + SFC[2]*(r-r0) + SFC[3]*z + SFC[4]*z^2 + SFC[5]*(r-r0)*z
sf_psi_dr(r, z, r0) = SFC[2] + SFC[5]*z
sf_psi_dz(r, z, r0) = SFC[3] + 2*SFC[4]*z + SFC[5]*(r-r0)

# Φ(t, r, θ, φ) — the ansatz of `mode_term`, folded. The dephasing term, the ×2, the phase and
# the θ-chain are NOT written here; χ(r) simply sits inside and the derivative is numerical.
function sf_Phi(orb, m, t, r, th, ph)
    z = cos(th)
    varphi = ph + kerr_varphi_offset(orb.a, r)
    term = sf_psi(r, z, orb.r0) * sinm_theta(m, z) / r * cis(m * (varphi - orb.Omega * t))
    return m == 0 ? real(term) : 2 * real(term)
end
sf_fd(f, x, h) = (f(x + h) - f(x - h)) / (2h)

# ---- ORACLE 2 (REGRESSION): spectre `ObserveSelfForce.cpp:52-77`, statement by statement ----
# `field` = Ψ^ℛ_m(x_p); `dfield` = ∂_{x0}Ψ^ℛ_m(x_p), x0 = r when `PenetratingHorizon: True`
# (the CRP setting and ours). Returns spectre's COMPLEX `self_force[0]` = our `Fr_complex`.
function sf_spectre_force(a, r0, M, m::Integer, field, dfield; penetrating = true)
    r_plus = M * (1 + sqrt(1 - (a / M)^2))
    r_minus = M * (1 - sqrt(1 - (a / M)^2))
    sf = dfield / r0
    if !penetrating                                     # ∂_{r*} → ∂_r, α = Δ/(r²+a²)
        sf /= 1 - 2 * M * r0 / (r0^2 + a^2)
    end
    sf -= field / r0^2
    if m > 0
        delta_phi = m * a / (r_plus - r_minus) * log((r0 - r_plus) / (r0 - r_minus))
        delta_phi_dr = m * a / (r0 - r_minus) / (r0 - r_plus)
        sf += complex(zero(delta_phi_dr), delta_phi_dr / r0) * field
        sf *= 2 * (cos(delta_phi) + im * sin(delta_phi))
    end
    return sf
end

# ---- ORACLE 2b (REGRESSION): the same algebra transcribed from the tex ---------------------
function sf_closed_form(a, r0, m::Integer, psi, dpsi_r, dpsi_z; phi_p = 0)
    rp = 1 + sqrt(1 - a^2)
    rm = 1 - sqrt(1 - a^2)
    Δ = r0^2 - 2 * r0 + a^2
    dphi = a == 0 ? zero(float(a)) : a / (rp - rm) * log((r0 - rp) / (r0 - rm))
    Ω = 1 / (a + sqrt(r0^3))
    br = dpsi_r / r0 - psi / r0^2 + im * m * (a / Δ) * psi / r0
    bth = -dpsi_z / r0
    if m == 0
        return (Fr = real(dpsi_r / r0 - psi / r0^2), Ft = zero(real(psi)),
                Fphi = zero(real(psi)), Ftheta = real(bth))
    end
    ph = exp(im * m * (phi_p + dphi))
    return (Fr = real(2 * ph * br), Ft = real(2 * ph * (-im * m * Ω) * psi / r0),
            Fphi = real(2 * ph * (im * m) * psi / r0), Ftheta = real(2 * ph * bth))
end

# ---- fill a ModeSolution's nodal vector with the collocation values of an analytic field ---
function sf_fill_collocation!(sol, f)
    mesh = sol.mesh
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        nx, ny = mesh.orders[e]
        em = element_matrices(Float64, mesh.blocks[mesh.blockindex[el.block]], el,
                              mesh.orders[e])
        off = mesh.offsets[e]
        for j in 1:ny, i in 1:nx
            sol.u[off + tensor_lin(i, j, nx)] = f(em.x[1][i], em.x[2][j])
        end
    end
    return sol
end

# =========== (a0) THE INDEPENDENT ORACLE: all four components against a FD of Φ ==============

@testset "(a0) all four components against a BigFloat finite difference of the full field" begin
    setprecision(BigFloat, 240) do
        worst = big(0.0)
        wtheta = big(0.0)
        for a in (big"0.0", big"0.5", big"-0.9"), r0 in (big"6.0", big"10.0"),
            m in (0, 1, 2, 3, 5), phi_p in (big"0.7", big"-1.3")
            orb = CircularOrbit(a, r0)
            h = BigFloat(10)^-25
            thp = big(pi) / 2
            got = mode_self_force(orb, m, sf_psi(orb.r0, big(0.0), orb.r0),
                                  sf_psi_dr(orb.r0, big(0.0), orb.r0),
                                  sf_psi_dz(orb.r0, big(0.0), orb.r0); phi_p = phi_p)
            Fr  = sf_fd(r  -> sf_Phi(orb, m, big(0.0), r, thp, phi_p), orb.r0, h)
            Ft  = sf_fd(t  -> sf_Phi(orb, m, t, orb.r0, thp, phi_p), big(0.0), h)
            Fth = sf_fd(th -> sf_Phi(orb, m, big(0.0), orb.r0, th, phi_p), thp, h)
            Fph = sf_fd(ph -> sf_Phi(orb, m, big(0.0), orb.r0, thp, ph), phi_p, h)
            rel(x, y) = abs(x - y) / max(abs(y), big(1e-40))
            worst = max(worst, rel(got.Fr, Fr), rel(got.Ft, Ft), rel(got.Fphi, Fph))
            wtheta = max(wtheta, rel(got.Ftheta, Fth))
            # F_θ is a REAL number here (the analytic Ψ has genuine z dependence), so this is
            # a live oracle for it, not a comparison of two zeros
            @test abs(Fth) > big(1e-3)
        end
        @printf("mode-self-force (a0) observed worst rel dev vs the FD oracle: %.3e (F_r/F_t/F_φ), \
                 %.3e (F_θ); gate %.1e\n", Float64(worst), Float64(wtheta), SF_FD_ORACLE_RTOL)
        @test worst <= SF_FD_ORACLE_RTOL
        @test wtheta <= SF_FD_ORACLE_RTOL
    end
end

@testset "(a0b) the FD oracle DISCRIMINATES each convention piece" begin
    setprecision(BigFloat, 240) do
        orb = CircularOrbit(big"0.5", big"10.0")
        h = BigFloat(10)^-25
        phi_p = big"0.7"
        m = 2
        Fr = sf_fd(r -> sf_Phi(orb, m, big(0.0), r, big(pi)/2, phi_p), orb.r0, h)
        p = sf_psi(orb.r0, big(0.0), orb.r0)
        dr = sf_psi_dr(orb.r0, big(0.0), orb.r0)
        ph = cis(m * (phi_p + kerr_varphi_offset(orb.a, orb.r0)))
        brad = dr / orb.r0 - p / orb.r0^2
        rel(x) = Float64(abs(x - Fr) / abs(Fr))
        nodeph = real(2 * ph * brad)                       # dephasing term dropped
        nox2   = real(ph * (brad + im*m*(orb.a/(orb.r0^2 - 2*orb.r0 + orb.a^2))*p/orb.r0))
        nophse = real(2 * (brad + im*m*(orb.a/(orb.r0^2 - 2*orb.r0 + orb.a^2))*p/orb.r0))
        @printf("mode-self-force (a0b) observed falsifier margins vs the FD oracle: \
                 no-dephasing %.3e, no-x2 %.3e, no-phase %.3e\n",
                rel(nodeph), rel(nox2), rel(nophse))
        @test rel(nodeph) > 1.0e-2
        @test rel(nox2) > 1.0e-1
        @test rel(nophse) > 1.0e-2
    end
end

# ================= (a1)-(a4) transcription regression checks and structure ===================

# spectre's non-penetrating branch is the one thing the port has that we deliberately do not:
# assert that it is a DIFFERENT number, so "we skip /= alpha" is a decision, not an omission.
function got_penetrating_only()
    a, r0, m = 0.5, 10.0, 2
    psi, dr = 0.31 - 0.47im, -0.13 + 0.29im
    pen = sf_spectre_force(a, r0, 1.0, m, psi, dr; penetrating = true)
    tor = sf_spectre_force(a, r0, 1.0, m, psi, dr; penetrating = false)
    return pen != tor && mode_self_force(CircularOrbit(a, r0), m, psi, dr, 0.0).Fr_complex ≈ pen
end

@testset "(a1) the formula against two TRANSCRIPTION regression checks" begin
    worst_c, worst_s = 0.0, 0.0
    for a in (0.0, 0.5, -0.9), m in (0, 1, 2, 5), r0 in (10.0, 6.0)
        orb = CircularOrbit(a, r0)
        psi, dr, dz = 0.31 - 0.47im, -0.13 + 0.29im, 0.0
        got = mode_self_force(orb, m, psi, dr, dz)
        ref = sf_closed_form(a, r0, m, psi, dr, dz)
        worst_c = max(worst_c, abs(got.Fr - ref.Fr) / abs(ref.Fr),
                      abs(got.Ft - ref.Ft) / max(abs(ref.Ft), 1e-300),
                      abs(got.Fphi - ref.Fphi) / max(abs(ref.Fphi), 1e-300))
        # spectre's COMPLEX self_force[0] IS our `Fr_complex` (the pair its rows print)
        s = sf_spectre_force(a, r0, 1.0, m, psi, dr)
        worst_s = max(worst_s, abs(got.Fr_complex - s) / abs(s))
    end
    @printf("mode-self-force (a1) observed worst rel dev vs the tex transcription %.3e ; \
             vs the spectre port %.3e ; gate %.1e\n", worst_c, worst_s, SF_FORMULA_RTOL)
    @test worst_c <= SF_FORMULA_RTOL
    @test worst_s <= SF_FORMULA_RTOL
    @test got_penetrating_only()
end

@testset "(a2) analytic field on a REAL ModeSolution mesh (collocation fill)" begin
    orb = CircularOrbit(0.5, 10.0)
    xp = puncture_chart_position(block_domain(orb, sf_cfg(; m = 1)), orb)
    @test xp == (10.0, 0.0)
    f(x, y)  = (0.7 - 0.3im) + (0.11 + 0.05im) * (x - 10) - (0.23 - 0.4im) * y +
               (0.017 + 0.009im) * (x - 10)^2 * y + (0.004 - 0.002im) * (x - 10)^3
    fx(x, y) = (0.11 + 0.05im) + 2 * (0.017 + 0.009im) * (x - 10) * y +
               3 * (0.004 - 0.002im) * (x - 10)^2
    psi_e, dr_e = f(xp[1], xp[2]), fx(xp[1], xp[2])
    worst = 0.0
    for m in (0, 1, 2, 5)
        solm = solve_mode(orb, sf_cfg(; m = m), m)
        @test solm.dom.chart === :z2
        sf_fill_collocation!(solm, f)
        got = mode_self_force(solm)                       # the MEAN reduction (default)
        ref = sf_closed_form(0.5, 10.0, m, psi_e, dr_e, 0.0)
        worst = max(worst, abs(got.Fr - ref.Fr) / abs(ref.Fr))
        @test got.reduction === :mean
        @test got.psi ≈ psi_e atol = 1e-13
        @test got.dpsi_r ≈ dr_e atol = 1e-13
        @test got.dpsi_z == 0                      # ∂_z = 2z∂_w with z = 0: EXACTLY zero
        @test got.Ftheta == 0                      # ... hence F_θ ≡ 0 bitwise
        @test got.m == m
    end
    @printf("mode-self-force (a2) observed worst rel dev on the mesh path: %.3e ; gate %.1e\n",
            worst, SF_FORMULA_RTOL)
    @test worst <= SF_FORMULA_RTOL
end

@testset "(a3) structural falsifiers on the raw-ingredient formula" begin
    orb = CircularOrbit(0.5, 10.0)
    psi, dr = 0.31 - 0.47im, -0.13 + 0.29im
    full = mode_self_force(orb, 2, psi, dr, 0.0)
    # the dephasing term vanishes IDENTICALLY at a = 0 (the term IS a/Δ)
    @test mode_self_force(CircularOrbit(0.0, 10.0), 2, psi, dr, 0.0).Fr_complex ==
          2 * (dr / 10 - psi / 100)
    # m = 0 is UNFOLDED: no ×2, no phase, no dephasing, no F_t / F_φ
    f0 = mode_self_force(orb, 0, psi, dr, 0.0)
    @test f0.Fr_complex == dr / 10 - psi / 100
    @test f0.Ft == 0 && f0.Fphi == 0
    @test f0.phase == 1
    # the phase is e^{+im(φ_p+χ)}, not its conjugate
    @test mode_self_force(orb, 2, psi, dr, 0.0; phi_p = 0.31).Fr !=
          mode_self_force(orb, 2, psi, dr, 0.0; phi_p = -0.31).Fr
    @test full.reduction === :direct && full.element == 0 && isempty(full.elements)
end

@testset "(a4) F_t, F_φ and F_θ structure" begin
    orb = CircularOrbit(0.5, 10.0)
    psi, dr = 0.31 - 0.47im, -0.13 + 0.29im
    for m in (0, 1, 2, 5)
        g = mode_self_force(orb, m, psi, dr, 0.0)
        @test g.Ft_complex == -orb.Omega * g.Fphi_complex     # by construction; F_t = −ΩF_φ
        @test g.Ft == -orb.Omega * g.Fphi
        @test g.Ftheta == 0                                   # dpsi_z = 0 ⇒ bitwise zero
        if m > 0
            @test g.Ft_complex ≈ 2 * cis(m * kerr_varphi(orb.a, orb.r0, 0)) *
                                 (-im * m * orb.Omega) * psi / orb.r0 rtol = 1e-14
        else
            @test g.Ft_complex == 0
        end
    end
    g = mode_self_force(orb, 1, psi, dr, 0.17 + 0.4im)
    @test g.Ftheta ≈ real(2 * cis(kerr_varphi(orb.a, orb.r0, 0)) *
                          (-(0.17 + 0.4im) / orb.r0)) rtol = 1e-14
end

# ================= the ONE production sweep, shared by gates (b), (c), (d) ===================
const SF_ORB = CircularOrbit(0.5, 10.0)
const SF_HIST = solve_mode_amr(SF_ORB, sf_cfg(; m = 1), 1; iterations = 8)
const SF_LEVELS = mode_self_force_history(SF_HIST)

# =========================== (b) THE SPECTRE REGRESSION ROWS =================================
# spectre's rows are its own coarse-mesh regression on a DIFFERENT lattice (that yaml is
# `ImposeEquatorialSymmetry: False`, a full-z 5×3 lattice on block bounds
# [r₊, 5, 6, 14, 25, 1e4], `Amr: Iterations: 2`, 480 → 975 grid points), and its 2e-7
# absolute tolerance is a build-to-build reproducibility bound on ITS numbers. Only the
# ITERATION-2 (final) row is used as an anchor: it is the only one of the two that its own
# sequence has converged towards, and GLOBAL_DAG row 44 quotes both rows at 2e-7 without
# saying so. At MATCHING AMR ITERATION the two sequences do not agree — the deviation is
# PRINTED below as documentation of the lattice difference, and deliberately NOT asserted: a
# "we disagree" assertion is green for any wrong code.
const SF_SPECTRE_ITER2 = (-1.365e-4, 9.138e-5, 0.0, 0.0)
const SF_SPECTRE_ITER1 = (-9.183e-5, 9.451e-5)

@testset "(b) the spectre regression rows at a = 0.5, r0 = 10, m = 1" begin
    sre, sim, sthr, sthi = SF_SPECTRE_ITER2
    mean = mode_self_force(SF_HIST.solution)                  # the DEFAULT: spectre parity
    @test mean.reduction === :mean
    dre = abs(real(mean.Fr_complex) - sre)
    dim = abs(imag(mean.Fr_complex) - sim)
    @printf("mode-self-force (b) observed converged F_r^{m=1} (mean of %d elements) = \
             %+.12e %+.12e im ; spectre iter-2 row %+.4e %+.4e ; |dRe| %.5e |dIm| %.5e \
             (their atol %.1e, our band %.1e)\n", length(mean.elements),
            real(mean.Fr_complex), imag(mean.Fr_complex), sre, sim, dre, dim,
            SF_SPECTRE_ATOL, SF_SPECTRE_BAND)
    @test dre <= SF_SPECTRE_ATOL                      # spectre's own tolerance
    @test dim <= SF_SPECTRE_ATOL
    @test dre <= SF_SPECTRE_BAND                      # OUR tighter regression band
    @test dim <= SF_SPECTRE_BAND
    @test mean.Ftheta == sthr && imag(mean.Ftheta_complex) == sthi   # their F_θ columns

    # the reduction IS spectre's: mean over every containing element. Averaging the
    # INGREDIENTS (what the implementation does) equals averaging the FORCES, by linearity.
    per_el = [mode_self_force(SF_HIST.solution; element = e).Fr for e in mean.elements]
    mean_of_forces = sum(per_el) / length(per_el)
    spread = maximum(per_el) - minimum(per_el)
    # `:mean` is the level route's DEFAULT now (spectre's reduction, as for the live
    # `mode_self_force(sol)`); a single-element trace is asked for explicitly.
    trace = mode_self_force(SF_HIST; level = lastindex(SF_HIST.levels),
                            reduction = :trace)
    @printf("mode-self-force (b) observed elements %d ; mean %+.12e ; first trace %+.12e ; \
             spread %.4e ; mean(ingredients)−mean(forces) %.3e\n", length(mean.elements),
            mean.Fr, trace.Fr, spread, abs(mean.Fr - mean_of_forces))
    @test length(mean.elements) >= 2
    @test abs(mean.Fr - mean_of_forces) <= SF_LINEARITY_RTOL * abs(mean.Fr)
    @test abs(mean.Fr - trace.Fr) <= SF_TRACE_SPREAD_ATOL      # tied to the DG trace spread
    @test spread <= SF_TRACE_SPREAD_ATOL
    @test trace.reduction === :trace && trace.element in mean.elements
    # `reduction = :mean` on the final level of the history IS the ModeSolution reduction
    @test mode_self_force(SF_HIST; reduction = :mean).Fr == mean.Fr
    # `reduction = :mean` is now available at EVERY level, not only where a `ModeSolution`
    # survives: `AmrLevel` stores one trace per containing element and their mean (row 43
    # seq 6), so the spectre reduction no longer needs the solution vector.
    @test mode_self_force(SF_HIST; level = 1, reduction = :mean).reduction === :mean
    @test mode_self_force(SF_HIST; level = 1, reduction = :mean).element == 0
    @test_throws ArgumentError mode_self_force(SF_HIST; level = 1, reduction = :median)

    # documentation only (printed, not asserted): the matching-iteration lattice difference
    @printf("mode-self-force (b) observed matching-iteration deviation, our level 1 vs \
             spectre iter 1: %.3e — a different lattice, not a disagreement about physics\n",
            abs(real(SF_LEVELS[2].Fr_complex) - SF_SPECTRE_ITER1[1]))
    # the sequence ENTERS the gate and stays there (a tail, not a golden level index)
    inside = [k - 1 for k in eachindex(SF_LEVELS)
              if abs(real(SF_LEVELS[k].Fr_complex) - sre) <= SF_SPECTRE_ATOL]
    @printf("mode-self-force (b) observed AMR levels inside the 2e-7 gate: %s\n", string(inside))
    @test !isempty(inside) && last(inside) == 8
    @test inside == collect(first(inside):8)
    @test length(inside) >= 4
end

# ============================== (c) AMR CONVERGENCE ==========================================

@testset "(c) two-level contraction envelope of F_r^m over the AMR sweep" begin
    cv = mode_force_convergence(SF_HIST)
    @test length(cv.values) == 9 && cv.steps == collect(0:8)
    @printf("mode-self-force (c) observed F_r^m per level: %s\n",
            join((@sprintf("%+.9e", v) for v in cv.values), " "))
    @printf("mode-self-force (c) observed diffs: %s\n",
            join((@sprintf("%.4e", d) for d in cv.diffs), " "))
    @printf("mode-self-force (c) observed 1-step ratios: %s (max %.4f, band %.1f)\n",
            join((@sprintf("%.4f", r) for r in cv.ratios), " "), maximum(cv.ratios),
            SF_RATIO_BAND)
    @printf("mode-self-force (c) observed 2-step ratios: %s (max %.4f, gate %.1f)\n",
            join((@sprintf("%.4f", r) for r in cv.two_step_ratios), " "),
            maximum(cv.two_step_ratios), SF_TWO_STEP_MAX)
    @printf("mode-self-force (c) observed resolution error: %.4e absolute, %.4e relative \
             (gates %.1e / %.1e)\n", cv.resolution_error, cv.relative_resolution_error,
            SF_RES_ERR_ATOL, SF_RES_ERR_RTOL)
    # THE contraction statement: every two-level ratio contracts. The one-step ratio of an
    # exponentially convergent sequence with plateaus is NOT bounded below 1 — it is held only
    # to a broad non-divergence band, and its observed maximum is historical, not a contract.
    @test all(<(SF_TWO_STEP_MAX), cv.two_step_ratios)
    @test maximum(cv.ratios) < SF_RATIO_BAND
    @test (cv.diffs[end] / cv.diffs[1])^(1 / (length(cv.diffs) - 1)) < SF_TWO_STEP_MAX
    @test cv.resolution_error <= SF_RES_ERR_ATOL
    @test cv.relative_resolution_error <= SF_RES_ERR_RTOL
    @test_throws ArgumentError mode_force_convergence(SF_HIST; component = :nope)
end

@testset "(c2) mode_force_convergence refuses a table too short to contract" begin
    short = solve_mode_amr(SF_ORB, sf_cfg(; m = 1), 1; iterations = 2)   # 3 levels
    @test length(short.levels) == 3
    @test_throws ArgumentError mode_force_convergence(short)
    @test length(mode_self_force_history(short)) == 3       # the sequence itself is fine
end

@testset "(d) F_θ, F_t and F_φ on the production sweep" begin
    conv = SF_LEVELS[end]
    @test conv.dpsi_z == 0                          # ∂_z = 2z∂_{z²}, z = 0 at the equator
    @test conv.Ftheta == 0                          # ... so F_θ is a BITWISE zero
    @test conv.Ftheta_complex == 0
    @test conv.Ft == -SF_ORB.Omega * conv.Fphi      # circular-orbit identity, per mode
    @printf("mode-self-force (d) observed F_r^{m=1} = %+.12e  F_t^{m=1} = %+.12e  \
             F_φ^{m=1} = %+.12e\n", conv.Fr, conv.Ft, conv.Fphi)
    @test conv.Ft_complex ≈ 2 * cis(kerr_varphi(SF_ORB.a, SF_ORB.r0, 0)) *
                            (-im * SF_ORB.Omega) * conv.psi / SF_ORB.r0 rtol = 1e-13
end

# rebuild a RunConfig with one nested field replaced (the config structs are immutable)
function sf_setfield(cfg::RunConfig{T}, group::Symbol, field::Symbol, value) where {T}
    g = getfield(cfg, group)
    G = typeof(g)
    g2 = G((f === field ? value : getfield(g, f) for f in fieldnames(G))...)
    return RunConfig{T}((f === group ? g2 : getfield(cfg, f)
                         for f in fieldnames(RunConfig{T}))...)
end

@testset "(d2) the radial chart contract is checked, not assumed" begin
    # spectre's `/= alpha` (ObserveSelfForce.cpp:61-65) converts ∂_{r_*} to ∂_r in its
    # NON-penetrating chart. Ours is horizon-penetrating BL r, so the factor is absent — and a
    # configuration that ever declared otherwise must THROW, not return a scaled number.
    @test sf_cfg(; m = 1).orbit.penetrating_horizon        # the production config is penetrating
    sol = SF_HIST.solution
    good = sol.config
    bad = sf_setfield(good, :orbit, :penetrating_horizon, false)
    @test !bad.orbit.penetrating_horizon
    sol.config = bad
    try
        @test_throws ArgumentError mode_self_force(sol)
    finally
        sol.config = good
    end
    # the AmrHistory carries its OWN config, so the guard is checked there too
    @test_throws ArgumentError ScalarSelfForce._sf_check_chart(bad)
    @test ScalarSelfForce._sf_check_chart(SF_HIST.config) === nothing
    @test mode_self_force(sol).reduction === :mean          # restored, and working again
end

# ============================== (e) THE m-SWEEP ==============================================
# m ∈ {0, 1, 2, 5, 10} at 6 AMR levels; after subtracting F^m_{r[4]} (through
# `mmode_reg_term_Fr`, which owns the m = 0 halving) the residual should fall like m^{-6}.
#
# ROW 45's GATE IS **NOT** DISCHARGED HERE and this file does not claim it is. What IS gated is
# the slope over the pair of modes whose residuals are RESOLVED at 6 levels — a point is
# resolved when |residual| is at least SF_SLOPE_HEADROOM times the level-to-level difference —
# together with the uncertainty that difference propagates into the slope. The largest-m point
# available here (m = 10) does NOT meet that bar at 6 levels; its slope is printed with its
# uncertainty and is not asserted. A discharge of row 45 needs the row-47 fan-out: 8 levels,
# m up to 20, and this same propagated uncertainty.

@testset "(e) m-sweep: the Heffernan-subtracted residual decay, on resolved modes only" begin
    ms = (0, 1, 2, 5, 10)
    frs, res, ldiff = Float64[], Float64[], Float64[]
    for m in ms
        h = m == 1 ? SF_HIST : solve_mode_amr(SF_ORB, sf_cfg(; m = m), m; iterations = 6)
        v = m == 1 ? SF_LEVELS : mode_self_force_history(h)
        k = m == 1 ? 7 : lastindex(v)                    # AMR level 6 in both cases
        f, d = v[k].Fr, abs(v[k].Fr - v[k-1].Fr)
        r = f - mmode_reg_term_Fr(4, SF_ORB, m)
        push!(frs, f); push!(res, r); push!(ldiff, d)
        @printf("mode-self-force (e) observed m = %2d : F_r^m %+.9e  term %+.9e  \
                 residual %+.9e  level-diff %.3e  headroom %.2f\n",
                m, f, mmode_reg_term_Fr(4, SF_ORB, m), r, d, abs(r) / d)
    end
    # `regularized_mode_force` IS that subtraction, and it goes through row 45's ONE helper
    @test regularized_mode_force(SF_HIST; level = 7) ≈ res[2]
    @test regularized_mode_force(SF_ORB, 0, 1.0) == 1.0 - mmode_reg_term_Fr(4, SF_ORB, 0)
    @test mmode_reg_term_Fr(4, SF_ORB, 0) == mmode_reg_Fr(4, SF_ORB, 0) / 2
    @test regularized_mode_force(SF_ORB, 3, 1.0) == 1.0 - mmode_reg_Fr(4, SF_ORB, 3)
    @test regularized_mode_force(SF_ORB, 3, 1.0; order = 2) == 1.0 - mmode_reg_Fr(2, SF_ORB, 3)

    headroom = [abs(res[k]) / ldiff[k] for k in eachindex(res)]
    # a two-point log-log slope, with the uncertainty the level differences propagate into it
    function slope2(i, j)
        s = log(abs(res[j]) / abs(res[i])) / log(ms[j] / ms[i])
        ds = (ldiff[i] / abs(res[i]) + ldiff[j] / abs(res[j])) / log(ms[j] / ms[i])
        (s, ds)
    end
    s25, d25 = slope2(3, 4)        # m = 2 → 5, both resolved
    s510, d510 = slope2(4, 5)      # m = 5 → 10, m = 10 not resolved at 6 levels
    lsq(idx) = begin
        xs = [log(Float64(ms[k])) for k in idx]; ys = [log(abs(res[k])) for k in idx]
        mx = sum(xs) / length(xs); my = sum(ys) / length(ys)
        sum((xs .- mx) .* (ys .- my)) / sum((xs .- mx) .^ 2)
    end
    @printf("mode-self-force (e) observed slopes: m∈{2,5} %+.4f ± %.4f (headroom %.0f / %.0f) ; \
             m∈{5,10} %+.4f ± %.4f (headroom %.0f / %.2f) ; least-squares {2,5,10} %+.4f ; \
             {1,2,5,10} %+.4f\n", s25, d25, headroom[3], headroom[4], s510, d510,
            headroom[4], headroom[5], lsq(3:5), lsq(2:5))
    # GATED: the resolved pair only
    @test headroom[3] >= SF_SLOPE_HEADROOM
    @test headroom[4] >= SF_SLOPE_HEADROOM
    @test d25 <= SF_SLOPE_SIGMA_MAX
    @test abs(s25 - SF_SLOPE_TARGET) <= SF_SLOPE_BAND
    # NOT gated, printed above: m = 10's residual and every fit that uses it. Recorded here as
    # the reason — it is a [FUTURE] obligation on row 45, not a passing gate.
    @test headroom[5] < SF_SLOPE_HEADROOM || d510 <= SF_SLOPE_SIGMA_MAX
    # the raw (unsubtracted) sequence is NOT monotone at small m — m = 0 and m = 1 dominate
    @test abs(frs[1]) > abs(frs[3]) && abs(frs[2]) > abs(frs[3])
end

# ============ (i) the `radial_bounds` override: spectre's own lattice, exactly ===============
# spectre's tests/InputFiles/SelfForce/ScalarCircularOrbit.yaml puts its radial block bounds at
# [r+, 5, 6, 14, 25, 1e4] with `HyperboloidalSlicingTransitions: [5., 5., 25., 25.]`. Those
# transitions are NOT the paper's Jinja2 rule (`r0 - 1.5 wr`, `r0 + 3 wr`), so no combination of
# `orbital_radius` and `worldtube_radius` can express that lattice; `radial_bounds` gives the six
# bounds directly and reads the transitions off them, exactly as the yaml's anchors do.
const SF_SPECTRE_BOUNDS = [1.8660254037844386, 5.0, 6.0, 14.0, 25.0, 1.0e4]

@testset "(i) radial_bounds overrides the derivation and nothing else" begin
    rp = SF_ORB.rp
    base = crp_paper_config(; spin = 0.5, orbital_radius = 10.0, m_mode = 1,
                            inner_radius = rp, P = 4, L = 0)
    # (1) THE DEFAULT IS UNCHANGED, BIT FOR BIT: passing the six bounds the Jinja2 rule derives
    # reproduces the omitted-keyword config exactly, field by field.
    derived = [rp, 10.0 - 4.8 * 1.5, 10.0 - 4.8, 10.0 + 4.8, 10.0 + 4.8 * 3, 1.0e4]
    same = crp_paper_config(; spin = 0.5, orbital_radius = 10.0, m_mode = 1,
                            inner_radius = rp, P = 4, L = 0, radial_bounds = derived)
    @test same.domain.block_bounds_r == base.domain.block_bounds_r
    @test same.orbit.hyperboloidal_slicing_transitions ==
          base.orbit.hyperboloidal_slicing_transitions
    @test same.observers.rcc_inner_radius == base.observers.rcc_inner_radius
    # (2) THE OVERRIDE REPRODUCES THE YAML: bounds and both transition pairs
    spec = crp_paper_config(; spin = 0.5, orbital_radius = 10.0, m_mode = 1,
                            inner_radius = rp, outer_radius = 1.0e4,
                            worldtube_radius = 4.0, P = 4, L = 0,
                            equatorial_symmetry = false, radial_bounds = SF_SPECTRE_BOUNDS)
    @test spec.domain.block_bounds_r == SF_SPECTRE_BOUNDS
    @test spec.orbit.hyperboloidal_slicing_transitions == (5.0, 5.0, 25.0, 25.0)
    # the yaml's angular bounds [-1, -0.5, 0.5, 1] are the SIGNED SQUARE ROOTS of the stored
    # z² rows, which `block_lattice_bounds` produces for the full-z chart (row-39 contract)
    @test spec.domain.block_bounds_z2 == [0.0, 0.25, 1.0]
    @test block_lattice_bounds(spec)[2] == [-1.0, -0.5, 0.5, 1.0]
    # (3) NOTHING IS SILENTLY IGNORED: every other keyword is checked against the bounds
    bad(; kw...) = crp_paper_config(; spin = 0.5, orbital_radius = 10.0, m_mode = 1,
                                  inner_radius = rp, P = 4, L = 0, kw...)
    @test_throws ArgumentError bad(radial_bounds = SF_SPECTRE_BOUNDS)          # wrong wr
    @test_throws ArgumentError bad(radial_bounds = [rp, 5.0, 6.0, 14.0, 25.0]) # 5 bounds
    @test_throws ArgumentError bad(radial_bounds = [rp, 6.0, 5.0, 14.0, 25.0, 1.0e4])
    @test_throws ArgumentError bad(worldtube_radius = 4.0,
                                   radial_bounds = [2.0, 5.0, 6.0, 14.0, 25.0, 1.0e4])
    @test_throws ArgumentError bad(worldtube_radius = 4.0,
                                   radial_bounds = [rp, 5.0, 6.0, 14.0, 25.0, 9.0e3])
end

# ============ (h) THE LIKE-FOR-LIKE LATTICE against spectre's regression rows ================
# With `equatorial_symmetry = false` AND `radial_bounds = SF_SPECTRE_BOUNDS` the lattice is
# spectre's ScalarCircularOrbit.yaml one: 5x3 blocks, per-row element pattern [1,1,4,2,2],
# 30 elements, 480 grid points, 4 elements touching the puncture (its NumContributingElements),
# angular bounds [-1, -0.5, 0.5, 1], P = 4, L = 0, the same Sommerfeld outer BC and the same
# DG/solver settings — the same discrete problem as far as the yaml specifies one, which is
# where the comparison is gated. It is called LIKE-FOR-LIKE and not "identical": the reference
# is four printed significant figures, so no comparison against it can resolve more than that.
#
# WHAT THE 480-POINT AGREEMENT DOES AND DOES NOT SAY. spectre publishes `-9.183e-5` and
# `9.451e-5`; the last printed digit is `1e-8`, so a round-to-nearest reading admits anything
# within a half-ULP of `5.0e-9`. Our observed deviations are `5.03505e-09` and `5.08179e-09` —
# AT that half-ULP, marginally OUTSIDE the printed interval (by `3.5e-11` and `8.2e-11`). The
# honest statement is therefore a BRACKET, not a margin: the true residual between the two
# codes is `> ~1e-10` and `<= 1e-8`, and it is limited by the printed reference, not by either
# solver. Quoting a "40x margin" against spectre's 2e-7 would be reading precision the
# reference does not carry. `SF_LIKE4LIKE_ATOL` is set at 10x the printed floor.
#
# WHICH VALUE: the MEAN over all 4 containing elements from the LIVE `ModeSolution` (spectre's
# reduction, ObserveSelfForce.hpp:89-107). Never the per-level one-sided trace.
#
# GRID CORRESPONDENCE: spectre's iteration 0 is its un-solved initial guess; its iteration 1 is
# the first solve on the 480-point mesh (our AMR level 0); its iteration 2 is the mesh after one
# AMR step (our AMR level 1) — where the two grids DIVERGE, 921 points against 975, so spectre's
# own tolerance is not the right criterion there and that grid is reported, not gated.
const SF_FULLZ_ITER1 = (-9.183e-5, 9.451e-5)      # spectre [1, 480, 4]
const SF_FULLZ_ITER2 = (-1.365e-4, 9.138e-5)      # spectre [2, 975, 4]
const SF_SPECTRE_PRINTED_ULP = 1.0e-8   # the last printed digit of the reference rows
const SF_LIKE4LIKE_ATOL = 5.0e-8   # OUR band, 10x the printed floor and 4x tighter than
                                   # spectre's own 2e-7 (which is also asserted)

sf_fullz(bounds) = crp_paper_config(; spin = 0.5, orbital_radius = 10.0, m_mode = 1,
                                    inner_radius = SF_ORB.rp, outer_radius = 1.0e4,
                                    worldtube_radius = bounds === nothing ? 4.8 : 4.0,
                                    P = 4, L = 0, equatorial_symmetry = false,
                                    radial_bounds = bounds)

@testset "(h) the like-for-like lattice against spectre's rows" begin
    cfg = sf_fullz(SF_SPECTRE_BOUNDS)
    dom = block_domain(SF_ORB, cfg)
    xp = puncture_chart_position(dom, SF_ORB)
    # (i) STRUCTURAL: spectre's lattice, element for element
    @test dom.chart === :z
    @test (dom.nr, dom.nz) == (5, 3)
    @test cfg.domain.block_bounds_r == SF_SPECTRE_BOUNDS
    @test length(dom.mesh.elements) == 30
    @test dom.mesh.ndofs == 480
    @test length(ScalarSelfForce._dr_containing(dom.mesh, xp[1], xp[2])) == 4
    per_block = Dict{Int,Int}()
    for el in dom.mesh.elements
        per_block[el.block] = get(per_block, el.block, 0) + 1
    end
    for j in 0:dom.nz-1                       # the yaml's 15 RefinedLevels patches, per row
        @test [get(per_block, i + dom.nr * j, 0) for i in 0:dom.nr-1] == [1, 1, 4, 2, 2]
    end

    # (ii) GRID 0 == spectre iteration 1, on the like-for-like lattice
    sol = solve_mode(SF_ORB, cfg, 1; dom = dom)
    f0 = mode_self_force(sol)
    @test f0.reduction === :mean && length(f0.elements) == 4
    d0re = abs(real(f0.Fr_complex) - SF_FULLZ_ITER1[1])
    d0im = abs(imag(f0.Fr_complex) - SF_FULLZ_ITER1[2])
    per = [mode_self_force(sol; element = e).Fr for e in f0.elements]
    @printf("mode-self-force (h) observed LIKE-FOR-LIKE grid 0 (%d pts, %d elements, mean of %d): \
             %+.12e %+.12e im ; spectre iter 1 %+.4e %+.4e ; dRe %.5e dIm %.5e \
             (our band %.1e, their atol %.1e)\n",
            sol.mesh.ndofs, length(sol.mesh.elements), length(f0.elements),
            real(f0.Fr_complex), imag(f0.Fr_complex), SF_FULLZ_ITER1..., d0re, d0im,
            SF_LIKE4LIKE_ATOL, SF_SPECTRE_ATOL)
    @printf("mode-self-force (h) observed per-element traces Re(F_r): %s ; spread %.4e — the \
             reduction is what makes the mean meet spectre, a single trace would not\n",
            join((@sprintf("%+.9e", v) for v in per), " "), maximum(per) - minimum(per))
    @test d0re <= SF_SPECTRE_ATOL && d0im <= SF_SPECTRE_ATOL   # spectre's own tolerance
    @test d0re <= SF_LIKE4LIKE_ATOL                            # and 10x the printed floor
    @test d0im <= SF_LIKE4LIKE_ATOL
    # the residual is BOUNDED BY THE PRINTED REFERENCE, not resolved by it: four significant
    # figures cannot say more than `<= 1e-8`, and our deviations sit at the half-ULP, a hair
    # outside the printed rounding interval. Assert the bracket, not a margin.
    @printf("mode-self-force (h) observed printed-precision bracket: dRe %.5e dIm %.5e against \
             a half-ULP of %.1e — the residual is > ~1e-10 and <= %.1e, limited by the four \
             significant figures spectre publishes\n",
            d0re, d0im, SF_SPECTRE_PRINTED_ULP / 2, SF_SPECTRE_PRINTED_ULP)
    @test d0re <= SF_SPECTRE_PRINTED_ULP && d0im <= SF_SPECTRE_PRINTED_ULP

    # F_theta AGAINST SPECTRE'S OWN PUBLISHED COLUMNS. This lattice does NOT impose equatorial
    # symmetry, so `dpsi_z` is a genuinely computed number rather than the z²-chart's
    # structural zero, and the comparison against the `0., 0.` of the regression row is a real
    # one — the only external check on F_theta that exists anywhere (spectre hard-codes the
    # value it prints, so this bounds OUR truncation error against ITS convention, and cannot
    # confirm a NON-ZERO F_theta).
    @printf("mode-self-force (h) observed F_theta on the full-z lattice: Re %+.5e Im %+.5e ; \
             spectre's columns 0., 0. at atol %.1e ; dpsi_z %+.5e (a computed number here, \
             not the z^2 chart's structural zero)\n", f0.Ftheta, imag(f0.Ftheta_complex),
            SF_SPECTRE_ATOL, real(f0.dpsi_z))
    @test f0.dpsi_z != 0
    @test abs(f0.Ftheta) <= SF_SPECTRE_ATOL
    @test abs(imag(f0.Ftheta_complex)) <= SF_SPECTRE_ATOL
    # the +z and -z traces agree TO 12 SIGNIFICANT DIGITS: the solution is equatorially
    # symmetric even though the domain does not impose it, so the 4-element mean is really a
    # 2-value mean. Bitwise equality is NOT expected and is not asserted — the mirrored
    # elements have different ids, offsets and hence different summation order.
    @test length(unique(round.(per; sigdigits = 12))) == 2

    # (iii) GRID 1: ours 921 points against spectre's 975 — the grids DIVERGE here, so this is
    # REPORTED, not gated. The cause is 2:1 BALANCE CLOSURE, not coarsening (neither criterion
    # ever emits Join): spectre's 975 decomposes as 20x(4,4) + 4x(4,5) + 23x(5,5)
    # = 320 + 80 + 575 over 47 elements, with block (1,1) split in BOTH dimensions and block
    # (3,1) split in r while carrying p+1 in theta. Ours performs only one of those three
    # closure splits, because our balancer checks the face-normal dimension alone — a row-43
    # defect, in flight there, not an extraction one.
    ScalarSelfForce.amr_step!(sol; flags = amr_flags(sol),
                              options = amr_level_options(cfg, 1))
    f1 = mode_self_force(sol)
    d1re = abs(real(f1.Fr_complex) - SF_FULLZ_ITER2[1])
    d1im = abs(imag(f1.Fr_complex) - SF_FULLZ_ITER2[2])
    @printf("mode-self-force (h) observed grid 1 (%d pts against spectre's 975, %d elements): \
             %+.12e %+.12e im ; spectre iter 2 %+.4e %+.4e ; dRe %.5e dIm %.5e — DIFFERENT \
             grids, reported not gated\n", sol.mesh.ndofs, length(sol.mesh.elements),
            real(f1.Fr_complex), imag(f1.Fr_complex), SF_FULLZ_ITER2..., d1re, d1im)
    @test sol.mesh.ndofs > 480 && f1.reduction === :mean

    # (iv) THE CONTRAST that isolates the cause: the same code on the CRP production bounds
    solc = solve_mode(SF_ORB, sf_fullz(nothing), 1)
    fc = mode_self_force(solc)
    @printf("mode-self-force (h) observed CRP-bounds grid 0 (%d pts): dRe %.5e dIm %.5e — only \
             the radial block bounds differ, and only the conservative part notices\n",
            solc.mesh.ndofs, abs(real(fc.Fr_complex) - SF_FULLZ_ITER1[1]),
            abs(imag(fc.Fr_complex) - SF_FULLZ_ITER1[2]))
    @test solc.mesh.ndofs == 480 && length(fc.elements) == 4
end

# =========================== API surface: guards and equivalences ============================

@testset "(g) API guards and cross-method equivalence" begin
    @test_throws ArgumentError mode_self_force(SF_ORB, -1, 1.0 + 0im, 1.0 + 0im)
    @test_throws ArgumentError mode_self_force(SF_HIST, 2)       # m ≠ hist.m
    @test_throws ArgumentError mode_self_force(SF_HIST.solution, 2)
    @test_throws BoundsError mode_self_force(SF_HIST; level = 99)
    @test_throws ArgumentError mode_self_force(SF_HIST; reduction = :median)
    # the history's :trace level and the ModeSolution's SAME element agree bitwise
    a = mode_self_force(SF_HIST; reduction = :trace)
    b = mode_self_force(SF_HIST.solution; element = SF_HIST.levels[end].element)
    @test (a.Fr, a.Ft, a.Ftheta, a.Fphi) == (b.Fr, b.Ft, b.Ftheta, b.Fphi)
    @test a.element == b.element == SF_HIST.levels[end].element
    @test a.step == 8 && a.ndofs == SF_HIST.levels[end].mesh.ndofs
    @test length(SF_LEVELS) == length(SF_HIST.levels)
    # a different element is a different one-sided trace
    other = first(e for e in mode_self_force(SF_HIST.solution).elements if e != a.element)
    @test mode_self_force(SF_HIST.solution; element = other).element == other
    @test mode_self_force(SF_HIST.solution; element = other).Fr != a.Fr
end

# ============================== (f) BigFloat =================================================

@testset "(f) BigFloat construction of the extraction formula" begin
    setprecision(BigFloat, 200) do
        orb = CircularOrbit(big"0.5", big"10.0")
        @test orb isa CircularOrbit{BigFloat}
        psi = Complex{BigFloat}(big"0.31", big"-0.47")
        dr = Complex{BigFloat}(big"-0.13", big"0.29")
        g = mode_self_force(orb, 2, psi, dr, zero(psi))
        @test g.Fr isa BigFloat
        @test g.Fr_complex isa Complex{BigFloat}
        ref = sf_closed_form(big"0.5", big"10.0", 2, psi, dr, zero(psi))
        @test abs(g.Fr - ref.Fr) / abs(ref.Fr) < big"1e-50"
        g64 = mode_self_force(CircularOrbit(0.5, 10.0), 2, ComplexF64(psi),
                              ComplexF64(dr), 0.0)
        @test abs(Float64(g.Fr) - g64.Fr) / abs(g64.Fr) < 1e-13
    end
end
