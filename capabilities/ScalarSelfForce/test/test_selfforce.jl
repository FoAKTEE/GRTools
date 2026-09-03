# test_selfforce.jl — ESSF node `mode-self-force-extraction` (GLOBAL_DAG row 44, wave 12,
# worker ESSF-W12-A), written BEFORE `src/Extraction/selfforce.jl` exists (TDD; `tests_red.log`
# captured against the missing file). THIS IS THE NUMBER THE PAPER IS ABOUT: F_r^m at the
# particle, out of a converged m-mode solution.
#
# WHAT THE NODE DELIVERS (`Extraction/selfforce.jl`):
#   * `mode_self_force(orbit, m, psi, dpsi_r, dpsi_z; phi_p)` — the pure Osburn–Nishimura
#     Eq. (3.10) formula on raw ingredients (type-generic, BigFloat-clean);
#   * `mode_self_force(sol::ModeSolution)` / `mode_self_force(hist::AmrHistory; level)` —
#     the same, read through the SEALED row-42 `field_at` / `field_grad_at` at the puncture;
#   * `regularized_mode_force(...; order = 4)` — `F_r^m − F^m_{r[order]}` with the sealed
#     m = 0 halving rule, the per-m term the mode sum (row 46) adds up;
#   * `mode_self_force_history(hist)` and `mode_force_convergence(hist)` — the per-AMR-level
#     sequence and its successive-difference table (rows 50 / 53's resolution floor).
#
# ORACLE COLUMN (the DAG row-44 gate column + the READ-ONLY mirrors, transcribed by hand):
#   * Osburn & Nishimura, arXiv:2206.07031, `src/main.tex:178-184` (Eq. 3.10):
#         F^m_α = q ∇_α(Ψ^R_m / r)                                            , m = 0
#         F^m_α = 2q Re[ ∇_α( Ψ^R_m e^{imΔφ} e^{im(φ−Ωt)} / r ) ]             , m > 0
#     with `Δφ(r) = a/(r₊−r₋) ln((r−r₊)/(r−r₋))` (`main.tex:94`) — so `dΔφ/dr = a/Δ` and the
#     radial derivative carries the DEPHASING term `+ i m (a/Δ) Ψ^R_m / r`.
#   * CRP arXiv:2606.04998 `main.tex:428-430` (`F^self_α = q ∇_α Φ^R|_{x_p}`) and
#     `main.tex:701-703` ("We extract the radial self-force contributions F_r^m from each
#     m-mode solution up to m = 20 following Eq. (3.10) in [Osburn_2022]").
#   * spectre `src/Elliptic/Systems/SelfForce/Scalar/Events/ObserveSelfForce.cpp:52-77` — the
#     reference implementation, ported LINE BY LINE below as `sf_spectre_force` (an
#     independent oracle, not a call into ours).
#   * spectre REGRESSION ROWS: `tests/InputFiles/SelfForce/ScalarCircularOrbit.yaml:19-23`
#     at `a = 0.5`, `r₀ = 10`, `m = 1`, `AbsoluteTolerance: 2.0e-7`; the column legend is
#     `tests/Unit/.../Events/Test_ObserveSelfForce.cpp:195-199` =
#     `IterationId | NumGridPoints | NumContributingElements | Re(F_r) | Im(F_r) | Re(F_θ) | Im(F_θ)`:
#         [0, 480, 4,  0.0,       0.0,      0, 0]
#         [1, 480, 4, -9.183e-5,  9.451e-5, 0, 0]
#         [2, 975, 4, -1.365e-4,  9.138e-5, 0, 0]
#     THEIR mesh is not ours (that yaml is `ImposeEquatorialSymmetry: False`, a 5×3 lattice on
#     block bounds `[r₊, 5, 6, 14, 25, 1e4]`, `Amr: Iterations: 2`); the honest comparison is
#     stated in gate (b) and measured, never fudged.
#   * `src/MModeOperator/ansatz.jl` `mode_term` — the LOCKED ansatz
#     `Φ = Σ_m Ψ_m sin^{|m|}θ / r · e^{+im(φ̄ − Ωt)}`, `φ̄ = φ + χ(r)`;
#     `src/Extraction/regularization.jl` `mmode_reg_factor_over_pi` — the LOCKED folding rule
#     (the ×2 is inside the parameter at every m, so an UNFOLDED m = 0 force term subtracts
#     `mmode_reg_Fr(order, orbit, 0)/2`).
#
# PREDECESSOR SCOPE (sealed; consumed, never re-derived here): `amr-driver` (`solve_mode_amr`,
# `AmrHistory`, `AmrLevel`, `puncture_chart_position`), `solver-stack-ssf` (`ModeSolution`,
# `solve_mode`, `field_at`, `field_grad_at`), `mode-ansatz-decomposition` (the conventions),
# `orbits-adapters` (`CircularOrbit`, `Omega`), `heffernan-regularization` (`mmode_reg_Fr`),
# Spacetime (`kerr_varphi`, `kerr_varphi_offset_dr = a/Δ`).
#
# NO mode sum (row 46), NO figures (rows 48-53), NO parallel fan-out (row 47).
#
# GATE CONSTANTS — every one calibrated on the throwaway prototypes
# `evidence/mode-self-force-extraction/proto{1,2,3}_w12a.jl` BEFORE this file was written,
# never tuned afterwards to pass. The measured value is quoted beside each gate.

using Test
using Printf
using ScalarSelfForce
using EllipticCore

# ---- configs (sf_-prefixed: sibling test files share the runtests.jl include scope) --------
sf_cfg(; spin = 0.5, r0 = 10.0, m = 1, L = 0, P = 4) =
    crp_paper_config(; spin = spin, orbital_radius = r0, m_mode = m,
                     inner_radius = CircularOrbit(spin, r0).rp, P = P, L = L)

# ---- ORACLE 1: spectre `ObserveSelfForce.cpp:52-77`, ported statement by statement ---------
# `field` = Ψ^R_m(x_p); `dfield` = ∂_{x0}Ψ^R_m(x_p) with x0 the block's radial coordinate
# (= r when `PenetratingHorizon: True`, which is the CRP setting and ours). Returns spectre's
# COMPLEX `self_force[0]`, i.e. the pair its regression rows print as Re/Im.
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

# ---- ORACLE 2: the hand-differentiated closed form, written from the tex, not from src -----
# Φ^R ⊃ Ψ^R_m(r,z) sin^{|m|}θ / r · e^{im(φ̄ − Ωt)}, φ̄ = φ + Δφ(r), evaluated at
# (t, r, θ, φ) = (0, r₀, π/2, φ_p). ∂_r hits Ψ^R_m, 1/r AND the phase; sin^{|m|}θ = 1 and its
# θ-derivative |m| sin^{|m|−1}θ cosθ vanishes at the equator, so ∂_θ hits only Ψ^R_m.
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

# ============================== (a) THE FORMULA, exactly =====================================

@testset "(a1) the extraction formula against two independent oracles" begin
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
    @test worst_c <= 1e-12          # MEASURED 0.0 (bitwise) on the proto1 grid
    @test worst_s <= 1e-12          # MEASURED 2.03e-16
    @printf("mode-self-force (a1) worst rel dev vs closed form %.3e ; vs spectre port %.3e\n",
            worst_c, worst_s)
end

@testset "(a2) analytic field on a REAL ModeSolution mesh (collocation fill)" begin
    cfg = sf_cfg(; m = 1)
    orb = CircularOrbit(0.5, 10.0)
    sol = solve_mode(orb, cfg, 1)
    @test sol.dom.chart === :z2
    xp = puncture_chart_position(sol.dom, orb)
    @test xp == (10.0, 0.0)
    # a manufactured polynomial of the element's own degree in the DISCRETIZATION chart
    # (x1 = r, x2 = w = z²), so interpolation and differentiation are exact to roundoff
    f(x, y)  = (0.7 - 0.3im) + (0.11 + 0.05im) * (x - 10) - (0.23 - 0.4im) * y +
               (0.017 + 0.009im) * (x - 10)^2 * y + (0.004 - 0.002im) * (x - 10)^3
    fx(x, y) = (0.11 + 0.05im) + 2 * (0.017 + 0.009im) * (x - 10) * y +
               3 * (0.004 - 0.002im) * (x - 10)^2
    sf_fill_collocation!(sol, f)
    psi_e, dr_e = f(xp[1], xp[2]), fx(xp[1], xp[2])

    worst = 0.0
    for m in (0, 1, 2, 5)
        cfgm = sf_cfg(; m = m)
        solm = solve_mode(CircularOrbit(0.5, 10.0), cfgm, m)
        sf_fill_collocation!(solm, f)
        got = mode_self_force(solm)
        ref = sf_closed_form(0.5, 10.0, m, psi_e, dr_e, 0.0)
        worst = max(worst, abs(got.Fr - ref.Fr) / abs(ref.Fr))
        @test got.psi ≈ psi_e atol = 1e-13
        @test got.dpsi_r ≈ dr_e atol = 1e-13
        @test got.dpsi_z == 0                      # ∂_z = 2z∂_w with z = 0: EXACTLY zero
        @test got.Ftheta == 0                      # ... hence F_θ ≡ 0 bitwise
        @test got.m == m
    end
    @test worst <= 1e-12                           # MEASURED ~4.4e-16 (proto1 (B))
    @printf("mode-self-force (a2) worst rel dev on the mesh path: %.3e\n", worst)
end

@testset "(a3) falsifiers: each convention piece changes the answer" begin
    orb = CircularOrbit(0.5, 10.0)
    psi, dr = 0.31 - 0.47im, -0.13 + 0.29im
    full = mode_self_force(orb, 2, psi, dr, 0.0)
    # (i) DROP the dephasing term +im(a/Δ)Ψ/r — an O(a) error
    nodeph = 2 * cis(2 * kerr_varphi(0.5, 10.0, 0)) * (dr / 10 - psi / 100)
    @test abs(real(nodeph) - full.Fr) / abs(full.Fr) > 1e-2      # MEASURED 9.6e-2 on proto1
    # ... and it vanishes identically at a = 0 (the term IS a/Δ)
    orb0 = CircularOrbit(0.0, 10.0)
    @test mode_self_force(orb0, 2, psi, dr, 0.0).Fr_complex ==
          2 * (dr / 10 - psi / 100)
    # (ii) DROP the ×2 real-part folding
    @test abs(full.Fr / 2 - full.Fr) / abs(full.Fr) ≈ 0.5
    # (iii) m = 0 must be UNFOLDED (no ×2, no phase, no dephasing)
    f0 = mode_self_force(orb, 0, psi, dr, 0.0)
    @test f0.Fr_complex == dr / 10 - psi / 100
    @test f0.Ft == 0 && f0.Fphi == 0
    @test abs(f0.Fr - full.Fr) / abs(full.Fr) > 1e-2             # MEASURED 5.5e-1
    # (iv) the phase is e^{+im(φ_p+χ)}, not its conjugate
    @test mode_self_force(orb, 2, psi, dr, 0.0; phi_p = 0.31).Fr !=
          mode_self_force(orb, 2, psi, dr, 0.0; phi_p = -0.31).Fr
end

@testset "(a4) F_t, F_φ and F_θ structure" begin
    orb = CircularOrbit(0.5, 10.0)
    psi, dr = 0.31 - 0.47im, -0.13 + 0.29im
    for m in (0, 1, 2, 5)
        g = mode_self_force(orb, m, psi, dr, 0.0)
        # the ∂_t → −imΩ / ∂_φ → +im rule of the ansatz phase, exactly
        @test g.Ft_complex == -orb.Omega * g.Fphi_complex
        @test g.Ft == -orb.Omega * g.Fphi                     # circular-orbit identity
        @test g.Ftheta == 0                                   # dpsi_z = 0 ⇒ bitwise zero
        if m > 0
            @test g.Ft_complex ≈ 2 * cis(m * kerr_varphi(orb.a, orb.r0, 0)) *
                                 (-im * m * orb.Omega) * psi / orb.r0 rtol = 1e-14
        else
            @test g.Ft_complex == 0
        end
    end
    # a NON-zero ∂_zΨ propagates into F_θ with the ∂_θ = −sinθ ∂_z sign
    g = mode_self_force(orb, 1, psi, dr, 0.17 + 0.4im)
    @test g.Ftheta ≈ real(2 * cis(kerr_varphi(orb.a, orb.r0, 0)) *
                          (-(0.17 + 0.4im) / orb.r0)) rtol = 1e-14
end

# ================= the ONE production sweep, shared by gates (b), (c), (d) ===================
# a = 0.5, r₀ = 10, m = 1 — spectre's own regression configuration — at the paper's 8 AMR
# iterations. MEASURED cost in a cold process: 227 s including JIT; ~20 s warm.
const SF_ORB = CircularOrbit(0.5, 10.0)
const SF_HIST = solve_mode_amr(SF_ORB, sf_cfg(; m = 1), 1; iterations = 8)
const SF_LEVELS = mode_self_force_history(SF_HIST)

# ============================ (b) THE SPECTRE REGRESSION ROWS ================================
# spectre's rows are its OWN coarse-mesh regression on a DIFFERENT lattice (that yaml is
# `ImposeEquatorialSymmetry: False`, a full-z 5×3 lattice on block bounds
# [r₊, 5, 6, 14, 25, 1e4], `Amr: Iterations: 2`, 480 → 975 grid points), and its 2e-7
# absolute tolerance is a build-to-build reproducibility bound on ITS numbers. So the
# comparison is made at three levels of strength, and the honest one is stated in each:
#   (i)   our ARITHMETIC agrees with spectre's to roundoff — gate (a1), 2.03e-16;
#   (ii)  at MATCHING AMR ITERATION the two sequences do NOT agree (different lattices):
#         our level 1 is 3.95e-5 from spectre's iteration-1 row. That is a FINDING, asserted
#         below so a later reader cannot mistake (iii) for a per-iteration match;
#   (iii) our CONVERGED value reproduces spectre's own final (iteration-2) row INSIDE the
#         DAG's 2e-7 absolute gate — MEASURED |ΔRe| 1.708e-8, |ΔIm| 8.665e-10, i.e. inside
#         the 4 significant figures spectre prints.
const SF_SPECTRE_ROWS = ((1, -9.183e-5, 9.451e-5, 0.0, 0.0),
                         (2, -1.365e-4, 9.138e-5, 0.0, 0.0))
const SF_SPECTRE_ATOL = 2.0e-7          # ScalarCircularOrbit.yaml:22

@testset "(b) the spectre regression rows at a = 0.5, r0 = 10, m = 1" begin
    conv = SF_LEVELS[end]
    _, sre, sim, sthr, sthi = SF_SPECTRE_ROWS[2]
    dre = abs(real(conv.Fr_complex) - sre)
    dim = abs(imag(conv.Fr_complex) - sim)
    @printf("mode-self-force (b) converged F_r^{m=1} = %+.12e %+.12e im ; spectre iter 2 row \
             %+.4e %+.4e ; |dRe| %.3e |dIm| %.3e (atol %.1e)\n",
            real(conv.Fr_complex), imag(conv.Fr_complex), sre, sim, dre, dim, SF_SPECTRE_ATOL)
    @test dre <= SF_SPECTRE_ATOL                    # MEASURED 1.708e-8
    @test dim <= SF_SPECTRE_ATOL                    # MEASURED 8.665e-10
    @test conv.Ftheta == sthr && imag(conv.Ftheta_complex) == sthi   # their F_θ columns

    # (ii) the honest non-match at matching AMR iteration — a lattice difference, not a bug
    _, s1re, s1im, _, _ = SF_SPECTRE_ROWS[1]
    d1 = abs(real(SF_LEVELS[2].Fr_complex) - s1re)
    @printf("mode-self-force (b) matching-iteration deviation, level 1 vs spectre iter 1: \
             %.3e (their lattice is a different one)\n", d1)
    @test d1 > SF_SPECTRE_ATOL                      # MEASURED 3.953e-5

    # (iii) the sequence ENTERS the gate and stays there — a contiguous tail ending at the
    # last level. Stated as a tail, not as a fixed first level: the AMR mesh sequence moves
    # slightly with EllipticCore's hp_refine order map, so where the crossing happens is a
    # re-measurement, not a contract. MEASURED here: levels 3-8 (3.5e-8 down to 1.7e-8),
    # level 2 still outside at 4.2e-7.
    inside = [k - 1 for k in eachindex(SF_LEVELS)
              if abs(real(SF_LEVELS[k].Fr_complex) - sre) <= SF_SPECTRE_ATOL]
    @printf("mode-self-force (b) AMR levels inside the 2e-7 gate: %s\n", string(inside))
    @test !isempty(inside) && last(inside) == 8
    @test inside == collect(first(inside):8)        # contiguous: it enters and never leaves
    @test length(inside) >= 4                       # MEASURED 6

    # the one-sided-trace choice is immaterial at production resolution: spectre AVERAGES the
    # extraction over its NumContributingElements (ObserveSelfForce.hpp:97-107, funcl::Plus
    # then funcl::Divides by datum 2); we take the first containing element. MEASURED spread.
    sol = SF_HIST.solution
    xp = SF_HIST.x_p
    es = [e for e in eachindex(sol.mesh.elements)
          if all(abs.(ScalarSelfForce._dr_logical(sol.mesh, e, xp[1], xp[2])) .<= 1 + 16eps())]
    fs = [mode_self_force(sol; element = e).Fr for e in es]
    @printf("mode-self-force (b) containing elements %d ; mean %+.12e ; first %+.12e ; \
             spread %.3e\n", length(es), sum(fs) / length(fs), conv.Fr,
            maximum(fs) - minimum(fs))
    @test length(es) >= 2
    @test abs(sum(fs) / length(fs) - conv.Fr) <= SF_SPECTRE_ATOL   # MEASURED ~1e-11
    # ... and the history's last level IS the final ModeSolution's extraction
    @test mode_self_force(sol).Fr == conv.Fr
end

# ============================== (c) AMR CONVERGENCE ==========================================

@testset "(c) F_r^m converges geometrically over the AMR sweep" begin
    cv = mode_force_convergence(SF_HIST)
    @test length(cv.values) == 9 && cv.steps == collect(0:8)
    @printf("mode-self-force (c) F_r^m per level: %s\n",
            join((@sprintf("%+.9e", v) for v in cv.values), " "))
    @printf("mode-self-force (c) diffs: %s\n",
            join((@sprintf("%.3e", d) for d in cv.diffs), " "))
    @printf("mode-self-force (c) 1-step ratios: %s\n",
            join((@sprintf("%.4f", r) for r in cv.ratios), " "))
    @printf("mode-self-force (c) 2-step ratios: %s\n",
            join((@sprintf("%.4f", r) for r in cv.two_step_ratios), " "))
    @printf("mode-self-force (c) resolution error: %.3e absolute, %.3e relative\n",
            cv.resolution_error, cv.relative_resolution_error)
    # The SINGLE-step ratio is NOT universally < 0.5 — the sealed row-43 contract says so and
    # MEASURED 1.3133 here at level 5. The honest geometric statement is the TWO-step ratio:
    @test all(<(0.5), cv.two_step_ratios[2:end])       # MEASURED max 0.2379
    # ... and the envelope over the whole sweep is geometric with mean ratio well under 0.5
    @test (cv.diffs[end] / cv.diffs[1])^(1 / (length(cv.diffs) - 1)) < 0.5   # MEASURED 0.168
    # final-level stability (the error budget's "numerical resolution floor", main.tex:723)
    @test cv.resolution_error <= 1.0e-9                # MEASURED 4.648e-11
    @test cv.relative_resolution_error <= 1.0e-5       # MEASURED 3.405e-7
    # a wrong component name is an error, not a silent zero
    @test_throws ArgumentError mode_force_convergence(SF_HIST; component = :nope)
end

@testset "(d) F_θ, F_t and F_φ on the production sweep" begin
    conv = SF_LEVELS[end]
    @test conv.dpsi_z == 0                          # ∂_z = 2z∂_{z²}, z = 0 at the equator
    @test conv.Ftheta == 0                          # ... so F_θ is a BITWISE zero
    @test conv.Ftheta_complex == 0
    @test conv.Ft == -SF_ORB.Omega * conv.Fphi      # circular-orbit identity, per mode
    @printf("mode-self-force (d) F_r^{m=1} = %+.12e  F_t^{m=1} = %+.12e  F_φ^{m=1} = %+.12e\n",
            conv.Fr, conv.Ft, conv.Fphi)
    # F_t follows from F_r's own ingredients through the ansatz phase rule ∂_t → −imΩ
    @test conv.Ft_complex ≈ 2 * cis(kerr_varphi(SF_ORB.a, SF_ORB.r0, 0)) *
                            (-im * SF_ORB.Omega) * conv.psi / SF_ORB.r0 rtol = 1e-13
end

# ============================== (e) THE m-SWEEP ==============================================
# m ∈ {0, 1, 2, 5, 10} at 6 AMR levels; after subtracting the Heffernan F^m_{r[4]} (HALVED at
# m = 0, the sealed folding rule of `mmode_reg_factor_over_pi`) the residual should fall like
# m^{-6} (row 45's [FUTURE] slope gate, discharged here on real data).

@testset "(e) m-sweep: the Heffernan-subtracted residual falls like m^-6" begin
    ms = (0, 1, 2, 5, 10)
    frs, res, ldiff = Float64[], Float64[], Float64[]
    for m in ms
        # m = 1 reuses the 8-level sweep TRUNCATED to level 6 — the same mesh sequence
        f, d = if m == 1
            (SF_LEVELS[7].Fr, abs(SF_LEVELS[7].Fr - SF_LEVELS[6].Fr))
        else
            h = solve_mode_amr(SF_ORB, sf_cfg(; m = m), m; iterations = 6)
            v = mode_self_force_history(h)
            (v[end].Fr, abs(v[end].Fr - v[end-1].Fr))
        end
        push!(frs, f); push!(ldiff, d)
        r = m == 1 ? regularized_mode_force(SF_HIST; level = 7) :
                     f - (m == 0 ? mmode_reg_Fr(4, SF_ORB, 0) / 2 : mmode_reg_Fr(4, SF_ORB, m))
        push!(res, r)
        @printf("mode-self-force (e) m = %2d : F_r^m %+.9e  F^m_{r[4]} %+.9e  residual %+.9e \
                 level-diff %.3e\n", m, f,
                m == 0 ? mmode_reg_Fr(4, SF_ORB, 0) / 2 : mmode_reg_Fr(4, SF_ORB, m), r, d)
    end
    slope(idx) = begin
        xs = [log(Float64(ms[k])) for k in idx]; ys = [log(abs(res[k])) for k in idx]
        mx = sum(xs) / length(xs); my = sum(ys) / length(ys)
        sum((xs .- mx) .* (ys .- my)) / sum((xs .- mx) .^ 2)
    end
    s_asym = slope(3:5)     # m = 2, 5, 10
    s_all = slope(2:5)      # m = 1, 2, 5, 10 — m = 1 is not in the asymptotic regime
    s_top = slope(4:5)      # m = 5, 10 — resolution-limited at 6 levels, see below
    @printf("mode-self-force (e) log-log slopes: m∈{2,5,10} %+.4f ; m∈{1,2,5,10} %+.4f ; \
             m∈{5,10} %+.4f\n", s_asym, s_all, s_top)
    @test -6.3 <= s_asym <= -5.7                    # MEASURED -6.1345 (row 45's -6 ± 0.3)
    @test s_all > s_asym                            # m = 1 flattens the fit: MEASURED -5.0607
    # HONEST BOUND on the largest-m point: at 6 AMR levels the m = 10 residual is only a few
    # times the level-to-level difference, so the {5,10} slope is resolution-limited.
    @printf("mode-self-force (e) m = 10 residual / level-diff = %.2f (resolution headroom)\n",
            abs(res[5]) / ldiff[5])
    @test abs(res[5]) > ldiff[5]                    # MEASURED 2.80x
    # the raw (unsubtracted) sequence is NOT monotone at small m — m = 0 and m = 1 dominate
    @test abs(frs[1]) > abs(frs[3]) && abs(frs[2]) > abs(frs[3])
    # `regularized_mode_force` IS the subtraction, with the m = 0 halving
    @test regularized_mode_force(SF_HIST) ≈ SF_LEVELS[end].Fr - mmode_reg_Fr(4, SF_ORB, 1)
    @test regularized_mode_force(SF_ORB, 0, 1.0) == 1.0 - mmode_reg_Fr(4, SF_ORB, 0) / 2
    @test regularized_mode_force(SF_ORB, 3, 1.0) == 1.0 - mmode_reg_Fr(4, SF_ORB, 3)
    @test regularized_mode_force(SF_ORB, 3, 1.0; order = 2) == 1.0 - mmode_reg_Fr(2, SF_ORB, 3)
end

# =========================== API surface: guards and equivalences ============================

@testset "(g) API guards and cross-method equivalence" begin
    @test_throws ArgumentError mode_self_force(SF_ORB, -1, 1.0 + 0im, 1.0 + 0im)
    @test_throws ArgumentError mode_self_force(SF_HIST, 2)       # m ≠ hist.m
    @test_throws ArgumentError mode_self_force(SF_HIST.solution, 2)
    @test_throws BoundsError mode_self_force(SF_HIST; level = 99)
    # the ModeSolution and AmrHistory paths agree bitwise on the final level
    a = mode_self_force(SF_HIST)
    b = mode_self_force(SF_HIST.solution)
    @test (a.Fr, a.Ft, a.Ftheta, a.Fphi) == (b.Fr, b.Ft, b.Ftheta, b.Fphi)
    @test a.element == b.element == SF_HIST.levels[end].element
    @test a.step == 8 && a.ndofs == SF_HIST.levels[end].mesh.ndofs
    @test length(SF_LEVELS) == length(SF_HIST.levels)
    # explicit `element` overrides the default choice
    other = first(e for e in eachindex(SF_HIST.solution.mesh.elements) if e != a.element &&
                  all(abs.(ScalarSelfForce._dr_logical(SF_HIST.solution.mesh, e,
                                                       SF_HIST.x_p[1], SF_HIST.x_p[2])) .<=
                      1 + 16eps()))
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
        # and it AGREES with the Float64 path to Float64 accuracy
        g64 = mode_self_force(CircularOrbit(0.5, 10.0), 2, ComplexF64(psi),
                              ComplexF64(dr), 0.0)
        @test abs(Float64(g.Fr) - g64.Fr) / abs(g64.Fr) < 1e-13
    end
end
