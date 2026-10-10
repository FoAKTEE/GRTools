# MR19/MR20 Newtonian planet-forced disc-mode solver (Miranda & Rafikov, arXiv:1811.09628 linear
# theory + arXiv:2007.13787 cooling extension) — CHARACTERIZATION gates.
#
# The model, its master-equation coefficients and its solver were migrated into `Modes` from the
# legacy `DiscModes` package (newtonian_disc.jl, master_equation_cool_newt.jl, solver_newt.jl,
# integrator.jl → NumericsCore.integrate_grid) together with the two solver variants that were
# FIGURE-LOCAL in the legacy MR19 reproduction scripts (fast_solver.jl → the `PotentialCache`
# throughput path; homog_solver.jl → `passive_mode`/`passive_sigma`).
#
# ORACLE: a one-shot capture of the legacy package on the MR19 fiducial disc
# (test/fixtures/mr19_newtonian.jl). The legacy code was itself verified against the authors'
# Python `planet-disk-modes-cooling` (grid rel 8.9e-16, δh max rel 4.5e-3 at 2000 points), so
# reproducing it to round-off reproduces the published pipeline.
#
# GATES:
#   (N1) background profiles + resonance denominator + WKB coefficients  — bitwise/1e-14 vs oracle
#   (N2) softened Laplace coefficients + the potential harmonics         — 1e-14 vs oracle
#   (N3) ε → 0 identity with the package's unsoftened `laplace_b`        — quadrature accuracy
#   (N4) master-equation coefficients C₁, C₀, Ψ at three (r, m)          — 1e-13 relative
#   (N5) forced mode solve δh, δh′, δΣ/Σ for m = 2, 5, 10 (adiabatic)    — 1e-9 relative
#        (m = 5, 10 only under GRTOOLS_SLOW_TESTS=1 — see `SLOW` below)
#   (N6) the cooling regime β = 0.1 (MR20 `parameters.py` value)         — 1e-9 relative (SLOW)
#   (N7) the passive/homogeneous problem of MR19 §passive (Figs. 10-12)  — 1e-9 relative
#   (N8) `PotentialCache` throughput path vs the exact solve             — 2e-3 relative (mode δh)
#        (an approximation: the table interpolation of Φ, Φ′, Φ″ is gated at 1e-7, 1e-6, 1e-6)

using Test
using AccretionDisc
const M = AccretionDisc.Modes

include(joinpath(@__DIR__, "..", "fixtures", "mr19_newtonian.jl"))

"""
Slow gates (the 2000-point forced solves that ground N5/N6) run only when
`GRTOOLS_SLOW_TESTS=1`. Cost at -t1, MEASURED: m = 2 → 162 s, m = 5 → 305 s, m = 10 → 733 s,
plus the cooling pair — ~30 min for the whole file, which is too slow for every-commit CI. The
DEFAULT run keeps the m = 2 end-to-end solve (the gate that actually exercises the shooting method,
the WKB boundaries and the refinement loop) and skips the higher-m and cooling repeats.
"""
const SLOW = get(ENV, "GRTOOLS_SLOW_TESTS", "0") == "1"

fid(beta) = NewtonianDisc(h0 = 0.1, q = 1.0, p = 1.0, gamma = 7 / 5, beta = beta,
                          soft = 0.6 * 0.1, indirect = 0, pole_disp = 1e-6)

relmax(a, b) = maximum(abs.(a .- b) ./ max.(abs.(b), eps()))

@testset "Newtonian planet-forced disc modes (MR19)" begin
    d = fid(80.0)                                  # adiabatic limit (MR19 has no cooling)
    rs = MR19.rs

    @testset "gate N1: background profiles" begin
        for (nm, f) in (("Sigma", M.Sigma), ("dSigmadr", M.dSigmadr), ("OmegaK", M.OmegaK),
                        ("Omega", M.Omega), ("dOmegadr", M.dOmegadr), ("kappa", M.kappa),
                        ("dkappadr", M.dkappadr), ("csiso", M.csiso), ("dcsisodr", M.dcsisodr),
                        ("csadi", M.csadi), ("Nr2", M.Nr2), ("dNr2dr", M.dNr2dr),
                        ("LSinv", M.LSinv), ("dLSinvdr", M.dLSinvdr), ("LTinv", M.LTinv),
                        ("dLTinvdr", M.dLTinvdr), ("tc", M.tc), ("dtcdr", M.dtcdr))
            got = [f(d, r) for r in rs]
            @test relmax(got, getfield(MR19, Symbol(nm))) < 1e-14
        end
        @test relmax([M.Dfreq(d, r, 5) for r in rs], MR19.Dfreq_m5) < 1e-14
        @test relmax([M.dDdr(d, r, 5) for r in rs], MR19.dDdr_m5) < 1e-14
        @test relmax([M.kwave(d, r, 5) for r in rs], MR19.kwave_m5) < 1e-14
        @test relmax([M.wkbC(d, r, 5) for r in rs], MR19.wkbC_m5) < 1e-14
        # the pattern frequency is the secondary's KEPLERIAN frequency ω_p = m·Ω_K(1) = m, so with a
        # pressure-supported disc (Ω(1) = sqrt(1 − h0²(p+q)) < 1) corotation is DISPLACED from r = 1:
        # D(1) = κ² − m²(Ω(1) − 1)² < κ². Both facts are gated here.
        for m in (1, 3, 7)
            @test M.omega_pat(d, m) == m
            @test isapprox(M.Dfreq(d, 1.0, m),
                           M.kappa(d, 1.0)^2 - m^2 * (M.Omega(d, 1.0) - 1)^2; rtol = 1e-14)
            @test M.Dfreq(d, 1.0, m) < M.kappa(d, 1.0)^2
        end
        @test M.Omega(d, 1.0) < 1                       # pressure support slows the disc
    end

    @testset "gate N2: softened Laplace coefficients / potential harmonics" begin
        @test relmax([M.Phi_m(d, r, 5) for r in rs], MR19.Phi_m5) < 1e-14
        @test relmax([M.dPhidr_m(d, r, 5) for r in rs], MR19.dPhidr_m5) < 1e-14
        @test relmax([M.d2Phidr2_m(d, r, 5) for r in rs], MR19.d2Phidr2_m5) < 1e-14
        # softening regularises the α → 1 divergence
        @test isfinite(laplace_b_soft(0.5, 5, 1.0, 0.06))
        @test laplace_b_soft(0.5, 5, 1.0, 0.06) < laplace_b_soft(0.5, 5, 1.0, 0.01)
        # the indirect term is m = 1 only
        di = NewtonianDisc(h0 = 0.1, q = 1.0, p = 1.0, gamma = 7 / 5, beta = 80.0,
                           soft = 0.06, indirect = 1, pole_disp = 1e-6)
        @test M.Phi_m(di, 0.7, 1) ≈ M.Phi_m(d, 0.7, 1) + 0.7
        @test M.Phi_m(di, 0.7, 2) == M.Phi_m(d, 0.7, 2)
    end

    @testset "gate N3: ε → 0 reduces to the unsoftened laplace_b" begin
        for (m, β) in ((0, 0.4), (1, 0.4), (2, 0.7), (5, 0.3))
            @test isapprox(laplace_b_soft(0.5, m, β, 0.0; nphi = 4000), laplace_b(0.5, m, β);
                           rtol = 1e-11)
        end
    end

    @testset "gate N4: master-equation coefficients" begin
        for (i, (r, m)) in enumerate(MR19.meq_cases)
            Om = M.Omega(d, r); dOm = M.dOmegadr(d, r); Sig = M.Sigma(d, r); dSig = M.dSigmadr(d, r)
            cs = M.csadi(d, r); ot = M.omega_pat(d, m) - m * Om
            Dv = M.Dfreq(d, r, m); dD = M.dDdr(d, r, m)
            LSi = M.LSinv(d, r); dLSi = M.dLSinvdr(d, r)
            LTi = M.LTinv(d, r); dLTi = M.dLTinvdr(d, r)
            N2 = M.Nr2(d, r); dN2 = M.dNr2dr(d, r); t = M.tc(d, r); dt = M.dtcdr(d, r)
            Ph = M.Phi_m(d, r, m); dPh = M.dPhidr_m(d, r, m); d2Ph = M.d2Phidr2_m(d, r, m)
            got = collect(master_equation_cool(r, d.gamma, cs, Om, dOm, Sig, dSig, LTi, dLTi,
                LSi, dLSi, ot, t, dt, Dv, dD, N2, dN2, d.pole_disp, Ph, dPh, d2Ph, m))
            @test relmax(got, MR19.meq[i]) < 1e-13
        end
    end

    x = loggrid(d, MR19.grid_rin, MR19.grid_rout, MR19.grid_n)
    idx = MR19.grid_idx

    @testset "gate N5: forced mode solve (adiabatic β = 80)" begin
        @test relmax(x[idx], MR19.grid_r) < 1e-15
        for m in (SLOW ? (2, 5, 10) : (2,))
            fm = solve_m_mode(d, m, x)
            @test relmax(fm[idx, 1], getfield(MR19, Symbol("dh_m$(m)"))) < 1e-9
            @test relmax(fm[idx, 2], getfield(MR19, Symbol("dhp_m$(m)"))) < 1e-9
            _, _, sig = vars_from_enthalpy(d, fm, x, m)
            frac = [sig[i] / M.Sigma(d, x[i]) for i in eachindex(x)]
            @test relmax(frac[idx], getfield(MR19, Symbol("frac_m$(m)"))) < 1e-9
        end
    end

    @testset "gate N6: cooling regime β = 0.1 (MR20 parameters.py)" begin
        if SLOW
            dc = fid(0.1)
            fm = solve_m_mode(dc, 5, x)
            @test relmax(fm[idx, 1], MR19.dh_cool_m5) < 1e-9
            # cooling DAMPS the outward-propagating wave: the adiabatic run carries more power at r ≈ 3
            j = findfirst(r -> r > 3.0, x)
            fa = solve_m_mode(d, 5, x)
            @test abs(fa[j, 1]) > abs(fm[j, 1])
        else
            @info "gate N6 skipped (two 2000-point forced solves, ~10 min); set GRTOOLS_SLOW_TESTS=1"
            @test_skip abs(solve_m_mode(fid(0.1), 5, x)[100, 1]) > 0
        end
    end

    @testset "gate N7: passive (homogeneous) problem, MR19 §passive" begin
        dh = NewtonianDisc(h0 = 0.05, q = 1.0, p = 1.0, gamma = 7 / 5, beta = 0.0,
                           soft = 0.06, indirect = 0, pole_disp = 1e-6)
        xh = loggrid(dh, MR19.homog_rin, MR19.homog_rout, MR19.homog_n)
        ih = MR19.homog_idx
        @test relmax(xh[ih], MR19.homog_r) < 1e-15
        for m in (5, 10)
            fh = passive_mode(dh, m, xh)
            @test fh[end, 1] ≈ 1.0                       # normalised at the outer edge
            @test relmax(fh[ih, 1], getfield(MR19, Symbol("homog_dh_m$(m)"))) < 1e-9
            sg = passive_sigma(dh, fh, xh, m)
            @test relmax(sg[ih], getfield(MR19, Symbol("homog_sig_m$(m)"))) < 1e-9
        end
    end

    @testset "gate N8: cached-potential path vs the exact solve" begin
        # The cache replaces the per-right-hand-side azimuthal quadrature by a 4-point local
        # interpolation of a tabulated potential, so it is an APPROXIMATION, not an identity. Gate it
        # (a) directly on the interpolation error and (b) on the resulting mode difference, and prove
        # the mode difference is interpolation-controlled by refining the table.
        cache = PotentialCache(d, 5, 0.25, 3.0; n = 3000)
        rs_c = 10 .^ range(log10(0.3), log10(2.5); length = 97)
        eP  = maximum(abs(M.potential_at(cache, r)[1] - M.Phi_m(d, r, 5)) / abs(M.Phi_m(d, r, 5)) for r in rs_c)
        edP = maximum(abs(M.potential_at(cache, r)[2] - M.dPhidr_m(d, r, 5)) / abs(M.dPhidr_m(d, r, 5)) for r in rs_c)
        ed2 = maximum(abs(M.potential_at(cache, r)[3] - M.d2Phidr2_m(d, r, 5)) / abs(M.d2Phidr2_m(d, r, 5)) for r in rs_c)
        @test eP < 1e-7        # MEASURED 4.2e-9
        @test edP < 1e-6       # MEASURED 1.7e-7
        @test ed2 < 1e-6       # MEASURED 1.8e-7

        xs = loggrid(d, 0.3, 2.5, 400)
        f_exact = solve_m_mode(d, 5, xs; rtol = 1e-8, atol = 1e-10)
        sel = findall(i -> abs(f_exact[i, 1]) > 1e-6 * maximum(abs.(f_exact[:, 1])), 1:length(xs))
        f_c3 = solve_m_mode(d, 5, xs, cache; rtol = 1e-8, atol = 1e-10)
        e3 = relmax(f_c3[sel, 1], f_exact[sel, 1])
        # MEASURED 2.0e-4: the shooting problem is a small residual of two large homogeneous
        # solutions, so it amplifies the 1.7e-7 source perturbation by ~1e3. Gate with 10x headroom.
        @test e3 < 2e-3
        # refining the table reduces the difference -> it is interpolation-controlled, not a defect
        f_c12 = solve_m_mode(d, 5, xs, PotentialCache(d, 5, 0.25, 3.0; n = 12_000);
                             rtol = 1e-8, atol = 1e-10)
        @test relmax(f_c12[sel, 1], f_exact[sel, 1]) < e3
    end
end
