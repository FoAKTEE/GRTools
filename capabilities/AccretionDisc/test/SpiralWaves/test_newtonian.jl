# test_newtonian.jl — EGD::sw-newtonian-limit: the Newtonian-limit pipeline of arXiv:2608.17003
# (src/SpiralWaves/newtonian.jl), i.e. the T_Newt denominator of Fig. 5. Gates, in order:
#   (a) the Newtonian background IS the M/r → 0 limit: Ω = r^{−3/2} and κ = Ω for 𝔭 + 𝔮 = 0 (and then
#       EXACTLY the relativistic χ = 0 rotation law), the closed form κ² = 1/r³ + (𝔭+𝔮)(𝔮−2)c_s²/r²,
#       Ω′ against a finite difference, r_c = a• for 𝔭 + 𝔮 = 0;
#   (b) the Plummer potential modes (60): against an INDEPENDENT dense quadrature written here, the
#       analytic Φ̂′ against a finite difference, and — the real gate — against the RELATIVISTIC
#       singular-field mode Φ = −ĥ_tt/2 as a• → ∞ (P_r, P_t → 1, K_tt → 2);
#   (c) the master system: `newtonian_master_matrices` (Newtonian background) vs
#       `master_matrices_newtonian` (same (50) algebra, relativistic background) vs the RELATIVISTIC
#       `master_matrices` — all three converge as M/a• → 0;
#   (d) THE BALANCE GATE: dF_J/dr = d𝔗/dr (71) for a SOLVED Newtonian mode away from corotation — the
#       one test that ties δΣ̂, δû_φ, the Fourier weight and the sign conventions together;
#   (e) the end-to-end limit: the per-mode torque of the Newtonian pipeline equals the RELATIVISTIC
#       one at a• = 2000M to O(M/a•);
#   (f) χ-independence (bitwise) — the property that makes Fig. 5 a measurement of relativity;
#   (g) exact scalings T ∝ m•², T ∝ ε₀;
#   (h) the survey front end (`model = :newtonian`) reproduces the stand-alone path, and its optimised
#       reduction is BITWISE the `corotation_split` it replaces (also for `model = :relativistic`).
using Test
using AccretionDisc
using AccretionDisc.SpiralWaves

_bgN(; χ = 0.0, aEMRI = 12.0, h0 = 0.05, 𝔭 = 0.0, 𝔮 = 0.0, ε0 = 1.0, b = 0.6) =
    DiscBackground(; χ = χ, aEMRI = aEMRI, h0 = h0, 𝔭 = 𝔭, 𝔮 = 𝔮, ε0 = ε0, b = b)

const GKN = (nL = 8, nR = 8, nC = 200, n_anc = 800)      # small but genuinely adaptive grids

@testset "SpiralWaves Newtonian limit: (50) on a Newtonian background (Fig. 5 denominator)" begin

    @testset "(a) Newtonian background = the M/r → 0 limit" begin
        bg = _bgN()
        for r in (7.0, 12.0, 40.0)
            rot = newtonian_rotation(bg, r)
            @test isapprox(rot.Ω, r^(-1.5); rtol = 1e-14)             # 𝔭 = 𝔮 = 0: Keplerian
            @test isapprox(rot.κ2, rot.Ω^2; rtol = 1e-12)             # κ = Ω
            @test isapprox(rot.Ω, Omega(bg, r); rtol = 1e-13)         # = the relativistic χ = 0 law
            # ω_d vanishes AT corotation (r = a• here), so compare on the scale of ω = mΩ•
            @test abs(newtonian_omega_d(bg, 3, r) - omega_d(bg, 3, r)) < 1e-13 * 3 * newtonian_Omega_emri(bg)
        end
        @test isapprox(newtonian_Omega_emri(bg), 12.0^(-1.5); rtol = 1e-15)
        @test isapprox(newtonian_corotation_radius(bg), 12.0; rtol = 1e-12)

        # pressure-supported disc: the closed forms and the derivative
        for (𝔭, 𝔮) in ((1.5, 1.0), (-1.5, 0.0), (2.0, 1.0))
            bgp = _bgN(; 𝔭 = 𝔭, 𝔮 = 𝔮)
            for r in (8.0, 12.0, 30.0)
                rot = newtonian_rotation(bgp, r)
                cs2 = newtonian_csound(bgp, r)^2
                @test isapprox(rot.Ω^2, 1 / r^3 - (𝔭 + 𝔮) * cs2 / r^2; rtol = 1e-13)
                @test isapprox(rot.κ2, 1 / r^3 + (𝔭 + 𝔮) * (𝔮 - 2) * cs2 / r^2; rtol = 1e-11)
                h = 1e-5 * r
                fd = (newtonian_rotation(bgp, r + h).Ω - newtonian_rotation(bgp, r - h).Ω) / (2h)
                @test isapprox(rot.dΩ, fd; rtol = 1e-8)
            end
            rc = newtonian_corotation_radius(bgp)
            @test isapprox(newtonian_rotation(bgp, rc).Ω, newtonian_Omega_emri(bgp); rtol = 1e-13)
            @test abs(rc / bgp.aEMRI - 1) < 0.01                       # a pressure correction, not a jump
            # The RELATIVISTIC rotation law differs from this one ONLY through the pressure term: the
            # linearised (11) carries Π(1/ν′ − 1/ψ′) = −(𝔭+𝔮)c_s²(r−3) where the Newtonian balance
            # carries −(𝔭+𝔮)c_s² r, so ΔΩ/Ω ≃ (3/2)(𝔭+𝔮)c_s² — O(h₀²M/r), and it must FALL LIKE 1/a•.
            rels = map((400.0, 4000.0)) do aa
                bgr = _bgN(; 𝔭 = 𝔭, 𝔮 = 𝔮, aEMRI = aa)
                rr = 1.05 * aa
                Ωn = newtonian_rotation(bgr, rr).Ω
                (abs(Ωn - Omega(bgr, rr)) / Ωn, 1.5 * abs(𝔭 + 𝔮) * newtonian_csound(bgr, rr)^2)
            end
            @test rels[1][1] < 3 * rels[1][2] && rels[2][1] < 3 * rels[2][2]   # the predicted size
            @test 5 < rels[1][1] / rels[2][1] < 20                             # ... and the 1/a• scaling
        end
    end

    @testset "(b) the Plummer potential modes (60)" begin
        bg = _bgN(; aEMRI = 12.0, h0 = 0.05, b = 0.6)
        a = bg.aEMRI; εs = softening_length(bg)
        # independent dense quadrature of (1/2π)∫cos(mΨ)/s_N dΨ (no shared code with the package)
        function Ibrute(m, r; N = 400_000)
            acc = 0.0
            for k in 0:(N - 1)
                Ψ = 2π * k / N
                acc += cos(m * Ψ) / sqrt(εs^2 + r^2 + a^2 - 2 * a * r * cos(Ψ))
            end
            acc / N
        end
        for m in (0, 1, 5, 20), r in (9.0, 12.0, 14.0)
            pm = newtonian_potential_mode(bg, m, r)
            @test isapprox(pm.Φ, -Ibrute(m, r); rtol = 1e-9, atol = 1e-14)
            @test pm.Φ ≤ 0                                            # attractive
        end
        # the analytic r-derivative
        for m in (1, 7), r in (10.0, 13.5)
            h = 1e-6 * r
            fd = (newtonian_potential_mode(bg, m, r + h).Φ - newtonian_potential_mode(bg, m, r - h).Φ) / (2h)
            @test isapprox(newtonian_potential_mode(bg, m, r).Φ_r, fd; rtol = 1e-6)
        end
        # spectral convergence of the trapezoid rule (the default N is already converged)
        pm = newtonian_potential_mode(bg, 10, 12.3)
        @test isapprox(newtonian_potential_mode(bg, 10, 12.3; N = 4 * pm.N).Φ, pm.Φ; rtol = 1e-14)

        # the relativistic (B1) quadratic form tends to the flat one, and h_tt → 2m•/s0
        for aa in (100.0, 1e4, 1e6)
            K = SpiralWaves.singular_field_coeffs(_bgN(; aEMRI = aa))
            @test isapprox(K.Pr, 1.0; atol = 40 / aa)
            @test isapprox(K.Pt, 1.0; atol = 40 / aa)
            @test isapprox(K.Ktt, 2.0; atol = 40 / aa)
        end
        # ... hence the relativistic Φ = −ĥ_tt/2 converges to the Newtonian potential mode
        for aa in (1e3, 1e5)
            bgf = _bgN(; aEMRI = aa)
            r = 1.03 * aa
            rel = abs(-singular_field_mode(bgf, 6, r).htt / 2 - newtonian_potential_mode(bgf, 6, r).Φ) /
                  abs(newtonian_potential_mode(bgf, 6, r).Φ)
            @test rel < 200 / aa
        end
        # threaded == serial, bitwise
        rs = collect(range(10.0, 14.0; length = 37))
        @test newtonian_potential_modes(bg, 4, rs; threads = true) == newtonian_potential_modes(bg, 4, rs)
    end

    @testset "(c) the master system: (50) on both backgrounds, and the relativistic limit" begin
        m = 5
        for aa in (1e3, 1e5)
            bg = _bgN(; aEMRI = aa, 𝔭 = 1.5, 𝔮 = 1.0)
            r = 1.04 * aa
            pm = newtonian_potential_mode(bg, m, r)
            CN, SN = newtonian_master_matrices(bg, m, r; hmode = pm, landau = false)
            # same (50) algebra on the RELATIVISTIC background: differs only through Ω, κ², Γ₁, c_s
            CH, SH = master_matrices_newtonian(bg, m, r, pm.Φ, pm.Φ_r; landau = false)
            @test maximum(abs, (CN .- CH) ./ abs.(CN)) < 300 / aa
            @test maximum(abs, (SN .- SH) ./ abs.(SN)) < 300 / aa
            # and the FULL relativistic system (47) with the singular-field source
            hm = singular_field_mode(bg, m, r)
            CR, SR = master_matrices(bg, m, r; hmodes = hm, landau = false)
            @test maximum(abs, (CN .- CR) ./ abs.(CR)) < 300 / aa
            @test maximum(abs, (SN .- SR) ./ abs.(SR)) < 300 / aa
        end
        # δû_φ: the covariant component is r × (50a)
        bg = _bgN()
        y = (0.3 + 0.1im, -0.2 + 0.4im)
        r = 12.7
        pm = newtonian_potential_mode(bg, 4, r)
        rc = newtonian_corotation_radius(bg)
        @test isapprox(newtonian_delta_u_phi(bg, 4, r, y, pm.Φ; rc = rc),
                       r * delta_u_phi_newtonian(bg, 4, r, y, pm.Φ; rc = rc); rtol = 1e-12)
    end

    @testset "(d) the balance law (71) for a solved Newtonian mode" begin
        bg = _bgN(; aEMRI = 15.0, h0 = 0.05, 𝔭 = 1.5, 𝔮 = 1.0)
        m = 4
        worst = Float64[]
        for gk in (GKN, (nL = 16, nR = 16, nC = 400, n_anc = 1600))
            g = newtonian_grid(bg, m; gk...)
            sol = newtonian_solve_mode(bg, m; grid = g)
            @test sol.residual < 1e-8
            rm, lhs, rhs, rel = newtonian_balance_residual(bg, sol)
            far = findall(r -> abs(r - sol.params.rc) > 3 * bg.σε, rm)
            @test !isempty(far)
            push!(worst, maximum(rel[far]))
        end
        @info "Newtonian balance (71) away from corotation: worst rel = $(worst) (coarse, fine)"
        @test worst[2] < 5e-3
        @test worst[2] < worst[1]                                     # converges with resolution
    end

    @testset "(e) the relativistic pipeline tends to the Newtonian one as M/a• → 0" begin
        m = 5
        aa = 2000.0
        bg = _bgN(; aEMRI = aa, h0 = 0.05, 𝔭 = 0.0, 𝔮 = 0.0)
        gk = (nL = 12, nR = 12, nC = 400, n_anc = 1200)
        solN = newtonian_solve_mode(bg, m; grid = newtonian_grid(bg, m; gk...))
        DN = newtonian_torque_density(bg, solN)
        TN = last(cumulative_torque(bg, solN; density = DN))
        solR = solve_mode(bg, m; grid = build_grid(bg, m; gk...))
        DR = torque_density(bg, solR)
        TR = last(cumulative_torque(bg, solR; density = DR))
        rel = abs(TN - TR) / abs(TR)
        @info "per-mode torque at a• = $(aa)M, m = $m: Newtonian $(TN) vs relativistic $(TR) (rel $(rel))"
        @test rel < 0.02
        @test sign(TN) == sign(TR)
    end

    @testset "(f) the Newtonian torque is χ-INDEPENDENT (bitwise)" begin
        m = 3
        ref = nothing
        for χ in (0.0, 0.9, -0.9)
            bg = _bgN(; χ = χ, aEMRI = 12.0, 𝔭 = 1.5, 𝔮 = 1.0)
            sol = newtonian_solve_mode(bg, m; grid = newtonian_grid(bg, m; GKN...))
            D = newtonian_torque_density(bg, sol)
            T = last(cumulative_torque(bg, sol; density = D))
            ref === nothing ? (ref = T) : (@test T === ref)
            @test newtonian_grid(bg, m; GKN...).r_in == 6.0            # (C1) with r_ISCO(0)
        end
    end

    @testset "(g) exact scalings" begin
        bg = _bgN(; aEMRI = 12.0, 𝔭 = 1.5, 𝔮 = 0.0)
        m = 3
        g = newtonian_grid(bg, m; GKN...)
        T(q, sol) = last(cumulative_torque(bg, sol; density = newtonian_torque_density(bg, sol)))
        s1 = newtonian_solve_mode(bg, m; grid = g, mEMRI = 1.0)
        s2 = newtonian_solve_mode(bg, m; grid = g, mEMRI = 2.0)
        @test isapprox(T(2.0, s2), 4 * T(1.0, s1); rtol = 1e-12)      # ∝ m•²
        bg2 = _bgN(; aEMRI = 12.0, 𝔭 = 1.5, 𝔮 = 0.0, ε0 = 3.0)
        s3 = newtonian_solve_mode(bg2, m; grid = newtonian_grid(bg2, m; GKN...))
        T3 = last(cumulative_torque(bg2, s3; density = newtonian_torque_density(bg2, s3)))
        @test isapprox(T3, 3 * T(1.0, s1); rtol = 1e-12)              # ∝ ε₀
    end

    @testset "(h) the survey front end and its optimised reduction" begin
        bg = _bgN(; aEMRI = 12.0, h0 = 0.05, 𝔭 = 1.5, 𝔮 = 1.0)
        ms = 1:4
        for model in (:relativistic, :newtonian)
            s = torque_survey(bg, ms; threads = false, grid_kwargs = GKN, model = model)
            sp_ = torque_survey(bg, ms; threads = true, grid_kwargs = GKN, model = model,
                                keep_profiles = true)
            for (a, b) in zip(s.modes, sp_.modes)                     # optimised == full-array path
                @test a.T === b.T && a.T_L === b.T_L && a.T_c === b.T_c
                @test a.F_in === b.F_in && a.F_out === b.F_out
            end
            # and the full-array path IS `corotation_split` of the stand-alone solve
            m = 3
            if model === :newtonian
                sol = newtonian_solve_mode(bg, m; grid = newtonian_grid(bg, m; GKN...))
                F = newtonian_amf(bg, sol); D = newtonian_torque_density(bg, sol)
            else
                sol = solve_mode(bg, m; grid = build_grid(bg, m; GKN...))
                F = amf(bg, sol); D = torque_density(bg, sol)
            end
            sp = corotation_split(bg, sol; flux = F, density = D)
            mt = s.modes[findfirst(==(m), s.ms)]
            @test mt.T === sp.T && mt.T_L === sp.T_L && mt.T_c === sp.T_c
            @test mt.F_in === F[1] && mt.F_out === F[end]
        end
        @test_throws ArgumentError torque_survey(bg, ms; model = :quantum, grid_kwargs = GKN)
        # the Newtonian survey is χ-independent; the relativistic one is not
        sN = [torque_survey(_bgN(; χ = χ, aEMRI = 12.0, 𝔭 = 1.5, 𝔮 = 1.0), ms;
                            threads = false, grid_kwargs = GKN, model = :newtonian).T_total
              for χ in (0.0, 0.9)]
        @test sN[1] === sN[2]
    end
end
