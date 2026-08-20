# test_boundary_response.jl — EGD::fig10-per-mode-reversal / EGD::sw-mode-solver: the INNER-BOUNDARY
# RESPONSE instruments added for the low-m torque-deficit investigation of arXiv:2608.17003.
#
# Three additive capabilities, every default unchanged:
#   (a) `wkb_amplitudes` / `mode_amplitudes` — split a solved mode into the two local WKB branches of
#       the master matrix C and report the complex reflection coefficient R = A_toward/A_away;
#   (b) `bc = :wkb1` — a FIRST-ORDER WKB (impedance) closure that replaces the zeroth-order row
#       "y ∥ Y_away" by "the incoming amplitude takes its locally slaved value", an INHOMOGENEOUS row;
#   (c) `build_grid(...; r_in = , r_out = )` — explicit domain-edge overrides, for displacing the
#       imposed boundary to r_ISCO(1+δ).
#
# Gates, in order:
#   1. `wkb_amplitudes` reconstructs y exactly (partition of unity of the two branches);
#   2. it agrees with `outgoing_branch` on which branch is "away" and returns the matching k;
#   3. at the inner boundary NODE of a `bc = :outgoing` solve the incoming amplitude is zero to
#      roundoff — the zeroth-order row is exactly "no incoming wave THERE";
#   4. the pure algebra `_impedance_row` reduces to the zeroth-order row when the eigenvector
#      derivatives and the edge source vanish (so `:wkb1` is a strict generalisation);
#   5. `bc = :outgoing` and `bc = :paper` are BITWISE unchanged by the refactor (compared against a
#      hand-assembled `solve_midpoint` with the historical rows);
#   6. `bc = :wkb1` solves, leaves the discrete residual at roundoff, satisfies its own inhomogeneous
#      boundary rows, and moves the mode by a finite (recorded, not asserted-small) amount;
#   7. `build_grid` edge overrides: `nothing` reproduces the adaptive grid BITWISE; explicit edges are
#      honoured; an override that swallows the corotation block is rejected.
using Test
using AccretionDisc
using AccretionDisc.SpiralWaves
using AccretionDisc.SpiralWaves: _impedance_row, _away_is_first

_bgB(::Type{T} = Float64; χ = 0, aEMRI = 15, h0 = 0.05, 𝔭 = 0, 𝔮 = 0) where {T} =
    DiscBackground{T}(; χ = T(χ), aEMRI = T(aEMRI), h0 = T(h0), 𝔭 = T(𝔭), 𝔮 = T(𝔮))
_gridB(bg, m; kw...) = build_grid(bg, m; nL = 30, nR = 30, nC = 600, n_anc = 2000, kw...)

@testset "SpiralWaves boundary response: WKB amplitudes, first-order closure, edge overrides" begin
    bg = _bgB(); m = 4
    g = _gridB(bg, m)

    # ---- (1)-(2) the two-branch split -----------------------------------------------------------
    @testset "wkb_amplitudes splits and reconstructs" begin
        for r in (g.r_in, g.r[end], (g.r_in + g.r_c) / 2)
            C, _ = master_matrices(bg, m, r; hmodes = singular_field_mode(bg, m, r))
            ωd = omega_d(bg, m, r)
            outward = r > g.r_c
            y = (0.3 + 0.7im, -1.1 + 0.2im)
            a = wkb_amplitudes(C, ωd, y; outward = outward)
            rec = (a.A_away .* a.Y_away) .+ (a.A_toward .* a.Y_toward)
            @test rec[1] ≈ y[1] rtol = 1e-12
            @test rec[2] ≈ y[2] rtol = 1e-12
            # the "away" branch is exactly the one outgoing_branch selects
            k, Y, _ = outgoing_branch(C, ωd; outward = outward)
            @test a.k_away == k
            @test a.Y_away == Y
            @test a.R ≈ a.A_toward / a.A_away
            # C·Y = kY for both branches
            for (kk, YY) in ((a.k_away, a.Y_away), (a.k_toward, a.Y_toward))
                @test C[1, 1] * YY[1] + C[1, 2] * YY[2] ≈ kk * YY[1] atol = 1e-9 * abs(kk)
                @test C[2, 1] * YY[1] + C[2, 2] * YY[2] ≈ kk * YY[2] atol = 1e-9 * abs(kk)
            end
        end
    end

    # ---- (3) the zeroth-order row is "no incoming wave at the edge node" -------------------------
    @testset "R = 0 at the boundary node of a bc = :outgoing solve" begin
        sol = solve_mode(bg, m; grid = g)
        for i in (1, length(sol.r))
            a = mode_amplitudes(bg, sol, i)
            @test abs(a.A_toward) / abs(a.A_away) < 1e-9
            @test abs(a.R) < 1e-9
        end
        # one node into the domain the split is a genuine mixture (the source and the mode coupling
        # feed the incoming branch): finite, non-zero, and much smaller than unity
        ain = mode_amplitudes(bg, sol, 6)
        @test isfinite(abs(ain.R))
        @test abs(ain.R) > 0
    end

    # ---- (4) the first-order row degenerates to the zeroth-order one -----------------------------
    @testset "_impedance_row reduces to the zeroth-order row" begin
        Ya = (0.6 + 0.1im, -0.2 + 0.4im)
        Yt = (0.3 - 0.5im, 0.9 + 0.0im)
        kt = -1.3 + 2.5im
        z2 = (0.0im, 0.0im)
        row, rhs = _impedance_row(Ya, Yt, kt, z2, z2, z2)
        @test row[1] ≈ -Ya[2]
        @test row[2] ≈ Ya[1]
        @test rhs == 0
        # a non-zero source alone makes the row inhomogeneous but leaves its direction untouched
        row2, rhs2 = _impedance_row(Ya, Yt, kt, z2, z2, (0.4 + 0.0im, -0.7im))
        @test row2[1] ≈ row[1]
        @test row2[2] ≈ row[2]
        @test rhs2 != 0
        # a non-zero dY_away tilts the row toward the incoming direction
        row3, _ = _impedance_row(Ya, Yt, kt, (0.05 + 0.0im, 0.02im), z2, z2)
        @test !(row3[1] ≈ row[1] && row3[2] ≈ row[2])
    end

    # ---- (5) the refactor left the landed branches bitwise unchanged -----------------------------
    @testset "bc = :outgoing / :paper unchanged (bitwise)" begin
        for bcv in (:outgoing, :paper)
            sol = solve_mode(bg, m; grid = g, bc = bcv)
            r = g.r; N = length(r)
            rmid = [(r[i] + r[i+1]) / 2 for i in 1:(N - 1)]
            hmid = source_modes(bg, m, rmid)
            CS = [master_matrices(bg, m, rmid[i]; hmodes = hmid[i], rc = g.r_c) for i in 1:(N - 1)]
            Cin, _ = master_matrices(bg, m, r[1]; hmodes = singular_field_mode(bg, m, r[1]), rc = g.r_c)
            Cout, _ = master_matrices(bg, m, r[N]; hmodes = singular_field_mode(bg, m, r[N]), rc = g.r_c)
            _, _, Win = outgoing_branch(Cin, omega_d(bg, m, r[1]); outward = false, bc = bcv)
            _, _, Wout = outgoing_branch(Cout, omega_d(bg, m, r[N]); outward = true, bc = bcv)
            yref = solve_midpoint(r, first.(CS), last.(CS), Win, 0.0im, Wout, 0.0im)
            @test sol.y == yref
        end
    end

    # ---- (6) the first-order closure end to end --------------------------------------------------
    @testset "bc = :wkb1 solves and closes its own rows" begin
        s0 = solve_mode(bg, m; grid = g, bc = :outgoing)
        s1 = solve_mode(bg, m; grid = g, bc = :wkb1)
        @test all(isfinite, s1.y)
        @test s1.residual < 1e-8
        @test s1.bc_residual < 1e-8 * maximum(abs, s1.y)
        @test s1.params.bc === :wkb1
        # it is a DIFFERENT closure: the boundary node no longer has exactly zero incoming amplitude
        a1 = mode_amplitudes(bg, s1, 1)
        @test abs(a1.R) > 1e-9
        # ... but the interior solution stays in the same ballpark (this is a boundary perturbation,
        # not a different problem): the mode torques agree to within a factor of two
        t0 = abs(lindblad_torque(bg, s0)); t1 = abs(lindblad_torque(bg, s1))
        @test 0.5 < t1 / t0 < 2.0
        @test_throws ArgumentError solve_mode(bg, m; grid = g, bc = :nonsense)
    end

    # ---- (7) explicit domain-edge overrides -------------------------------------------------------
    @testset "build_grid edge overrides" begin
        gd = _gridB(bg, m)
        @test gd.r == g.r                                     # nothing ⇒ unchanged, bitwise
        δ = 1e-3
        rin2 = r_isco(bg) * (1 + δ)
        g2 = _gridB(bg, m; r_in = rin2)
        @test g2.r_in == rin2
        @test g2.r[1] == rin2
        @test g2.r_out == g.r_out
        g3 = _gridB(bg, m; r_out = 4 * bg.aEMRI)
        @test g3.r_out == 4 * bg.aEMRI
        @test g3.r[end] == 4 * bg.aEMRI
        @test_throws ErrorException _gridB(bg, m; r_in = g.r_c)   # swallows the corotation block
    end
end
