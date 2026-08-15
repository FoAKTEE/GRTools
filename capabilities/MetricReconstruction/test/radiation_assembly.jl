# Grounding for Radiation.Assembly — the CCK mode-sum surface: sum coordinate h_μν over an (l,m) mode list
# (PARALLEL over modes) and evaluate the linearized-Einstein residual over an (r,θ) grid (PARALLEL over points).
#
# (1) CHARACTERIZATION vs legacy MetricReconstruction.h_coord_modesum (HEAVY oracle → baked one-shot literals,
#     NOT a live dep). Captured with the legacy package standalone at (a=0.5,ω=0.2,r=8,θ=π/3, gauge=:ORG,
#     coords=:BL, whichweyl=:unit):
#       cd packages/MetricReconstruction && julia --project=. -e 'h_coord_modesum(:ORG,:BL,0.5,0.2,modes,8.0,π/3)'
#     The TIGHT (rel<1e-7) fixture uses two distinct l=2 modes [(2,2),(2,1)] (both at the ~1e-9 cross-solver
#     floor demonstrated by the single-mode Projection oracle). The 3-mode list [(2,2),(3,2),(4,2)] adds l=3,4
#     whose GrTools-vs-legacy agreement is set by the cross-Teukolsky-solver accuracy floor (l=3 ~3e-6, l=4 ~1e-2),
#     so those higher-l modes are used ONLY for the assembly-logic gates (parallel==serial, vacuum) below, not the
#     tight value-characterization (each higher-l GrTools mode is itself vacuum to machine precision — test (2)).
# (2) VACUUM GATE (reuses green Einstein): the linearized-Einstein operator is LINEAR, so G⁽¹⁾[Σ_modes h] =
#     Σ_modes G⁽¹⁾[h]. Each mode is a vacuum CCK reconstruction ⇒ einstein_grid returns ‖G⁽¹⁾‖/‖h‖<1e-10 at
#     every grid point for every mode ⇒ the assembled mode sum is vacuum to machine precision.
# (3) PARALLEL==SERIAL bit-identical: the threaded h_coord_modesum equals a serial reference sum EXACTLY (the
#     per-mode slots are summed in fixed order k=1..n).
# (4) a→0 regression: h_coord_modesum at a→0 stays finite and equals its serial reference; einstein_grid vacuum.

using Test
import MetricReconstruction.Radiation as R
import Base.Threads

@testset "Radiation.Assembly (CCK mode-sum + Einstein grid)" begin
    a = 0.5; ω = 0.2; r = 8.0; θ = π/3
    modes = [(2,2),(3,2),(4,2)]                 # assembly-logic gates (parallel==serial, vacuum)
    modes2 = [(2,2),(2,1)]                       # tight value-characterization (both l=2, ~1e-9 cross-solver floor)

    # ---- (1) CHARACTERIZATION vs legacy MR h_coord_modesum (baked 2-mode literals, rel<1e-7) ---------------
    ref2 = ComplexF64[0.00037945038383895263 + 2.6469779601696886e-22im 0.00047645321722370875 + 0.0im 0.001677757637972141 + 0.0im -0.0027816388612175195 + 0.0im;
                      0.00047645321722370875 + 0.0im 0.0005935206521723547 + 0.0im 0.002194030030838203 + 0.0im -0.00394949547861435 + 0.0im;
                      0.001677757637972141 + 1.6263032587282567e-19im 0.002194030030838203 + 0.0im 0.03429145194069208 + 0.0im -0.0038679585035335095 + 0.0im;
                      -0.0027816388612175195 - 4.0657581468206416e-20im -0.00394949547861435 + 0.0im -0.0038679585035335095 + 0.0im -0.02368572001983328 + 0.0im]
    @testset "(1) legacy MR h_coord_modesum fixture rel<1e-7 (two l=2 modes)" begin
        hs = R.h_coord_modesum(:ORG, :BL, a, ω, modes2, r, θ; whichweyl=:unit)
        @test maximum(abs, hs .- ref2) < 1e-7 * maximum(abs, ref2)
    end

    # ---- (2) VACUUM GATE: einstein_grid ‖G⁽¹⁾‖/‖h‖ → 0 at every grid point, every mode --------------------
    @testset "(2) vacuum gate ‖G⁽¹⁾‖/‖h‖<1e-10 on grid (⇒ assembled sum vacuum by linearity)" begin
        pts = [(8.0, π/3), (10.0, π/4), (12.0, 2π/5)]
        for (l, m) in modes
            w = R.weyl_amplitudes(:unit, ω, l, m, a)
            Gs = R.einstein_grid(:ORG, a, ω, l, m, pts; weyl=w)
            for (k, (rr, θθ)) in enumerate(pts)
                H = R.h_coordinate(:ORG, :BL, a, ω, l, m, rr, θθ, 0.0, 0.0, w)
                @test maximum(abs, Gs[k]) < 1e-10 * maximum(abs, H)
            end
        end
    end

    # ---- (3) PARALLEL==SERIAL bit-identical --------------------------------------------------------------
    @testset "(3) parallel==serial bit-identical (nthreads=$(Threads.nthreads()))" begin
        hs = R.h_coord_modesum(:ORG, :BL, a, ω, modes, r, θ; whichweyl=:unit)
        acc = zeros(ComplexF64, 4, 4)
        for (l, m) in modes
            w = R.weyl_amplitudes(:unit, ω, l, m, a)
            acc .+= R.h_coordinate(:ORG, :BL, a, ω, l, m, r, θ, 0.0, 0.0, w)
        end
        @test hs == acc
    end

    # ---- (4) a→0 regression -----------------------------------------------------------------------------
    @testset "(4) a→0 regression: finite + serial-identical + vacuum" begin
        aa = 1e-6
        hs = R.h_coord_modesum(:ORG, :BL, aa, ω, modes, r, θ; whichweyl=:unit)
        @test all(isfinite, hs)
        acc = zeros(ComplexF64, 4, 4)
        for (l, m) in modes
            w = R.weyl_amplitudes(:unit, ω, l, m, aa)
            acc .+= R.h_coordinate(:ORG, :BL, aa, ω, l, m, r, θ, 0.0, 0.0, w)
        end
        @test hs == acc
        w22 = R.weyl_amplitudes(:unit, ω, 2, 2, aa)
        G = R.einstein_grid(:ORG, aa, ω, 2, 2, [(8.0, π/3)]; weyl=w22)[1]
        H = R.h_coordinate(:ORG, :BL, aa, ω, 2, 2, 8.0, π/3, 0.0, 0.0, w22)
        @test maximum(abs, G) < 1e-10 * maximum(abs, H)
    end
end
