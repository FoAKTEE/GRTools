# Grounding for Radiation.API — the user-facing single-mode surface (metric_tetrad/h_tetrad/h_coord) + the
# unified reconstruct/psi entry points.
#
# (1) CHARACTERIZATION vs legacy MetricReconstruction.{metric_tetrad,h_tetrad,h_coord} (HEAVY oracle → baked
#     one-shot literals, NOT a live dep). Captured with the legacy package standalone at (a=0.5,ω=0.2,l=2,m=2,
#     r=8,θ=π/3), default :unit weyl, IRG+ORG:
#       cd packages/MetricReconstruction && julia --project=. capture.jl
#     rel<1e-7 (l=2 cross-Teukolsky-solver floor).
# (2) WRAPPER-IDENTITY self-consistency (machine tol 1e-12): reconstruct == h_coord, psi == weyl_scalar (same
#     args), and metric_tetrad default-weyl == explicit weyl_amplitudes(:unit) — the AHA consolidation is a
#     no-semantics-change re-expression of the green primitives.

using Test
import MetricReconstruction.Radiation as R

@testset "Radiation.API (user-facing single-mode surface + reconstruct/psi)" begin
    a = 0.5; ω = 0.2; l = 2; m = 2; r = 8.0; θ = π/3

    # ---- (1) CHARACTERIZATION vs legacy MR api.jl single-mode outputs (baked literals) --------------------
    mt_IRG = Dict(:nn => 5.564805156407711 + 0.0im,
                  :nm => -2.70375834375853 - 2.0289142972271748im,
                  :nmbar => -2.70375834375853 + 2.0289142972271748im,
                  :mm => -0.510796474507245 + 4.756503171521537im,
                  :mbarmbar => -0.510796474507245 - 4.756503171521537im)
    mt_ORG = Dict(:ll => 0.0020029077162305427 + 0.0im,
                  :lm => 6.596770619609661e-5 - 0.0016867282292346473im,
                  :lmbar => 6.596770619609661e-5 + 0.0016867282292346473im,
                  :mm => -0.0013081301249652346 + 0.00019999746486349065im,
                  :mbarmbar => -0.0013081301249652346 - 0.00019999746486349065im)
    @testset "(1a) metric_tetrad fixture rel<1e-7" begin
        for (comp, ref) in mt_IRG
            @test abs(R.metric_tetrad(:IRG, comp, a, ω, l, m, r, θ) - ref) / abs(ref) < 1e-7
        end
        for (comp, ref) in mt_ORG
            @test abs(R.metric_tetrad(:ORG, comp, a, ω, l, m, r, θ) - ref) / abs(ref) < 1e-7
        end
    end

    ht_IRG = ComplexF64[0.0 0.0 0.0 0.0;
                        0.0 5.564805156407711+0.0im -2.70375834375853-2.0289142972271748im -2.70375834375853+2.0289142972271748im;
                        0.0 -2.70375834375853-2.0289142972271748im -0.510796474507245+4.756503171521537im 0.0;
                        0.0 -2.70375834375853+2.0289142972271748im 0.0 -0.510796474507245-4.756503171521537im]
    hc_ORG = ComplexF64[0.0003849338921832856 - 5.293955920339377e-23im 0.0004415476795752422 + 0.0im 0.0004545719144789502 + 0.0im -0.006854654066541326 + 0.0im;
                        0.0004415476795752422 + 0.0im 0.0005007269290576357 + 0.0im 0.000671343941352303 + 0.0im -0.008418728171356776 + 0.0im;
                        0.00045457191447895016 + 1.6263032587282567e-19im 0.000671343941352303 + 0.0im -0.08443855972441863 + 0.0im 0.0063721993299521405 + 0.0im;
                        -0.006854654066541326 + 0.0im -0.008418728171356776 + 0.0im 0.0063721993299521405 + 0.0im 0.06841577901463165 + 0.0im]
    @testset "(1b) h_tetrad / h_coord fixture rel<1e-7" begin
        Ht = R.h_tetrad(:IRG, a, ω, l, m, r, θ)
        @test maximum(abs, Ht .- ht_IRG) < 1e-7 * maximum(abs, ht_IRG)
        Hc = R.h_coord(:ORG, :BL, a, ω, l, m, r, θ)
        @test maximum(abs, Hc .- hc_ORG) < 1e-7 * maximum(abs, hc_ORG)
    end

    # ---- (2) WRAPPER-IDENTITY self-consistency (machine tol) ---------------------------------------------
    @testset "(2) reconstruct==h_coord, psi==weyl_scalar, default==:unit" begin
        w = R.weyl_amplitudes(:unit, ω, l, m, a)
        for g in (:IRG, :ORG)
            rc = R.reconstruct(g, a, ω, l, m, r, θ; coords=:BL, weyl=w)
            hc = R.h_coord(g, :BL, a, ω, l, m, r, θ; weyl=w)
            @test maximum(abs, rc .- hc) < 1e-12 * maximum(abs, hc)
            @test rc == hc                      # exact passthrough
        end
        for (which, g) in ((:psi0, :IRG), (:psi4, :ORG))
            @test R.psi(which, g, a, ω, l, m, r, θ) == R.weyl_scalar(which, g, a, ω, l, m, r, θ)
        end
        # default-weyl == explicit weyl_amplitudes(:unit)
        @test R.metric_tetrad(:IRG, :nn, a, ω, l, m, r, θ) ==
              R.metric_tetrad(:IRG, :nn, a, ω, l, m, r, θ; weyl=w)
        @test R.h_coord(:ORG, :BL, a, ω, l, m, r, θ) ==
              R.h_coord(:ORG, :BL, a, ω, l, m, r, θ; weyl=w)
    end
end
