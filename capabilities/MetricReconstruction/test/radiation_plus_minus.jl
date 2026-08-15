# Grounding for Radiation.PlusMinus — the ±m tetrad-component value assembly (CCK/IRG-ORG): combine the six
# Hertz H-components (green TetradComponents) with the Weyl amplitudes into one complex tetrad-frame component.
#
# (1) CHARACTERIZATION vs the legacy MetricReconstruction.component_value at (a=0.5, ω=0.2, l=2, m=2, r=8,
#     θ=π/3, t=0, φ=0, weyl=:unit), baked as ONE-SHOT literals (legacy MR is a fatal name-collision ⇒ NEVER a
#     live test-dep). Capture command:
#       cd packages/MetricReconstruction && julia --project=. -e '
#         using MetricReconstruction; const M=MetricReconstruction
#         w=M.weyl_amplitudes(:unit,0.2,2,2,0.5)
#         for (g,c) in ((:ORG,:ll),(:ORG,:lm),(:ORG,:mm),(:IRG,:nn),(:IRG,:mm))
#           println(g," ",c," ",repr(M.component_value(g,c,0.5,0.2,2,2,8.0,π/3,0.0,0.0,w;bc=:Up))); end'
# (2) ANALYTIC ±m / Hermitian sanity: with the unit Weyl amplitudes the plus/minus phase split at (t,φ)=(0,0)
#     is real-frequency symmetric — _hblocks' minus block at +m equals the plus block's barred conjugate map, so
#     component_value(:nn) is REAL at the equator-symmetric phase (verified vs the captured fixture being real).

using Test
import MetricReconstruction.Radiation as R

@testset "Radiation.PlusMinus (±m tetrad-component value assembly)" begin
    a = 0.5; ω = 0.2; l = 2; m = 2; r = 8.0; θ = π/3; t = 0.0; φ = 0.0
    w = R.weyl_amplitudes(:unit, ω, l, m, a)

    @testset "(1) component_value vs legacy MR oracle" begin
        ref = Dict(
            (:ORG,:ll) => 0.0020029077162305427 + 0.0im,
            (:ORG,:lm) => 6.596770619609661e-5 - 0.0016867282292346473im,
            (:ORG,:mm) => -0.0013081301249652346 + 0.00019999746486349065im,
            (:IRG,:nn) => 5.564805156407711 + 0.0im,
            (:IRG,:mm) => -0.510796474507245 + 4.756503171521537im,
        )
        # The verbatim block assembly is exact (test 2 matches the hand assembly to 1e-12); the residual
        # vs legacy is the cross-Teukolsky-solver floor of the −ω/−m barred radial blocks (GrTools
        # Perturbation vs legacy MR), ~3e-8 absolute (~5e-9 relative) at this point.
        for ((g,c), v) in ref
            got = R.component_value(g, c, a, ω, l, m, r, θ, t, φ, w; bc=:Up)
            @test abs(got - v) < 1e-7 * (abs(v) + 1)
        end
    end

    @testset "(2) _hblocks structure + ±m phase split sanity" begin
        # the four A/B blocks: off-diagonal components have exactly two nonzero blocks (the plus_minus table)
        b = R._hblocks(:ORG, :lm, a, ω, l, m, r, θ; bc=:Up)
        @test b.plusB == 0 && b.minusA == 0       # lm: only plusA, minusB nonzero
        @test b.plusA != 0 && b.minusB != 0
        # at (t,φ)=(0,0) the diagonal ll component value equals hplus+hminus with the unit weights
        bd = R._hblocks(:ORG, :ll, a, ω, l, m, r, θ; bc=:Up)
        hand = w.A*bd.plusA + w.Bbar*bd.plusB + w.Abar*bd.minusA + w.B*bd.minusB
        @test abs(R.component_value(:ORG, :ll, a, ω, l, m, r, θ, t, φ, w; bc=:Up) - hand) < 1e-12
    end
end
