# Grounding for Radiation.Projection — the tetrad→coordinate projection (CCK): build the 4×4 tetrad-frame
# H-matrix from component_value, form the dual co-tetrad legs (e^l=−n_down,…) from Spacetime kinnersley/metric,
# and contract to the coordinate h_μν.
#
# (1) CHARACTERIZATION vs the legacy MetricReconstruction.h_coordinate at (gauge=:ORG, coords=:BL, a=0.5, ω=0.2,
#     l=2, m=2, r=8, θ=π/3, t=0, φ=0, weyl=:unit), baked as ONE-SHOT literals (capture cmd below). Tolerance is
#     the cross-Teukolsky-solver floor inherited through the H-blocks.
#       cd packages/MetricReconstruction && julia --project=. -e '
#         using MetricReconstruction; const M=MetricReconstruction; w=M.weyl_amplitudes(:unit,0.2,2,2,0.5)
#         H=M.h_coordinate(:ORG,:BL,0.5,0.2,2,2,8.0,π/3,0.0,0.0,w;bc=:Up); for i=1:4,j=1:4; println(repr(H[i,j])); end'
# (2) ANALYTIC Kinnersley duality identity (corpus-free, live via Spacetime): the covariant dual legs D[μ,a]=
#     (e^a)_μ contracted with the contravariant Kinnersley legs E[μ,b]={l,n,m,m̄} reproduce δ — i.e.
#     Σ_μ D[μ,a] E[μ,b] = δ_{ab} (the l·n=−1, m·m̄=+1 normalization), the dual-basis property the projection rests on.

using Test
import MetricReconstruction.Radiation as R
import Spacetime

@testset "Radiation.Projection (tetrad→coordinate CCK projection)" begin
    a = 0.5; ω = 0.2; l = 2; m = 2; r = 8.0; θ = π/3; t = 0.0; φ = 0.0
    w = R.weyl_amplitudes(:unit, ω, l, m, a)

    @testset "(1) h_coordinate(:ORG,:BL) vs legacy MR oracle" begin
        Href = [
            0.0003849338921832856      0.0004415476795752422  0.0004545719144789502   -0.006854654066541326;
            0.0004415476795752422      0.0005007269290576357  0.000671343941352303    -0.008418728171356776;
            0.00045457191447895016     0.000671343941352303  -0.08443855972441863      0.0063721993299521405;
           -0.006854654066541326      -0.008418728171356776   0.0063721993299521405    0.06841577901463165 ]
        H = R.h_coordinate(:ORG, :BL, a, ω, l, m, r, θ, t, φ, w; bc=:Up)
        for i in 1:4, j in 1:4
            @test abs(real(H[i,j]) - Href[i,j]) < 1e-7 * (abs(Href[i,j]) + 1)
            @test abs(imag(H[i,j])) < 1e-7      # the static (t=φ=0) ORG h is real to solver floor
        end
    end

    @testset "(2) Kinnersley duality: D' E = δ (live Spacetime, corpus-free)" begin
        D = R.covariant_legs(:BL, a, r, θ)
        lup, nup, mup, mbarup = Spacetime.kinnersley_tetrad(a, r, θ; coords=:BL)
        E = hcat(lup, nup, mup, mbarup)         # contravariant legs as columns (l,n,m,m̄)
        G = transpose(D) * E                    # G[a,b] = Σ_μ D[μ,a] E[μ,b]
        δ = [i == j ? one(eltype(G)) : zero(eltype(G)) for i in 1:4, j in 1:4]
        @test maximum(abs, G - δ) < 1e-12
    end
end
