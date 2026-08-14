# Characterization of HybridC module 2 (wasabic_codegen): the machine-generated 9-arg SF-PN
# closed forms (δx + ∂δx/∂χ1), exercised through the public `wc_δx` / `wc_dδxdχ` wrappers.
#
# GROUNDING (WaSABI NOT a live dep — baked WL-oracle literals captured from the legacy oracle16a
# run, the "DX" rows of `wolframscript -file oracle16a.wl` over the upstream paclet; reproduced
# verbatim from packages/WaSABI/test/inspiral_WaSABI-C_deps_test.jl). The codegen file is a
# deterministic WL InputForm translation artifact (codegen_wasabic.wl); the emitted values — NOT
# the new code — are re-validated against these WL literals.

using Test
using Hybrid.HybridC

# Baked fixture copy of the upstream WaSABI data (provenance: test/fixtures/wasabi/README.md).
const _CG_BASE = joinpath(@__DIR__, "fixtures", "wasabi")

@testset "HybridC.wasabic_codegen" begin
    mdl = wasabic_model(joinpath(_CG_BASE, "InspiralModels", "sf_data"),
                        joinpath(_CG_BASE, "AmplitudeModels", "sf_amp_data"))

    @testset "δx, dδx/dχ vs WL oracle (5 (ω,χ1) points)" begin
        # WL: δx[ω,0,0.1,1,χ1,0.3,0,0,0] and dδxdχ[...] (oracle16a.wl, "DX" rows).
        cases = [
            (0.02, 0.5, -0.0011703049117363458, 0.0002806120176106345),
            (0.03, 0.9, -0.0021729525904774, 0.000675073226653616),
            (0.015, -0.5, -0.0007940112450649676, 0.00019830999046972536),
            (0.025, 0.0, -0.0021738625892952983, 0.0009301346198665411),
            (0.04, 0.7, -0.0034262225831149205, 0.0017410546136173775),
        ]
        for (ω, χ1, dxwl, ddxwl) in cases
            dx = wc_δx(mdl, ω, 0.0, 0.1, 1, χ1, 0.3, 0.0, 0.0, 0.0)
            ddx = wc_dδxdχ(mdl, ω, 0.0, 0.1, 1, χ1, 0.3, 0.0, 0.0, 0.0)
            @test isapprox(dx, dxwl; rtol = 1e-9)
            @test isapprox(ddx, ddxwl; rtol = 1e-8)
        end
    end
end
