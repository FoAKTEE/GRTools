# Grounding for Radiation TetradComponents — the six Hertz/CCK tetrad components (IRG s=−2 / ORG s=+2).
#
# Two grounding axes:
#   (1) SELF-CONSISTENCY (analytic, no dep): scalar HIRG_*/HORG_* == .v of the corresponding RTJet
#       (_jet) form at (r,θ) — the two implementations run the same operator chain, one collapsing
#       with value(), the other carrying the full RTJet. (legacy cov_components_api 'scalar==jet.v').
#   (2) ABSOLUTE VALUES vs a one-shot fixture captured from the legacy MetricReconstruction package
#       (HIRG_nn/.../HORG_mm) at (a=0.1, ω=0.3, l=2, m=2, r=6.0, θ=1.0). Legacy MR carries the
#       Teukolsky closure (heavy) + shares this package's NAME, so it is NEVER a live test-dep —
#       captured once via:  julia --project=packages/MetricReconstruction -e '... MR.HIRG_nn(...) ...'.
#       The GrTools radial_jet (MST) was checked to match the legacy (HeunC) mode to ~2e-11 rel, so
#       the linear operator chains reproduce the legacy H-values to ~1e-13 (observed); tol 1e-10.

using Test
using MetricReconstruction.Radiation

@testset "Radiation.TetradComponents (Hertz/CCK H-blocks)" begin
    a, ω, l, m, r, θ = 0.1, 0.3, 2, 2, 6.0, 1.0

    # captured legacy MetricReconstruction fixtures
    LEG = Dict(
        :HIRG_nn =>  1.3289661882711823 + 0.6122327105121582im,
        :HIRG_nm => -2.826204514253861  + 1.2889708801813613im,
        :HIRG_mm => -0.34899696332780117 - 3.730529954438228im,
        :HORG_ll =>  0.0017763536830153856 - 0.0011346769931373658im,
        :HORG_lm => -5.132645279836483e-5  - 0.000995522769560826im,
        :HORG_mm => -0.0003076187765635874 - 4.5482848483232326e-5im,
    )

    pairs = ((:HIRG_nn, HIRG_nn, HIRG_nn_jet), (:HIRG_nm, HIRG_nm, HIRG_nm_jet),
             (:HIRG_mm, HIRG_mm, HIRG_mm_jet), (:HORG_ll, HORG_ll, HORG_ll_jet),
             (:HORG_lm, HORG_lm, HORG_lm_jet), (:HORG_mm, HORG_mm, HORG_mm_jet))

    @testset "absolute values vs legacy fixture" begin
        for (name, f, _) in pairs
            @test f(a, ω, l, m, r, θ) ≈ LEG[name] rtol=1e-10
        end
    end

    @testset "self-consistency scalar == jet.v" begin
        for (name, f, fj) in pairs
            @test f(a, ω, l, m, r, θ) ≈ fj(a, ω, l, m, r, θ).v rtol=1e-12
        end
    end

    @testset "barred variants run + scalar==jet.v" begin
        for (_, f, fj) in pairs
            vb  = f(a, ω, l, m, r, θ; barred=true)
            vbj = fj(a, ω, l, m, r, θ; barred=true).v
            @test isfinite(real(vb)) && isfinite(imag(vb))
            @test vb ≈ vbj rtol=1e-12
        end
    end

    @testset "BigFloat inheritance (Complex{BigFloat})" begin
        setprecision(BigFloat, 256) do
            ab = big"0.1"; ωb = big"0.3"; rb = big"6.0"; θb = big"1.0"
            v = HIRG_mm(ab, ωb, l, m, rb, θb)
            @test v isa Complex{BigFloat}
            # BigFloat path agrees with the Float64 legacy fixture (radial/angular carry precision)
            @test Float64(real(v)) ≈ real(LEG[:HIRG_mm]) rtol=1e-8
            @test Float64(imag(v)) ≈ imag(LEG[:HIRG_mm]) rtol=1e-8
        end
    end
end
