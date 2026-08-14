# Characterization of `kummer_hyp` (HypGeo-delegated confluent ₁F₁). Expected values are CLOSED
# FORMS (M(a;a;z)=eᶻ, M(1;2;z)=(eᶻ−1)/z) + a cross-check vs NumericsCore's self-contained series
# `kummer_1F1` (the distinct decision) — never == new-code.

@testset "kummer" begin
    # ---- closed forms: M(a;a;z) = eᶻ ;  M(1;2;z) = (eᶻ − 1)/z -----------------------------
    for z in (0.7, -1.3, 1.0 + 0.5im, -2.0 + 1.0im)
        @test kummer_hyp(2.3, 2.3, z) ≈ exp(z) atol=1e-12
        @test kummer_hyp(1.0, 2.0, z) ≈ (exp(z) - 1) / z atol=1e-12
    end

    # ---- cross-check vs NumericsCore.kummer_1F1 (vendored series — distinct primitive) ----
    for (a, b, z) in ((0.5, 1.5, 0.8), (1.2, 2.7, -0.6), (0.5, 1.5, 0.4 + 0.3im))
        @test kummer_hyp(a, b, z) ≈ NumericsCore.kummer_1F1(a, b, z) atol=1e-11
    end

    # ---- type-genericity: Complex{BigFloat} flows through (closed form at full precision) ----
    setprecision(BigFloat, 256) do
        zb = Complex{BigFloat}(big"0.7", big"0.3")
        v = kummer_hyp(Complex{BigFloat}(big"2.3"), Complex{BigFloat}(big"2.3"), zb)
        @test v isa Complex{BigFloat}
        @test abs(v - exp(zb)) < 1e-60                     # M(a;a;z)=eᶻ, BigFloat precision
        v2 = kummer_hyp(Complex{BigFloat}(big"1.0"), Complex{BigFloat}(big"2.0"), zb)
        @test abs(v2 - (exp(zb) - 1) / zb) < 1e-60
    end
end
