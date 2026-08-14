# Characterization of the self-contained ₂F₁ stratum (consolidates Teukolsky `_f21_maclaurin`/
# `_f21_zinf` + Punctures `_f21_term`). Expected values are CLOSED FORMS (Gauss special cases,
# DLMF connection), the `HypergeometricFunctions` oracle, or the LIGHT live Teukolsky internals.

@testset "hyp2f1" begin
    # ---- f21_maclaurin: closed forms 2F1(1,1;2;z)=-log(1-z)/z, 2F1(a,b;b;z)=(1-z)^{-a} ----
    let z = 0.3 + 0.1im
        @test f21_maclaurin(1.0, 1.0, 2.0, z) ≈ -log(1 - z) / z atol=1e-14
        @test f21_maclaurin(2.0, 3.0, 3.0, z) ≈ (1 - z)^(-2) atol=1e-14   # b=c degeneracy
        @test f21_maclaurin(0.7, 1.3, 2.2, z) ≈ HypergeometricFunctions._₂F₁(0.7, 1.3, 2.2, z) atol=1e-13
    end
    # BigFloat Maclaurin vs -log(1-z)/z (cov_internals mac_err < 1e-70 at 256-bit)
    setprecision(BigFloat, 256) do
        zm = Complex{BigFloat}(big"0.3", big"0.1")
        @test abs(f21_maclaurin(Complex{BigFloat}(1), Complex{BigFloat}(1), Complex{BigFloat}(2), zm)
                  - (-log(1 - zm) / zm)) < 1e-70
        # BAKED Teukolsky oracle: bit-identical (same incremental-ratio series), captured once
        #   Teukolsky._f21_maclaurin(big"0.7", big"1.3", big"2.2", zm)
        @test f21_maclaurin(Complex{BigFloat}(big"0.7"), Complex{BigFloat}(big"1.3"),
                            Complex{BigFloat}(big"2.2"), zm) ==
              Complex{BigFloat}(parse(BigFloat, "1.147837651962715510344857385556806640224883509644200625775181054482067342870217"), parse(BigFloat, "0.06293996589476951626874441144032223525013659058504810844562028683795678253078431"))
    end

    # ---- Gauss z=1 summation 2F1(a,b;c;1)=Γ(c)Γ(c-a-b)/(Γ(c-a)Γ(c-b)) via gamma_big ----
    let a = 0.5, b = 0.5, c = 2.0
        gauss = gamma_big(complex(c)) * gamma_big(complex(c - a - b)) /
                (gamma_big(complex(c - a)) * gamma_big(complex(c - b)))
        @test gauss ≈ HypergeometricFunctions._₂F₁(a, b, c, 1.0) atol=1e-13
    end

    # ---- f21_zinf: z→∞ connection matches HypGeo at large |z| (cov zinf_err < 1e-13) -----
    setprecision(BigFloat, 256) do
        a, b, c = 0.7, 1.3, 2.2                              # a−b = −0.6 (non-integer)
        zL = -30.0 + 5.0im
        ref = HypergeometricFunctions._₂F₁(a, b, c, zL)
        zb = Complex{BigFloat}(big"-30.0", big"5.0")
        val = f21_zinf(Complex{BigFloat}(a), Complex{BigFloat}(b), Complex{BigFloat}(c), zb)
        @test abs(ComplexF64(val) - ref) / abs(ref) < 1e-13
        # BAKED Teukolsky oracle: bit-identical (same DLMF 15.8.2 + gamma_big), captured once
        #   Teukolsky._f21_zinf(Complex{BigFloat}(0.7), Complex{BigFloat}(1.3), Complex{BigFloat}(2.2), zb)
        @test val == Complex{BigFloat}(parse(BigFloat, "0.1590808150237040589763233267185029994234123800001219874185703102187751868318813"), parse(BigFloat, "0.01592703086206977066133857640284767828460583180348097940226054467961637610356034"))
    end

    # ---- safe_2f1: byte-identical pass-through for Float64/ComplexF64 (cov f2f1_log/pow) --
    let z = 0.3 + 0.1im
        @test safe_2f1(1.0, 1.0, 2.0, z) === HypergeometricFunctions._₂F₁(1.0, 1.0, 2.0, z)
        @test abs(safe_2f1(1.0, 1.0, 2.0, z) - (-log(1 - z) / z)) < 1e-14
        @test abs(safe_2f1(2.0, 3.0, 3.0, z) - (1 - z)^(-2)) < 1e-14
        @test safe_2f1(1.0 + 0im, 0.5 + 0im, 2.0 + 0im, -3.5 + 0im) ===
              HypergeometricFunctions._₂F₁(1.0 + 0im, 0.5 + 0im, 2.0 + 0im, -3.5 + 0im)
    end
    # BigFloat: correct value (closed form) + matches Teukolsky._safe_2f1 (same try-path / Pfaff)
    setprecision(BigFloat, 256) do
        zb = Complex{BigFloat}(big"0.3", big"0.1")
        @test abs(safe_2f1(Complex{BigFloat}(2), Complex{BigFloat}(3), Complex{BigFloat}(3), zb)
                  - (1 - zb)^(-2)) < 1e-68
        # BAKED Teukolsky oracle (same try-path / Pfaff), captured once
        #   Teukolsky._safe_2f1(big"0.7", big"1.3", big"2.2", zb)
        @test safe_2f1(Complex{BigFloat}(big"0.7"), Complex{BigFloat}(big"1.3"),
                       Complex{BigFloat}(big"2.2"), zb) ==
              Complex{BigFloat}(parse(BigFloat, "1.147837651962715510344857385556806640224883509644200625775181054482067342870286"), parse(BigFloat, "0.06293996589476951626874441144032223525013659058504810844562028683795678253081885"))
    end

    # ---- f21_terminating: finite == series, + explicit l=1 closed form 1 − 2z/c ----------
    let c = 2.5, z = 0.4
        @test f21_terminating(1, c, z) ≈ 1 - 2*z/c atol=1e-14          # 2F1(-1,2;c;z)
        for l in 0:5
            @test f21_terminating(l, c, z) ≈ f21_maclaurin(-l, 1 + l, c, z) atol=1e-12
        end
    end
    setprecision(BigFloat, 256) do
        cb = Complex{BigFloat}(big"2.5"); zb = Complex{BigFloat}(big"0.4")
        @test abs(f21_terminating(3, cb, zb)
                  - f21_maclaurin(Complex{BigFloat}(-3), Complex{BigFloat}(4), cb, zb)) < 1e-66
    end
end
