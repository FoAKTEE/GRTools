# Module 4 (Scattering): MST scattering-phase pipeline (numerical branch).
# GROUNDING: (a) LIVE ScatteringPhaseMST oracle at its frozen tolerances (starobinsky_constant
# rtol 1e-11, the brefl/mst_coefficient/scattering_phase_numerical accessors to machine precision —
# same MST core); (b) the analytic Starobinsky–Chandrasekhar identities baked as literals
# (s=0 ⇒ C*=1; s=−1,a=0 ⇒ Re C*=l(l+1); s=−2,a=0 ⇒ Im C*=12Pω); (c) parallel==serial table.
# Legacy oracle RHS BAKED once into fixtures/leg_scattering.jl (Phase 7-B); tolerances unchanged.
include("fixtures/leg_scattering.jl")  # REF_STAR / REF_BREFL / REF_MSTC / REF_SPN

@testset "Scattering :: starobinsky_constant vs ScatteringPhaseMST oracle (baked, captured once)" begin
    for (s, l, m, q, ε, P) in [(0, 0, 0, 0.0, 0.1, +1), (-1, 2, 1, 0.5, 0.3, +1),
                               (-2, 3, -2, 0.9, 0.2, -1), (-2, 2, 2, 0.5, 0.25, +1),
                               (-1, 3, -2, 0.7, 0.2, +1)]
        @test isapprox(starobinsky_constant(s, l, m, q, ε; P=P),
                       REF_STAR[(s, l, m, q, ε, P)]; rtol=1e-11)
    end
end

@testset "Scattering :: Starobinsky–Chandrasekhar analytic identities (a=0)" begin
    @test starobinsky_constant(0, 2, 2, 0.0, 0.2) == complex(1.0)       # scalar normalization
    for (l, m) in [(1, 0), (2, 1), (3, -2)], ε in (0.1, 0.4)
        c = starobinsky_constant(-1, l, m, 0.0, ε)                      # EM: |C*|² = λ = l(l+1)
        @test real(c) ≈ l * (l + 1) rtol = 1e-12
    end
    for (l, m) in [(2, 0), (3, 1)], ε in (0.1, 0.3), P in (+1, -1)
        c = starobinsky_constant(-2, l, m, 0.0, ε; P=P)                 # grav: Im C* = 12 P ω
        @test imag(c) ≈ 12 * P * (ε / 2) rtol = 1e-12
    end
end

@testset "Scattering :: brefl_over_binc / mst_coefficient vs oracle" begin
    for (s, l, m, q, ε) in [(-2, 2, 2, 0.5, 0.25), (0, 2, 1, 0.6, 0.2), (-1, 3, -2, 0.7, 0.2)]
        @test isapprox(Perturbation._brefl_over_binc(s, l, m, q, ε; n=10, nUp=5),
                       REF_BREFL[(s, l, m, q, ε)]; rtol=1e-9)
        for nn in (-2, 0, 1, 3)
            @test isapprox(mst_coefficient(nn, s, l, m, q, ε),
                           REF_MSTC[(nn, s, l, m, q, ε)]; rtol=1e-9)
        end
    end
end

@testset "Scattering :: scattering_phase_numerical vs oracle + dispatcher" begin
    for (s, l, m, q, ε, P) in [(-2, 2, 0, 0.0, 0.2, -1), (0, 0, 0, 0.5, 0.15, +1),
                               (-1, 2, 1, 0.6, 0.2, +1)]
        v = scattering_phase_numerical(s, l, m, q, ε; n=40, nUp=5, P=P)
        @test isapprox(v, REF_SPN[(s, l, m, q, ε, P)]; rtol=1e-9)
        @test scattering_phase(s, l, m, q, ε; mode=:numerical, n=40, nUp=5, P=P) == v
    end
    @test_throws ArgumentError scattering_phase(0, 2, 2, 0.0, 0.1; mode=:bogus)
end

@testset "Scattering :: scattering_phase_analytic closed-form ε-series" begin
    # GROUNDING (a): the published phase{slm}Expansion series transcribed VERBATIM from
    # NumericalMST.nb cells (In[10]/In[8]/In[9]/In[5]) — the SAME closures the legacy test
    # packages/ScatteringPhaseMST/test bakes. The new closed-form Julia transcription must
    # reproduce them. ζ(3) and γE are published math constants.
    γE = Float64(Base.MathConstants.eulergamma)
    ζ3 = 1.2020569031595942854
    phase000(q, ε) = 0.5*(-1 + 2γE + 2log(2) + 2log(ε))*ε +
        (6im - 6im*q^2 + 6im*sqrt(1-q^2) + 11π*sqrt(1-q^2))*ε^2/(12*sqrt(1-q^2)) +
        (1/(36*sqrt(1-q^2)))*(-33 + 18im*π + 33q^2 - 18im*π*q^2 + 51*sqrt(1-q^2) -
            36γE*sqrt(1-q^2) + 18im*π*sqrt(1-q^2) + 11π^2*sqrt(1-q^2) + 3q^2*sqrt(1-q^2) -
            36*sqrt(1-q^2)*log(2) - 36*sqrt(1-q^2)*log(sqrt(1-q^2)*ε) - 12*sqrt(1-q^2)*ζ3)*ε^3
    phase010(q, ε) = 0.5*(-3 + 2γE + 2log(2) + 2log(ε))*ε + 19π*ε^2/60 +
        (19π^2 - 9q^2 - 60ζ3)*ε^3/180
    phase110(q, ε) = 0.25*(-5 + 4γE + 4log(2) + 4log(ε))*ε + (94π - 15im*q^2)*ε^2/240 +
        (105 + 188π^2 - 180q^2 - 480ζ3)*ε^3/1440
    phase220(q, ε) = (1/3)*(-5 + 3γE + 3log(2) + 3log(ε))*ε + 107π*ε^2/420 +
        (-1015 + 1926π^2 - 945q^2 - 7560ζ3)*ε^3/22680

    for q in (0.0, 0.5, 0.9), ε in (0.1, 0.2)
        @test isapprox(scattering_phase_analytic(0,  0, 0, q, ε), phase000(q, ε); rtol=1e-13)
        @test isapprox(scattering_phase_analytic(0,  1, 0, q, ε), phase010(q, ε); rtol=1e-13)
        @test isapprox(scattering_phase_analytic(-1, 1, 0, q, ε), phase110(q, ε); rtol=1e-13)
        @test isapprox(scattering_phase_analytic(-2, 2, 0, q, ε), phase220(q, ε); rtol=1e-13)
        # dispatcher mode=:analytic routes to it (returns δ, NOT e^{2iδ})
        @test scattering_phase(0, 0, 0, q, ε; mode=:analytic) == scattering_phase_analytic(0, 0, 0, q, ε)
    end

    # GROUNDING (b): independent agreement with the already-green numerical branch via
    # δ = log(scattering_phase_numerical(...))/2i at small ε, with the series' own O(ε⁴)
    # truncation rate (ratio ≈ 2⁴ = 16 between ε=0.1 and ε=0.2). The P-sign is load-bearing
    # for s=-2 (the C* term 12·P·i·ω).
    δnum(s, l, m, q, ε; P=1) = log(scattering_phase_numerical(s, l, m, q, ε; n=40, nUp=5, P=P)) / (2im)
    q = 0.5
    for (s, l, m, P, ref) in ((0, 0, 0, +1, phase000), (0, 1, 0, +1, phase010),
                              (-1, 1, 0, +1, phase110), (-2, 2, 0, -1, phase220))
        d1 = abs(δnum(s, l, m, q, 0.1; P=P) - scattering_phase_analytic(s, l, m, q, 0.1))
        d2 = abs(δnum(s, l, m, q, 0.2; P=P) - scattering_phase_analytic(s, l, m, q, 0.2))
        @test d1 < 40 * 0.1^4
        @test d2 < 40 * 0.2^4
        @test 12 < d2 / d1 < 22                          # O(ε⁴) convergence rate
    end
    # opposite P does NOT match the s=-2 analytic branch (P load-bearing sanity)
    @test abs(δnum(-2, 2, 0, q, 0.1; P=+1) - scattering_phase_analytic(-2, 2, 0, q, 0.1)) > 1e-2

    # unsupported mode is an honest [HOLE] boundary (no published series)
    @test_throws ArgumentError scattering_phase_analytic(-2, 3, 1, 0.5, 0.1)

    # firewall: the analytic branch pulled in NO symbolic dependency
    @test !("Symbolics" in Set(string(mod) for mod in values(Base.loaded_modules)))

    # BigFloat precision is inherited end-to-end (type-generic closed form)
    setprecision(BigFloat, 256) do
        vb = scattering_phase_analytic(-2, 2, 0, big"0.5", big"0.1")
        @test vb isa Complex{BigFloat}
        @test isapprox(ComplexF64(vb), scattering_phase_analytic(-2, 2, 0, 0.5, 0.1); rtol=1e-12)
    end
end

@testset "Scattering :: scattering_phase_table parallel == serial" begin
    εvals = (0.1, 0.15, 0.2, 0.25, 0.3)
    tab = scattering_phase_table(-2, 2, 0, 0.5, εvals; n=20, nUp=5, P=-1)
    ser = [scattering_phase_numerical(-2, 2, 0, 0.5, ε; n=20, nUp=5, P=-1) for ε in εvals]
    @test tab == ser
end

@testset "Scattering :: mst_convergence diagnostic" begin
    d = mst_convergence(-2, 2, 2, 0.5, 0.2; N=20)
    @test d.tail < 1e-6
    @test d.ratio < 1.0
    @test d.monotone
end
