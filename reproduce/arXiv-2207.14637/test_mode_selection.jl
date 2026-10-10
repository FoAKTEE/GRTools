# Run directly with --project=reproduce; no shared test registry is needed.
using Test
using GrTools.AccretionDisc.Modes: EccentricOrbit, lindblad_locations
include("_fr22_setup.jl")

@testset "FR22 positive-frequency mode mask and nested truncations" begin
    setup = FR22Setup()
    modes = fr22_modes(setup)
    @test length(modes) == length(Set(modes)) == 2514
    # Independent rectangular-domain membership oracle includes forbidden sectors.
    expected = Set((m,l) for m in 0:151 for l in -9:159
                   if 1 <= m <= 150 && l >= 1 && abs(l-m) <= 8)
    @test Set(modes) == expected
    @test issorted(modes)
    @test (1,0) ∉ modes
    @test all(((m,l),) -> m > 0 && l > 0, modes)
    @test isempty(fr22_modes(setup; m_max=0))
    @test fr22_modes(setup; lm_max=0) == [(m,m) for m in 1:150]
    @test_throws ArgumentError fr22_modes(setup; m_max=-1)
    @test_throws ArgumentError fr22_modes(setup; lm_max=-1)
    truncations = fr22_truncations(setup)
    @test first.(truncations.by_m) == collect(0:10:150)
    @test first.(truncations.by_offset) == collect(0:8)
    for family in (truncations.by_m, truncations.by_offset)
        @test last(family).second == modes
        @test all(issubset(Set(family[i].second), Set(family[i+1].second)) for i in 1:length(family)-1)
    end
    @test length.(last.(truncations.by_m)) == [0; [17m-36 for m in 10:10:150]]
    @test length.(last.(truncations.by_offset)) == [150*(2k+1)-k*(k+1)÷2 for k in 0:8]
    println("[PASS] modes=", length(modes), "; m counts=", length.(last.(truncations.by_m)),
            "; offset counts=", length.(last.(truncations.by_offset)))
end

@testset "FR22 pressureless resonance domain" begin
    for T in (Float64, BigFloat)
        setprecision(192) do
            setup = FR22Setup(T)
            orbit = EccentricOrbit(one(T), T(1)/4, one(T))
            radii = Tuple{Int,Int,Symbol,T}[]
            for (m,l) in fr22_modes(setup)
                ri, ro = lindblad_locations(orbit,m,l)
                @test (ri === nothing) == (m == 1)
                for (sign,radius) in ((-1,ri),(1,ro))
                    radius === nothing && continue
                    @test radius isa T
                    @test setup.values.r_min <= radius <= setup.values.r_max
                    # Unsquared equation independently checks sign and radius.
                    Om = orbit.n / (radius/orbit.a)^(T(3)/2)
                    @test isapprox(l*orbit.n-m*Om,sign*Om; rtol=2048eps(T),atol=2048eps(T))
                    push!(radii,(m,l,sign == -1 ? :inner : :outer,radius))
                end
            end
            @test length(radii) == 5019 # nine m=1 modes have only an outer branch
            lo, hi = radii[argmin(last.(radii))], radii[argmax(last.(radii))]
            @test lo[1:3] == (2,10,:inner)
            @test hi[1:3] == (9,1,:outer)
            @test isapprox(lo[4]^3,T(1)/100;rtol=16eps(T))
            @test isapprox(hi[4]^3,T(100);rtol=16eps(T))
            @test lindblad_locations(EccentricOrbit(T(2),T(1)/2,T(7)),2,10) ==
                  2 .* lindblad_locations(orbit,2,10)
            @test_throws ArgumentError lindblad_locations(orbit,1,0)
            @test_throws ArgumentError lindblad_locations(orbit,0,1)
            @test_throws ArgumentError lindblad_locations(orbit,1,-1)
            println("[PASS] ",T," resonances=",length(radii)," excluded=0 min=",lo," max=",hi)
        end
    end
end
