# Grounding for rhs_amplitudes (RHSamplitudescalculation, Sourceincoordinates.nb).
#
# STRONGEST closed loop, no external oracle: rhs_amplitudes IS the generator of the shipped
# `*.dat` source-mode tables (it sums source_amplitudes over l=max(|m|,2)..10 at the ring
# radius). We regenerate every row of i1.0_p12.0_e0_r20.0.dat (circular a=0, p=12, x=1,
# ringrad=20, ringmass=1e-4; 10 modes m=1..10, n=k=0) and assert the recomputed Sl/Sn/Sm/Sm†
# and ω match the shipped columns. The file is read back with the (already-ported)
# read_source_table, so the literal file columns ARE the reference — no value is taken from new
# code output. Agreement is MST-truncation limited (~2e-6 relative, dominated by the small
# imaginary parts of the m=1 mode); tolerances frozen verbatim from the legacy
# rhs_amplitudes_test.jl (ω rtol=1e-6, legs rtol=1e-5, antisymmetry atol=1e-14).
#
# Cross-check: the Sm = -Sm† antisymmetry visible in the file (cols 9-10 == -cols 11-12) is
# reproduced to machine precision by the summed amplitudes.

using ModifiedGravity, Test

@testset "rhs_amplitudes" begin
    datafile = joinpath(@__DIR__, "data", "i1.0_p12.0_e0_r20.0.dat")
    modes = read_source_table(datafile; radius=20.0, ringmass=1e-4)
    @test length(modes) == 10

    @testset "regenerate row (m=$(md.m),n=$(md.n),k=$(md.k))" for md in modes
        r = rhs_amplitudes(md.m, md.n, md.k, 0.0, 12.0, 0.0, 1.0, 20.0, 1e-4)
        @test r.m == md.m && r.n == md.n && r.k == md.k
        @test r.radius == 20.0 && r.ringmass == 1e-4
        @test r.ω ≈ md.ω rtol = 1e-6           # file ω column carries 8 sig figs
        @test r.Sources.Sl    ≈ md.Sl    rtol = 1e-5
        @test r.Sources.Sn    ≈ md.Sn    rtol = 1e-5
        @test r.Sources.Sm    ≈ md.Sm    rtol = 1e-5
        @test r.Sources.Smdag ≈ md.Smdag rtol = 1e-5
        # Sm = -Sm† antisymmetry (cols 9-10 == -cols 11-12 in the shipped file).
        @test r.Sources.Sm ≈ -r.Sources.Smdag atol = 1e-14
        @test md.Sm ≈ -md.Smdag atol = 1e-14
    end


    @testset "argument validation" begin
        args = (1, 0, 0, 0.0, 12.0, 0.0, 1.0, 20.0, 1e-4)
        @test_throws ArgumentError rhs_amplitudes(args...; lmax=1)
        @test_throws ArgumentError rhs_amplitudes(5, args[2:end]...; lmax=4)
        @test_throws ArgumentError rhs_amplitudes(args...; Nmax=0)
        @test_throws ArgumentError rhs_amplitudes(args...; Nmax=-1)
    end
end
