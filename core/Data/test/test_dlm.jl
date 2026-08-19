# Module 2: dlm — DelimitedFiles readers. GROUNDING:
#  * read_dat: the checksum-verified on-disk FEW `.dat` files — raw `readdlm` is the reference for the
#    shape; spot literal values transcribed from the file head/tail.
#  * read_csv: round-trip a small Float64 matrix through DelimitedFiles.
#  * _mma_sci: hand-known Mathematica sci-notation strings (Printf closed form is the oracle).
#  * read_source_table/_parse_source_line: a hand-authored 12-col tab record (the format is the oracle).

using Test
using Data
using DelimitedFiles: readdlm, writedlm
using TestKit: oracle_root, oracle_or_skip

const _FEWDATA = oracle_root("GRTOOLS_FEWDATA_DIR",
    "/data/haiyangw/claude/emri-disk/ref-code/FastEMRIWaveforms/.fewdata")

@testset "dlm" begin
    @testset "read_dat vs raw readdlm (checksum-verified files)" begin
        normf = joinpath(_FEWDATA, "AmplitudeVectorNorm.dat")
        fluxf = joinpath(_FEWDATA, "FluxNewMinusPNScaled_fixed_y_order.dat")
        if oracle_or_skip(isfile(normf) && isfile(fluxf),
                          "read_dat vs raw readdlm on the FEW .dat files"; env = "GRTOOLS_FEWDATA_DIR")
            a = Data.read_dat(normf)
            @test a == readdlm(normf)                  # faithful to raw readdlm
            @test size(a) == (1650, 3)
            @test a[1, :] == [1.368639426, 0.0, 1.494163808]      # file head
            @test a[end, :] == [3.818639426, 0.8, 0.08806320917]  # file tail
            fl = Data.read_dat(fluxf)
            @test size(fl) == (1650, 4)
            @test fl[1, :] == [1.3686394259, 0.0, 5.3489022, 5.3489022]
        end
    end

    @testset "read_csv round-trip" begin
        M = [1.0 2.5 -3.0; 4.0 5.5 6.0]
        f = joinpath(mktempdir(), "m.csv")
        writedlm(f, M, ',')
        @test Data.read_csv(f) == M
    end

    @testset "_mma_sci closed form" begin
        # Printf-formatted Mathematica sci notation; trailing-zero mantissa trimmed.
        @test Data._mma_sci(0.024056261; sigfigs = 8) == "2.4056261e-2"
        @test Data._mma_sci(6.0; sigfigs = 8) == "6e0"
        @test Data._mma_sci(1234.5; sigfigs = 8) == "1.2345e3"
        @test Data._mma_sci(0.0; sigfigs = 8) == "0e0"
        @test_throws ErrorException Data._mma_sci(Inf)
    end

    @testset "read_source_table 12-col contract" begin
        # One hand-authored mode record (m n k ω, then Re/Im of Sl,Sn,Sm,Sm†).
        line = "2\t0\t1\t0.15\t1.0\t2.0\t3.0\t4.0\t5.0\t6.0\t7.0\t8.0"
        modes = Data.read_source_table(IOBuffer(line * "\n\n"); radius = 10.0, ringmass = 1.0)
        @test length(modes) == 1                       # blank line skipped
        s = modes[1]
        @test (s.m, s.n, s.k) == (2, 0, 1)
        @test s.ω == 0.15
        @test s.radius == 10.0 && s.ringmass == 1.0
        @test s.Sl == 1.0 + 2.0im
        @test s.Sn == 3.0 + 4.0im
        @test s.Sm == 5.0 + 6.0im
        @test s.Smdag == 7.0 + 8.0im
        @test s.Snorm == 0.0
        @test_throws ErrorException Data._parse_source_line("1\t2\t3", 0.0, 0.0)  # wrong column count
    end
end
