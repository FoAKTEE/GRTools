# Grounding for SourceTables: parse/format of the 12-column .dat source-mode table.
#
# ORACLE: the shipped *.dat files in test/data/ ARE the oracle — they are the literal output of
# the upstream `mywrite` (Sourceincoordinates.nb, Polcar & Witzany 2025). Expected field values
# are transcribed directly from those files; the column→field mapping and Snorm=0 are read from
# `Sourcefiletotable` in the same notebook. The read path is reused from Data (DRY); this module
# adds the WRITE side, grounded by the read∘write round-trip identity + byte-exact ω column.
# (Transcribed verbatim from packages/ModifiedTeukolsky/test/source_table_io_test.jl.)

using ModifiedGravity, Test

const DATADIR = joinpath(@__DIR__, "data")
const F_I1 = joinpath(DATADIR, "i1.0_p12.0_e0_r20.0.dat")
const F_I05_R20 = joinpath(DATADIR, "i0.5_p12.0_e0_r20.0.dat")
const F_I05_R400 = joinpath(DATADIR, "i0.5_p12.0_e0_r400.0.dat")

rawtokens(line) = split(line, '\t')

@testset "SourceTables" begin

    @testset "parse: literal column values (i1.0 row 1)" begin
        modes = read_source_table(F_I1; radius=20.0, ringmass=1e-4)
        @test length(modes) == 10
        m1 = modes[1]
        @test m1.m == 1
        @test m1.n == 0
        @test m1.k == 0
        @test m1.ω == parse(Float64, "2.4056261e-2")
        @test m1.Sl == complex(parse(Float64, "4.41089376980272e-4"),
                               parse(Float64, "-3.532461734331e-6"))
        @test m1.Sn == complex(parse(Float64, "-1.924506402659909e-4"),
                               parse(Float64, "-1.78535749256213e-5"))
        @test m1.Sm == complex(parse(Float64, "-1.8600511940717e-5"),
                               parse(Float64, "-8.106123685776e-6"))
        @test m1.Smdag == complex(parse(Float64, "1.8600511940717e-5"),
                                  parse(Float64, "8.106123685776e-6"))
        @test m1.radius == 20.0
        @test m1.ringmass == 1e-4
    end

    @testset "parse: last row (i1.0 row 10)" begin
        modes = read_source_table(F_I1; radius=20.0, ringmass=1e-4)
        m10 = modes[10]
        @test (m10.m, m10.n, m10.k) == (10, 0, 0)
        @test m10.ω == parse(Float64, "2.4056261e-1")
        @test m10.Sl == complex(parse(Float64, "3.86316616741e-5"),
                                parse(Float64, "3.1365730552e-6"))
    end

    @testset "Snorm = 0 on read (Sourcefiletotable convention)" begin
        modes = read_source_table(F_I1; radius=20.0, ringmass=1e-4)
        @test all(m -> m.Snorm == 0, modes)
    end

    @testset "Sm = -Sm† antisymmetry (i1.0 circular-equatorial)" begin
        modes = read_source_table(F_I1; radius=20.0, ringmass=1e-4)
        @test all(m -> m.Sm == -m.Smdag, modes)
    end

    @testset "schema: 12 tab columns in mywrite order" begin
        for line in eachline(F_I1)
            isempty(strip(line)) && continue
            @test length(rawtokens(line)) == 12
        end
        modes = read_source_table(F_I1; radius=20.0, ringmass=1e-4)
        line = source_mode_line(modes[1])
        @test endswith(line, '\n')
        @test length(rawtokens(chomp(line))) == 12
    end

    @testset "write: omega reproduced byte-exact (8 sig figs)" begin
        modes = read_source_table(F_I1; radius=20.0, ringmass=1e-4)
        rawlines = [l for l in eachline(F_I1) if !isempty(strip(l))]
        for (mode, raw) in zip(modes, rawlines)
            orig_omega = rawtokens(raw)[4]
            written_omega = rawtokens(chomp(source_mode_line(mode)))[4]
            @test written_omega == orig_omega
        end
    end

    @testset "round-trip: read∘write is identity on parsed values" begin
        for f in (F_I1, F_I05_R20, F_I05_R400)
            modes = read_source_table(f; radius=1.0, ringmass=1.0)
            io = IOBuffer()
            for mode in modes
                write_source_mode(io, mode)
            end
            roundtrip = read_source_table(IOBuffer(take!(io)); radius=1.0, ringmass=1.0)
            @test length(roundtrip) == length(modes)
            for (a, b) in zip(modes, roundtrip)
                @test (a.m, a.n, a.k) == (b.m, b.n, b.k)
                @test a.ω == b.ω
                @test a.Sl == b.Sl
                @test a.Sn == b.Sn
                @test a.Sm == b.Sm
                @test a.Smdag == b.Smdag
            end
        end
    end

    @testset "write_source_table writes all rows" begin
        modes = read_source_table(F_I1; radius=1.0, ringmass=1.0)
        io = IOBuffer()
        write_source_table(io, modes)
        roundtrip = read_source_table(IOBuffer(take!(io)); radius=1.0, ringmass=1.0)
        @test length(roundtrip) == length(modes)
    end

    @testset "parse: special zero token 0.e-N (i0.5)" begin
        modes = read_source_table(F_I05_R20; radius=20.0, ringmass=1e-4)
        @test modes[1].Sl == 0.0 + 0.0im
    end
end
