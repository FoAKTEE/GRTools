# Module 5: ampgrid — de-duped Teukolsky amplitude-grid loader over AmpGridSchema. GROUNDING:
#  * schwarzecc_amp_modes: closed-form mode-discovery (storage order, counts) — the FEW formula is the oracle.
#  * read_ampgrid: the checksum-verified FEW file Teuk_amps_a0.0_lmax_10_nmax_30_new.h5 — RAW HDF5 is the
#    oracle: grid shape, and the C-order+u-reversal reshape `idx=(i-1)*num_u+(num_u-j)+1` verified at nodes.
#  * AmpGridSchema round-trip: construct a small schema, validate.

using Test
using Data
using HDF5   # fires DataHDF5Ext — read_ampgrid body activates

const _AMPFILE = joinpath(@__DIR__, "..", "..", "..", "..",
    "ref-code", "FastEMRIWaveforms", ".fewdata", "Teuk_amps_a0.0_lmax_10_nmax_30_new.h5")

@testset "ampgrid" begin
    @testset "schwarzecc_amp_modes closed form" begin
        modes = Data.schwarzecc_amp_modes()
        @test modes[1] == (2, 0, -30)               # m=0 block first, n ascends from -nmax
        @test (2, 2, 0) in modes
        @test (10, 10, 30) in modes
        @test length(modes) == 3843                 # sum_{l=2}^{10}(l+1)*61 = 63*61
        @test count(t -> t[2] == 0, modes) == 549   # (lmax-1)*(2nmax+1) = 9*61
        @test all(t -> t[2] == 0, modes[1:549])     # the leading block is exactly the m=0 modes
    end

    @testset "AmpGridSchema round-trip" begin
        s = AmpGridSchema{Float64}([0.0, 0.5, 1.0], [0.0, 1.0, 2.0], [(2, 2, 0)])
        @test validate_schema(s)
        @test schema_version(s) isa Integer
        @test s.modes == [(2, 2, 0)]
    end

    @testset "read_ampgrid vs raw HDF5 (checksum-verified file)" begin
        if isfile(_AMPFILE)
            # raw oracle
            valr, vali = h5open(_AMPFILE, "r") do f
                @test size(f["grid"]) == (4, 1650)        # Julia-order read of the (1650,4) dataset
                d = read(f["l2m2/n0k0"])                  # (2, N)
                (d[1, :], d[2, :])
            end

            schema, re, im = Data.read_ampgrid(_AMPFILE; modes = [(2, 2, 0)])
            @test schema.modes == [(2, 2, 0)]
            @test length(schema.w) == 33                  # unique e
            @test length(schema.u) == 50                  # unique y(p,e)
            @test validate_schema(schema)
            num_u = length(schema.u); num_w = length(schema.w)
            @test size(re[1]) == (num_w, num_u)

            # C-order + u-reversal reshape: idx=(i-1)*num_u+(num_u-j)+1, verified at corner/edge nodes.
            for (i, j) in ((1, 1), (1, 50), (2, 50), (33, 1), (17, 25))
                idx = (i - 1) * num_u + (num_u - j) + 1
                @test re[1][i, j] == valr[idx]
                @test im[1][i, j] == vali[idx]
            end
        else
            @info "FEW amplitude grid absent — skipping read_ampgrid on-disk grounding"
        end
    end
end
