# Module: fluxgrid — de-duped Kerr eccentric-equatorial flux-grid loader over FluxGridSchema.
# GROUNDING: characterization vs a DIRECT in-test raw HDF5 read of the checksum-verified FEW file
#   KerrEccEqFluxData.h5. The region axis literals (NU/NW/NZ) and the HDF5 dim-reversal +
#   permutedims(3,2,1) stride are read verbatim from packages/FastEMRIWaveforms/src/trajectory/
#   kerrflux.jl (_region_grid / _read_region) — NO heavy FastEMRIWaveforms dep pulled; HDF5 is the
#   live light test-dep that fires DataHDF5Ext.read_fluxgrid.

using Test
using Data
using HDF5   # fires DataHDF5Ext — read_fluxgrid body activates

const _FLUXFILE = joinpath(@__DIR__, "..", "..", "..", "..",
    "ref-code", "FastEMRIWaveforms", ".fewdata", "KerrEccEqFluxData.h5")

@testset "fluxgrid" begin
    @testset "FluxGridSchema round-trip" begin
        s = FluxGridSchema{Float64}(:A, [0.0, 0.5, 1.0], [0.0, 0.5, 1.0], [0.0, 0.5, 1.0])
        @test validate_schema(s)
        @test schema_version(s) isa Integer
        @test s.region == :A
    end

    @testset "read_fluxgrid vs raw HDF5 (checksum-verified file)" begin
        if isfile(_FLUXFILE)
            for (region, name, NU, NW, NZ) in (
                (:A, "regionA", 129, 65, 65),
                (:B, "regionB", 65, 33, 33),
            )
                # raw oracle: read Edot/Ldot exactly as the FEW loader does.
                Eraw, Lraw = h5open(_FLUXFILE, "r") do f
                    (read(f["$name/Edot"]), read(f["$name/Ldot"]))
                end
                Eref = permutedims(Eraw, (3, 2, 1))
                Lref = permutedims(Lraw, (3, 2, 1))

                schema, ch = Data.read_fluxgrid(_FLUXFILE, region)
                @test schema isa FluxGridSchema{Float64}
                @test schema.region == region
                # axes are the canonical [0,1] tricubic ranges (kerrflux._region_grid literals).
                @test schema.u == collect(range(0.0, 1.0; length = NU))
                @test schema.w == collect(range(0.0, 1.0; length = NW))
                @test schema.z == collect(range(0.0, 1.0; length = NZ))
                @test validate_schema(schema)
                # channels are byte-identical to the permuted raw read.
                @test ch[:Edot] == Eref
                @test ch[:Ldot] == Lref
                @test size(ch[:Edot]) == (NU, NW, NZ)
            end
        else
            @info "FEW flux grid absent — skipping read_fluxgrid on-disk grounding"
        end
    end

    @testset "read_fluxgrid downsample stride" begin
        if isfile(_FLUXFILE)
            ds = (2, 2, 2)
            schema, ch = Data.read_fluxgrid(_FLUXFILE, :A; downsample = ds)
            full_u = collect(range(0.0, 1.0; length = 129))
            @test schema.u == full_u[1:2:end]
            @test length(schema.w) == length(collect(range(0.0, 1.0; length = 65))[1:2:end])
            # strided channel matches the strided permuted oracle.
            Eraw = h5open(_FLUXFILE, "r") do f
                read(f["regionA/Edot"])
            end
            Eref = permutedims(Eraw, (3, 2, 1))[1:2:end, 1:2:end, 1:2:end]
            @test ch[:Edot] == Eref
            @test size(ch[:Edot]) == (length(schema.u), length(schema.w), length(schema.z))
        end
    end

    @testset "read_fluxgrid invalid region" begin
        @test_throws ArgumentError Data.read_fluxgrid(_FLUXFILE, :C)
    end
end
