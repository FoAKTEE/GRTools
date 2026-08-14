# Module: kerrampgrid + roman NN — S4-b: the two EMRIWaveforms in-src HDF5 reader bodies moved into
# DataHDF5Ext (design §7.2/§10: ALL grid reader/writer BODIES live in Data's HDF5 weakdep-ext).
#
# CHARACTERIZATION grounding (tests-first): fixtures/hdf5_readers_ref.jl was baked from the COMMITTED
# EMRIWaveforms loaders (_load_roman_nn @ src/AmplitudeModels.jl:55, KerrAmpSpline.read_kerrampgrid @
# src/KerrAmpSpline.jl:178, kerr_mode_table @ :137, repo commit f409fc7) BEFORE any code move — the
# moved Data bodies must reproduce it BIT-EXACTLY (pure Float64 file reads, no arithmetic).
# The on-disk grids are gitignored (ZNAmps is 4.74GB): absent files skip with @info, exactly like the
# existing ampgrid/fluxgrid on-disk groundings.

using Test
using Data
using HDF5   # fires DataHDF5Ext — read_kerrampgrid / read_roman_nn bodies activate

include(joinpath(@__DIR__, "fixtures", "hdf5_readers_ref.jl"))

const _FEWDATA_KAG = joinpath(@__DIR__, "..", "..", "..", "..",
    "ref-code", "FastEMRIWaveforms", ".fewdata")
const _ZNAMPS_KAG = joinpath(_FEWDATA_KAG, "ZNAmps_l10_m10_n55_DS2Outer.h5")
const _NN_KAG = joinpath(_FEWDATA_KAG, "SchwarzschildEccentricInput.hdf5")

@testset "kerrampgrid + roman NN (S4-b moved bodies)" begin

    @testset "kerr_mode_table (pre-change fixture)" begin
        modes, index = Data.kerr_mode_table(; lmax = 10, nmax = 55)
        @test length(modes) == KERR_MODE_TABLE_REF.n_modes
        @test modes[1] == KERR_MODE_TABLE_REF.first
        @test modes[550] == KERR_MODE_TABLE_REF.m550
        @test modes[end] == KERR_MODE_TABLE_REF.last
        @test index[(2, 2, 0)] == KERR_MODE_TABLE_REF.idx_220
        @test index[(3, 2, 1)] == KERR_MODE_TABLE_REF.idx_321
        @test index[(4, 0, -3)] == KERR_MODE_TABLE_REF.idx_40m3
        @test index[(10, 10, 55)] == KERR_MODE_TABLE_REF.idx_101055
    end

    @testset "read_kerrampgrid vs pre-change fixture" begin
        if isfile(_ZNAMPS_KAG)
            g = Data.read_kerrampgrid(_ZNAMPS_KAG; modes = KERRAMP_REF.modes)
            @test g isa Data.KerrAmpGrid{Float64}
            @test g.modes == KERRAMP_REF.modes
            @test (g.lmax, g.mmax, g.nmax) == (KERRAMP_REF.lmax, KERRAMP_REF.mmax, KERRAMP_REF.nmax)
            @test g.num_total_modes == KERRAMP_REF.num_total_modes
            @test g.z_knots == KERRAMP_REF.z_knots           # bit-exact (pure file reads)
            @test g.twA == KERRAMP_REF.twA
            @test g.tuA == KERRAMP_REF.tuA
            @test g.twB == KERRAMP_REF.twB
            @test g.tuB == KERRAMP_REF.tuB
            @test size(g.coeffsA) == KERRAMP_REF.sizeA
            @test size(g.coeffsB) == KERRAMP_REF.sizeB
            for (i1, i2, i3, i4, v) in KERRAMP_REF.coeffsA_probes
                @test g.coeffsA[i1, i2, i3, i4] == v
            end
            for (i1, i2, i3, i4, v) in KERRAMP_REF.coeffsB_probes
                @test g.coeffsB[i1, i2, i3, i4] == v
            end
        else
            @info "ZNAmps grid absent (gitignored 4.74GB) — skipping read_kerrampgrid fixture grounding"
        end
    end

    @testset "read_roman_nn vs pre-change fixture" begin
        if isfile(_NN_KAG)
            nn = Data.read_roman_nn(_NN_KAG)
            @test length(nn.weights) == ROMAN_REF.num_layers
            @test length(nn.biases) == ROMAN_REF.num_layers
            @test [size(w) for w in nn.weights] == ROMAN_REF.weight_sizes
            @test [length(b) for b in nn.biases] == ROMAN_REF.bias_lengths
            @test size(nn.reduced_basis) == ROMAN_REF.rb_size
            @test nn.transform_factor_inv == ROMAN_REF.transform_factor_inv
            @test nn.break_index == ROMAN_REF.break_index
            @test nn.num_teuk_modes == ROMAN_REF.num_teuk_modes
            for (L, i, j, v) in ROMAN_REF.weight_probes
                @test nn.weights[L][i, j] == v
            end
            for (L, i, v) in ROMAN_REF.bias_probes
                @test nn.biases[L][i] == v
            end
            for (i, j, v) in ROMAN_REF.rb_probes
                @test nn.reduced_basis[i, j] == v
            end
            # positional destructuring matches the legacy _load_roman_nn 6-tuple contract
            w, b, rb, tfi, bi, ntm = nn
            @test w === nn.weights && b === nn.biases && rb === nn.reduced_basis
            @test (tfi, bi, ntm) ==
                  (ROMAN_REF.transform_factor_inv, ROMAN_REF.break_index, ROMAN_REF.num_teuk_modes)
        else
            @info "SchwarzschildEccentricInput.hdf5 absent — skipping read_roman_nn fixture grounding"
        end
    end
end
