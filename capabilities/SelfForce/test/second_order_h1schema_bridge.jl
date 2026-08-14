# Production H1Schema (4 rows/field) -> SecondOrder radial-jet bridge.
#
# The live Lorenz writer stores raw trace-reversed Barack-Lousto R_i and dR_i/dr.  SecondOrder's
# kernels consume the physical, non-trace-reversed h_i radial jet, including d²h_i/dr².  These tests
# pin the adapter at both sides of that convention boundary and retain the legacy 6-row tests in
# second_order_readcorpus.jl as an independent compatibility gate.

using Test
using SelfForce
import Data
import HDF5
using NumericsCore: a_il

const _H1SCHEMA_FIX = joinpath(@__DIR__, "fixtures", "h1lorenz")
const _H1LEGACY_FIX = joinpath(@__DIR__, "fixtures", "h1")
const _BRIDGE_SO = SelfForce.SecondOrder
const _BRIDGE_PF = SelfForce.Puncture
const _BRIDGE_LF = SelfForce.LorenzField

function _bridge_one_mode(f, l, m)
    mktempdir() do dir
        cp(joinpath(_H1SCHEMA_FIX, "h1-l$(l)m$(m).h5"),
           joinpath(dir, "h1-l$(l)m$(m).h5"))
        return f(dir)
    end
end

# Expand one H1Schema split table into the duplicated-r0 raw trace-reversed field arrays.
function _bridge_raw_bar(s, l, m)
    fields = Data.h1_parity_fields(l, m)
    N_left = size(s.inhom_left, 2)
    N = N_left + size(s.inhom_right, 2)
    hbar = Dict(i => zeros(ComplexF64, N) for i in fields)
    dhbar = Dict(i => zeros(ComplexF64, N) for i in fields)
    for (k, i) in enumerate(fields)
        c = 4 * (k - 1)
        for j in axes(s.inhom_left, 2)
            hbar[i][j] = complex(s.inhom_left[c + 1, j], s.inhom_left[c + 2, j])
            dhbar[i][j] = complex(s.inhom_left[c + 3, j], s.inhom_left[c + 4, j])
        end
        for j in axes(s.inhom_right, 2)
            q = N_left + j
            hbar[i][q] = complex(s.inhom_right[c + 1, j], s.inhom_right[c + 2, j])
            dhbar[i][q] = complex(s.inhom_right[c + 3, j], s.inhom_right[c + 4, j])
        end
    end
    return fields, hbar, dhbar
end

@testset "SecondOrder production 4-row H1Schema bridge" begin
    @testset "real nonstatic fixture: typed jet, conversion, reflection" begin
        _bridge_one_mode(2, 2) do dir
            jets = _BRIDGE_SO.read_h1(dir; l_max = 2)
            @test jets isa _BRIDGE_SO.H1RadialJetCorpus{Float64}
            @test length(jets.r) == 8
            @test all(v -> all(isfinite, v), values(jets.h))
            @test all(v -> all(isfinite, v), values(jets.dh))
            @test all(v -> all(isfinite, v), values(jets.ddh))

            # Raw storage i=3 becomes physical i=6; d2 is reconstructed before this conversion.
            path = joinpath(dir, "h1-l2m2.h5")
            s = Data.read_h1(path)
            fields, hbar, dhbar = _bridge_raw_bar(s, 2, 2)
            j = 2
            expected_h6 = a_il(6, 2) * hbar[3][j] / jets.r[j]
            expected_dh6 = a_il(6, 2) * (dhbar[3][j] - hbar[3][j] / jets.r[j]) / jets.r[j]
            @test jets.h[(6, 2, 2)][j] ≈ expected_h6 rtol = 2e-15
            @test jets.dh[(6, 2, 2)][j] ≈ expected_dh6 rtol = 2e-15

            for i in fields
                ih = i == 3 ? 6 : (i == 6 ? 3 : i)
                @test jets.h[(ih, 2, -2)] ≈ conj.(jets.h[(ih, 2, 2)]) rtol = 2e-15
                @test jets.dh[(ih, 2, -2)] ≈ conj.(jets.dh[(ih, 2, 2)]) rtol = 2e-15
                @test jets.ddh[(ih, 2, -2)] ≈ conj.(jets.ddh[(ih, 2, 2)]) rtol = 2e-15
            end
        end
    end

    @testset "ddbar agrees with independent Puncture field equation" begin
        _bridge_one_mode(2, 2) do dir
            jets = _BRIDGE_SO.read_h1(dir; l_max = 2)
            s = Data.read_h1(joinpath(dir, "h1-l2m2.h5"))
            fields, hbar, dhbar = _bridge_raw_bar(s, 2, 2)
            residual = 0.0
            for j in eachindex(jets.r)
                hv = zeros(ComplexF64, 10)
                dv = zeros(ComplexF64, 10)
                for i in fields
                    hv[i] = hbar[i][j]
                    dv[i] = dhbar[i][j]
                end
                for i in fields
                    ih = i == 3 ? 6 : (i == 6 ? 3 : i)
                    ddbar = _BRIDGE_PF.ddhret(i, 2, 2, jets.r[j], jets.r0, 1.0, hv, dv)
                    expected = a_il(ih, 2) *
                               (ddbar - 2 * (dv[i] - hv[i] / jets.r[j]) / jets.r[j]) /
                               jets.r[j]
                    residual = max(residual, abs(jets.ddh[(ih, 2, 2)][j] - expected))
                end
            end
            @test residual < 2e-16
        end
    end

    @testset "static routing is explicit" begin
        # The analytic monopole is covered by d2hdr2's static-mass identity and is safe to close.
        _bridge_one_mode(0, 0) do dir
            jets = _BRIDGE_SO.read_h1(dir; l_max = 0)
            @test all(v -> all(isfinite, v), values(jets.ddh))
        end

        # l>=2 even-static R6/R7 require solver-provided/validated second derivatives.
        _bridge_one_mode(2, 0) do dir
            err = try
                _BRIDGE_SO.read_h1(dir; l_max = 2)
                nothing
            catch ex
                ex
            end
            @test err isa ArgumentError
            msg = sprint(showerror, err)
            @test occursin("even-static", msg)
            @test occursin("legacy 6-row radial-jet", msg)
        end
    end

    @testset "malformed row layout is rejected" begin
        s = Data.read_h1(joinpath(_H1SCHEMA_FIX, "h1-l2m2.h5"))
        mktempdir() do dir
            path = joinpath(dir, "h1-l2m2.h5")
            HDF5.h5open(path, "w") do f
                f["grid"] = s.grid
                f["inhom_left"] = s.inhom_left[1:end-1, :]
                f["inhom_right"] = s.inhom_right[1:end-1, :]
            end
            @test_throws DimensionMismatch _BRIDGE_SO.read_h1(dir; l_max = 2)
        end
    end

    @testset "all selected modes must share the exact radial grid" begin
        mktempdir() do dir
            cp(joinpath(_H1SCHEMA_FIX, "h1-l2m2.h5"), joinpath(dir, "h1-l2m2.h5"))
            s = Data.read_h1(joinpath(_H1SCHEMA_FIX, "h1-l3m3.h5"))
            drifted_grid = copy(s.grid)
            drifted_grid[2] = nextfloat(drifted_grid[2])
            drifted = Data.H1Schema{Float64}(
                drifted_grid, s.inhom_left, s.inhom_right, s.C_inf, s.C_horiz,
                s.r0index, s.ra_index, s.rb_index)
            Data.write_mode_h5(joinpath(dir, "h1-l3m3.h5"), drifted)

            err = try
                _BRIDGE_SO.read_h1(dir; l_max = 3)
                nothing
            catch ex
                ex
            end
            @test err isa ArgumentError
            msg = sprint(showerror, err)
            @test occursin("radial grid mismatch for (l=3,m=3) at index 2", msg)
            @test occursin("first selected mode (l=2,m=2)", msg)
        end
    end

    @testset "mixed four/six-row corpus and BigFloat four-row path" begin
        mktempdir() do dir
            cp(joinpath(_H1SCHEMA_FIX, "h1-l2m2.h5"), joinpath(dir, "h1-l2m2.h5"))
            cp(joinpath(_H1LEGACY_FIX, "h1-l3m3.h5"), joinpath(dir, "h1-l3m3.h5"))
            jets = _BRIDGE_SO.read_h1(dir; l_max = 3)
            @test jets isa _BRIDGE_SO.H1RadialJetCorpus{Float64}
            @test haskey(jets.h, (1, 2, 2))
            @test haskey(jets.h, (1, 3, 3))
            @test all(v -> all(isfinite, v), values(jets.ddh))
        end

        _bridge_one_mode(2, 2) do dir
            jets = _BRIDGE_SO.read_h1(BigFloat, dir; l_max = 2)
            @test jets isa _BRIDGE_SO.H1RadialJetCorpus{BigFloat}
            @test eltype(jets.r) == BigFloat
            @test eltype(first(values(jets.ddh))) == Complex{BigFloat}
            @test all(v -> all(isfinite, v), values(jets.ddh))
        end
    end

    @testset "LorenzField write -> radial jet -> finite R2 mode" begin
        grid = [4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0]
        nmode, orbit = _BRIDGE_LF.compute_one_mode(8.0, 2, 2;
                                                   grid = grid, indices = (0, 3, 6))
        mktempdir() do dir
            _BRIDGE_LF.write_h1_mode(joinpath(dir, "h1-l2m2.h5"), nmode, orbit)
            jets = _BRIDGE_SO.read_h1(dir; l_max = 2)
            # The sole available pair (2,2)x(2,2) couples to target (l,m)=(4,4).
            mode = _BRIDGE_SO.second_order_ricci(jets, 1, 4, 4;
                                                  h1lmax = 2, dlmax = 2)
            @test all(isfinite, mode)
            @test maximum(abs, mode) > 1e-6
        end
    end
end
