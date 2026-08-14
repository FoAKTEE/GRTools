# Module 3: wl — Mathematica InputForm numeric reader. GROUNDING: closed-form fixtures — small
# hand-authored `.wl`/`.m` strings exercising each grammar feature; the InputForm grammar is the
# oracle, so every expected value is a hand-known literal. WaSABI is HEAVY → NO live dep.

using Test
using Data

# Write `s` to a temp file with extension `ext`, return the path.
function _wlfile(s::AbstractString; ext = ".wl")
    p = joinpath(mktempdir(), "fixture" * ext)
    write(p, s)
    return p
end

@testset "wl" begin
    @testset "read_mma_inputform grammar features" begin
        # Nested list.
        @test Data.read_mma_inputform(_wlfile("{{1.0, 2.0}, {3.0, 4.0}}")) ==
              Any[Any[1.0, 2.0], Any[3.0, 4.0]]
        # Association <|"k"->v|>.
        d = Data.read_mma_inputform(_wlfile("<|\"k\" -> 1.5, \"j\" -> 2.0|>"))
        @test d isa Dict && d["k"] == 1.5 && d["j"] == 2.0
        # Around[c, e] -> c (Around-stripping).
        @test Data.read_mma_inputform(_wlfile("Around[6.5, 0.1]")) == 6.5
        # *^ exponent.
        @test Data.read_mma_inputform(_wlfile("2.486*^-7")) == 2.486e-7
        # `precision backtick + *^ exponent.
        @test Data.read_mma_inputform(_wlfile("5.469`96.37*^6")) == 5.469e6
        # (* nested (* comment *) *) stripping.
        @test Data.read_mma_inputform(_wlfile("(* outer (* inner *) *) 42.0")) == 42.0
    end

    @testset "read_wl_table flat pairs" begin
        @test Data.read_wl_table(_wlfile("{{1.0, 2.0}, {3.0, 4.5}}")) ==
              [(1.0, 2.0), (3.0, 4.5)]
        # *^ exponents honoured in the table.
        @test Data.read_wl_table(_wlfile("{{1.0, 2.486*^-7}}"))[1][2] ≈ 2.486e-7
    end

    @testset "read_wl_data BigFloat full precision" begin
        digs = "3.141592653589793238462643383279502884197"   # 40 sig figs > Float64
        v = Data.read_wl_data(_wlfile("{<|\"x\" -> $(digs)`50, \"y\" -> Around[6.5, 0.1]|>}"))
        @test v isa Vector{Dict{String,Any}}
        @test length(v) == 1
        @test v[1]["x"] isa BigFloat
        @test v[1]["x"] == parse(BigFloat, digs)               # all file digits preserved
        @test v[1]["x"] != BigFloat(Float64(parse(BigFloat, digs)))  # truly beyond Float64
        @test v[1]["y"] == BigFloat(6.5)                       # Around-stripped
    end

    @testset "read_wl_grid2d / grid3d" begin
        @test Data.read_wl_grid2d(_wlfile("{{{1.0, 2.0}, 3.0}}")) == [((1.0, 2.0), 3.0)]
        @test Data.read_wl_grid3d(_wlfile("{{{1.0, 2.0, 3.0}, 4.0}}")) ==
              [((1.0, 2.0, 3.0), 4.0)]
    end

    @testset "numeric-dialect ownership (symbolic Num walker delegated to PostNewtonian)" begin
        # Format separation (design §10): Data.read_mma_inputform owns the NUMERIC InputForm dialect
        # ONLY — a record of nested lists/associations/numbers parses into plain Julia values, never a
        # symbolic `Symbolics.Num`. The SYMBOLIC `.m` SeriesData walker is PostNewtonian.MMAInputForm's,
        # so NO symbolic `series_data` stub exists in Data to mis-mark as a [HOLE].
        v = Data.read_mma_inputform(_wlfile(
            "{<|\"p\" -> 12.0, \"flux\" -> 2.486*^-7|>, Around[6.5, 0.1]}"))
        @test v isa Vector{Any}
        @test v[1] isa Dict{String,Any}
        @test v[1]["p"] === 12.0 && v[1]["flux"] ≈ 2.486e-7
        @test v[2] === 6.5
        # Every leaf is a concrete Julia number — and Data declares NO Symbolics dep (the symbolic
        # walker's home, PostNewtonian, owns that quarantine), so the numeric path stays Symbolics-free.
        @test v[1]["p"] isa Float64
        @test Base.identify_package(Data, "Symbolics") === nothing
    end

    @testset "read_wl_table_complex / _wl_complex" begin
        out = Data.read_wl_table_complex(_wlfile("{{1.0, 2.0 + 3.0*I}, {2.0, 4.0 - 5.0*I}}"))
        @test out == [(1.0, 2.0 + 3.0im), (2.0, 4.0 - 5.0im)]
        @test Data._wl_complex("5.0") == 5.0 + 0.0im
        @test Data._wl_complex("4.0*I") == 0.0 + 4.0im
        @test Data._wl_complex("I") == 0.0 + 1.0im
        @test Data._wl_complex("-I") == 0.0 - 1.0im
        @test Data._wl_complex("2.0 - 3.0*I") == 2.0 - 3.0im
        @test Data._wl_complex("1.5*^-2*I") == 0.0 + 0.015im     # *^ exponent in imag part
    end
end
