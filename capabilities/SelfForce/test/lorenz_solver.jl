# Solver — the VoP/two-point-shoot integrator that turns the green inputs into the h1 corpus.
#
# Grounding:
#   • PRIMARY characterization vs the CAPTURED genuine legacy h1 oracle fixtures under
#     test/fixtures/h1lorenz/h1-l{l}m{m}.h5 (the H1Lorenz 4-part writer layout — distinct from the
#     6-part SecondOrderRicci corpus under fixtures/h1/ that SecondOrder.read_h1_corpus ingests). These were
#     captured ONE-SHOT, offline, from the legacy packages/H1Lorenz solver (HEAVY → never a live test
#     dep) on the small 7-point grid grid=[4,5,6,8,10,12,14], r0=8 (r0index=3), ra=0, rb=6, via
#       Pkg.develop(H1Lorenz); compute_one_mode(8.0,l,m; grid, indices=(0,3,6)); write_mode_h5(...)
#     (scratchpad/capture_h1*.jl). The test reads GRID + r0/ra/rb FROM each
#     fixture, solves that mode, builds inhom_left/right in the loader layout (h1_schema), and asserts
#     maxabs/rel vs the fixture's inhom_left/inhom_right at the LEGACY tolerance (5e-8 for l≤2, 1e-7
#     high-l — verbatim from H1Lorenz/test/ref_oracle.jl). Representative spread (one per coupled set):
#     (2,2)=even, (2,1)=odd, (1,1)=even_dipole, (1,0)=odd_dipole, (0,0)=monopole, (3,3).
#   • h1-corpus WRITE → Data.read_h1 round-trip (HDF5 live dep) reproduces the written schema.
#   • static cross-check: compute_mode(0,0).inhom == construct_static_monopole (analytic identity).
#   • BigFloat precision inheritance of the analytic monopole (same field, widened precision).

using Test
using SelfForce
import Data
using HDF5

const LF = SelfForce.LorenzField

# maxabs + masked-relative error, matching H1Lorenz/test/ref_oracle.jl _maxabs_relmask.
function _maxabs_relmask(A, B; thresh=1e-3)
    s = abs.(A .- B); mab = maximum(s)
    mask = abs.(B) .> thresh
    mr = any(mask) ? maximum((s ./ max.(abs.(B), 1e-300))[mask]) : 0.0
    (mab, mr)
end

_fixture(name) = joinpath(@__DIR__, "fixtures", "h1lorenz", name)

@testset "LorenzField.Solver" begin

    @testset "h1 oracle characterization (genuine legacy fixtures)" begin
        cases = [(2,2),(2,1),(1,1),(1,0),(0,0),(3,3)]
        for (l, m) in cases
            path = _fixture("h1-l$(l)m$(m).h5")
            @test isfile(path)
            grid, r0i, ra, rb, Lo, Ro = h5open(path, "r") do f
                g = f["grid"]; a = attrs(g)
                (read(g), Int(a["r0index"]), Int(a["raindex"]), Int(a["rbindex"]),
                 read(f["inhom_left"]), read(f["inhom_right"]))
            end
            any(isnan, Lo) && continue   # NaN even_static oracle (not in this spread)
            p = grid[r0i+1]
            nmode, orbit = LF.compute_one_mode(p, l, m; grid=grid, indices=(ra, r0i, rb))
            s = LF.h1_schema(nmode, orbit)
            @test size(s.inhom_left) == size(Lo)
            @test size(s.inhom_right) == size(Ro)
            mabL, mrL = _maxabs_relmask(s.inhom_left, Lo)
            mabR, mrR = _maxabs_relmask(s.inhom_right, Ro)
            tol = l >= 4 ? 1e-7 : 5e-8
            @test mabL < tol
            @test mabR < tol
            @test mrL  < tol
            @test mrR  < tol
        end
    end

    @testset "h1-corpus write → Data.read_h1 round-trip" begin
        path = _fixture("h1-l2m2.h5")
        grid, r0i, ra, rb = h5open(path, "r") do f
            g = f["grid"]; a = attrs(g)
            (read(g), Int(a["r0index"]), Int(a["raindex"]), Int(a["rbindex"]))
        end
        p = grid[r0i+1]
        nmode, orbit = LF.compute_one_mode(p, 2, 2; grid=grid, indices=(ra, r0i, rb))
        s = LF.h1_schema(nmode, orbit)
        mktempdir() do dir
            out = LF.write_h1_mode(joinpath(dir, "h1-l2m2.h5"), nmode, orbit)
            back = Data.read_h1(out)
            @test back.r0index == r0i
            @test back.ra_index == ra
            @test back.rb_index == rb
            @test back.grid ≈ Float64.(grid)
            @test back.inhom_left  == s.inhom_left      # exact Float64 round-trip
            @test back.inhom_right == s.inhom_right
            @test size(back.inhom_right, 1) ÷ 4 == nmode.nf   # nf == 7 for even (5 coupled + 2 gauge)
        end
    end

    @testset "static l=m=0 == construct_static_monopole (analytic identity)" begin
        grid = [4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0]; ra, r0i, rb = 0, 3, 6
        p = grid[r0i+1]
        orbit = LF.OrbitalParams(p, grid, ra, r0i, rb)
        nmode = LF.compute_mode(orbit, 0, 0)
        mono = LF.construct_static_monopole(orbit)
        @test nmode.cset_kind === :monopole
        @test nmode.inhom == mono
        @test nmode.omega == 0
    end

    @testset "BigFloat precision inheritance (analytic monopole)" begin
        grid = [4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0]; ra, r0i, rb = 0, 3, 6
        nf64, ob64 = LF.compute_one_mode(8.0, 0, 0; grid=grid, indices=(ra, r0i, rb))
        nfbig, obbig = LF.compute_one_mode(8.0, 0, 0; grid=grid, indices=(ra, r0i, rb), bits=128)
        @test eltype(nfbig.inhom) == BigFloat
        @test maximum(abs.(Float64.(nfbig.inhom) .- nf64.inhom)) < 1e-12
    end
end
