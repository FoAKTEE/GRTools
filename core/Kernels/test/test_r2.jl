# Characterization of the 10 quadratic 2nd-order-Ricci kernels (Kernels.r2_term, n=1..10).
#
# GROUNDING — BAKED g++/GSL UPSTREAM-ORACLE FIXTURE (no heavy live dep; SecondOrderRicci is FORBIDDEN as a
# test dep). The fixture test/fixtures/r2_oracle_ref.jl holds the verbatim ComplexF64 outputs of the
# GENUINE upstream kernels (ref-code/SecondOrderRicci/src/R2_<n>.cc + Coupling.cc + utils.cc, compiled via
# test/oracle/r2_oracle.cpp with -DKERNEL=R2_<n> against GSL+Boost) on the deterministic seeded fields
# below — captured one-shot by scratchpad/capture_r2.jl (the `capture_r2_oracle` harness mirrored here for
# reproducibility). Branch coverage IS the merge gate: every reachable (l1,l2) low-l branch
# (≥2,≥2),(0,0),(1,0),(0,1),(1,1),(0,≥2),(≥2,0),(1,≥2),(≥2,1) is exercised for all 10 n. The (l1=0,l2=0)
# branch of R2_4/R2_5 is identically zero upstream — exact-zero reproduction is part of the gate. Tested
# element-wise at the legacy tolerance rtol=1e-11 atol=1e-12, PLUS the parallel===serial bit-identity gate.
#
# COUPLING PURITY (design §5.5/§7.7.1): r2_term takes the PRECOMPUTED coupling array `c` — Wigner is NOT
# a Kernels dep, so these tests feed the BAKED coupling fixture test/fixtures/r2_coupling_ref.jl (the
# c-vectors the pre-hoist committed kernels evaluated internally via Wigner.coupling; the coupling-
# evaluating host wrapper now lives in SelfForce.SecondOrder). R2_COUPLING_SPEC_REF gates
# Kernels.R2_COUPLING_SPEC against transcription drift.

using Test
using Kernels

# Fixture-regeneration inputs only (capture_r2_oracle below; NOT used by the CI tests).
# Override via ENV at a site that has the upstream SecondOrderRicci mirror; the recorded
# defaults are the capture site's paths.
const SOR_ROOT = get(ENV, "GRTOOLS_SECONDORDERRICCI_ROOT",
                     "/data/haiyangw/claude/emri-disk/ref-code/SecondOrderRicci")
const SOR_INC = joinpath(SOR_ROOT, "include")
const SOR_SRC = joinpath(SOR_ROOT, "src")
const R2_ORACLE_SRC = get(ENV, "GRTOOLS_R2_ORACLE_CPP",
                          "/data/haiyangw/claude/emri-disk/packages/SecondOrderRicci/test/oracle/r2_oracle.cpp")

# (l1,l2) representatives hitting each reachable branch, with (l3,m3,m1,m2) giving a nonzero c1.
const R2_BRANCHES = [
    (l1=2, l2=2, l3=2, m3=2, m1=1, m2=1),   # l1>=2 && l2>=2
    (l1=0, l2=0, l3=0, m3=0, m1=0, m2=0),   # l1==0 && l2==0
    (l1=1, l2=0, l3=1, m3=1, m1=1, m2=0),   # l1==1 && l2==0
    (l1=0, l2=1, l3=1, m3=1, m1=0, m2=1),   # l1==0 && l2==1
    (l1=1, l2=1, l3=2, m3=2, m1=1, m2=1),   # l1==1 && l2==1
    (l1=0, l2=2, l3=2, m3=1, m1=0, m2=1),   # l1==0 (l2>=2)
    (l1=2, l2=0, l3=2, m3=1, m1=1, m2=0),   # l2==0 (l1>=2)
    (l1=1, l2=2, l3=3, m3=2, m1=1, m2=1),   # l1==1 (l2>=2)
    (l1=2, l2=1, l3=3, m3=2, m1=1, m2=1),   # l2==1 (l1>=2)
]
const RGRID = [4.0, 6.0, 8.0, 10.0, 12.0]
const M_ = 1.0
const R0_ = 8.0

# Deterministic seeded complex field value for (tag,i,j): distinct, O(1), nonzero — identical to the
# legacy ground-truth test so the baked oracle and the live Kernels call see the SAME inputs.
seedc(tag, i, j) = complex(0.3 + 0.11 * tag + 0.07 * i + 0.05 * j,
                           -0.2 + 0.09 * tag - 0.04 * i + 0.03 * j)

function make_fields(b, N)
    h1  = Dict((i, b.l1, b.m1) => [seedc(1, i, j) for j in 1:N] for i in 1:10)
    dh1 = Dict((i, b.l1, b.m1) => [seedc(2, i, j) for j in 1:N] for i in 1:10)
    ddh1 = Dict((i, b.l1, b.m1) => [seedc(3, i, j) for j in 1:N] for i in 1:10)
    h2  = Dict((i, b.l2, b.m2) => [seedc(4, i, j) for j in 1:N] for i in 1:10)
    dh2 = Dict((i, b.l2, b.m2) => [seedc(5, i, j) for j in 1:N] for i in 1:10)
    ddh2 = Dict((i, b.l2, b.m2) => [seedc(6, i, j) for j in 1:N] for i in 1:10)
    return h1, dh1, ddh1, h2, dh2, ddh2
end

# capture_r2_oracle (NEW test helper, mirrored from scratchpad/capture_r2.jl): the one-shot harness that
# REGENERATES the baked fixture by compiling the genuine upstream g++/GSL oracle and running it on the
# seeded inputs. Not exercised in CI (it needs g++/GSL + ref-code); kept here so the fixture is
# reproducible from source. Returns R2_ORACLE_REF[n][branch] :: Vector{ComplexF64}.
function capture_r2_oracle(; ns = 1:10)
    f = [1 - 2M_ / r for r in RGRID]
    fp = [2M_ / r^2 for r in RGRID]
    gsl_cflags = split(strip(read(`gsl-config --cflags`, String)))
    gsl_libs = split(strip(read(`gsl-config --libs`, String)))
    ref = Vector{Vector{ComplexF64}}[]
    for n in ns
        bin = tempname()
        srcs = [R2_ORACLE_SRC, joinpath(SOR_SRC, "R2_$(n).cc"),
                joinpath(SOR_SRC, "Coupling.cc"), joinpath(SOR_SRC, "utils.cc")]
        run(`g++ -O2 -std=c++14 -DKERNEL=R2_$(n) -o $bin $srcs
             -I$SOR_INC -I/usr/include $gsl_cflags $gsl_libs`)
        branches = Vector{ComplexF64}[]
        for b in R2_BRANCHES
            h1, dh1, ddh1, h2, dh2, ddh2 = make_fields(b, length(RGRID))
            io = IOBuffer()
            print(io, "$(b.l3) $(b.m3) $(b.l1) $(b.m1) $(b.l2) $(b.m2) $(length(RGRID)) $M_ $R0_\n")
            print(io, join(RGRID, " "), "\n"); print(io, join(f, " "), "\n"); print(io, join(fp, " "), "\n")
            for (fld, l, m) in ((h1, b.l1, b.m1), (dh1, b.l1, b.m1), (ddh1, b.l1, b.m1),
                                (h2, b.l2, b.m2), (dh2, b.l2, b.m2), (ddh2, b.l2, b.m2))
                for i in 1:10, j in 1:length(RGRID)
                    c = fld[(i, l, m)][j]
                    print(io, real(c), " ", imag(c), " ")
                end
                print(io, "\n")
            end
            outs = read(pipeline(IOBuffer(String(take!(io))), `$bin`), String)
            rows = [parse.(Float64, split(line)) for line in split(strip(outs), "\n")]
            push!(branches, [complex(row[1], row[2]) for row in rows])
        end
        rm(bin; force = true)
        push!(ref, branches)
    end
    return ref
end

include(joinpath(@__DIR__, "fixtures", "r2_oracle_ref.jl"))    # const R2_ORACLE_REF[n][branch]
include(joinpath(@__DIR__, "fixtures", "r2_coupling_ref.jl"))  # const R2_COUPLING_{SPEC_REF,REF[n][branch]}

@testset "r2" begin
    f = [1 - 2M_ / r for r in RGRID]
    fp = [2M_ / r^2 for r in RGRID]

    @testset "R2_COUPLING_SPEC matches the pre-hoist committed coupling calls" begin
        @test Kernels.R2_COUPLING_SPEC == R2_COUPLING_SPEC_REF
    end

    @testset "R2_$n vs baked g++/GSL oracle, all (l1,l2) branches + parallel===serial" for n in 1:10
        any_nonzero = false
        for (bi, b) in enumerate(R2_BRANCHES)
            h1, dh1, ddh1, h2, dh2, ddh2 = make_fields(b, length(RGRID))
            oracle = R2_ORACLE_REF[n][bi]
            c = R2_COUPLING_REF[n][bi]
            jl = r2_term(n, M_, R0_, RGRID, f, fp, h1, dh1, ddh1, h2, dh2, ddh2,
                         b.l3, b.m3, b.l1, b.m1, b.l2, b.m2, c)
            @test length(jl) == length(oracle)
            @test all(isfinite, oracle)            # guard: branch is not NaN-poisoned (mirrors _ssqrt use)
            any_nonzero |= maximum(abs, oracle) > 0
            for j in 1:length(oracle)
                @test jl[j] ≈ oracle[j] rtol = 1e-11 atol = 1e-12
            end
            # parallel === serial bit-identity (matches test_aak / test_amp_gemm gate)
            jt = r2_term(n, M_, R0_, RGRID, f, fp, h1, dh1, ddh1, h2, dh2, ddh2,
                         b.l3, b.m3, b.l1, b.m1, b.l2, b.m2, c; threaded = true)
            @test all(jt[j] === jl[j] for j in eachindex(jl))
        end
        @test any_nonzero
    end

    @testset "r2_term coupling-array length validation" begin
        b = R2_BRANCHES[1]
        h1, dh1, ddh1, h2, dh2, ddh2 = make_fields(b, length(RGRID))
        @test_throws ArgumentError r2_term(1, M_, R0_, RGRID, f, fp, h1, dh1, ddh1, h2, dh2, ddh2,
                                           b.l3, b.m3, b.l1, b.m1, b.l2, b.m2,
                                           zeros(4))   # R2_1 needs length(R2_COUPLING_SPEC[1]) == 5
    end


    @testset "absent whole parent modes are exact bilinear zeros" begin
        b = R2_BRANCHES[1]
        h1, dh1, ddh1, h2, dh2, ddh2 = make_fields(b, length(RGRID))
        empty_mode = Dict{NTuple{3,Int},Vector{ComplexF64}}()
        z = zeros(ComplexF64, length(RGRID))
        for n in 1:10
            c = R2_COUPLING_REF[n][1]
            left = r2_term(n, M_, R0_, RGRID, f, fp,
                           empty_mode, empty_mode, empty_mode, h2, dh2, ddh2,
                           b.l3, b.m3, b.l1, b.m1, b.l2, b.m2, c)
            right = r2_term(n, M_, R0_, RGRID, f, fp,
                            h1, dh1, ddh1, empty_mode, empty_mode, empty_mode,
                            b.l3, b.m3, b.l1, b.m1, b.l2, b.m2, c;
                            threaded = true)
            @test left == z
            @test right == z
        end
    end

    @testset "component zero-padding and length validation precede @inbounds kernels" begin
        b = R2_BRANCHES[1]
        h1, dh1, ddh1, h2, dh2, ddh2 = make_fields(b, length(RGRID))
        missing_dh = empty(dh1)
        partial = r2_term(
            1, M_, R0_, RGRID, f, fp, h1, missing_dh, ddh1, h2, dh2, ddh2,
            b.l3, b.m3, b.l1, b.m1, b.l2, b.m2, R2_COUPLING_REF[1][1])
        @test all(isfinite, partial)

        short_h = copy(h1)
        short_h[(1, b.l1, b.m1)] = short_h[(1, b.l1, b.m1)][1:end-1]
        @test_throws DimensionMismatch r2_term(
            1, M_, R0_, RGRID, f, fp, short_h, dh1, ddh1, h2, dh2, ddh2,
            b.l3, b.m3, b.l1, b.m1, b.l2, b.m2, R2_COUPLING_REF[1][1])
    end
end
