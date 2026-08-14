# Grounded tests for SelfForce.SecondOrder end-to-end Ricci-source assembly + writer.
#
# GROUNDING — captured g++/GSL oracle (exactly as legacy SecondOrderRicci grounded vs its
# r2_oracle.cpp / Ricci.cc capture) + HDF5 round-trip, no heavy live dep:
#
# (1) ORACLE: test/fixtures/src_self.h5 is the captured genuine compiled-C++ second-order Ricci
#     source (the upstream Ricci.cc driver, GSL coupling + R2_1..R2_10) over the captured h1 corpus,
#     invoked `./Ricci dir dir 2 2` — i.e. the SELF-coupling with corpus h1lmax=3 but targets
#     lmax=2, dlmax=2 (so l=3 first-order modes DO couple into the l≤2 targets via dlmax). Each
#     "src i={i} l={l} m={m}" dataset is on-disk shape [N][2]=(re,im) (Julia (2,N)). We read the FULL
#     l_max=3 corpus via SecondOrder.read_h1, assemble ricci_source(; lmax=2, dlmax=2), and compare
#     EVERY dataset to the oracle at rtol 1e-11. (Confirmed: the legacy SecondOrderRicci.ricci_source
#     reproduces this exact oracle at maxrel 1.5e-14 with these same args; Kernels.r2_term is itself
#     grounded vs the compiled R2 oracle in the Kernels core.)
# (2) The assembled mode set equals selection_modes(2) (the driver target set for lmax=2).
# (3) write_src → read HDF5 round-trip reproduces every dataset bit-for-bit, plus the "r" grid.
# (4) ricci_source_mass / ricci_source_spin: grounded vs captured oracle fixtures src_mass.h5 /
#     src_spin.h5 at rtol 1e-11. These were captured one-shot from the legacy SecondOrderRicci
#     ricci_source_mass/spin over the same l_max=3 corpus (lmax=2) — itself grounded element-wise vs
#     the genuine compiled Ricci_M.cc / Ricci_S.cc drivers at rtol 1e-10 (legacy ref_assembly_test).
#     Capture: scratchpad/capture_ms.jl (Pkg.activate packages/SecondOrderRicci → write_src).

using Test
using SelfForce.SecondOrder: read_h1, second_order_ricci, ricci_source, ricci_source_mass,
                             ricci_source_spin, write_src, selection_modes
import HDF5

const _FIXH1 = joinpath(@__DIR__, "fixtures", "h1")
const _ORACLE = joinpath(@__DIR__, "fixtures", "src_self.h5")

# Read an oracle "src i l m" dataset [2,N] back into a complex vector.
_oracle_vec(file, i, l, m) = (D = read(file["src i=$i l=$l m=$m"]); complex.(D[1, :], D[2, :]))

@testset "ricci_source vs captured g++/GSL oracle (rtol 1e-11)" begin
    fields = read_h1(_FIXH1; l_max = 3)       # full corpus h1lmax=3
    src = ricci_source(fields; lmax = 2, dlmax = 2)   # ./Ricci dir dir 2 2 invocation
    @test Set(keys(src)) == Set(selection_modes(2))
    HDF5.h5open(_ORACLE, "r") do file
        oracle_keys = [k for k in keys(file) if startswith(k, "src ")]
        @test length(oracle_keys) == length(selection_modes(2))
        for (i, l, m) in selection_modes(2)
            got = src[(i, l, m)]
            want = _oracle_vec(file, i, l, m)
            @test length(got) == length(want)
            @test got ≈ want rtol = 1e-11
        end
    end
end

@testset "negative dlmax preserves the legacy truncation sentinel" begin
    full = read_h1(_FIXH1; l_max = 3)
    truncated = read_h1(_FIXH1; l_max = 2)

    # Ricci.cc interprets dlmax=-1 at h1lmax=3 as: lower the usable H1 cutoff to 2, use dlmax=2,
    # then cap the requested target lmax at 2.  It must not become an empty bandwidth/all-zero source.
    sentinel = ricci_source(full; lmax = 3, dlmax = -1)
    explicit = ricci_source(truncated; lmax = 3, dlmax = 2)

    @test Set(keys(sentinel)) == Set(selection_modes(2))
    @test sentinel == explicit
    @test any(v -> any(x -> !iszero(x), v), values(sentinel))

    direct_sentinel = second_order_ricci(full, 1, 0, 0; h1lmax = 3, dlmax = -1)
    direct_explicit = second_order_ricci(full, 1, 0, 0; h1lmax = 2, dlmax = 2)
    @test direct_sentinel == direct_explicit
    @test any(x -> !iszero(x), direct_sentinel)
end

@testset "write_src → read HDF5 round-trip" begin
    fields = read_h1(_FIXH1; l_max = 3)
    src = ricci_source(fields; lmax = 2, dlmax = 2)
    tmp = joinpath(mktempdir(), "src_roundtrip.h5")
    write_src(tmp, src; r = fields.r)
    HDF5.h5open(tmp, "r") do file
        @test read(file["r"]) ≈ collect(float.(fields.r)) rtol = 1e-15
        for (i, l, m) in selection_modes(2)
            D = read(file["src i=$i l=$l m=$m"])
            @test size(D) == (2, length(fields.r))
            rt = complex.(D[1, :], D[2, :])
            @test rt == src[(i, l, m)]        # exact: write→read is lossless Float64
        end
    end
end

@testset "ricci_source_mass vs captured Ricci_M oracle (rtol 1e-11)" begin
    fields = read_h1(_FIXH1; l_max = 3)
    sm = ricci_source_mass(fields; lmax = 2)
    @test Set(keys(sm)) == Set(selection_modes(2))
    HDF5.h5open(joinpath(@__DIR__, "fixtures", "src_mass.h5"), "r") do file
        for (i, l, m) in selection_modes(2)
            @test sm[(i, l, m)] ≈ _oracle_vec(file, i, l, m) rtol = 1e-11
        end
    end
end

@testset "ricci_source_spin vs captured Ricci_S oracle (rtol 1e-11)" begin
    fields = read_h1(_FIXH1; l_max = 3)
    ss = ricci_source_spin(fields; lmax = 2)
    @test Set(keys(ss)) == Set(selection_modes(2))
    HDF5.h5open(joinpath(@__DIR__, "fixtures", "src_spin.h5"), "r") do file
        for (i, l, m) in selection_modes(2)
            @test ss[(i, l, m)] ≈ _oracle_vec(file, i, l, m) rtol = 1e-11
        end
    end
end
