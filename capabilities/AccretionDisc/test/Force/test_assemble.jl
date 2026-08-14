# Module 6 — Force/Assemble.jl. Parity-aware BLS field → metric perturbation h_μν(r,θ) (DDKW dual-leg
# assembly, the iter-245 normalization-fix form). GROUNDING:
#   • kerr_parity truth table; _f/_fK closed forms (_f(2)=0, _f(r)=1−2/r, _fK(r,0)=_f(r)); ladder
#     identities (_λ̂(l)²=l(l+1), _Λ̂(2)=√24); _sY at θ=π/2 vs the captured legacy SWSH literal (the
#     Harmonics oracle reproduces the legacy spin-weighted spherical-Y).
#   • the FULL hmunu / tetmodes_from_bls tensor — CHARACTERIZED vs a CAPTURED one-shot
#     DiscEMRIForce.hmunu fixture (a synthetic H1Mode with baked Complex literals at representative
#     (l,m,r,θ), a∈{0,0.5}, even + odd parity). NOT a live DiscEMRIForce dep (HDF5-heavy). Captured
#     once by packages/DiscEMRIForce hmunu/tetmodes_from_bls; legacy reconstruction tolerance 1e-9.

using Test
using AccretionDisc
using AccretionDisc.Force
using Data
const F = AccretionDisc.Force

# Synthetic H1Mode (single radial index k=1 at r), identical layout to the capture script.
function _synth(::Type{T}, nf; r = 8.0) where {T}
    vals = [[Complex{T}(T(0.1) * j + T(1)/100, -T(0.05) * j + T(1)/50)] for j in 1:nf]
    drs  = [[Complex{T}(T(0.001) * j, T(0.002) * j)] for j in 1:nf]
    par = nf == 7 ? :even : :odd
    m = par === :even ? 2 : 1
    F.H1Mode{T}(2, m, T[r], vals, drs, T[r], vals, drs, T(r), T[], nf, par)
end

function _assemble_schema(nf; gs = 6, r0index = 2)
    grid = collect(10.0:(10.0 + gs - 1))
    nL = r0index + 1; nR = gs - r0index
    left = Float64[row * 1000 + col for row in 1:(4nf), col in 1:nL]
    right = Float64[5_000_000 + row * 1000 + col for row in 1:(4nf), col in 1:nR]
    Data.H1Schema{Float64}(grid, left, right, collect(1.0:nf), collect(101.0:(100.0 + nf)),
                           r0index, 0, gs - 1)
end

@testset "Force.Assemble" begin

    @testset "kerr_parity / _f / _fK / ladder closed forms" begin
        @test F.kerr_parity(2, 2) === :even
        @test F.kerr_parity(3, 1) === :even
        @test F.kerr_parity(2, 1) === :odd
        @test F.kerr_parity(3, 2) === :odd
        @test F._f(2.0) == 0.0
        @test F._f(8.0) ≈ 1 - 2/8.0 rtol = 1e-15
        @test F._fK(8.0, 0.0) ≈ F._f(8.0) rtol = 1e-15
        @test F._fK(8.0, 0.5) ≈ (8.0^2 - 16 + 0.25)/8.0^2 rtol = 1e-15
        @test F._λ̂(2)^2 ≈ 2 * 3 rtol = 1e-14
        @test F._λ̂(5)^2 ≈ 5 * 6 rtol = 1e-14
        @test F._Λ̂(2) ≈ sqrt(24.0) rtol = 1e-14
        @test F._Λ̂(2) ≈ 4.898979485566356 rtol = 1e-13      # captured seed
        @test F._λ̂(2) ≈ 2.449489742783178 rtol = 1e-13
    end

    @testset "_sY vs captured legacy SWSH literals (Harmonics oracle)" begin
        @test F._sY(2, 2, 2, π/2)  ≈ 0.15769578262625997 rtol = 1e-12
        @test F._sY(-2, 2, 2, π/2) ≈ 0.15769578262626005 rtol = 1e-12
        @test F._sY(1, 2, 1, π/3)  ≈ -0.31539156525252 rtol = 1e-12
        @test F._sY(2, 1, 1, π/2) == 0.0 + 0.0im            # l<|s| guard
    end

    @testset "tetmodes_from_bls vs captured DiscEMRIForce (even l=2)" begin
        me = _synth(Float64, 7)
        p0 = F.tetmodes_from_bls(me, 1, 2; a = 0.0)
        @test p0.lplp ≈ (0.16000000000000003 - 0.0688888888888889im) rtol = 1e-12
        @test p0.mpmp ≈ (0.8328265125462807 - 0.37558842722675406im) rtol = 1e-12
        @test p0.htr  ≈ (-0.025 + 0.012500000000000002im) rtol = 1e-12
        p5 = F.tetmodes_from_bls(me, 1, 2; a = 0.5)
        @test p5.lplp ≈ (0.15834626432924376 - 0.0681768638084244im) rtol = 1e-12
    end

    @testset "hmunu vs captured DiscEMRIForce — even l=2 m=2 (θ=π/2)" begin
        me = _synth(Float64, 7)
        H0 = F.hmunu(me, 1, 2, 2, π/2; a = 0.0)
        @test H0[1,1] ≈ (0.010079342459042606 - 0.0039834527083641436im) rtol = 1e-9
        @test H0[1,4] ≈ (0.021245081111275432 + 0.045709113906077445im) rtol = 1e-9
        @test H0[3,3] ≈ (0.3901369440434214 - 0.15322210013586524im) rtol = 1e-9
        @test H0[4,4] ≈ (0.258803715355537 - 0.09399338915897615im) rtol = 1e-9
        @test H0[2,2] ≈ (-0.008476572766619991 + 0.004506532356937211im) rtol = 1e-9
        # hmunu_eq is the θ=π/2 entry point
        @test F.hmunu_eq(me, 1, 2, 2) ≈ H0 rtol = 1e-12
        H5 = F.hmunu(me, 1, 2, 2, π/2; a = 0.5)
        @test H5[1,1] ≈ (0.00983122365159975 - 0.004733269971610854im) rtol = 1e-9
        @test H5[1,4] ≈ (0.01427883886024403 + 0.04882325662888178im) rtol = 1e-9
        @test H5[3,3] ≈ (0.3938139231175199 - 0.1549097013127491im) rtol = 1e-9
        @test H5[2,2] ≈ (-0.008508274170698073 + 0.004512335267231415im) rtol = 1e-9
    end

    @testset "hmunu vs captured DiscEMRIForce — odd l=2 m=1 (θ=π/3)" begin
        mo = _synth(Float64, 3)
        H0 = F.hmunu(mo, 1, 2, 1, π/3; a = 0.0)
        @test H0[1,3] ≈ (-0.0041846371885845566 - 0.009978750218932403im) rtol = 1e-9
        @test H0[3,4] ≈ (-0.020279395606217464 + 0.007725484040463794im) rtol = 1e-9
        @test H0[1,1] == 0.0 + 0.0im                          # odd parity: no h̄¹ content
        H5 = F.hmunu(mo, 1, 2, 1, π/3; a = 0.5)
        @test H5[1,3] ≈ (-0.004128633480384186 - 0.009989466628025537im) rtol = 1e-9
        @test H5[3,4] ≈ (-0.01873115805107339 + 0.01147153402597337im) rtol = 1e-9
        @test H5[1,1] ≈ (0.0001334118296429038 - 6.0092056627549e-5im) rtol = 1e-9
    end

    @testset "l<2 guard (no spin-2 content)" begin
        me = _synth(Float64, 7)
        p = F.tetmodes_from_bls(me, 1, 1; a = 0.0)            # l=1 → Λ̂=0 branch
        @test p.mpmp == 0.0 + 0.0im
        @test p.mmmm == 0.0 + 0.0im
    end

    @testset "canonical low-multipole contractions are finite" begin
        dir = mktempdir()
        for (l, m) in ((0, 0), (1, 0), (1, 1), (2, 1), (2, 2))
            fields = Data.h1_parity_fields(l, m)
            Data.write_mode_h5(joinpath(dir, "h1-l$(l)m$(m).h5"), _assemble_schema(length(fields)))
            mode = F.load_h1_mode(l, m; dir = dir)
            p = F.tetmodes_from_bls(mode, 1, l; a = 0.0)
            @test all(isfinite, values(p))
            if l < 2
                @test p.mpmp == 0.0 + 0.0im
                @test p.mmmm == 0.0 + 0.0im
            end
            if l == 0
                @test p.ρlpmp == 0.0 + 0.0im
                @test p.ρclpmm == 0.0 + 0.0im
                @test p.ρclmmp == 0.0 + 0.0im
                @test p.ρlmmm == 0.0 + 0.0im
            end
            H = F.hmunu(mode, 1, l, m, π/3; a = 0.0)
            @test all(isfinite, H)
        end
    end

    @testset "BigFloat inheritance (radial T, angular Float64 cap promoted)" begin
        meBF = _synth(BigFloat, 7)
        HBF = F.hmunu(meBF, 1, 2, 2, BigFloat(π)/2; a = big"0.0")
        @test HBF isa Matrix{Complex{BigFloat}}
        me = _synth(Float64, 7)
        H0 = F.hmunu(me, 1, 2, 2, π/2; a = 0.0)
        @test ComplexF64(HBF[3,3]) ≈ H0[3,3] rtol = 1e-9
    end

end
