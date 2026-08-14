# Grounding for Lorenz.SchwarzschildComponents — the a=0 ten-tetrad-mode → BLS h̄^i assembly.
#
# CHARACTERIZATION vs the legacy KerrLorenzCirc.to_BLS oracle (the modern port is a line-for-line
# consolidation of components.jl, so it must reproduce the legacy h̄^i to machine precision). The
# expected values are a ONE-SHOT fixture captured from the legacy package (corpus-free — build_recon/
# to_BLS need no external h1Lorenz HDF5 corpus):
#
#   cd packages/KerrLorenzCirc && julia --project=. -e '
#     using KerrLorenzCirc; const K=KerrLorenzCirc; r0=8.0; m=2
#     for l in (2,3); R=K.build_recon(r0,l;m=m)
#       for r in (r0+1.0,r0+2.0,r0+4.0); println(l," ",r," ",K.to_BLS(R,r)); end; end'
#
# The legacy package itself is validated against the independent h1Lorenz Lorenz-gauge C code at
# tol ~3.1e-3 (gate G-K8, ref_h1lorenz.jl). Here we ground the PORT against the legacy oracle output
# at the parity structure (even l: h̄1..h̄7 carry; odd l: h̄8..h̄10 carry) and machine tolerance.

using Test
import MetricReconstruction.Lorenz as L

# Captured legacy KerrLorenzCirc.to_BLS(build_recon(8.0,l;m=2), r) — fixture (see header for the command).
const REF_BLS = Dict(
    (2, 9.0)  => (2.509875668416624 - 0.2721935425094275im, -0.22739232566363576 + 1.028081119823104im,
                  4.027356016366492 + 0.5332843145734976im, -0.8425797398220529 + 7.331432226509658im,
                  -2.1419354693466106 - 2.3640124313170676im, 1.085541996500769 + 0.5170383244270548im,
                  -4.357827881226153 - 4.806646791794709im, 0.0+0im, 0.0+0im, 0.0+0im),
    (2, 10.0) => (2.0653884801597093 - 0.30186371403517576im, -0.2814411363682439 + 1.0676046569065194im,
                  3.3002060973302854 + 0.581132521577077im, -1.0199008069454705 + 6.499058307843765im,
                  -2.0001539151579473 - 2.607246391178892im, 0.9265689405510784 + 0.5586452858252525im,
                  -3.852296025796448 - 5.275132018805271im, 0.0+0im, 0.0+0im, 0.0+0im),
    (2, 12.0) => (1.539493414063801 - 0.343049758481299im, -0.3988277256843809 + 1.0738508872993844im,
                  2.366084158568481 + 0.6655142397093373im, -1.4009702532945727 + 5.397871970824054im,
                  -1.5502563566231315 - 3.00332414674522im, 0.6785077209945956 + 0.6265977687305763im,
                  -2.747028980156574 - 6.032113903686006im, 0.0+0im, 0.0+0im, 0.0+0im),
    (3, 9.0)  => (0.0+0im, 0.0+0im, 0.0+0im, 0.0+0im, 0.0+0im, 0.0+0im, 0.0+0im,
                  -4.115142042961451 - 0.0005206856512402074im, 0.0033805859111949604 - 0.3667706032848982im,
                  0.018512551068730218 - 4.857454217974467im),
    (3, 10.0) => (0.0+0im, 0.0+0im, 0.0+0im, 0.0+0im, 0.0+0im, 0.0+0im, 0.0+0im,
                  -2.9659509642700583 - 0.0007954103071233527im, 0.0047749549212570816 - 0.3489019852728182im,
                  0.02588613830067799 - 3.74142910193913im),
    (3, 12.0) => (0.0+0im, 0.0+0im, 0.0+0im, 0.0+0im, 0.0+0im, 0.0+0im, 0.0+0im,
                  -1.717614823599953 - 0.0016388826796135878im, 0.008508320017446879 - 0.29848537966132194im,
                  0.04547240794031405 - 2.490897710698977im),
)

@testset "Lorenz.SchwarzschildComponents (a=0 to_BLS vs legacy KLC oracle)" begin
    r0 = 8.0; m = 2
    for l in (2, 3)
        R = L.build_recon(r0, l; m=m)
        for r in (r0+1.0, r0+2.0, r0+4.0)
            bls = L.to_BLS(R, r)
            ref = REF_BLS[(l, r)]
            for i in 1:10
                if abs(ref[i]) > 1e-6
                    @test abs(bls[i] - ref[i]) / abs(ref[i]) < 1e-8     # machine-level port reproduction
                else
                    @test abs(bls[i]) < 1e-5                            # parity-zero component stays negligible
                end
            end
        end
    end

    @testset "tetrad_modes diagonal lplp/mpmp/h (l=2, r0+2) vs oracle" begin
        R = L.build_recon(8.0, 2; m=2)
        M = L.tetrad_modes(R, 10.0)
        @test abs(M.lplp - (0.2787417724674164 + 0.11964702232364743im)) / abs(0.2787417724674164 + 0.11964702232364743im) < 1e-8
        @test abs(M.mpmp - (-7.863466334460668 - 10.767818143242193im)) / abs(-7.863466334460668 - 10.767818143242193im) < 1e-8
        @test abs(M.h - (0.2373637156779207 + 0.0022487235751824375im)) / abs(0.2373637156779207 + 0.0022487235751824375im) < 1e-8
    end

    # F2 DEDUP GATE: the length-4 tuple jet kernel j4mul/j4add/j4scale/j4inv is now a set of thin ALIASES
    # (SchwarzschildComponents.jl) to the canonical NumericsCore jet algebra (NTuple jmul/jinv + untyped
    # jadd/jscale-on-tuples). Gate the aliases against the FORMER hand-unrolled bodies (VERBATIM oracle),
    # BIT-EXACT (isequal: ±0-distinguishing, value-based on the mutable BigFloat), Float64/BigFloat-256 ×
    # Real/Complex. (Deep coverage lives in core/NumericsCore/test/test_jets.jl; this guards the alias wiring.)
    @testset "F2 j4* aliases bit-exact to legacy length-4 tuple kernel" begin
        @inline _leg_j4mul(f, g) = (f[1]*g[1], f[1]*g[2]+f[2]*g[1],
                                    f[1]*g[3]+2f[2]*g[2]+f[3]*g[1],
                                    f[1]*g[4]+3f[2]*g[3]+3f[3]*g[2]+f[4]*g[1])
        _leg_j4add(f, g)   = (f[1]+g[1], f[2]+g[2], f[3]+g[3], f[4]+g[4])
        _leg_j4scale(s, f) = (s*f[1], s*f[2], s*f[3], s*f[4])
        function _leg_j4inv(f)
            g1 = inv(f[1]); g2 = -(1*f[2]*g1)*g1; g3 = -(2*f[2]*g2 + 1*f[3]*g1)*g1
            g4 = -(3*f[2]*g3 + 3*f[3]*g2 + 1*f[4]*g1)*g1; (g1, g2, g3, g4)
        end
        _beq(x::Real, y::Real) = isequal(x, y)
        _beq(x::Complex, y::Complex) = isequal(real(x), real(y)) && isequal(imag(x), imag(y))
        tbeq(a, b) = all(_beq(x, y) for (x, y) in zip(a, b))
        # dependency-free SplitMix64 (MR test target = Test/TestKit/Data; no Random) → reproducible tuples
        seed = Ref(0x0000000012345678)
        function nxt()
            seed[] += 0x9e3779b97f4a7c15; z = seed[]
            z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9; z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb; z ⊻ (z >> 31)
        end
        rr(::Type{Float64}) = (Float64(nxt() >> 11) / 9007199254740992.0 * 2 - 1) * 3
        function rr(::Type{BigFloat})
            hi = BigFloat(nxt() >> 11); lo = BigFloat(nxt() >> 11)
            (2 * ((hi + lo / (BigFloat(2)^53)) / (BigFloat(2)^53)) - 1) * 3
        end
        rel(::Type{T}) where {T<:Real} = rr(T)
        rel(::Type{Complex{T}}) where {T} = complex(rr(T), rr(T))
        rtup(::Type{T}) where {T} = ntuple(_ -> rel(T), 4)
        setprecision(BigFloat, 256) do
            for T in (Float64, ComplexF64, BigFloat, Complex{BigFloat})
                for _ in 1:500
                    f = rtup(T); g = rtup(T)
                    while abs(f[1]) < 0.1; f = rtup(T); end
                    s = rel(T)
                    @test tbeq(L.j4mul(f, g),   _leg_j4mul(f, g))
                    @test tbeq(L.j4add(f, g),   _leg_j4add(f, g))
                    @test tbeq(L.j4scale(s, f), _leg_j4scale(s, f))
                    @test tbeq(L.j4inv(f),      _leg_j4inv(f))
                    @test L.j4mul(f, g) isa NTuple{4}      # stays a tuple (not a Vector)
                end
            end
        end
        # captured PRE-CHANGE literals (independent anchor: fixture from the former inline kernel)
        ft = (ComplexF64(1.7,-0.4), ComplexF64(-0.3,0.9), ComplexF64(0.11,0.02), ComplexF64(-0.02,0.05))
        gt = (ComplexF64(0.9,0.2),  ComplexF64(0.4,-0.1), ComplexF64(-0.05,0.3), ComplexF64(0.007,-0.11))
        @test tbeq(L.j4mul(ft, gt), (1.61 - 0.020000000000000018im, 0.18999999999999995 + 0.41999999999999993im,
                                     0.07000000000000002 + 1.35im, -0.6871 - 0.5628im))
        @test tbeq(L.j4inv(ft), (0.5573770491803279 + 0.13114754098360656im, 0.21961838215533463 - 0.22026337006181143im,
                                 -0.10772055810838799 - 0.3506665315599105im, -0.6209016602169567 - 0.14050066355665944im))
    end
end
