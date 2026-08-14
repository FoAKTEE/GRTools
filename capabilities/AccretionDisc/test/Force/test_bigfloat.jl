# Force BigFloat inheritance — the assembly/force/source arithmetic carries BigFloat precision
# (HARD REQUIREMENT, design §9). The a=0 corpus fixture is Float64-accurate (the oracle ceiling),
# but loading widened to T=BigFloat must produce BigFloat-typed metric/force/source, and the
# FD-FREE stages self-converge across precisions. GROUNDING (legacy ref_bigfloat.jl gates, frozen):
#   • gbf: hmunu_eq / source_F carry Complex{BigFloat} at the requested precision; the BigFloat
#     source_F agrees with the Float64 one to <1e-3 (the Float64-corpus accuracy ceiling).
#   • gbf_converge: the FD-free metric assembly (hmunu_eq) and the reconstruction δp=(ε+p)𝔥
#     identity at 256 vs 512 bits self-converge to >200 bits (precision-inheritance, NOT capped at
#     Float64's 53 bits). Fixtures only (no packages/* reach).
using Test
using AccretionDisc
using AccretionDisc.Force

const FIXDIR_BF = joinpath(@__DIR__, "..", "fixtures")

# representative FD-FREE calculation at the given precision: equatorial metric assembly hmunu_eq
# (l=2, m=2, r≈11) + the reconstruction pressure identity δp=(ε+p)𝔥. Pure precision-inheriting
# algebra ⇒ self-converges across precisions.
function _bf_repr_force(bits)
    setprecision(BigFloat, bits) do
        d = disc_a0(BigFloat); m = 2; r_p = BigFloat(8)
        modes = Dict(l => load_h1_mode(l, m; dir = FIXDIR_BF, T = BigFloat) for l in 2:4)
        g = modes[2].r_right; k = argmin(abs.(g .- BigFloat(11)))
        ω0 = m * disc_OmegaK(d, r_p)
        H = hmunu_eq(modes[2], k, 2, 2)
        hval = BigFloat(1)/1000 + (BigFloat(2)/1000)*im
        R = reconstruct(d, modes, k, m, r_p, ω0, hval,
                        -BigFloat(3)/10000 + (BigFloat(1)/10000)*im; lmax = 4)
        (H = H, δp = R.δp)
    end
end

@testset "Force BigFloat-native (a=0 fixtures widened)" begin

    @testset "type smoke: BigFloat carried through hmunu_eq & source_F" begin
        setprecision(BigFloat, 256) do
            d64 = disc_a0(Float64); dbf = disc_a0(BigFloat); m = 2
            mdbf = load_h1_mode(2, 2; dir = FIXDIR_BF, T = BigFloat)
            md64 = load_h1_mode(2, 2; dir = FIXDIR_BF, T = Float64)
            g = md64.r_right; k = argmin(abs.(g .- 11.0))
            Hbf = hmunu_eq(mdbf, k, 2, 2)
            @test eltype(Hbf) === Complex{BigFloat}
            @test precision(real(Hbf[1, 1])) == 256
            modes64 = Dict(l => load_h1_mode(l, m; dir = FIXDIR_BF, T = Float64)  for l in 2:4)
            modesbf = Dict(l => load_h1_mode(l, m; dir = FIXDIR_BF, T = BigFloat) for l in 2:4)
            F64 = source_F(d64, modes64, k, m, 8.0, m * disc_OmegaK(d64, 8.0); lmax = 4)
            Fbf = source_F(dbf, modesbf, k, m, BigFloat(8), m * disc_OmegaK(dbf, BigFloat(8)); lmax = 4)
            @test typeof(Fbf) === Complex{BigFloat}
            @test abs(Complex{Float64}(Fbf) - F64) / (abs(F64) + 1e-300) < 1e-3   # Float64-corpus ceiling
        end
    end

    @testset "self-convergence 256 vs 512 bits (FD-free assembly + δp identity)" begin
        r256 = _bf_repr_force(256); r512 = _bf_repr_force(512)
        @test eltype(r256.H) === Complex{BigFloat} && precision(real(r256.H[1, 1])) == 256
        @test eltype(r512.H) === Complex{BigFloat} && precision(real(r512.H[1, 1])) == 512
        @test typeof(r256.δp) === Complex{BigFloat} && typeof(r512.δp) === Complex{BigFloat}
        relH  = maximum(abs.(Complex{BigFloat}.(r256.H) .- r512.H)) / maximum(abs.(r512.H))
        reldp = abs(Complex{BigFloat}(r256.δp) - r512.δp) / abs(r512.δp)
        @test Float64(-log2(relH))  > 200                  # digits ~double, ≫ Float64's 53 bits
        @test Float64(-log2(reldp)) > 200
    end
end
