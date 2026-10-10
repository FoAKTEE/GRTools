# Run in the existing derivation_fr environment, which contains AccretionDisc.
using AccretionDisc
using Test
include(joinpath(@__DIR__, "..", "src", "R3Derivation.jl"))
const EC = R3Derivation.EccentricConventions
const SW = AccretionDisc.SpiralWaves

@testset "C-NORM circular implementation crosscheck" begin
    for bits in (128,256)
        setprecision(BigFloat,bits) do
            for T in (Float64,BigFloat)
                for m in (1,3)
                    @test EC.reality_pair_weight(m,0,T) == SW._fourier_weight(T,m)
                end
                # Circular m=0 denotes a single static coefficient; an eccentric
                # k>0 axisymmetric reality pair has twice that weight.
                @test EC.reality_pair_weight(0,1,T) == 2SW._fourier_weight(T,0)
                bg = SW.DiscBackground(; χ=T(0),aEMRI=T(12),h0=T(1)/20,
                    𝔭=T(0),𝔮=T(0),ε0=T(1),b=T(3)/5)
                K = SW.newtonian_potential_coeffs(bg)
                @test K.Ktt == T(2)
                for m in (0,1,3)
                    potential = SW.newtonian_potential_mode(bg,m,T(16); N=128,mEMRI=T(1)/7)
                    metric = SW._singular_field_mode(K,m,T(16),128,T(1)/7)
                    @test metric.htt == -2potential.Φ
                    @test metric.htt_r == -2potential.Φ_r
                    @test metric.htt isa T
                    # For m>0, Phi_FR=2Phi_hat; m=0 here is static and has no
                    # independent dynamic pair in the circular calculation.
                    if m > 0
                        Phi_FR = 2potential.Φ
                        @test metric.htt == -Phi_FR
                    end
                end
                mode = AccretionDisc.ForcingModes.NewtonianMode(2,-1,T(7)/3)
                @test (mode.m,mode.k,mode.omega) == (2,-3,-T(7)/3)
            end
        end
    end
end
println("[PASS] C-NORM circular code: weights, Newtonian metric/potential amplitudes and radial derivatives")
