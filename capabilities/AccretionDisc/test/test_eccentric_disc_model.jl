isdefined(@__MODULE__, :eccentric_legacy_sample) || include("eccentric_test_support.jl")
using AccretionDisc.SpiralWaves: DiscModel, Omega_ref

@testset "ECC disc model" begin
    @test_throws UndefKeywordError DiscModel(χ=0,r_ref=20,h_ref=0.05,evolution=:fixed)
    @test_throws UndefKeywordError DiscModel(χ=0,r_ref=20,h_ref=0.05,closure=:adiabatic)
    @test_throws ArgumentError DiscModel(χ=0,r_ref=20,h_ref=0.05,closure=:cooling,evolution=:fixed)
    @test_throws ArgumentError DiscModel(χ=0,r_ref=20,h_ref=0.05,closure=:adiabatic,evolution=:implicit)
    for T in (Float64,BigFloat), χ in (-0.7,0,0.6), closure in (:adiabatic,:locally_isothermal)
        bg = DiscBackground{T}(χ=χ,aEMRI=20,h0=0.06,𝔭=1.5,𝔮=0.75,ε0=2)
        d = DiscModel(bg;closure,evolution=:fixed)
        @test d isa DiscModel{T}
        @test d.r_ref == 20 && d.h_ref == bg.h0 && d.epsilon_ref == 2
        @test Omega_ref(d) === Omega_emri(bg) || Omega_ref(d) == Omega_emri(bg)
        @test F_J0(d;q=T(0.01)) == d.epsilon_ref*d.r_ref^4*Omega_ref(d)^2*T(0.01)^2/d.h_ref^3
        @test F_J1(d) == F_J0(d)/d.r_ref
        for r in T.((9,17,23,61))
            @test csound(d,r) == bg.h0*bg.aEMRI*Omega_emri(bg)*(r/bg.aEMRI)^(-bg.𝔮/2)
            @test eps_profile(d,r) == d.epsilon_ref*(r/d.r_ref)^(-d.pSigma)
            @test rotation_profile(d,r) == rotation_profile(bg,r)
            @test dOmega_dr(d,r) == dOmega_dr(bg,r)
            pot = kerr_equatorial_potentials(T(χ),r)
            dr = kerr_equatorial_potentials_dr(T(χ),r)
            v = rotation_profile(d,r).v
            Π = -(d.pSigma+d.qTemp)*pressure(d,r)/(r*(eps_profile(d,r)+pressure(d,r)))
            residual = Π+(dr.dν-v^2*dr.dψ+exp(pot.ψ-pot.ν)*v*dr.dω̃)/(1-v^2)
            @test abs(residual) < 100eps(T)
        end
        moved = DiscModel(bg;closure,evolution=:moving)
        @test moved.evolution === :moving && d.evolution === :fixed
        @test Omega(moved,T(15)) == Omega(d,T(15))
        # Lander review R1: the keyword constructor (the public entry point) maps every physical
        # parameter onto the same frozen evaluator as the explicit-background constructor.
        kw = DiscModel(;χ=T(χ),r_ref=20,h_ref=0.06,pSigma=1.5,qTemp=0.75,epsilon_ref=2,
                       closure,evolution=:fixed)
        @test kw isa DiscModel{T}
        @test all(getproperty(kw,n) == getproperty(d,n) for n in propertynames(d) if n !== :reference)
        @test pressure(kw,T(13)) == pressure(d,T(13)) && rotation_profile(kw,T(13)) == rotation_profile(d,T(13))
    end
end
println("[PASS] disc profiles, reference units, hydrostatic balance, explicit closure/evolution")
