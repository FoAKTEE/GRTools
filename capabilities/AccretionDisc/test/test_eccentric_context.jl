isdefined(@__MODULE__, :eccentric_legacy_sample) || include("eccentric_test_support.jl")
using AccretionDisc.SpiralWaves: DiscModel, Omega_ref, EccentricPerturber, ModeContext,
    EccentricMasterInputs, reconstruct_fluid
import AccretionDisc.SpiralWaves as ESW

@testset "ECC arbitrary frequency context" begin
    for T in (Float64,BigFloat)
        d = DiscModel(χ=T(-0.4),r_ref=T(20),h_ref=T(0.06),pSigma=T(1.5),qTemp=T(0.75),
                      closure=:adiabatic,evolution=:fixed)
        for ω in (T(-0.1),zero(T),T(0.02),T(0.4))
            c = ModeContext(d,HarmonicMode(3,1,ω);closure=:adiabatic)
            @test omega_d(c,T(13)) == ω-3Omega(d,T(13))
            @test c.closure == Adiabatic
            if c.r_c === nothing
                @test landau_shift(c,T(13)) == 0
            else
                @test abs(omega_d(c,c.r_c)) < 100eps(T)
                @test landau_shift(c,c.r_c) > 0
            end
            inp = master_inputs(c,T(13);hmodes=eccentric_test_source(T(13)))
            @test inp isa EccentricMasterInputs{T}
            @test !iszero(imag(inp.htr)) && !iszero(imag(inp.hrphi))
            C,S = master_matrices(c,T(13);hmodes=eccentric_test_source(T(13)))
            @test eltype(C) == eltype(S) == Complex{T}
            @test all(isfinite,C) && all(isfinite,S)
            h = eccentric_test_source(T(13)); a = Complex{T}(2,3)
            ha = map(x -> a*x,h)
            C2,S2 = master_matrices(c,T(13);hmodes=ha)
            @test C2 == C
            @test S2 ≈ a*S rtol=10000eps(T)
            for field in (:htr,:hrphi)
                hz = merge(ESW._zero_source(c),NamedTuple{(field,)}((h[field],)))
                _,Sodd = master_matrices(c,T(13);hmodes=hz)
                @test maximum(abs,Sodd) > 0
            end
            f = reconstruct_fluid(c,T(13),(a,-a);hmodes=h)
            @test all(isfinite,values(f)) && f.dut isa Complex{T}
        end
        @test ModeContext(d,HarmonicMode(0,1,T(0.03));closure=:adiabatic).r_c === nothing
        @test_throws ArgumentError ModeContext(d,HarmonicMode(0,0,zero(T));closure=:adiabatic)
        @test_throws ArgumentError ModeContext(d,HarmonicMode(2,0,T(0.02));closure=:locally_isothermal)
        @test_throws UndefKeywordError ModeContext(d,HarmonicMode(2,0,T(0.02)))
        @test_throws ArgumentError ModeContext(d,HarmonicMode(2,0,T(0.02));closure=:adiabatic,domain=(2,10))
        orbit = EquatorialOrbit(T(-0.4),T(20),T(0.1))
        p = EccentricPerturber(orbit;mass=T(1e-5),softening=FixedSoftening(T(0.3)))
        c = ModeContext(d,p,3,-1;closure=:adiabatic)
        @test c.mode.omega == orbit_frequency(orbit,3,-1)
        @test_throws ArgumentError ModeContext(d,p,HarmonicMode(3,-1,T(0.1));closure=:adiabatic)
    end
    d = DiscModel(χ=0,r_ref=20,h_ref=0.05,closure=:adiabatic,evolution=:fixed)
    @test_throws ArgumentError ModeContext(d,HarmonicMode(1,0,big"0.02");closure=:adiabatic)
    c = ModeContext(d,HarmonicMode(2,0,2Omega(d,20.0));closure=:adiabatic)
    @test_throws DomainError master_inputs(c,20.0;hmodes=ESW._zero_source(c),landau=false)
    @test master_inputs(c,20.0;hmodes=ESW._zero_source(c)).omega_d != 0
    @test !isdefined(ESW,:Symbolics)
    @test all(m -> nameof(m) != :Symbolics,values(Base.loaded_modules))
end
println("[PASS] frequency/type/closure validation, complex htr/hrphi, root/no-root, runtime isolation")

@testset "ECC e to zero through emitted eccentric operator" begin
    for χ in (-0.6,0.0,0.6), q in (0.0,0.75), m in (1,3)
        bg = DiscBackground(χ=χ,aEMRI=20,h0=0.06,𝔭=1.5,𝔮=q)
        d = DiscModel(bg;closure=:adiabatic,evolution=:fixed)
        orbit = EquatorialOrbit(χ,20.0,0.0)
        p = EccentricPerturber(orbit;softening=FixedSoftening(0.3))
        c = ModeContext(d,p,m,0;closure=:adiabatic)
        for r in (10.0,17.0,23.0,45.0)
            h = singular_field_mode(bg,m,r;N=128)
            he = merge(h,(htr=0im,hrphi=0im))
            A,b = master_matrices(bg,m,r;hmodes=h,landau=false)
            E,f = master_matrices(c,r;hmodes=he,landau=false)
            @test E ≈ A rtol=2e-11
            @test f ≈ b rtol=2e-11
            y = (0.03+0.02im,-0.01+0.04im)
            @test delta_u_phi(c,r,y;hmodes=he,landau=false) ≈
                  delta_u_phi(bg,m,r,y;hmodes=h,landau=false) rtol=2e-11
        end
    end
    d = DiscModel(χ=0.4,r_ref=20,h_ref=0.06,closure=:adiabatic,evolution=:fixed)
    function coefficients(e)
        p = EccentricPerturber(EquatorialOrbit(0.4,20.0,e);softening=FixedSoftening(0.3))
        c = ModeContext(d,p,3,0;closure=:adiabatic)
        first(master_matrices(c,13.0;hmodes=ESW._zero_source(c),landau=false))
    end
    A0,A1,A2 = coefficients(0.0),coefficients(1e-3),coefficients(1e-4)
    @test maximum(abs,A1-A0)/maximum(abs,A2-A0) > 80
end
println("[PASS] e→0: eccentric generated operator/source/reconstruction parity and O(e²) approach")
