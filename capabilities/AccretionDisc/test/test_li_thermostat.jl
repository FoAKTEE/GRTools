using Test, AccretionDisc
import AccretionDisc.SpiralWaves as TSW
isdefined(TSW,:ThermostattedLI) || Base.include(TSW,joinpath(@__DIR__,"../src/SpiralWaves/li_thermostat.jl"))

@testset "thermostat distinct rest-density law and compression" begin
    tag=TSW.ThermostattedLI()
    for temperature in (0.0,0.001,0.03,0.3), rho in (0.5,2.0)
        st=TSW.thermostat_state(tag,rho,temperature)
        @test st.pressure == rho*temperature
        @test st.energy ≈ rho+st.pressure/(st.gamma_gas-1)
        @test st.energy isa Float64
        @test st.sound_speed2 ≈ st.theta/(1+st.theta)
        @test 0<=st.sound_speed2<1
        drho=0.02*rho
        dp=TSW.thermostat_lagrangian(tag,st,drho,0.0,0.0)
        @test dp.pressure ≈ temperature*drho
        @test dp.energy ≈ st.energy/rho*drho
        @test st.gamma1_eff==1
        if temperature>0
            @test st.gamma1_eff != 1+st.theta
            @test dp.pressure/drho ≈ st.pressure/rho
            @test !isapprox(dp.pressure/drho,(1+st.theta)*st.pressure/rho;rtol=1e-8)
        end
        xi=0.1; tr=-0.75temperature/7
        lag=TSW.thermostat_lagrangian(tag,st,drho,xi,tr)
        @test lag.pressure ≈ temperature*drho+rho*xi*tr
        @test lag.energy ≈ st.energy/rho*drho+rho*st.heat_capacity*xi*tr
    end
    for args in ((0.0,0.1),(-1.0,0.1),(1.0,-0.1),(Inf,0.1),(1.0,NaN))
        @test_throws DomainError TSW.thermostat_state(tag,args...)
    end
    @test_throws DomainError TSW.thermostat_state(tag,1.0,0.1;gamma_gas=1.0)
    @test_throws DomainError TSW.thermostat_state(tag,1.0,0.1;gamma_gas=3.0)
end

@testset "bath heat, force and Killing work" begin
    tag=TSW.ThermostattedLI(); st=TSW.thermostat_state(tag,2.0,0.03)
    # Nonzero azimuthal velocity makes reservoir angular momentum observable.
    u=(1.25,0.0,0.75,0.0); low=(-1.25,0.0,0.75,0.0)
    for theta in (-0.2,0.0,0.2), dt in (-0.01,0.0,0.01)
        b=TSW.thermostat_bath(tag,st,u,low;expansion=theta,temperature_advection=dt)
        @test b.heat ≈ st.pressure*theta+st.rho*st.heat_capacity*dt
        @test -sum(low[i]*b.force_up[i] for i in 1:4) ≈ b.heat
        @test all(abs(b.force_up[i]-b.heat*u[i])<=16eps(Float64) for i in 1:4)
        @test b.killing_energy ≈ 1.25b.heat
        @test b.killing_angular_momentum ≈ 0.75b.heat
        # Independent time derivative of epsilon(rho,T), with baryon continuity.
        dr=-st.rho*theta
        de=(1+st.heat_capacity*st.temperature)*dr+st.rho*st.heat_capacity*dt
        @test de+(st.energy+st.pressure)*theta ≈ b.heat atol=32eps(Float64)
    end
    @test TSW.thermostat_bath(tag,st,u,low;expansion=0.0,temperature_advection=0.0).heat==0
    @test_throws ArgumentError TSW.thermostat_bath(tag,st,u,u;expansion=0.0,temperature_advection=0.0)
    @test_throws DimensionMismatch TSW.thermostat_bath(tag,st,(1.0,),(-1.0,);expansion=0.0,temperature_advection=0.0)
end

@testset "sensitivity cannot select primary LI dispatch" begin
    tag=TSW.ThermostattedLI()
    @test !(tag isa AccretionDisc.ForcingModes.ClosureTag)
    @test TSW._closure_tag(:locally_isothermal) == AccretionDisc.ForcingModes.LocallyIsothermalOptionA
    for bad in (tag,:thermostatted_li,:li_thermostat,:thermostat)
        @test_throws ArgumentError TSW._closure_tag(bad)
        @test_throws ArgumentError TSW.DiscModel(χ=0.0,r_ref=20.0,h_ref=0.06,closure=bad,evolution=:fixed)
    end
    d=TSW.DiscModel(χ=0.0,r_ref=20.0,h_ref=0.06,closure=:locally_isothermal,evolution=:fixed)
    mode=AccretionDisc.ForcingModes.HarmonicMode(2,1,0.03)
    @test_throws ArgumentError TSW.ModeContext(d,mode;closure=tag)
    @test_throws MethodError TSW.thermostat_state(AccretionDisc.ForcingModes.LocallyIsothermalOptionA,1.0,0.1)
end
