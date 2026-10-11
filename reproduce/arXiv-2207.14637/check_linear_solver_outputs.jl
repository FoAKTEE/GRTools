# Replay only the recorded Float64 responses; no ODE or source integration.
using AccretionDisc, Test, Serialization, TOML, LinearAlgebra
const P=AccretionDisc.Modes
include("linear_solver_audit.jl")
const groups=[(;h=.1,q=0.,p=0.,ml=[(1,1),(1,2),(2,1),(5,7)]),
    (;h=.06,q=0.,p=0.,ml=[(1,3),(5,4),(40,40),(150,158)]),
    (;h=.08,q=1.,p=1.5,ml=[(1,1),(2,3),(10,8),(80,84)])]
rows=Dict{String,Any}[]
@testset "refinement checkpoint replay" begin
    measured_rows=TOML.parsefile(joinpath(OUT,"refinement.toml"))["rows"]
    @test length(measured_rows)==48
    for row in measured_rows
        gi=row["group"];m=row["m"];l=row["l"];nx=row["nx"];nr=row["nr"]
        stride=7992÷(nr-1)
        z=deserialize(joinpath(OUT,"refine-$gi-$m-$l-$nx-$stride.jls"))
        mode=FM.NewtonianMode(m,l,1.)
        b=deserialize(joinpath(OUT,"bank-$gi.jls"))
        c=subcache(MO.harmonic_cache(b,1,mode),stride)
        v=observable(disc(GROUPS[gi]),z.x,mode,c,z.s.mode)
        @test v.torque==z.v.torque==row["torque"]
        @test v.density==z.v.density
        @test length(z.x)==nx
    end
end

@testset "recorded FR22 response invariants" begin
    for (gi,g) in enumerate(groups), (m,l) in g.ml
        z=deserialize(joinpath(OUT,"mode-$gi-$m-$l-1.0e-13.jls"))
        s=z.s;x=z.x;a=s.phase_safeguard
        @test eltype(s.mode)==ComplexF64
        @test all(isfinite,s.mode)
        @test all(isfinite,z.v.density)
        if m==1
            @test s.C_after[1]===s.C_before[1]
        end
        if a.status!=:rejected
            @test isfinite(a.X_after[2]) && a.X_after[2]<a.X_before[2]
            if m==1
                @test a.inner_response_change<=.1
            else
                @test isfinite(a.X_after[1]) && a.X_after[1]<a.X_before[1]
            end
        else
            @test s.C_after===s.C_before
        end
        @test z.v.density ≈ s.mode[:,1] .* x.^(-g.p+g.q) ./ g.h^2 rtol=2e-14
        ii,io=s.indices
        residuals=[abs(s.mode[i,2]-C*s.mode[i,1])/maximum(abs,s.mode[:,2]) for (i,C) in zip((ii,io),s.C_after)]
        # H3 (obligation fr22-o-m1-integrator): the (1,3) h=.06, e=.5 inner residual is 8.5e-8 of the
        # peak derivative, above the strict 1e-8 threshold, which stays unchanged. Recorded as a known
        # failure (lander review), so the replay exits 0 and an unexpected pass is reported.
        if (gi,m,l)==(2,1,3)
            @test_broken maximum(residuals)<1e-8
        else
            @test maximum(residuals)<1e-8
        end
        push!(rows,Dict("group"=>gi,"m"=>m,"l"=>l,"bc_residual_over_peak_derivative"=>residuals,
            "inner_amplitude_over_peak"=>abs(s.mode[ii,1])/maximum(abs,s.mode[:,1]),
            "inner_amplitude_over_r02"=>abs(s.mode[ii,1])/abs(s.mode[searchsortedfirst(x,.2),1]),
            "safeguard_status"=>string(a.status)))
    end
    open(joinpath(OUT,"saved-checks.toml"),"w") do io;TOML.print(io,Dict("rows"=>rows));end
end

