# Run under the Spacetime project: source algebra has no Symbolics dependency.
using Spacetime, LinearAlgebra, Test
include(joinpath(@__DIR__,"..","src","eccentric_source.jl"))
const ES = EccentricSource
const GEO = Spacetime

@testset "C-SOURCE Spacetime rates, charts, odd components and remainder" begin
    for T in (Float64,BigFloat)
        tol = 2000eps(T)
        flat = Matrix(Diagonal(T[-1,1,1,1]))
        flat[3,3] = T(7)^2
        newt = ES.field_jet(flat,T[1,0,0,0],T(7),T(9),T(1)/3,T(1)/4,T(1)/5)
        plummer = -T(1)/5/sqrt(T(1)/16+T(9)^2+T(7)^2-2T(9)*T(7)*cos(T(1)/3))
        @test isapprox(newt.h[1,1],-2plummer;rtol=tol)
        for spin in (-T(3)/5,zero(T),T(3)/5)
            circular = GEO.EquatorialOrbit(spin,T(20),zero(T))
            sc = GEO.worldline_rates(circular,T(1))
            gcbl = Spacetime.kerr_metric(spin,sc.r,T(pi)/2)
            gc = Matrix(gcbl[[1,2,4,3],[1,2,4,3]])
            gc[4,4] /= sc.r^2
            hc = ES.field_jet(gc,T[sc.ut,sc.ur,sc.uphi,0],sc.r,T(21),T(1)/3,T(1)/5,T(1)/13)
            @test iszero(hc.h[1,2]) && iszero(hc.h[2,3])
            o = GEO.EquatorialOrbit(spin,T(20),T(1)/5)
            # Consume the canonical owner at its documented positive-spin chart;
            # source construction needs rates only, never element derivatives.
            tp = GEO.turning_point_jacobian(abs(spin),o.p,o.e,spin<0 ? -one(T) : one(T))
            @test all(isfinite,tp)
            for q in (T(1),2T(pi)-T(1))
                s = GEO.worldline_rates(o,q)
                R = s.r
                gbl = Spacetime.kerr_metric(spin,R,T(pi)/2)
                g = Matrix(gbl[[1,2,4,3],[1,2,4,3]])
                g[4,4] /= R^2
                u = T[s.ut,s.ur,s.uphi,0]
                ul = g*u
                @test isapprox(dot(u,ul),-one(T);atol=tol)
                @test isapprox(ul[1],-o.E;rtol=tol)
                @test isapprox(ul[3],o.L;rtol=tol)
                @test sign(s.ur) == (q<T(pi) ? 1 : -1)
                phase = T(2)/7
                r, psi, soft, mass = R+T(1)/3,T(1)/11,T(1)/5,T(1)/13
                f = ES.field_jet(g,u,R,r,psi,soft,mass)
                @test f.s2 >= soft^2
                @test f.h[1,2] != 0 && f.h[2,3] != 0
                @test isapprox(f.h[1,2],-4mass*o.E*ul[2]/sqrt(f.s2);rtol=tol)
                @test isapprox(f.h[2,3],4mass*o.L*ul[2]/sqrt(f.s2);rtol=tol)
                @test isapprox(sum(inv(g).*f.h),4mass/sqrt(f.s2);rtol=tol)
                # Full Cartesian transform at arbitrary worldline azimuth,
                # with nonzero Delta t, then restrict time to zero.
                dx,dy = r*cos(phase+psi)-R*cos(phase),r*sin(phase+psi)-R*sin(phase)
                P = ES.chart_projector(g,u,R,phase)
                delta = T[0,dx,dy,0]
                @test isapprox(dot(delta,P*delta)+soft^2,f.s2;rtol=tol)
                delta[1] = T(1)/17
                direct = ES.chart_jacobian(R,phase)*delta
                @test isapprox(dot(delta,P*delta),dot(direct,g*direct)+dot(ul,direct)^2;rtol=tol)
                # Radial derivatives tested numerically in addition to the
                # symbolic chain rule; independent symmetric field evaluations.
                step = eps(T)^(one(T)/4)
                fp = ES.field_jet(g,u,R,r+step,psi,soft,mass)
                fm = ES.field_jet(g,u,R,r-step,psi,soft,mass)
                @test isapprox((fp.h-fm.h)/(2step),f.h_r;rtol=100sqrt(eps(T)))
                @test isapprox((fp.h_r-fm.h_r)/(2step),f.h_rr;rtol=100sqrt(eps(T)))
                twist = -ul[2]*o.L/(g[3,3]+o.L^2)
                alpha20 = g[2,2]+ul[2]^2-(g[3,3]+o.L^2)*twist^2
                beta = 4*(g[3,3]+o.L^2)
                B = ES.chart_remainder_bound(g,u,R,twist)
                for scale in (T(1)/20,T(1)/40,T(1)/80), radial in (-1,1), azimuth in (-1,1)
                    dr, dp = radial*scale,azimuth*scale
                    cart = ES.distance_jet(g,u,R,R+dr,dp,zero(T)).s2
                    ref = alpha20*dr^2+beta*sin((dp-twist*dr)/2)^2
                    @test abs(cart-ref) <= B*scale^3
                end
                # Complex Fourier weighting must preserve the odd components.
                weighted = f.h*cis(T(2)/3)
                @test !iszero(imag(weighted[1,2])) && !iszero(imag(weighted[2,3]))
                @test real(weighted[1,2]) isa T
            end
        end
    end
end
println("[PASS] C-SOURCE runtime: signed Spacetime rates, trace, h_tr/h_rphi, rotated chart, derivatives, cubic bound, Float64/BigFloat")
