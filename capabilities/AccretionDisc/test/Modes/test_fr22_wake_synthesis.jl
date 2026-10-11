using Test, AccretionDisc

@testset "FR22 wake: unit weight, interpolation, frames and mass" begin
    M, F = AccretionDisc.Modes, AccretionDisc.ForcingModes
    d = M.NewtonianDisc(h0=.1, p=.7, q=0., beta=0.)
    o = M.EccentricOrbit(1., .5, 1.)
    r = [.3,.55,.9,1.4,2.2]
    modes = [F.NewtonianMode(2,3,1.), F.NewtonianMode(2,1,1.)]
    a(r) = (1+2r-r^3) + im*(2-r+r^2)
    b(r) = (-2+r^2) + im*(1+3r)
    wake = M.FR22Wake(r,modes,hcat(a.(r),b.(r)))
    rq, ph, t = [.4,.8,1.2,2.], [-.4,.2,1.3,4.], .73
    expected = [d.h0^3/M.Sigma(d,x)*real(a(x)*cis(2y-3t)+b(x)*cis(2y-t)) for x in rq,y in ph]
    @test M.fr22_wake(wake,d,o,rq,ph,t;frame=:inertial) ≈ expected rtol=3e-14
    for (frame,shift) in ((:guiding,o.n*t),(:planet,M.planet_orbit(o.n*t,o).phi_p))
        @test M.fr22_wake(wake,d,o,rq,ph,t;frame) ≈
              M.fr22_wake(wake,d,o,rq,ph .+ shift,t;frame=:inertial) rtol=3e-14
    end
    thermal = M.FR22Wake(r,modes,d.h0^3 .* hcat(a.(r),b.(r)))
    @test M.fr22_wake(thermal,d,o,rq,ph,t;Gmu=d.h0^3) ≈ M.fr22_wake(wake,d,o,rq,ph,t)
    circular = M.EccentricOrbit(1.,0.,1.)
    diagonal = [F.NewtonianMode(1,1,1.),F.NewtonianMode(3,3,1.)]
    wc = M.FR22Wake(r,diagonal,hcat(a.(r),b.(r)))
    stationary = [d.h0^3/M.Sigma(d,x)*real(a(x)*cis(y)+b(x)*cis(3y)) for x in rq,y in ph]
    for phase in (0.,.3,2.,9.)
        @test M.fr22_wake(wc,d,circular,rq,ph,phase;frame=:planet) ≈ stationary rtol=3e-14
        @test M.fr22_wake(wc,d,circular,rq,ph,phase;frame=:guiding) ≈ stationary rtol=3e-14
    end
    # Total frequency and a non-unit mean motion; no accidental extra n in phase.
    n=1.7^(-1.5)
    os=M.EccentricOrbit(1.7,.2,n)
    ms=[F.NewtonianMode(2,3,n),F.NewtonianMode(2,1,n)]
    ws=M.FR22Wake(r,ms,hcat(a.(r),b.(r)))
    hp=d.h0*os.a^((1-d.q)/2)
    analytic=[hp^3/M.Sigma(d,x)*real(a(x)*cis(2y-3n*t)+b(x)*cis(2y-n*t)) for x in rq,y in ph]
    @test M.fr22_wake(ws,d,os,rq,ph,t;frame=:inertial) ≈ analytic rtol=3e-14
    # Lander test-only addition (L-G44): q != 0, so h_p = h0*a^((1-q)/2) is exercised beyond q = 0.
    dq=M.NewtonianDisc(h0=.1,p=.7,q=.5,beta=0.)
    hq=.1*1.7^.25
    analytic_q=[hq^3/M.Sigma(dq,x)*real(a(x)*cis(2y-3n*t)+b(x)*cis(2y-n*t)) for x in rq,y in ph]
    @test M.fr22_wake(ws,dq,os,rq,ph,t;frame=:inertial) ≈ analytic_q rtol=3e-14
    @test_throws ArgumentError M.fr22_wake(wake,d,o,[.1],ph,t)
    @test_throws ArgumentError M.fr22_wake(wake,d,o,rq,ph,t;frame=:bad)
    @test_throws DomainError M.fr22_wake(wake,d,o,rq,ph,t;Gmu=0)
    @test_throws ArgumentError M.FR22Wake(r,[modes[1],modes[1]],hcat(a.(r),b.(r)))
    @test_throws DimensionMismatch M.FR22Wake(r,modes,zeros(4,2))
end
