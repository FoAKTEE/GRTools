using Test, AccretionDisc

@testset "Sole-owner softening policies: FR22 / FR24 / FD25" begin
    M, F = AccretionDisc.Modes, AccretionDisc.ForcingModes
    h, b, a = .06, .3, 1.
    fixed, moving, localH = F.FixedSoftening(b*h*a), F.InstantaneousScaleHeight(b), F.LocalScaleHeight(b)
    for q in (0.,.5,1.), e in (0.,.12,.5), anomaly in (0.,pi/2,pi)
        orbit = M.EccentricOrbit(a,e,1.)
        rp = M.planet_orbit(anomaly,orbit).r_p
        power = (3-q)/2
        @test M._bank_softlength(fixed,h,power,a,rp,.7) == b*h*a
        @test M._bank_softlength(moving,h,power,a,rp,.7) ≈ b*h*a*(rp/a)^power
        @test M._bank_softlength(localH,h,power,a,rp,.7) ≈ b*h*a*(.7/a)^power
    end
    # FD25 2_linear_theory.tex:189 and 5_discussion.tex:40-51, transcribed by hand.
    xs = sqrt(2e-6/(3fixed.eps))
    @test xs ≈ .006085806194501846
    @test 8pi/(3xs)/(2pi) ≈ 219.08902300206645
    @test (1-.5)^1.5 ≈ .3535533905932738
    @test (1+.5)^1.5 ≈ 1.8371173070873836
    d=M.NewtonianDisc(h0=h,p=0.,q=0.,beta=0.,soft=b*h)
    modes=[F.NewtonianMode(2,2,1.),F.NewtonianMode(2,3,1.)]
    r=[.8,1.,1.2]
    circular=M.EccentricOrbit(1.,0.,1.)
    c=M.potential_bank(d,[circular],r,modes;softening=fixed,Nt=256,Nphi=256)
    s=M.potential_bank(d,[circular],r,modes;softening=moving,Nt=256,Nphi=256)
    @test maximum(abs.(c.caches[1,1].P-s.caches[1,1].P)) < 1e-13
    @test maximum(abs.(c.caches[1,2].P)) < 1e-13
    ecc=M.EccentricOrbit(1.,.12,1.)
    c=M.potential_bank(d,[ecc],r,modes;softening=fixed,Nt=256,Nphi=256)
    s=M.potential_bank(d,[ecc],r,modes;softening=moving,Nt=256,Nphi=256)
    @test maximum(abs.(c.caches[1,2].P-s.caches[1,2].P)) > .001
    # Lander test-only addition (L-G44): the bank's "fixed" FR22 reading is SemimajorScaleHeight, b*H(a) =
    # b*h0*a^((3-q)/2) at every time; check it against the constant length on an eccentric orbit, a != 1, q != 0.
    dq=M.NewtonianDisc(h0=h,p=0.,q=.5,beta=0.,soft=b*h)
    nq=1.3^-1.5
    oq=M.EccentricOrbit(1.3,.3,nq)
    mq=[F.NewtonianMode(2,2,nq),F.NewtonianMode(2,3,nq)]
    semimajor=M.potential_bank(dq,[oq],r,mq;softening=F.SemimajorScaleHeight(b),Nt=256,Nphi=256)
    constant=M.potential_bank(dq,[oq],r,mq;softening=F.FixedSoftening(b*h*1.3^1.25),Nt=256,Nphi=256)
    @test all(j -> maximum(abs.(semimajor.caches[1,j].P-constant.caches[1,j].P)) <= 1e-14, 1:2)
end
