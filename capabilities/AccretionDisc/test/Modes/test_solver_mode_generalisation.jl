using Test
using AccretionDisc
using AccretionDisc.ForcingModes: NewtonianMode, HarmonicMode

@testset "frequency-aware Newtonian solver" begin
    P = AccretionDisc.Modes
    bits(a, b) = reinterpret(UInt64, vec(a)) == reinterpret(UInt64, vec(b))
    synthetic(d, x) = PotentialCache(log10.(x), x.^2, 2 .* x, fill(eltype(x)(2), length(x)))

    @testset "circular mode and Phi_m cache: bit identity" begin
        d = NewtonianDisc(beta=80., lap_nphi=128)
        x = loggrid(d, .3, 3., 120)
        for m in (1, 2, 5)
            mode = NewtonianMode(m,m,1.)
            cache = PotentialCache(d,m,.3,3.; n=320)
            for fun in (P.Dfreq,P.dDdr,P.kwave,P.wkbC)
                @test bits([fun(d,r,m) for r in x], [fun(d,r,mode) for r in x])
            end
            for extra in ((), (cache,))
                a = solve_m_mode(d,m,x,extra...; rtol=1e-8,atol=1e-10)
                b = solve_m_mode(d,mode,x,extra...; rtol=1e-8,atol=1e-10)
                @test bits(a,b)
                for (u,v) in zip(vars_from_enthalpy(d,a,x,m,extra...),
                                 vars_from_enthalpy(d,b,x,mode,extra...))
                    @test bits(u,v)
                end
            end
            @test bits(solve_m_mode(d,m,x,cache),P.solve_lm_mode(d,mode,x,cache))
            @test P.corotation_index(d,x,m) == P.corotation_index(d,x,mode)
            @test bits(solve_m_mode(d,m,x,cache; phgrad_steps=1),
                       P.solve_lm_mode(d,mode,x,cache; phgrad_steps=1))
        end
    end

    @testset "off-diagonal corotation and Lindblad roots" begin
        d = NewtonianDisc(p=0.,q=0.,beta=0.)
        for (m,l,a) in ((2,1,1.),(2,3,1.),(3,5,1.7),(5,2,.8),(1,2,1.))
            n = a^(-3/2); mode = NewtonianMode(m,l,n)
            rc = (m/l)^(2/3)*a
            rin = ((m-1)/l)^(2/3)*a
            rout = ((m+1)/l)^(2/3)*a
            # Include the analytic root exactly so the gate is independent of mesh error.
            x = sort!(unique!([collect(loggrid(d,.15a,4a,201)); rc]))
            ic = P.corotation_index(d,x,mode)
            @test x[ic] == rc
            @test abs(m*P.Omega(d,rc)-l*n) < 32eps()*max(1.,l*n)
            for r in (m == 1 ? (rout,) : (rin,rout))
                @test abs(P.Dfreq(d,r,mode)) < 64eps()*P.kappa(d,r)^2
                dr = 1e-6*r
                @test P.Dfreq(d,r-dr,mode)*P.Dfreq(d,r+dr,mode) < 0
                @test P.dDdr(d,r,mode) ≈ (P.Dfreq(d,r+dr,mode)-P.Dfreq(d,r-dr,mode))/(2dr) rtol=1e-8
            end
            c = synthetic(d,x)
            bc = (.02+.03im,.04-.01im)
            for extra in ((),(c,))
                f = P.integrate_master(d,x,bc,mode,0,extra...; rtol=1e-8,atol=1e-10)
                @test f[ic,:] == collect(bc)
            end
            @test P._mode_anchor(d,x,mode,:auto) == ic
            @test P._mode_anchor(d,x,mode,:corotation) == ic
            println("RESONANCE m=$m l=$l a=$a rc=$(x[ic]) rminus=$rin rplus=$rout")
        end
    end

    @testset "FR22 (16), (19)-(21): arbitrary q, l != m" begin
        # Hand transcription of ref-paper/arXiv-2207.14637/src/Sections/
        # 2_framework.tex:140-145,165-169. The independent oracle constructs
        # omega, D and their derivatives without the frequency-aware helpers.
        for T in (Float64,BigFloat)
            setprecision(192) do
                worst = zero(T)
                for q in (0//1,3//4,3//2), (m,l) in ((2,1),(2,3),(1,2)), r0 in (3//5,23//10)
                    r=T(r0); n=T(4)/5; mode=NewtonianMode(m,l,n)
                    d=NewtonianDisc{T}(h0=T(1)/10,p=T(1)/2,q=T(q),gamma=T(7)/5,
                        beta=zero(T),pole_disp=zero(T),lap_nphi=64)
                    x=T[r/2,r,2r,3r]; c=synthetic(d,x)
                    h=T(2)/7+im*T(3)/11; h1=-T(5)/13+im*T(7)/17
                    y=T[real(h),imag(h),real(h1),imag(h1)]
                    O=P.Omega(d,r); O1=P.dOmegadr(d,r); K=P.kappa(d,r)^2
                    w=l*n-m*O; D=K-w^2
                    D1=2P.kappa(d,r)*P.dkappadr(d,r)+2m*w*O1
                    S=P.Sigma(d,r); cs2=P.csiso(d,r)^2; L=-d.q/r
                    A=(1-d.p)/r-D1/D; B=-d.p/r+O1/O-D1/D
                    C1=A-L
                    C0=-(2m*O/(r*w)*(L+B)+L*A+d.q/r^2+m^2/r^2+D/cs2)
                    for extra in ((),(c,))
                        phi,phi1,phi2=isempty(extra) ? (P.Phi_m(d,r,m),P.dPhidr_m(d,r,m),P.d2Phidr2_m(d,r,m)) : P.potential_at(c,r)
                        psi=-phi2-A*phi1+(2m*O/(r*w)*B+m^2/r^2)*phi
                        rhs=P.master_rhs(d,r,y,mode,1,extra...)
                        expected=psi-C1*h1-C0*h
                        err=abs(complex(rhs[3],rhs[4])-expected)/max(one(T),abs(expected))
                        worst=max(worst,err)
                        @test err <= 1024eps(T)
                        g=isempty(extra) ? P.MasterRHS(d,mode,1) : P.CachedMasterRHS(d,mode,1,c)
                        out=zero(y); g(out,r,y)
                        @test out == collect(rhs) == g(r,y)
                        ur=im/D*(w*(h1+phi1)-2m*O/r*(h+phi)-L*w*h)
                        up=(K/(2O)*(h1+phi1)-m*w/r*(h+phi)-L*K/(2O)*h)/D
                        sig=S*h/cs2
                        got=vars_from_enthalpy(d,reshape([h,h1],1,2),[r],mode,extra...)
                        for (v,oracle) in zip(got,(ur,up,sig))
                            @test eltype(v)==Complex{T}
                            @test abs(v[1]-oracle) <= 1024eps(T)*max(one(T),abs(oracle))
                        end
                    end
                end
                println("ASSEMBLY type=$T max_scaled_error=$worst")
            end
        end
    end

    @testset "pole displacement and cooling assembly" begin
        for beta in (0.,.1,80.), (m,l) in ((1,2),(3,2))
            d=NewtonianDisc(q=.75,p=.5,beta=beta,lap_nphi=64)
            r=.7; mode=NewtonianMode(m,l,.8); y=[.2,.3,-.4,.7]
            x=[.4,.7,1.2,3.]; c=synthetic(d,x)
            O=P.Omega(d,r); w=l*.8-m*O; K=P.kappa(d,r)^2; D=K-w^2
            D1=2P.kappa(d,r)*P.dkappadr(d,r)+2m*w*P.dOmegadr(d,r)
            phi,phi1,phi2=P.potential_at(c,r)
            coeff=P.master_equation_cool(r,d.gamma,P.csadi(d,r),O,P.dOmegadr(d,r),
                P.Sigma(d,r),P.dSigmadr(d,r),P.LTinv(d,r),P.dLTinvdr(d,r),
                P.LSinv(d,r),P.dLSinvdr(d,r),w,beta/O,-beta/O^2*P.dOmegadr(d,r),
                D,D1,P.Nr2(d,r),P.dNr2dr(d,r),d.pole_disp,phi,phi1,phi2,m)
            c1=complex(coeff[1],coeff[2]); c0=complex(coeff[3],coeff[4]); psi=complex(coeff[5],coeff[6])
            expected=psi-c1*complex(y[3],y[4])-c0*complex(y[1],y[2])
            actual=P.master_rhs(d,r,y,mode,1,c)
            @test complex(actual[3],actual[4]) ≈ expected rtol=1e-13
        end
    end

    @testset "m=1 inner evanescence and phase outer-only update" begin
        d=NewtonianDisc(p=0.,q=0.,beta=0.)
        for l in (1,2,3)
            mode=NewtonianMode(1,l,1.)
            x=loggrid(d,.1,3.,600)
            ci,co=P.wkb_boundary_coeffs(d,x,1,length(x),mode)
            r=x[1]; O=r^(-3/2); w=l-O; D=r^(-3)-w^2
            D1=-3r^(-4)-3w*r^(-5/2)
            # The decaying inward branch has +sqrt(D)/cs in h'/h.
            @test D > 0
            expected=sqrt(D)/d.h0+D1/(4D)-1/(2r)
            @test ci ≈ expected rtol=2e-13
            @test co == P.wkbC(d,x[end],mode)
            c=synthetic(d,x)
            sol=P.solve_lm_mode(d,mode,x,c; diagnostics=true)
            ii,io=sol.indices
            @test abs(sol.mode[ii,2]-sol.C_after[1]*sol.mode[ii,1]) < 1e-7*max(1,abs(sol.mode[ii,2]))
            @test sol.C_after[1] ≈ P.wkb_boundary_coeffs(d,x,ii,io,mode)[1]
            phase=P.solve_lm_mode(d,mode,x,c; phgrad_steps=1,reint_steps=0,diagnostics=true)
            @test phase.C_after[1] == phase.C_before[1]
            @test all(isfinite,phase.mode)
            println("M1 l=$l inner_coefficient=$ci residual=$(abs(sol.mode[ii,2]-sol.C_after[1]*sol.mode[ii,1]))")
        end
    end

    @testset "potential bank consumption, domains and inherited precision" begin
        d=NewtonianDisc(beta=0.,p=.5,q=.75,lap_nphi=64)
        mode=NewtonianMode(2,3,1.)
        x=loggrid(d,.3,3.,80)
        bank=P.potential_bank(d,[P.EccentricOrbit(1.,.1,1.)],x,[mode]; Nt=16,Nphi=32)
        cache=P.harmonic_cache(bank,1,mode)
        f=P.solve_lm_mode(d,mode,x,cache; reint_steps=0)
        @test all(isfinite,f)
        @test maximum(abs,f)>0
        doubled=PotentialCache(cache.lg,2cache.P,2cache.dP,2cache.d2P)
        # Reconstruction must consume the supplied cache, never silently Phi_m.
        z=zero(f)
        @test vars_from_enthalpy(d,z,x,mode,doubled)[1] ≈ 2vars_from_enthalpy(d,z,x,mode,cache)[1]
        @test_throws ArgumentError P.solve_lm_mode(d,NewtonianMode(0,1,1.),x,cache)
        @test_throws ArgumentError P.solve_lm_mode(d,NewtonianMode(2,0,1.),x,cache)
        @test_throws ArgumentError P.solve_lm_mode(d,NewtonianMode(2,-1,1.),x,cache)
        @test_throws ArgumentError P.solve_lm_mode(d,NewtonianMode(2,100,1.),x,cache)
        @test_throws ArgumentError P.solve_lm_mode(d,mode,reverse(x),cache)
        @test_throws ArgumentError P.solve_lm_mode(d,mode,x,cache; anchor=:unknown)
        short=PotentialCache(cache.lg[2:end],cache.P[2:end],cache.dP[2:end],cache.d2P[2:end])
        @test_throws ArgumentError P.solve_lm_mode(d,mode,x,short)
        for prec in (128,192)
            setprecision(prec) do
                db=NewtonianDisc{BigFloat}(h0=big(1)/10,p=big(0),q=big(0),beta=big(0))
                mb=NewtonianMode(2,3,big(1)); xb=loggrid(db,big(7)/10,big(9)/10,12)
                cb=synthetic(db,xb)
                fb=P.integrate_master(db,xb,(big(.1)+0im,big(.2)+0im),mb,1,cb; rtol=big(1e-8),atol=big(1e-10))
                @test eltype(fb)==Complex{BigFloat}
                @test precision(real(fb[1,1]))==prec
                @test precision(real(P.wkbC(db,xb[1],mb)))==prec
                @test precision(BigFloat)==prec
            end
        end
    end
end

@testset "direct circular forcing refuses non-circular harmonic solves (lander review)" begin
    # The direct path forces with the circular Laplace Phi_m. A HarmonicMode is a circular-planet
    # mode there only for k == 0 and omega == m; any other mode must bring its own potential.
    d = NewtonianDisc(beta=80., lap_nphi=64)
    x = loggrid(d, .3, 3., 40)
    @test_throws ArgumentError solve_m_mode(d, NewtonianMode(2,3,1.), x)
    @test_throws ArgumentError solve_m_mode(d, NewtonianMode(2,2,.9), x)
    @test_throws ArgumentError solve_m_mode(d, HarmonicMode(2,1,2.), x)  # omega == m, but k != 0
    cache = PotentialCache(d, 2, .3, 3.; n=80)
    @test all(isfinite, solve_m_mode(d, NewtonianMode(2,3,1.), x, cache))  # the cached path stays open
end

@testset "phase-step descent safeguard is forwarded by solve_m_mode (lander review)" begin
    # Hand-off of ECC::fr22-phase-gradient: phgrad_safeguard reaches phase_gradient_refine through the
    # production solver. MR beta = 80, m = 1 (cached) is the case whose unguarded step ascends at both
    # edges; its adaptive boundaries stay at the grid ends.
    bits(a, b) = reinterpret(UInt64, vec(a)) == reinterpret(UInt64, vec(b))
    P = AccretionDisc.Modes
    d = NewtonianDisc(beta=80., lap_nphi=800)
    x = loggrid(d, .05, 5., 2000)
    n = length(x)
    cache = PotentialCache(d, 1, .05, 5.; n)
    plain = solve_m_mode(d, 1, x, cache; phgrad_steps=1, diagnostics=true)
    off = solve_m_mode(d, 1, x, cache; phgrad_steps=1, diagnostics=true, phgrad_safeguard=false)
    @test !haskey(plain, :phase_safeguard) && !haskey(off, :phase_safeguard)
    @test bits(plain.mode, off.mode) && plain.C_after === off.C_after
    @test all(plain.X_after .> plain.X_before)  # the unguarded MR step ascends here
    on = solve_m_mode(d, 1, x, cache; phgrad_steps=1, diagnostics=true, phgrad_safeguard=true)
    a = on.phase_safeguard
    @test on.indices == plain.indices == (1, n)
    @test a.status == :backtracked && a.backtracks >= 1 && a.fraction == 2.0^-a.backtracks
    @test a.X_before == on.X_before && a.X_after == on.X_after
    @test all(on.X_after .< on.X_before)
    @test on.C_before === plain.C_before
    @test on.C_after[1] === on.C_before[1]  # m = 1: only C_out moves
    @test on.C_after[2] == on.C_before[2] + a.fraction * (plain.C_after[2] - on.C_before[2])
    @test bits(solve_m_mode(d, 1, x, cache; phgrad_steps=1, phgrad_safeguard=true), on.mode)
    rejected = solve_m_mode(d, 1, x, cache; phgrad_steps=1, diagnostics=true,
                            phgrad_safeguard=true, phgrad_max_backtracks=0)
    @test rejected.phase_safeguard.status == :rejected
    @test rejected.C_after === rejected.C_before
    @test bits(rejected.mode, solve_m_mode(d, 1, x, cache))  # the unrefined default solution
    # MR beta = 80, m = 10 (cached): the full step descends at both edges and is accepted unchanged.
    c10 = PotentialCache(d, 10, .05, 5.; n)
    full = solve_m_mode(d, 10, x, c10; phgrad_steps=1, diagnostics=true)
    acc = solve_m_mode(d, 10, x, c10; phgrad_steps=1, diagnostics=true, phgrad_safeguard=true)
    @test acc.indices == (1, n)
    @test acc.phase_safeguard.status == :accepted && acc.phase_safeguard.backtracks == 0
    @test acc.phase_safeguard.fraction == 1.0
    @test bits(acc.mode, full.mode) && acc.C_after === full.C_after
    @test all(acc.X_after .< acc.X_before)
    # The safeguard is integer-only in PhaseGradient.jl: a harmonic request must refuse, never skip.
    @test_throws ArgumentError P.solve_lm_mode(d, NewtonianMode(1,1,1.), x, cache;
                                               phgrad_steps=1, phgrad_safeguard=true)
end

@testset "adaptive-boundary anchor rule (lander review)" begin
    # :auto keeps the legacy r = 1 anchor for every circular mode (integer m, or k == 0 with omega == m),
    # even where corotation of the pressure-supported disc is a different grid point, and uses corotation
    # for any other frequency. MR beta = 80 (h0 = 0.1, p = q = 1): Omega(1) < 1, corotation at r ~ 0.993.
    P = AccretionDisc.Modes
    d = NewtonianDisc(beta=80.)
    x = loggrid(d, .05, 5., 2000)
    i1 = argmin(abs.(x .- 1))
    for m in (1, 2, 10)
        ic = P.corotation_index(d, x, m)
        @test ic != i1
        @test P._mode_anchor(d, x, m, :auto) == i1
        @test P._mode_anchor(d, x, NewtonianMode(m, m, 1.), :auto) == i1
        @test P._mode_anchor(d, x, m, :corotation) == ic
        @test P._mode_anchor(d, x, NewtonianMode(m, m, 1.), :legacy) == i1
    end
    mode = NewtonianMode(2, 3, 1.)
    @test P._mode_anchor(d, x, mode, :auto) == P.corotation_index(d, x, mode) != i1
    @test P._mode_anchor(d, x, mode, :legacy) == i1
end

@testset "cooled WKB wavenumber and outgoing coefficient use the mode frequency (lander review)" begin
    # With beta > 0 the cooling terms of kwave/wkbC (MR20 section 3) depend on omega_tilde = omega - m Omega
    # directly, not only through D: independent assembly at l != m.
    P = AccretionDisc.Modes
    for beta in (.1, 1., 80.), (m, l, n) in ((2, 3, 1.), (3, 2, .8), (1, 2, 1.), (5, 4, 1.3))
        d = NewtonianDisc(h0=.07, p=.5, q=.75, beta=beta)
        mode = NewtonianMode(m, l, n)
        g = d.gamma
        for r in (.4, .9, 1.7, 3.)
            O = P.Omega(d, r); w = l * n - m * O; t = beta / O
            D = P.kappa(d, r)^2 - w^2
            cs = ((1 + (g * w * t)^2) / (1 + (w * t)^2))^(1 / 4) * P.csiso(d, r)
            k = sqrt(complex(-(D + (g * w * t)^2 / (1 + (g * w * t)^2) * P.Nr2(d, r)))) / cs
            @test P.kwave(d, r, mode) ≈ k rtol=1e-12
            D1 = 2P.kappa(d, r) * P.dkappadr(d, r) + 2m * w * P.dOmegadr(d, r)
            C = im * (k + im * (g - 1) / 2 * w * t / (1 + g * (w * t)^2) * k) + D1 / (4D) -
                P.dSigmadr(d, r) / (2P.Sigma(d, r)) + P.dcsisodr(d, r) / (2P.csiso(d, r)) - 1 / (2r) +
                (1 + im * g * w * t) / 2 / (1 + (g * w * t)^2) * P.LTinv(d, r)
            @test P.wkbC(d, r, mode) ≈ C rtol=1e-12
        end
    end
end

@testset "m = 1 inner boundary is evanescent and inward-decaying for every admissible (1,l) (lander review)" begin
    # FR22 setup (p = q = 0, h_p = 0.1, beta = 0) on its grid [0.1, 5] a, n = 1. Keplerian
    # D = kappa^2 - (Omega - l)^2 = l (2 Omega - l) > 0 inside r_+ = (2/l)^(2/3) > r_c = l^(-2/3); so whenever
    # corotation lies on the grid (l <= 31) the inner boundary is evanescent, and the m = 1 correction must
    # select the branch with h'/h = +sqrt(D)/c_s, which decays toward small r.
    P = AccretionDisc.Modes
    d = NewtonianDisc(p=0., q=0., h0=.1, beta=0.)
    x = loggrid(d, .1, 5., 1000)
    r = x[1]; O = r^(-3/2)
    for l in 1:31
        mode = NewtonianMode(1, l, 1.)
        w = l - O; D = l * (2O - l)
        @test D > 0
        @test P.Dfreq(d, r, mode) ≈ D rtol=1e-12
        k = P.kwave(d, r, mode)
        @test real(k) == 0 && imag(k) > 0
        ci, co = P.wkb_boundary_coeffs(d, x, 1, length(x), mode)
        @test ci ≈ sqrt(D) / d.h0 + (-3r^(-4) - 3w * r^(-5/2)) / (4D) - 1 / (2r) rtol=1e-12
        @test real(ci) > 0
        @test co == P.wkbC(d, x[end], mode)  # the outer boundary is never corrected
    end
    c = PotentialCache(log10.(x), x.^2, 2 .* x, fill(2., length(x)))
    @test_throws ArgumentError P.solve_lm_mode(d, NewtonianMode(1, 32, 1.), x, c)  # corotation 32^(-2/3) < 0.1
end
