using Test, AccretionDisc, NumericsCore

# Exercise the opt-in crossing keyword whenever this registered FR22 suite runs.
include("test_fr22_lindblad_crossing.jl")

@testset "FR22 harmonic safeguarded response" begin
    P = AccretionDisc.Modes
    d = NewtonianDisc(h0=.1, q=0., p=0., beta=0., gamma=1., soft=.03, lap_nphi=256)
    x = loggrid(d,.1,5.,1000)
    c = PotentialCache(d,1,.1,5.;n=1000)
    mode = AccretionDisc.ForcingModes.NewtonianMode(1,1,1.)
    opts = (;phgrad_steps=1, phgrad_safeguard=true, diagnostics=true, anchor=:corotation)
    s = solve_lm_mode(d,mode,x,c;opts...)
    circ = solve_m_mode(d,1,x,c;opts...)
    @test isequal(s,circ)
    @test s.C_before[1] == s.C_after[1]
    @test all(isfinite,s.mode)
    a = s.phase_safeguard
    if a.status != :rejected
        @test isfinite(a.X_after[2]) && a.X_after[2] < a.X_before[2]
        @test a.inner_response_change <= .1
    end
end

@testset "grid-capped integration reaches actual target states" begin
    P = AccretionDisc.Modes
    NC = NumericsCore
    d = NewtonianDisc(h0=.1,p=0.,q=0.,beta=0.,gamma=1.,lap_nphi=64)
    x = loggrid(d,.2,3.,101)
    mode = AccretionDisc.ForcingModes.NewtonianMode(2,3,1.)
    c = PotentialCache(d,2,.2,3.;n=100)
    bc = (1.0+.2im, .3-.1im)
    y0 = [real(bc[1]),imag(bc[1]),real(bc[2]),imag(bc[2])]
    ico = P.corotation_index(d,x,mode)
    rhs = P.CachedMasterRHS(d,mode,1,c)
    f = P.integrate_master(d,x,bc,mode,1,c;rtol=1e-7,atol=1e-9,dense=false)
    for dir in (-1,1)
        Y = NC.integrate_grid(rhs,y0,x,ico,dir;rtol=1e-7,atol=1e-9,dense=false)
        inds = dir == -1 ? (1:ico) : (ico:length(x))
        @test f[inds,1] == complex.(Y[inds,1],Y[inds,2])
        @test f[inds,2] == complex.(Y[inds,3],Y[inds,4])
    end
    @test P.integrate_master(d,x,bc,mode,1,c) == P.integrate_master(d,x,bc,mode,1,c;dense=true)
end

@testset "off-diagonal FR22 m=1 and m>=2 branches" begin
    P=AccretionDisc.Modes; F=AccretionDisc.ForcingModes
    d=NewtonianDisc(h0=.1,p=0.,q=0.,beta=0.,gamma=1.,soft=.03)
    x=loggrid(d,.1,5.,1000)
    modes=[F.NewtonianMode(m,l,1.) for (m,l) in ((1,2),(1,3),(2,3))]
    b=P.potential_bank(d,[P.EccentricOrbit(1.,.1,1.)],loggrid(d,.1,5.,256),modes;
        softening=F.FixedSoftening(.03),Nt=128,Nphi=128)
    for mode in modes
        c=P.harmonic_cache(b,1,mode)
        s=solve_lm_mode(d,mode,x,c;phgrad_steps=1,phgrad_safeguard=true,
            reint_steps=2,anchor=:corotation,diagnostics=true,rtol=1e-9,atol=1e-11,dense=false)
        a=s.phase_safeguard
        @test all(isfinite,s.mode)
        @test a.status in (:accepted,:backtracked,:rejected)
        @test P.mr19_boundary_coeffs(d,x,1,length(x),F.NewtonianMode(1,1,1.)) ==
              P.mr19_boundary_coeffs(d,x,1,length(x),1)
        if mode.m==1
            @test s.C_before[1] === s.C_after[1]
        end
        if a.status==:rejected
            @test s.C_before === s.C_after
        elseif mode.m==1
            @test a.X_after[2] < a.X_before[2]
            @test a.inner_response_change <= .1
        else
            @test all(a.X_after .< a.X_before)
        end
        _,_,sig=vars_from_enthalpy(d,s.mode,x,mode,c)
        @test sig ≈ P.Sigma.(Ref(d),x) .* s.mode[:,1] ./ P.csiso.(Ref(d),x).^2 rtol=5e-15
    end
end

@testset "FR22 rtol, atol and dense reach every integration (lander review)" begin
    # solve_lm_mode must equal a replay of the solver's own steps (two re-seeds, adaptive boundaries,
    # phase step off) assembled from integrate_master with the SAME settings, and must differ from the
    # replay with the dense default or with the cached default tolerances.
    P = AccretionDisc.Modes; F = AccretionDisc.ForcingModes
    d = NewtonianDisc(h0=.1, p=0., q=0., beta=0., gamma=1., lap_nphi=64)
    x = loggrid(d, .2, 3., 301); n = length(x)
    mode = F.NewtonianMode(2, 3, 1.); c = PotentialCache(d, 2, .2, 3.; n=300)
    function replay(; rtol, atol, dense)
        kw = (; rtol, atol, dense); z, o = zero(ComplexF64), one(ComplexF64)
        h1 = P.integrate_master(d, x, (o, z), mode, 0, c; kw...)
        h2 = P.integrate_master(d, x, (z, o), mode, 0, c; kw...)
        ih = P.integrate_master(d, x, (z, z), mode, 1, c; kw...)
        i1 = P.corotation_index(d, x, mode)
        Ci, Co = P.wkb_boundary_coeffs(d, x, 1, n, mode)
        a1, a2 = P.homogeneous_coefficients(h1, h2, ih, Ci, Co, 1, n)
        f = a1 .* h1 .+ a2 .* h2 .+ ih
        ii, io = P.mode_boundaries(f, i1); Ci, Co = P.wkb_boundary_coeffs(d, x, ii, io, mode)
        for _ in 1:2
            ico = P.corotation_index(d, x, mode)
            ih = P.integrate_master(d, x, (f[ico, 1], f[ico, 2]), mode, 1, c; kw...)
            a1, a2 = P.homogeneous_coefficients(h1, h2, ih, Ci, Co, ii, io)
            f = a1 .* h1 .+ a2 .* h2 .+ ih
            ii, io = P.mode_boundaries(f, i1); Ci, Co = P.wkb_boundary_coeffs(d, x, ii, io, mode)
            a1, a2 = P.homogeneous_coefficients(h1, h2, ih, Ci, Co, ii, io)
            f = a1 .* h1 .+ a2 .* h2 .+ ih
        end
        f[1:ii-1, :] .= 0; f[io+1:end, :] .= 0
        f
    end
    bits(a, b) = reinterpret(UInt64, vec(a)) == reinterpret(UInt64, vec(b))
    s = solve_lm_mode(d, mode, x, c; rtol=1e-9, atol=1e-11, dense=false,
                      reint_steps=2, phgrad_steps=0, anchor=:corotation)
    @test bits(s, replay(rtol=1e-9, atol=1e-11, dense=false))
    @test !bits(s, replay(rtol=1e-9, atol=1e-11, dense=true))
    @test !bits(s, replay(rtol=1e-8, atol=1e-10, dense=false))
end
