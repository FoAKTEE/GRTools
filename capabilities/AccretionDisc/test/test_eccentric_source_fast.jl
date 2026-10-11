using Test, Printf, LinearAlgebra
using AccretionDisc, Spacetime, NumericsCore
using AccretionDisc.ForcingModes
import AccretionDisc.SpiralWaves as ESW

function fast_source_fixture(spin=0.5,p=20.0,e=0.3;soft=0.36)
    orbit = EquatorialOrbit(spin,p,e)
    disc = ESW.DiscModel(χ=spin,r_ref=p,h_ref=0.06,closure=:adiabatic,evolution=:fixed)
    perturber = ESW.EccentricPerturber(orbit;softening=FixedSoftening(soft))
    disc,perturber
end

fast_source_modes(p,pairs) = [HarmonicMode(m,k,orbit_frequency(p.orbit,m,k)) for (m,k) in pairs]

function fast_source_error(a,b,p;derivatives=true)
    # Independent statement of the preregistered dimensional budget.
    r=p.orbit.p
    scales=(1.0,r,r*r,1.0,1.0,1.0,r)
    maximum(eachindex(a)) do i
        maximum(eachindex(values(a[i]))) do j
            unit=p.mass*scales[derivatives ? (j+1)÷2 : j]/r/(derivatives && iseven(j) ? r : 1)
            abs(a[i][j]-b[i][j])/(1e-8*abs(b[i][j])+1e-12*unit)
        end
    end
end

const FAST_SOURCE_STAGE = get(ENV,"FAST_SOURCE_STAGE","all")
if FAST_SOURCE_STAGE in ("all","fft")
@testset "Eccentric source fast phase and compensated fallback" begin
    for (spin,e) in ((-0.7,0.3),(0.5,0.6),(0.9,0.0))
        d,p=fast_source_fixture(spin,20.0,e)
        b=ESW.source_azimuth_bank(d,p,[9.0,20.0,70.0];Nq=128,mmax=3)
        modes=fast_source_modes(p,[(m,k) for m in (0,1,3) for k in -12:12])
        reference=ESW.source_modes(b,modes)
        work=ESW.SourceProjectionWork()
        fast=ESW.source_modes(b,modes;method=:fft,work)
        error=fast_source_error(fast,reference,p)
        @test error<=1
        @printf("FAST_PARITY,%.1f,%.1f,128,%.9g,%d,%d,%d\n",spin,e,error,
                work.spread_terms,work.fft_butterflies,work.direct_terms)
        flush(stdout)
        @test ESW.source_modes(b,modes;method=:direct)==reference
        nod=ESW.source_modes(b,modes;method=:fft,derivatives=false)
        @test fast_source_error(nod,ESW.source_modes(b,modes;derivatives=false),p;derivatives=false)<=1
        @test all(nod[i].htt==fast[i].htt && nod[i].hrphi==fast[i].hrphi for i in eachindex(fast))
        selected=modes[[end,2,2,1]]
        @test ESW.source_modes(b,selected;method=:fft)==fast[:,[end,2,2,1]]
        @test size(ESW.source_modes(b,HarmonicMode{Float64}[];method=:fft))==(3,0)
        @test_throws ArgumentError ESW.source_modes(b,modes;method=:invalid)
        @test_throws ArgumentError ESW.source_modes(b,[HarmonicMode(4,0,4p.orbit.Ωφ)];method=:fft)
        @test_throws ArgumentError ESW.source_modes(b,[HarmonicMode(1,0,1.0)];method=:fft)
        # Sharing reduces phase-dependent work even when a signed band grows.
        @test work.spread_terms<14*length(b.q)*length(b.rs)*length(modes)
        @test work.direct_terms<14*length(b.q)*length(b.rs)*length(modes)÷2
        # The omitted Omega_r*Delta_t term must be observable at nonzero e.
        if e>0
            m,k=3,2
            naive=sum(b.time_weights[i]*b.events[i].tensor[1]*b.kernels[m+1,i,1,1]*
                cis(k*b.q[i]+m*p.orbit.Ωφ*b.delta_t[i]-m*b.delta_phi[i]) for i in eachindex(b.q))
            j=findfirst(mode->mode.m==m && mode.k==k,modes)
            @test abs(naive-reference[1,j].htt)>1e-5abs(reference[1,j].htt)
        end
    end
    # Exact R1, including a power-of-two bank that really enters the FFT.
    _,p=fast_source_fixture()
    event=ESW.EccentricSourceEvent(20.0,0.36,1.0,1.0,0.0,0.0,ntuple(_->1.0,7))
    for n in (3,4)
        kernels=zeros(ComplexF64,1,n,1,2)
        for derivative in 1:2
            kernels[1,1:3,1,derivative]=[complex(1e16,-1e16),complex(1,1),complex(-1e16,1e16)]
        end
        b=ESW.SourceAzimuthBank(p,[20.0],zeros(n),ones(n),zeros(n),zeros(n),fill(event,n),
                               kernels,[256],[0.1],1,0)
        work=ESW.SourceProjectionWork()
        h=only(ESW.source_modes(b,[HarmonicMode(0,0,0.0)];method=:fft,work))
        @test all(==(1+1im),values(h))
        @test work.direct_terms==14n
    end
    # Lander addition: cancellation ACROSS lattice rows. Uniform nodes take the plain-FFT plan,
    # whose first butterfly forms 1e16+1 before -1e16 arrives; only the roundoff screen and its
    # compensated fallback recover the exact 1+1im (a disabled screen returns 0 here).
    q4=[2pi*(i-1)/4 for i in 1:4]
    kernels=zeros(ComplexF64,1,4,1,2)
    for derivative in 1:2
        kernels[1,:,1,derivative]=[complex(1e16,-1e16),complex(-1e16,1e16),complex(1,1),0]
    end
    b=ESW.SourceAzimuthBank(p,[20.0],q4,ones(4),zeros(4),zeros(4),fill(event,4),kernels,[256],[0.1],1,0)
    @test ESW._source_k_plan(b,0).order==1
    h=only(ESW.source_modes(b,[HarmonicMode(0,0,0.0)];method=:fft))
    @test all(==(1+1im),values(h))
    # High |k| and non-power-of-two Nq are legal; no wrapping to another mode.
    d,p=fast_source_fixture()
    b=ESW.source_azimuth_bank(d,p,[70.0];Nq=30,mmax=1)
    modes=fast_source_modes(p,[(1,-81),(1,0),(1,65)])
    @test fast_source_error(ESW.source_modes(b,modes;method=:fft),ESW.source_modes(b,modes),p)<=1
    roomy=ESW.source_azimuth_bank(d,p,[70.0];Nq=16,mmax=3)
    tight=ESW.source_azimuth_bank(d,p,[70.0];Nq=16,mmax=3,max_memory_bytes=30000)
    @test tight.workspace_bound_bytes<=30000<roomy.workspace_bound_bytes
    @test tight.kernels==roomy.kernels
    uniform=ESW.SourceAzimuthBank(p,roomy.rs,roomy.q,roomy.time_weights,zeros(16),
        roomy.delta_phi,roomy.events,roomy.kernels,roomy.nphi,roomy.strips,roomy.blocksize,roomy.workspace_bound_bytes)
    @test ESW._source_k_plan(uniform,4).order==1
    uniform_modes=fast_source_modes(p,[(3,k) for k in -4:4])
    @test fast_source_error(ESW.source_modes(uniform,uniform_modes;method=:fft),
        ESW.source_modes(uniform,uniform_modes),p)<=1
end
end

if FAST_SOURCE_STAGE in ("all","radial")
@testset "Source adapted Chebyshev radial bank" begin
    d,p=fast_source_fixture()
    modes=fast_source_modes(p,[(m,k) for m in (0,1,3) for k in -4:4])
    # Cover both turning points, the interior, and both geometric tails.
    for center in (9.0,p.orbit.p/(1+p.orbit.e),20.0,p.orbit.p/(1-p.orbit.e),70.0)
        domain=(center-0.12,center+0.12)
        table=ESW.source_radial_bank(d,p,domain,modes;Nq=128)
        radii=collect(range(domain...;length=9)) .+ 0.001
        radii[end]=domain[2]
        direct=ESW.source_modes(ESW.source_azimuth_bank(d,p,radii;Nq=128,mmax=3),modes)
        interpolated=ESW.source_modes(table,modes,radii)
        err=fast_source_error(interpolated,direct,p)
        @test err<=1
        @printf("RADIAL_PARITY,%.9g,%d,%.9g\n",center,length(table.grid.rs),err)
        flush(stdout)
        s=ESW.mode_source(table,modes[1])
        @test s(radii[1])==interpolated[1,1]
        @test_throws DomainError s(prevfloat(domain[1]))
        @test_throws DomainError s(nextfloat(domain[2]))
        @test_throws ArgumentError ESW.mode_source(table,HarmonicMode(1,0,1.0))
        @test_throws ArgumentError ESW.mode_source(table,HarmonicMode{Int}(0,0,0))
        @test_throws ArgumentError ESW.mode_source(table,only(fast_source_modes(p,[(2,0)])))
        @test maximum(maximum(panel.error_bound ./ ESW._source_units(p,true)) for panel in table.grid.panels)<=2.5e-13
        # Reusing a callable table on different BVP grids performs zero orbit/
        # azimuth work. Allocation grows with output size, not Nq or mmax.
        s(center)
        @test @allocated(s(center)) < 4096
    end
    for (spin,porb,e) in ((-0.7,20.0,0.6),(0.9,20.0,0.0),(0.5,1e5,0.5))
        scale=porb/20
        d,p=fast_source_fixture(spin,porb,e;soft=0.36scale)
        modes=fast_source_modes(p,[(m,k) for m in (0,1,3) for k in -4:4])
        center=porb/(1+e)
        domain=(center-0.15scale,center+0.15scale)
        table=ESW.source_radial_bank(d,p,domain,modes;Nq=128)
        radii=collect(range(domain...;length=11))
        direct=ESW.source_modes(ESW.source_azimuth_bank(d,p,radii;Nq=128,mmax=3),modes)
        err=fast_source_error(ESW.source_modes(table,modes,radii),direct,p)
        @test err<=1
        @printf("RADIAL_REGIME,%.1f,%.0f,%.1f,%d,%.9g\n",spin,porb,e,length(table.grid.rs),err)
        flush(stdout)
    end
    d,p=fast_source_fixture()
    @test_throws ArgumentError ESW.source_radial_bank(d,p,(20.0,19.0),modes;Nq=16)
    @test_throws ArgumentError ESW.source_radial_bank(d,p,(9.0,70.0),modes;Nq=16,max_memory_bytes=1000)
    # Lander addition: multi-panel windows across each turning point and into its tail, so the
    # softening-scale panel adaptivity is exercised by parity (the narrow windows above fit in
    # one or two panels whatever the width rule).
    d,p=fast_source_fixture()
    modes=fast_source_modes(p,[(m,k) for m in (0,2) for k in (-3,0,2)])
    for domain in ((13.0,17.0),(26.0,31.0))
        table=ESW.source_radial_bank(d,p,domain,modes;Nq=16)
        @test length(table.grid.panels)>=3
        radii=[domain[1]+(domain[2]-domain[1])*(j+0.37)/17 for j in 0:16]
        direct=ESW.source_modes(ESW.source_azimuth_bank(d,p,radii;Nq=16,mmax=2),modes)
        err=fast_source_error(ESW.source_modes(table,modes,radii),direct,p)
        @test err<=1
        @printf("RADIAL_WIDE,%.1f,%.1f,%d,%d,%.9g\n",domain...,length(table.grid.panels),length(table.grid.rs),err)
        flush(stdout)
    end
    # Lander addition: restate the analytic panel bound independently (2x2 eigenvalues by
    # LinearAlgebra, the ellipse of each panel's own lo/hi) with the ACTUAL node count, so a
    # stored bound cannot drift from the interpolant that is really used.
    for (spin,e) in ((-0.7,0.3),(0.9,0.6))
        d,p=fast_source_fixture(spin,20.0,e)
        data=ESW._source_orbit_data(p,256)
        grid=ESW._source_radial_grid(p,(6.0,100.0),data)
        Q(K)=Symmetric([K.grr+K.vr^2 K.vr*K.vtan; K.vr*K.vtan K.gtan+K.vtan^2])
        lo,hi=minimum(K->eigmin(Q(K)),data.events),maximum(K->eigmax(Q(K)),data.events)
        peri,apo=extrema(K.rp for K in data.events)
        soft=p.softening.eps
        amp=[sum(abs(data.tw[i]*data.events[i].tensor[c]) for i in eachindex(data.q)) for c in 1:7]
        r=p.orbit.p
        scales=(1.0,r,r*r,1.0,1.0,1.0,r)
        unit=[p.mass*scales[(j+1)÷2]/r/(iseven(j) ? r : 1) for j in 1:14]
        @test first(grid.edges)==6.0 && last(grid.edges)==100.0
        @test grid.rs==reduce(vcat,[panel.nodes for panel in grid.panels])
        worst=0.0
        for panel in grid.panels
            h,c=(panel.hi-panel.lo)/2,(panel.hi+panel.lo)/2
            B=soft^2+lo*max(peri-(c+1.25h),(c-1.25h)-apo,0.0)^2-hi*(0.75h)^2
            @test B>0
            n=length(panel.nodes)-1
            for j in 1:14
                M=amp[(j+1)÷2]*(isodd(j) ? 1/sqrt(B) : sqrt(hi)*sqrt(1+2hi*(0.75h)^2/B)/B)
                worst=max(worst,4M*2.0^(-n)/(2.5e-13*unit[j]))
            end
        end
        @test worst<=1+1e-10
        @printf("RADIAL_BOUND,%.1f,%.1f,%d,%d,%.9g\n",spin,e,length(grid.panels),length(grid.rs),worst)
        flush(stdout)
    end
end
end

if FAST_SOURCE_STAGE in ("all","solver")
@testset "Shared source on different grids and closures" begin
    d,p=fast_source_fixture(0.4,20.0,0.2;soft=2.0)
    domain=(28.0,32.0)
    modes=fast_source_modes(p,[(2,1),(3,-1)])
    table=ESW.source_radial_bank(d,p,domain,modes;Nq=32)
    for closure in (Adiabatic,LocallyIsothermalOptionA), (j,mode) in enumerate(modes)
        disc=ESW.DiscModel(;χ=0.4,r_ref=20.0,h_ref=0.06,qTemp=0.75,closure,evolution=:fixed)
        c=ESW.ModeContext(disc,p,mode;closure,domain)
        radii=collect(range(domain...;length=33+16j))
        grid=ESW.EccentricGrid(radii,c.r_c,1:0,ESW.mode_cache_key(c),(;source_fast=true))
        source=ESW.mode_source(table,mode)
        before=deepcopy(table.jets)
        fast=ESW.solve_mode(c;source,grid,boundary_rows=ESW.mode_boundary_rows(c))
        direct=ESW.solve_mode(c,p;grid,bank_options=(;Nq=32))
        relative=norm(fast.y-direct.y)/norm(direct.y)
        @test relative<1e-8
        @test table.jets==before
        @printf("SHARED_SOLVER,%s,%d,%d,%d,%.9g\n",string(closure),mode.m,mode.k,length(radii),relative)
        flush(stdout)
    end
end
end
