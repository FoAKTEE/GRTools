# ECC::shared-newtonian-projector: the all-mode Newtonian source bank `potential_bank` and its
# per-mode `harmonic_cache` views (DESIGN v3 §3, Modes/EccentricPlanet.jl). The two API test sets
# always run. The calibrated holdout, independent-oracle, circular, full FR22 mask and three-level
# convergence sets (~3 min at -t1) run under GRTOOLS_SLOW_TESTS=1.
#
# Gate (registered from self-convergence before any parity was inspected; emri-disk
# results/emri-eccentric/paper_EMRI-ECC/evidence/shared-newtonian-projector/GATES.json, sha256
# 88ec0f316cf006f2d403391cf23db3df4529379211fbd90f2530bc746f2e0d5f): per component j = Φ, Φ′, Φ″,
# |A-B| <= 4(δA+δB+F_j), F_j = PROJECTOR_FLOOR_EPS_UNITS[j]*eps*S_j, with δ the final self-refinement
# of each method, δ <= 1e-7 S_j, and S_j = max(max |amplitude| over the case, Gmu/a^j). No pointwise
# relative gate. Coefficients with max(|A|,|B|) <= F_j are noise-floor rows under the same gate.
module NewtonianProjectorTests

using Test
using AccretionDisc
using AccretionDisc.ForcingModes
using NumericsCore: gl_nodes_weights

const MO = AccretionDisc.Modes
const SLOW = get(ENV, "GRTOOLS_SLOW_TESTS", "0") == "1"
const PROJECTOR_FLOOR_EPS_UNITS = [8192, 4096, 65536]

bank_values(b) = [getfield(b.caches[e,k],(:P,:dP,:d2P)[j])[i]
                 for e in axes(b.caches,1), i in eachindex(b.caches[1,1].P),
                     k in axes(b.caches,2), j in 1:3]

component_scales(values, a=1.0, Gmu=1.0) =
    [max(maximum(abs,selectdim(values,4,j)),Gmu/a^j) for j in 1:3]

component_deltas(a,b) = [maximum(abs,selectdim(a.-b,4,j)) for j in 1:3]

# Test-only oracles. The runtime has one sampler; these do not call it.
# Factorised orbit integral: midpoint mean-anomaly nodes (offset from the runtime's endpoint
# nodes) over the softened Laplace coefficients, with analytic jets of the softening length.
function laplace_jet(d, o, r, mode, soft; N=1024, Gmu=1.0)
    l = mode.m+mode.k
    power = (3-d.q)/2
    out = zeros(3)
    for i in 0:N-1
        M = 2pi*(i+0.5)/N
        p = MO.planet_orbit(M,o)
        rp = p.r_p
        ep = MO._bank_softlength(soft,d.h0,power,o.a,rp,r)
        dep = soft isa LocalScaleHeight ? power*ep/r : 0.0
        d2ep = soft isa LocalScaleHeight ? power*(power-1)*ep/r^2 : 0.0
        alpha,epsilon = r/rp,ep/rp
        b = MO.laplace_b_soft(0.5,mode.m,alpha,epsilon;nphi=N)
        db = MO.laplace_b_soft_deriv(0.5,mode.m,alpha,epsilon;nphi=N)
        d2b = MO.laplace_b_soft_2deriv(0.5,mode.m,alpha,epsilon;nphi=N)
        b3 = iszero(dep) ? 0.0 : MO.laplace_b_soft(1.5,mode.m,alpha,epsilon;nphi=N)
        db3 = iszero(dep) ? 0.0 : MO.laplace_b_soft_deriv(1.5,mode.m,alpha,epsilon;nphi=N)
        b5 = iszero(dep) ? 0.0 : MO.laplace_b_soft(2.5,mode.m,alpha,epsilon;nphi=N)
        w = -Gmu*cos(mode.m*p.phi_p-l*M)/N
        out[1] += w*b/rp
        out[2] += w*(db/rp^2-ep*dep*b3/rp^3)
        out[3] += w*(d2b/rp^3-2ep*dep*db3/rp^4-
                     (dep^2+ep*d2ep)*b3/rp^3+3ep^2*dep^2*b5/rp^5)
    end
    out
end

# Composite Gauss in eccentric anomaly and ABSOLUTE azimuth: Cartesian separation, Kepler time
# Jacobian (1-e cos E)/n and the physical n/(2pi^2) of FR22 (E-ECC-2). Independent of the runtime
# mean-anomaly grid, Kepler solver, conjunction-stable distance and Laplace factorisation; it shares
# only the softening-length helper, whose four policies are checked against hand formulas below.
function double_jet(d, o, r, mode, soft; panels=16, order=32, Gmu=1.0)
    x,w = gl_nodes_weights(Float64,order)
    angles = [(j+(v+1)/2)*2pi/panels for j in 0:panels-1 for v in x]
    weights = [v*pi/panels for _ in 1:panels for v in w]
    c,s = cos.(angles),sin.(angles)
    cm,sm = cos.(mode.m.*angles),sin.(mode.m.*angles)
    out = zeros(3)
    power=(3-d.q)/2
    for i in eachindex(angles)
        E=angles[i]; M=E-o.e*s[i]
        xp,yp=o.a*(c[i]-o.e),o.a*sqrt(1-o.e^2)*s[i]
        rp=hypot(xp,yp)
        ep=MO._bank_softlength(soft,d.h0,power,o.a,rp,r)
        dep=soft isa LocalScaleHeight ? power*ep/r : 0.0
        d2ep=soft isa LocalScaleHeight ? power*(power-1)*ep/r^2 : 0.0
        cl,sl=cos((mode.m+mode.k)*M),sin((mode.m+mode.k)*M)
        dt=weights[i]*(1-o.e*c[i])/o.n
        for j in eachindex(angles)
            dx,dy=r*c[j]-xp,r*s[j]-yp
            S=dx^2+dy^2+ep^2
            U=dx*c[j]+dy*s[j]+ep*dep
            v=inv(sqrt(S))
            weight=Gmu*o.n/(2pi^2)*dt*weights[j]*(cm[j]*cl+sm[j]*sl)
            out[1] -= weight*v
            out[2] += weight*U*v^3
            out[3] += weight*((1+dep^2+ep*d2ep)*v^3-3U^2*v^5)
        end
    end
    out
end

function check_pair(name, a, b, da, db, scale)
    err = abs.(a.-b)
    floor = eps().*PROJECTOR_FLOOR_EPS_UNITS.*scale
    gate = 4 .* (da.+db.+floor)
    for j in 1:3
        @test da[j] <= 1e-7*scale[j]
        @test db[j] <= 1e-7*scale[j]
        @test err[j] <= gate[j]
    end
    maximum(err./gate)
end

@testset "Newtonian bank contract and inherited type" begin
    d=MO.NewtonianDisc{Float64}()
    o=MO.EccentricOrbit(1,0,1)
    modes=[NewtonianMode(1,1,1),NewtonianMode(1,2,1),NewtonianMode(0,1,1)]
    rs=[.7,.8,1.2,1.3]
    b=MO.potential_bank(d,[o],rs,modes;Nt=64,Nphi=64)
    c=MO.harmonic_cache(b,1,modes[1])
    @test c isa MO.PotentialCache{Float64}
    @test c === MO.harmonic_cache(b,1,modes[1])
    @test c.lg === b.caches[1,2].lg
    @test MO.potential_at(c,rs[2]) == (c.P[2],c.dP[2],c.d2P[2])
    @test length(b.mode_indices)==3
    @test all(isfinite,bank_values(b))
    @test maximum(abs,b.caches[1,2].P)<1e-13
    @test maximum(abs,b.caches[1,3].P)<1e-13
    b2=MO.potential_bank(d,[o],rs,modes;Nt=64,Nphi=64,softening=FixedSoftening(.12))
    @test b.key != b2.key
    @test b.caches[1,1].P != b2.caches[1,1].P
    @test MO.potential_bank(d,[o],rs,HarmonicMode{Float64}[];Nt=64,Nphi=64).caches |> isempty
    @test_throws KeyError MO.harmonic_cache(b,1,NewtonianMode(1,1,2))
    @test_throws ArgumentError MO.potential_bank(d,[o],reverse(rs),modes)
    @test_throws ArgumentError MO.potential_bank(d,[o],rs,[modes[1],modes[1]])
    @test_throws ArgumentError MO.potential_bank(d,[o],rs,[NewtonianMode(1,1,2)])
    @test_throws ArgumentError MO.potential_bank(d,[o],rs,[NewtonianMode(32,1,1)];Nt=64,Nphi=64)
    @test_throws ArgumentError MO.potential_bank(d,[o],rs,[NewtonianMode(1,-32,1)];Nt=64,Nphi=64)
    @test_throws ArgumentError MO.potential_bank(d,[o],rs,modes;Nt=63,method=:fft)
    @test_throws ArgumentError MO.potential_bank(d,[o],rs,modes;derivatives=1)
    @test_throws ArgumentError MO.potential_bank(MO.NewtonianDisc{Float64}(indirect=1),[o],rs,modes)
    @test_throws DomainError MO.potential_bank(d,[o],[1.0],modes;softening=FixedSoftening(0),Nt=64,Nphi=64)
    # Precision inheritance is a small API unit test, not a multiprecision audit.
    # Do not change global precision: inherit whichever the test caller supplies.
    bits=precision(BigFloat)
    for method in (:direct,:fft)
        bb=MO.potential_bank(d,[MO.EccentricOrbit(big(1),big(1)/4,big(1))],
            [big(2)],modes;softening=FixedSoftening(big(1)/10),Nt=8,Nphi=8,method)
        @test bb.caches[1,1] isa MO.PotentialCache{BigFloat}
        @test precision(bb.caches[1,1].P[1])==bits
        @test precision(BigFloat)==bits
    end
end

@testset "Softening locations, local derivatives and rectangular nodes" begin
    h,power,a,rp,r=.1,1.25,1.4,.8,2.0
    @test MO._bank_softlength(FixedSoftening(.04),h,power,a,rp,r)==.04
    @test MO._bank_softlength(SemimajorScaleHeight(.3),h,power,a,rp,r)==.3h*a^power
    @test MO._bank_softlength(InstantaneousScaleHeight(.3),h,power,a,rp,r)==.3h*rp^power
    @test MO._bank_softlength(LocalScaleHeight(.3),h,power,a,rp,r)==.3h*r^power
    d=MO.NewtonianDisc{Float64}(h0=h,q=.5)
    o=MO.EccentricOrbit(a,.25,7/3)
    mode=NewtonianMode(3,2,o.n)
    step=.001
    rs=r.+(-2:2).*step
    c=MO.harmonic_cache(MO.potential_bank(d,[o],rs,[mode];softening=LocalScaleHeight(.3),
        Nt=256,Nphi=512,method=:direct),1,mode)
    fd1=(c.P[1]-8c.P[2]+8c.P[4]-c.P[5])/(12step)
    fd2=(-c.P[1]+16c.P[2]-30c.P[3]+16c.P[4]-c.P[5])/(12step^2)
    @test abs(fd1-c.dP[3])<=1e-8
    @test abs(fd2-c.d2P[3])<=1e-8
    # The absolute 2Re convention includes the zero bin (twice the mean).
    zero_mode=NewtonianMode(0,0,1)
    circ=MO.EccentricOrbit(1,0,1)
    z=MO.potential_bank(d,[circ],[2.0],[zero_mode];Nt=63,Nphi=65)
    @test abs(z.caches[1,1].P[1]+MO.laplace_b_soft(.5,0,2.0,d.soft;nphi=512))<1e-13
    bf=MO.potential_bank(d,[o],[2.0],[mode];Nt=128,Nphi=256,method=:fft)
    bd=MO.potential_bank(d,[o],[2.0],[mode];Nt=127,Nphi=255,method=:direct)
    @test maximum(abs.(bank_values(bf).-bank_values(bd)))<1e-12
    # Semimajor and moving lengths agree at e=0, including nonunit a.
    oc=MO.EccentricOrbit(a,0,o.n)
    bc=MO.potential_bank(d,[oc],[2.0],[mode];softening=FixedSoftening(.3h*a^power),Nt=128,Nphi=256)
    for soft in (SemimajorScaleHeight(.3),InstantaneousScaleHeight(.3))
        bs=MO.potential_bank(d,[oc],[2.0],[mode];softening=soft,Nt=128,Nphi=256)
        @test bank_values(bc)==bank_values(bs)
    end
end

# Lander addition (review finding F6): without it the always-on tier had no independent oracle, so
# a convention error shared by both backends (a mirrored orbit, Gmu dropped from one jet) passed the
# routine suite and was caught only under GRTOOLS_SLOW_TESTS=1. Smooth cases off the orbit reach
# roundoff at 128 nodes; the gate is the registered absolute gate of the slow sets.
@testset "Routine-tier independent oracle off the orbit" begin
    d=MO.NewtonianDisc{Float64}(h0=.1,q=.5)
    a,n,Gmu=1.4,7/3,3/7
    o=MO.EccentricOrbit(a,.3,n)
    modes=[NewtonianMode(m,l,n) for (m,l) in ((1,-1),(0,1),(2,2),(3,2),(2,1))]
    for r in (.45a,1.8a), soft in (FixedSoftening(.12a),LocalScaleHeight(.3)), method in (:direct,:fft)
        b1,b2=[bank_values(MO.potential_bank(d,[o],[r],modes;softening=soft,Nt=N,Nphi=N,method,Gmu))
               for N in (128,256)]
        sc=component_scales(b2,a,Gmu)
        db=component_deltas(b1,b2)
        for (k,mode) in enumerate(modes)
            r1,r2=[double_jet(d,o,r,mode,soft;panels=P,order=32,Gmu) for P in (4,8)]
            check_pair("routine/$method/r=$(r/a)a/k=$k",b2[1,1,k,:],r2,db,abs.(r2.-r1),sc)
        end
    end
end

if SLOW
@testset "Calibrated holdout and independent oracles" begin
    d=MO.NewtonianDisc{Float64}(h0=.1,q=.5)
    a,n,Gmu=1.4,7/3,3/7
    modes=[NewtonianMode(m,l,n) for (m,l) in ((1,-1),(0,1),(2,2),(3,2),(8,10))]
    rs=a.*[.73,1.0,1.36]
    worst=Dict{String,Float64}()
    for (pn,soft) in (("fixed",FixedSoftening(.03a)),("semimajor",SemimajorScaleHeight(.3)),
                      ("instantaneous",InstantaneousScaleHeight(.3)),("local",LocalScaleHeight(.3))), e in (.25,.6)
        o=MO.EccentricOrbit(a,e,n)
        bs=Dict()
        for method in (:direct,:fft), N in (512,1024,2048)
            bs[(method,N)]=bank_values(MO.potential_bank(d,[o],rs,modes;
                softening=soft,Nt=N,Nphi=N,method,Gmu))
        end
        sc=component_scales(bs[(:fft,2048)],a,Gmu)
        floor=eps().*PROJECTOR_FLOOR_EPS_UNITS.*sc
        deltas=Dict(method=>component_deltas(bs[(method,1024)],bs[(method,2048)]) for method in (:fft,:direct))
        for method in (:fft,:direct)
            prev=component_deltas(bs[(method,512)],bs[(method,1024)])
            @test all(deltas[method] .<= max.(prev./4,floor))
        end
        for i in eachindex(rs), k in eachindex(modes)
            av,bv=bs[(:fft,2048)][1,i,k,:],bs[(:direct,2048)][1,i,k,:]
            worst["backends"]=max(get(worst,"backends",0.0),
                check_pair("backends/$pn/e=$e/r=$i/k=$k",av,bv,deltas[:fft],deltas[:direct],sc))
        end
        # Two independent integrations at three resolutions, including a signed
        # frequency, m=0, and off-diagonal modes. Holdout r=a crosses the orbit.
        k=e==.25 ? (pn=="fixed" ? 1 : pn=="local" ? 2 : 4) : 5
        r=rs[2]; mode=modes[k]
        fs=[laplace_jet(d,o,r,mode,soft;N,Gmu) for N in (512,1024,2048)]
        ds=[double_jet(d,o,r,mode,soft;panels=P,order=32,Gmu) for P in (16,32,64)]
        df,dd=abs.(fs[3].-fs[2]),abs.(ds[3].-ds[2])
        @test all(df .<= max.(abs.(fs[2].-fs[1])./4,floor))
        @test all(dd .<= max.(abs.(ds[2].-ds[1])./4,floor))
        ref=bs[(:fft,2048)][1,2,k,:]
        worst["factorised"]=max(get(worst,"factorised",0.0),
            check_pair("factorised/$pn/e=$e",ref,fs[3],deltas[:fft],df,sc))
        worst["double"]=max(get(worst,"double",0.0),
            check_pair("double/$pn/e=$e",ref,ds[3],deltas[:fft],dd,sc))
        # Negative controls: dropping the mean motion or doubling the projection is resolved.
        if maximum(abs,ref) > 100maximum(floor)
            @test maximum(abs.(ref.-ds[3]./n))>100maximum(floor)
            @test maximum(abs.(ref.-2ds[3]))>100maximum(floor)
        end
        println("ORACLES $pn e=$e scale=$sc factor_delta=$df double_delta=$dd");flush(stdout)
    end
    println("ORACLES max error/gate ",sort(collect(worst)));flush(stdout)
end

@testset "Circular analytic diagonal and vanishing sidebands" begin
    d=MO.NewtonianDisc{Float64}(soft=.03)
    o=MO.EccentricOrbit(1,0,1)
    modes=[NewtonianMode(m,l,1) for m in (1,3,8,50,150) for l in (m-1,m,m+1)]
    rs=[.6,1.0,1.4]
    worst=0.0
    for method in (:fft,:direct)
        b1=bank_values(MO.potential_bank(d,[o],rs,modes;Nt=1024,Nphi=1024,method))
        b2=bank_values(MO.potential_bank(d,[o],rs,modes;Nt=2048,Nphi=2048,method))
        sc=component_scales(b2)
        db=component_deltas(b1,b2)
        for (i,r) in enumerate(rs), (k,m) in enumerate(modes)
            refs=[]
            for N in (1024,2048)
                dd=MO.NewtonianDisc{Float64}(soft=.03,lap_nphi=N)
                push!(refs,m.k==0 ? [MO.Phi_m(dd,r,m.m),MO.dPhidr_m(dd,r,m.m),MO.d2Phidr2_m(dd,r,m.m)] : zeros(3))
            end
            worst=max(worst,check_pair("circular/$method/r=$r/m=$(m.m)/l=$(m.m+m.k)",
                b2[1,i,k,:],refs[2],db,abs.(refs[2].-refs[1]),sc))
        end
    end
    println("CIRCULAR max error/gate $worst");flush(stdout)
end

@testset "FR22 all 2514 modes and nested truncations" begin
    # Landed enumeration: reproduce/arXiv-2207.14637/_fr22_setup.jl:19-24.
    modes=[NewtonianMode(m,l,1) for m in 1:150 for l in max(1,m-8):m+8]
    @test length(modes)==2514
    d=MO.NewtonianDisc{Float64}(soft=.03)
    os=[MO.EccentricOrbit(1,e,1) for e in (0.0,.1,.5)]
    rs=[.7,1.0,1.3]
    coarse=MO.potential_bank(d,os,rs,modes;Nt=1024,Nphi=1024,method=:direct)
    fine=MO.potential_bank(d,os,rs,modes;Nt=2048,Nphi=2048,method=:direct)
    ff=MO.potential_bank(d,os,rs,modes;Nt=2048,Nphi=2048,method=:fft)
    c,b,f=bank_values(coarse),bank_values(fine),bank_values(ff)
    for ei in eachindex(os)
        sc=component_scales(b[ei:ei,:,:,:])
        delta=component_deltas(c[ei:ei,:,:,:],b[ei:ei,:,:,:])
        err=component_deltas(b[ei:ei,:,:,:],f[ei:ei,:,:,:])
        check_pair("FR22/e=$(os[ei].e)/full_max",err,zeros(3),delta,zeros(3),sc)
        println("FR22 e=$(os[ei].e) backend_max=$err refinement_max=$delta scale=$sc");flush(stdout)
    end
    counts_m=[count(x->x.m<=m,modes) for m in 0:10:150]
    counts_l=[count(x->abs(x.k)<=k,modes) for k in 0:8]
    @test counts_m==[0,134,304,474,644,814,984,1154,1324,1494,1664,1834,2004,2174,2344,2514]
    @test counts_l==[150,449,747,1044,1340,1635,1929,2222,2514]
    for mcut in 0:10:150, mode in filter(m->m.m<=mcut,modes)
        @test MO.harmonic_cache(fine,2,mode) === fine.caches[2,fine.mode_indices[mode]]
    end
    for kcut in 0:8, mode in filter(m->abs(m.k)<=kcut,modes)
        @test MO.harmonic_cache(fine,3,mode) === fine.caches[3,fine.mode_indices[mode]]
    end
    # The independently projected sparse set must match its bank handles.
    subset=modes[[1,17,149,700,2000,2514]]
    sb=MO.potential_bank(d,[os[2]],rs,subset;Nt=2048,Nphi=2048,method=:direct)
    for mode in subset
        cc=MO.harmonic_cache(sb,1,mode); cf=MO.harmonic_cache(fine,2,mode)
        @test maximum(abs.(cc.P.-cf.P)) <= eps()*PROJECTOR_FLOOR_EPS_UNITS[1]*2
        @test maximum(abs.(cc.dP.-cf.dP)) <= eps()*PROJECTOR_FLOOR_EPS_UNITS[2]*10
        @test maximum(abs.(cc.d2P.-cf.d2P)) <= eps()*PROJECTOR_FLOOR_EPS_UNITS[3]*1000
    end
end

@testset "FR22 three-level convergence cap and decay" begin
    ms=[NewtonianMode(m,l,1) for m in 1:150 for l in max(1,m-8):m+8]
    d=MO.NewtonianDisc{Float64}(soft=.03)
    for e in (0.0,.1,.5)
        vals=[bank_values(MO.potential_bank(d,[MO.EccentricOrbit(1,e,1)],[.7,1.0,1.3],ms;
            Nt=N,Nphi=N,method=:direct)) for N in (512,1024,2048)]
        sc=component_scales(vals[3]); prev=component_deltas(vals[1],vals[2]); last=component_deltas(vals[2],vals[3])
        for j in 1:3
            @test last[j]<=1e-7*sc[j]
            @test last[j]<=max(prev[j]/4,eps()*PROJECTOR_FLOOR_EPS_UNITS[j]*sc[j])
        end
        println("FR22 SELF e=$e previous=$prev last=$last scale=$sc");flush(stdout)
    end
end
else
    @info "shared-newtonian-projector oracle, circular, full-mask and convergence sets skipped (~3 min); set GRTOOLS_SLOW_TESTS=1"
end

end # module NewtonianProjectorTests
