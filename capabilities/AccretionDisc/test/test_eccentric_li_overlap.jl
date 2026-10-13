# ECC::gr-li-circular-overlap (C-LI-0, runtime gate): at e = 0, qTemp = 0 the OPTION A and adiabatic
# closures coincide (symbolic gate: derivation/scripts/certify_li_overlap.jl), so both eccentric paths
# must give the same solution and the same diagnostics. No circular kernel is used; the receipt driver
# also runs this file with the circular kernels poisoned (test/eccentric_modes_isolation.jl).
#
# Owner decision 2026-10-13 (gr-li-o-amf-definition): subtract the conserved baryon
# current h u_phi Mdot. Both accessors then give (72a), including the pressure current.
# The formerly broken parity/balance gates retain their original thresholds.
module LICircularOverlapTests
using Test, Printf
using AccretionDisc, Spacetime
using AccretionDisc.ForcingModes
const SW = AccretionDisc.SpiralWaves
# Preregistered: accumulated Float64 rounding, fixed before the run.
# Symmetric maximum norm, with no unit-sized floor hiding small amplitudes.
const OVERLAP_TOL = 65536eps(Float64)
function overlap_error(a, b)
    scale = max(maximum(abs,a), maximum(abs,b))
    iszero(scale) ? 0.0 : maximum(abs,a-b)/scale
end

function shared_limit(spin)
    li = SW.DiscModel(;χ=spin,r_ref=20.0,h_ref=0.06,pSigma=1.5,qTemp=0.0,
                       εp=1e-3,closure=LocallyIsothermalOptionA,evolution=:fixed)
    ad = SW.DiscModel(li.reference;closure=Adiabatic,evolution=:fixed)
    pert = SW.EccentricPerturber(EquatorialOrbit(spin,20.0,0.0);mass=1e-5,
                                 softening=FixedSoftening(0.3*0.06*20))
    (;li,ad,pert)
end

# Both closures on one grid, one source bank and the certified LI edge rows.
function solve_pair(s,ms,n)
    lc = [SW.ModeContext(s.li,s.pert,m,0;closure=LocallyIsothermalOptionA,domain=(10.0,40.0)) for m in ms]
    ac = [SW.ModeContext(s.ad,s.pert,c.mode;closure=Adiabatic,domain=c.domain) for c in lc]
    grids = [sort!(unique!(vcat(collect(range(10.0,40.0;length=n+1)),c.r_c))) for c in lc]
    samples = sort!(unique!(vcat([vcat(r,(r[1:end-1]+r[2:end])/2) for r in grids]...)))
    bank = SW.source_azimuth_bank(s.li,s.pert,samples;Nq=16,Nphi=1024,mmax=maximum(ms),blocksize=1)
    map(zip(lc,ac,grids)) do (lctx,actx,radii)
        @test lctx.mode.k == 0 && lctx.mode == actx.mode
        @test lctx.r_c == actx.r_c
        # The two default exterior approximations differ. Match them using
        # the independently certified LI impedance including transport.
        boundary = SW.mode_boundary_rows(lctx)
        for (i,outward) in ((1,false),(2,true))
            @test boundary[i] == (SW.li_boundary(lctx,lctx.domain[i];outward).row,0)
        end
        ls,as = map((lctx,actx)) do ctx
            grid = SW.EccentricGrid(copy(radii),ctx.r_c,1:0,SW.mode_cache_key(ctx),(;li_overlap=true))
            only(SW.solve_modes(_ -> ctx,s.pert,[ctx.mode];bank,
                grid_builder=_ -> grid,boundary_builder=_ -> boundary,landau=true,linsolve=:banded))
        end
        @test ls.closure == LocallyIsothermalOptionA && as.closure == Adiabatic
        @test ls.r == as.r && ls.source.jets == as.source.jets
        @test ls.params.rows == as.params.rows == boundary
        @test ls.params.landau && as.params.landau && ls.params.linsolve == as.params.linsolve == :banded
        @test max(ls.residual,as.residual) < 1e-8 && max(ls.bc_residual,as.bc_residual) < 1e-12
        @test maximum(abs,ls.y) > 0 && all(isfinite,ls.y) && all(isfinite,as.y)
        (;ls,as)
    end
end

@testset "C-LI-0 Float64 e=0 q=0 end-to-end eccentric parity" begin
    rows = 0
    for spin in (-0.9,0.9), pair in solve_pair(shared_limit(spin),2:10,256)
        ls,as = pair.ls,pair.as
        solutions = (ls,as)
        state = maximum(overlap_error(ls.y[i,:],as.y[i,:]) for i in 1:2)
        lf,af = SW.reconstruct_fluid.(solutions)
        fluid = maximum(overlap_error(getproperty.(lf,n),getproperty.(af,n)) for n in keys(first(lf)))
        lb,ab = SW.balance_diagnostics.(solutions)
        @test first(lb.emitted).tag == "li" && first(ab.emitted).tag == "adiabatic"
        @test eltype(lb.amf) == eltype(ab.amf) == Float64
        torque = overlap_error([lb.torque],[ab.torque])
        power = overlap_error([lb.power],[ab.power])
        td = overlap_error(lb.torque_density,ab.torque_density)
        pd = overlap_error(lb.power_density,ab.power_density)
        flux = overlap_error(lb.amf,ab.amf)
        @printf("OVERLAP,%+.1f,%d,%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                spin,ls.context.mode.m,length(ls.r),state,fluid,torque,power,td,pd,flux,lb.torque,ab.torque)
        for error in (state,fluid,torque,power,td,pd)
            @test error <= OVERLAP_TOL
        end
        @test flux <= OVERLAP_TOL
        rows += 1
    end
    @test rows == 18
end
println("[PASS] C-LI-0: 18 eccentric mode pairs; states, fluid, torque and power within ",OVERLAP_TOL,
        "; conserved AMF parity at the same tolerance")

# The balance dF/dr = sqrt(-g)<T> cell by cell in the midpoint scheme's own form (as the circular
# balance_residual), outside 3 sigma_eps of corotation, on two uniform levels. The adiabatic (72a)
# flux and the conserved OPTION A flux both converge at second order.
function balance_defects(sol)
    c,r = sol.context,sol.r
    F = SW.amf(sol)
    defect = Float64[]; density = Float64[]
    for i in 1:length(r)-1
        rm = (r[i]+r[i+1])/2
        abs(rm-c.r_c) < 3c.disc.σε && continue
        ym = (sol.y[:,i]+sol.y[:,i+1])/2
        t = SW.torque_density_mode(c,rm,ym;hmodes=sol.source(rm),landau=sol.params.landau)
        push!(defect,(F[i+1]-F[i])/(r[i+1]-r[i])-t); push!(density,t)
    end
    maximum(abs,defect)/maximum(abs,density)
end

@testset "C-LI-0 shared-limit angular-momentum balance of both AMF accessors" begin
    for spin in (-0.9,0.9)
        s = shared_limit(spin)
        levels = [only(solve_pair(s,(3,),n)) for n in (1024,2048)]
        ad = [balance_defects(x.as) for x in levels]
        li = [balance_defects(x.ls) for x in levels]
        @printf("AMF_BALANCE,%+.1f,m=3,adiabatic=%.6e,%.6e,ratio=%.3f,li=%.6e,%.6e,ratio=%.3f\n",
                spin,ad...,ad[1]/ad[2],li...,li[1]/li[2])
        @test 3.5 < ad[1]/ad[2] < 4.5 && ad[2] < 1e-3
        @test 3.5 < li[1]/li[2] < 4.5 && li[2] < 1e-3
    end
end
println("[PASS] C-LI-0 balance: both conserved AMF accessors converge at second order")
end
