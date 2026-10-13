# Owner decision 2026-10-13, gr-li-o-amf-definition. General closure: no pressure
# law or Gamma1 restriction. Run in derivation/; ARGS[1] is the evidence directory.
# Text sources: ref-paper/arXiv-2608.17003/src/main.tex:260-268 (19),
# :784-812 (72a); no mirrored content is executed. R3 thermal_balance.jl certifies
# baryon conservation + energy conservation -> (19), with zero advected constants.
using Symbolics, Test
include(joinpath(@__DIR__,"..","src","R3Derivation.jl"))
const R3=R3Derivation
const OUT=isempty(ARGS) ? mktempdir() : mkpath(ARGS[1])
@variables lam rho rho_r eps_r hspec measure u_phi ur2 ut wd
const enth=R3.eps+R3.p
# Independent barred amplitudes: no closure or reality specialization in the proof.
const de=R3.y0
const dp=R3.htt
const ur=im*R3.y1
const uphi=R3.hrphi
const logmeasure=R3.hrr
const xi=R3.htphi
const gradient=eps_r-hspec*rho_r
const drho=(de+xi*eps_r)/hspec-xi*rho_r
const H=Dict(hspec=>enth/rho)
# Since ur=i*y1, (16) gives xi=-y1/(ut*wd). Keep substitutions in
# the real symbolic alphabet; Complex{Num}(x,0) nests an opaque complex node.
const KIN=Dict(xi=>-R3.y1/(ut*wd),R3.cbar(xi)=>-R3.y1c/(ut*wd))

# Real field represented by its two Fourier partners. Four equally spaced phases
# integrate the quadratic coefficient exactly (only harmonics -2:2 are present).
field(a,z)=a*z+R3.cbar(a)/z
coefficient2(e)=R3.subs(R3.dsym(R3.dsym(e,lam),lam),Dict(lam=>0))/2
average2(f)=sum(coefficient2(f(z)) for z in (1,im,-1,-im))/4
stress(z)=measure*(1+lam*field(logmeasure,z))*(enth+lam*field(de+dp,z))*
    (lam*field(ur,z)+lam^2*ur2)*(u_phi+lam*field(uphi,z))
baryons(z)=measure*(1+lam*field(logmeasure,z))*(rho+lam*field(drho,z))*
    (lam*field(ur,z)+lam^2*ur2)
enthalpy_current(z)=measure*(1+lam*field(logmeasure,z))*(enth+lam*field(de+dp,z))*
    (lam*field(ur,z)+lam^2*ur2)
flux(z)=measure*(u_phi*field(dp,z)*field(ur,z)+enth*field(ur,z)*field(uphi,z))

@testset "Conserved baryon AMF, general closure, exact symbolic" begin
    # All averages here are normalized by 2pi; multiply every current by 2pi
    # for the package convention (4pi per dynamic reality pair).
    total=average2(stress)
    mdot=average2(baryons)
    f72=2measure*R3.re_part(u_phi*dp*R3.cbar(ur)+enth*ur*R3.cbar(uphi))
    gap=total-hspec*u_phi*mdot-f72
    displacement=2measure*u_phi*gradient*R3.re_part(xi*R3.cbar(ur))
    residuals=[
        "relation19" => (drho+xi*rho_r)/rho-(de+xi*eps_r)/enth,
        "quadratic_stress_expansion" => total-measure*enth*u_phi*ur2-
            2measure*R3.re_part(enth*ur*R3.cbar(uphi)+u_phi*(de+dp+enth*logmeasure)*R3.cbar(ur)),
        "baryon_current_expansion" => mdot-measure*rho*ur2-
            2measure*R3.re_part((drho+rho*logmeasure)*R3.cbar(ur)),
        "pre_average_gradient_term" => gap+displacement,
        "xi_ur_average" => R3.subs(R3.re_part(xi*R3.cbar(ur)),KIN),
        "stress_equals_h_u_phi_Mdot_plus_72a" => R3.subs(gap,KIN),
        "mean_baryon_reduces_to_deps" => R3.subs(hspec*u_phi*mdot-
            measure*enth*u_phi*ur2-2measure*u_phi*R3.re_part((de+enth*logmeasure)*R3.cbar(ur)),KIN),
        "unknown_second_order_mean_cancels" => R3.dsym(gap,ur2),
        "enthalpy_subtraction_gap_is_minus_pressure_current" =>
            total-u_phi*average2(enthalpy_current)-f72+2measure*u_phi*R3.re_part(dp*R3.cbar(ur))]
    open(joinpath(OUT,"symbolic_residuals.tsv"),"w") do io
        println(io,"identity\treal_residual\timag_residual")
        for (name,ex) in residuals
            ex=R3.subs(ex,H)
            ok=R3._zero_num(Num(real(ex))) && R3._zero_num(Num(imag(ex)))
            @test ok
            println(io,name,'\t',ok ? "0\t0" : string(ex))
            println("RESIDUAL ",name," = ",ok ? "0" : ex)
        end
    end
    # Two different erroneous omissions must not be conflated: deleting the
    # xi-gradient from (19) violates that relation; omitting PHASE AVERAGING
    # leaves a nonzero pointwise term. Once correctly averaged, that term is
    # identically zero, so dropping it cannot be killed by an averaged-flux test.
    mutants=[
        "enthalpy_instead_of_baryon_current" => total-u_phi*average2(enthalpy_current)-f72,
        "drop_xi_gradient_from_relation19" => (de/hspec+xi*rho_r)/rho-(de+xi*eps_r)/enth,
        "drop_phase_average" => coefficient2(stress(1))-hspec*u_phi*coefficient2(baryons(1))-flux(1)]
    point=Dict{Any,Any}(R3.eps=>5//3,R3.p=>2//7,rho=>11//5,rho_r=>3//13,
        eps_r=>-7//11,measure=>17//3,u_phi=>5//7,ur2=>2//9,ut=>7//5,wd=>3//11)
    for (a,b,z) in ((R3.y0,R3.y0c,2//3+im*4//7),(R3.y1,R3.y1c,5//7+im*2//9),
                    (R3.htt,R3.httc,3//5+im*1//11),(R3.hrphi,R3.hrphic,2//11-im*3//7),
                    (R3.hrr,R3.hrrc,1//13+im*2//17))
        point[a]=z; point[b]=conj(z)
    end
    open(joinpath(OUT,"mutant_residuals.tsv"),"w") do io
        println(io,"mutant\texact_nonzero_witness")
        for (name,ex) in mutants
            value=R3.exact_eval(R3.subs(R3.subs(ex,KIN),H),point)
            @test !iszero(value)
            println(io,name,'\t',value)
            println("MUTANT ",name," REJECTED: ",value)
        end
    end
end
println("[PASS] general closure baryon-current identity: 9 zero residuals; 3 nonzero mutants")
println("Assumptions: u^r_0=0; stationary axisymmetric background; conserved baryons; relation (19);")
println("real nonzero u^t omega_d; periodic averaging; arbitrary complex amplitudes and pressure closure.")
println("Landau continuation is a runtime convention, not literal conjugation of complex omega_d.")
