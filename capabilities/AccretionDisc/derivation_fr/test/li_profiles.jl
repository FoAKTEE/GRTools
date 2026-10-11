@testset "Num and inherited precision through actual Newtonian profiles" begin
    @variables x h p q
    F=FRDerivation; M=F.Modes
    d=F.certificate_disc(h,p,q)
    omega=M.Omega(d,x); op=M.dOmegadr(d,x); kappa=M.kappa(d,x)
    defect=kappa^2-4omega^2-2x*omega*op
    @test F.certify_symbolic("Num kappa identity","NewtonianPlanet.jl:53–64",defect).pass
    @test F.certify_symbolic("Num differentiated kappa identity","NewtonianPlanet.jl:53–64",
        Symbolics.derivative(Symbolics.simplify(defect;expand=true),x)).pass
    # Build numeric evaluators of the unreduced derivative expressions; this
    # separately exercises both sides without Num retaining unevaluated roots.
    kd=Symbolics.build_function(Symbolics.derivative(kappa^2,x),x,h,p,q;expression=Val(false))
    od=Symbolics.build_function(10omega*op+2x*(op^2+omega*Symbolics.derivative(op,x)),x,h,p,q;expression=Val(false))
    for bits in (128,256)
        setprecision(bits) do
            for temp in (0,3//4,1), radius in (3//5,7//3)
                r=BigFloat(radius); disc=F.certificate_disc(BigFloat("0.06"),BigFloat("1.5"),BigFloat(temp))
                o=M.Omega(disc,r); a=M.dOmegadr(disc,r); k=M.kappa(disc,r)
                @test o>0 && k>0 && precision(o)==bits && precision(k)==bits
                @test isapprox(k^2,4o^2+2r*o*a;rtol=128eps(BigFloat))
                # Differentiated identity uses the actual analytic profile
                # derivative, evaluated through Num then at BigFloat values.
                lhs=kd(r,disc.h0,disc.p,disc.q); rhs=od(r,disc.h0,disc.p,disc.q)
                @test isapprox(lhs,rhs;rtol=512eps(BigFloat),atol=512eps(BigFloat)*abs(k^2/r))
                println("PROFILE bits=$bits q=$temp r=$radius kappa_defect=",k^2-4o^2-2r*o*a," derivative_defect=",lhs-rhs)
            end
        end
    end
end
