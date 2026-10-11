@testset "C-LI-N coefficients and source basis" begin
    @test length(CERTS)==46
    @test count(c->c.modality==:symbolic,CERTS)==18
    @test all(c->c.pass,CERTS)
    # Negative control independent of a sampled response: source curvature has
    # unit coefficient on the left, so reversing its sign must fail.
    e=FRDerivation.fr22_master_residuals()
    @test !FRDerivation.certify_symbolic("wrong Phi curvature","negative control",
        Symbolics.derivative(e.master,e.symbols.P2)+1;verbose=false).pass
end
@testset "C-LI-N non-scalar metric components are post-Newtonian" begin
    # Scope of the scalar source: every other metric component, at its PN order
    # (spatial s^2, mixed s^3, as in R3 newtonian_env) and with a random rational
    # coefficient, must enter each scaled R3 row beyond that row's Newtonian order.
    R=LINewtonian.R3; ex=LINewtonian.expressions(); rng=LINewtonian.MersenneTwister(2608191)
    comps=((R.hrr,2),(R.hrr_r,2),(R.hzz,2),(R.hzz_r,2),(R.hphiphi,2),(R.hphiphi_r,2),
           (R.hrphi,2),(R.hrphi_r,2),(R.htphi,3),(R.htphi_r,3),(R.htr,3),(R.htr_r,3))
    nonzero()=R.RatB(rand(rng,1:97)*rand(rng,(-1,1)),rand(rng,1:89))
    for (row,order) in zip(R.subs.(ex.raw,Ref(ex.scaling)),(2,1,0,1,1)), j in 1:3, (x,power) in comps
        env=Dict{Any,R.LSeries}(Symbolics.value(v)=>R.series_const(R.RatB(rand(rng,1:997),rand(rng,1:991)))
            for v in ex.independent)
        env[Symbolics.value(ex.mN)]=R.series_const(R.RatB(rand(rng,1:9)))
        env[Symbolics.value(ex.sN)]=R.series_mono(R.RatB(1),1)
        for (y,_) in comps; env[Symbolics.value(y)]=R.series_zero(); end
        for y in (R.y0,R.y1,R.htt,R.htt_r); env[Symbolics.value(y)]=R.series_zero(); end
        env[Symbolics.value(x)]=R.series_mono(nonzero(),power)
        d=R.eval_series(row,env)
        @test all(a->R._lim(a)>order && (R.lorder(a)===nothing || R.lorder(a)>order),
                  d isa R.LSeries ? (d,) : d)
    end
end
include("li_profiles.jl")
