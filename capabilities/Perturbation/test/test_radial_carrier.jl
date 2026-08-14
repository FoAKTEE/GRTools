# Callable-dispatch grounding for the TeukolskyRadialFunction carrier R(r). The two SUPPORTED
# backends are real + grounded (no "method not built" error): :Static routes to _static_radial
# (already vs REF_RADSTATIC) and :MST routes to _mst_radial_in/_up (already vs REF_RADIAL). This
# file locks the public callable surface: R(r) == the internal it dispatches to, and (BHPToolkit
# rescaled) the oracle table. The :HeunC/:NumericalIntegration/:SasakiNakamura backends remain an
# honest [HOLE] (firewall forbids an OrdinaryDiffEq edge, no oracle this round) — @test_throws locks
# that boundary so it is never faked green.

@isdefined(REF_RADIAL) || include("ref_teuk_radial.jl")
@isdefined(REF_RADSTATIC) || include("ref_teuk_radial_static.jl")

@testset "Radial: carrier R(r) :MST dispatch == internal _mst_radial + BHPToolkit oracle" begin
    for r in REF_RADIAL
        R = Perturbation.TeukolskyRadial(r.s, r.l, r.m, r.a, r.w)
        @test R.In.method === :MST && R.Up.method === :MST
        for (j, rr) in enumerate(r.r)
            # callable dispatches exactly to the internal MST evaluator (norm = R.norm)
            @test R.In(rr) === Perturbation._mst_radial_in(R.In.s, R.In.m, R.In.a, R.In.omega,
                                  R.In.nu, R.In.lambda, rr; norm=R.In.norm)
            @test R.Up(rr) === Perturbation._mst_radial_up(R.Up.s, R.Up.m, R.Up.a, R.Up.omega,
                                  R.Up.nu, R.Up.lambda, rr; norm=R.Up.norm)
            # rescale by the transmission amplitude → the norm=1 BHPToolkit table value
            @test isapprox(real(R.In(rr) * R.In.norm), r.inre[j]; rtol=1e-8, atol=1e-9)
            @test isapprox(imag(R.In(rr) * R.In.norm), r.inim[j]; rtol=1e-8, atol=1e-9)
            @test isapprox(real(R.Up(rr) * R.Up.norm), r.upre[j]; rtol=1e-8, atol=1e-9)
            @test isapprox(imag(R.Up(rr) * R.Up.norm), r.upim[j]; rtol=1e-8, atol=1e-9)
        end
    end
end

@testset "Radial: carrier R(r) :Static dispatch == internal _static_radial" begin
    for r in REF_RADSTATIC
        r.ok == 1 || continue
        R = Perturbation.TeukolskyRadial(r.s, r.l, r.m, r.a, 0.0)
        @test R.In.method === :Static && R.Up.method === :Static
        for rr in r.r
            @test R.In(rr) === Perturbation._static_radial(R.In.s, R.In.l, R.In.m, R.In.a,
                                  R.In.bc, rr; norm=R.In.norm)
            @test R.Up(rr) === Perturbation._static_radial(R.Up.s, R.Up.l, R.Up.m, R.Up.a,
                                  R.Up.bc, rr; norm=R.Up.norm)
        end
    end
end

@testset "Radial: carrier R(r) unsupported backend is an honest [HOLE] throw" begin
    # Build a carrier flagged with an unbuilt backend; the callable must throw (NOT silently fake).
    # :HeunC (test_radial_heunc.jl) and :SasakiNakamura (test_radial_sn.jl) are real backends now.
    # :NumericalIntegration remains an honest [HOLE] (firewall forbids an OrdinaryDiffEq edge; no
    # in-firewall replacement built this round).
    base = Perturbation.TeukolskyRadial(-2, 2, 2, 0.5, 0.3).In
    for meth in (:NumericalIntegration,)
        Rh = Perturbation.TeukolskyRadialFunction(base.s, base.l, base.m, base.a, base.omega,
                 base.nu, base.lambda, base.bc, base.norm, base.amplitudes, meth)
        @test Rh.method === meth
        @test_throws ErrorException Rh(6.0)
    end
end
