using Test
using AccretionDisc

@testset "FR22 opt-in phase-gradient descent" begin
    P = AccretionDisc.Modes
    d = NewtonianDisc(beta=80.0)
    x = loggrid(d, .05, 5., 1200)
    h = exp.(im .* (2 .* x .+ x.^2 ./ 2))
    h1 = hcat(h, im .* (2 .+ x) .* h)
    h2 = conj.(h1)
    ih = hcat(.03 .* x.^2 .+ .02im .* x, .06 .* x .+ .02im)
    Ci, Co = -.4+2.1im, .2+7.1im
    a1, a2 = P.homogeneous_coefficients(h1,h2,ih,Ci,Co,1,length(x))
    baseline = a1 .* h1 .+ a2 .* h2 .+ ih
    saved = copy(baseline)

    for m in (1, 10)
        args = (d,x,m,h1,h2,ih,Ci,Co,1,length(x))
        legacy = phase_gradient_refine(args...; baseline)
        off = phase_gradient_refine(args...; baseline, safeguard=false)
        @test keys(off) == (:mode, :C, :X, :jacobian)
        @test isequal(legacy, off)
        @test reinterpret(UInt64, vec(legacy.mode)) == reinterpret(UInt64, vec(off.mode))
        for max_backtracks in (0, 6)
            guarded = phase_gradient_refine(args...; baseline, safeguard=true, max_backtracks)
            audit = guarded.safeguard
            @test guarded.X == legacy.X
            @test guarded.jacobian == legacy.jacobian
            @test audit.X_before == legacy.X
            @test audit.X_after == phase_error_diagnostic(d,guarded.mode,x,m,1,length(x))
            @test audit.status in (:accepted, :backtracked, :rejected)
            if audit.status == :rejected
                @test guarded.C === (Ci, Co)
                @test isequal(guarded.mode, baseline)
                @test audit.fraction == 0
                @test audit.backtracks == max_backtracks
            else
                @test all(isfinite, audit.X_after)
                @test all(audit.X_after .< audit.X_before)
                @test audit.fraction == 2.0^(-audit.backtracks)
                @test (audit.status == :accepted) == (audit.backtracks == 0)
            end
            m == 1 && (@test guarded.C[1] === Ci)
        end
        @test_throws ArgumentError phase_gradient_refine(args...; baseline,
                                                         safeguard=true, max_backtracks=-1)
    end
    @test baseline == saved

    # A zero peak-to-peak baseline cannot strictly improve. A zero secant must
    # therefore be rejected, with the original coefficients and their combination.
    flat = hcat(ones(ComplexF64,length(x)), im .* real.([P.kwave(d,r,1) for r in x]))
    rejected = phase_gradient_refine(d,x,1,h1,h2,ih,Ci,Co,1,length(x);
                                     baseline=flat, safeguard=true)
    @test rejected.X == (0., 0.)
    @test rejected.safeguard.status == :rejected
    @test rejected.C === (Ci, Co)
    @test isequal(rejected.mode, baseline)
    @test rejected.safeguard.X_after == phase_error_diagnostic(d,baseline,x,1,1,length(x))

    invalid = phase_gradient_refine(d,x,1,h1,h2,ih,Ci,Co,1,length(x);
                                    baseline=zeros(ComplexF64,size(baseline)), safeguard=true)
    @test all(isnan, invalid.X)
    @test invalid.safeguard.status == :rejected
    @test invalid.C === (Ci, Co)
    @test isequal(invalid.mode, baseline)
    mode = AccretionDisc.ForcingModes.NewtonianMode(1,1,1.0)
    @test_throws ArgumentError phase_gradient_refine(d,x,mode,h1,h2,ih,Ci,Co,1,length(x);
                                                     baseline, safeguard=true)
end

@testset "a zero secant leaves X unchanged and is rejected (lander review R2)" begin
    # The zero-baseline cases above are rejected because the coefficients' own combination has a
    # larger X than the zero baseline, not by equality. Here the particular solution satisfies both
    # boundary rows and has epsilon == 0, so a1 = a2 = 0, X == (0, 0) and the secant step is
    # exactly zero: every trial leaves X bit-identical, and only strict descent rejects it.
    P = AccretionDisc.Modes
    d = NewtonianDisc(beta=80.0)
    x = loggrid(d, .05, 5., 1200)
    h = exp.(im .* (2 .* x .+ x.^2 ./ 2))
    h1 = hcat(h, im .* (2 .+ x) .* h)
    h2 = conj.(h1)
    for m in (1, 10)
        flat = hcat(ones(ComplexF64, length(x)), im .* real.([P.kwave(d, r, m) for r in x]))
        Cf = (flat[1, 2] / flat[1, 1], flat[end, 2] / flat[end, 1])
        unguarded = phase_gradient_refine(d, x, m, h1, h2, flat, Cf..., 1, length(x); baseline=flat)
        @test unguarded.X == (0., 0.)
        @test unguarded.C == Cf
        still = phase_gradient_refine(d, x, m, h1, h2, flat, Cf..., 1, length(x); baseline=flat, safeguard=true)
        @test still.safeguard.X_after == (0., 0.)
        @test still.safeguard.status == :rejected
        @test still.C === Cf
        @test isequal(still.mode, flat)
    end
end

@testset "safeguard on real shooting bases: accepted, backtracked, rejected (lander review R1)" begin
    # The synthetic case above reaches only :rejected. These three solver modes reach every branch.
    # Their adaptive boundaries stay at the grid ends (asserted), so the basis is rebuilt from the
    # solver's own steps, and the unguarded step must equal solve_m_mode(...; phgrad_steps=1)
    # bitwise: the test cannot drift from the production basis.
    P = AccretionDisc.Modes
    function shooting_basis(d, m, x, cache)
        kw = (; rtol=1e-8, atol=1e-10)  # the cached solve_m_mode defaults
        z, o = zero(ComplexF64), one(ComplexF64)
        h1 = P.integrate_master(d, x, (o, z), m, 0, cache; kw...)
        h2 = P.integrate_master(d, x, (z, o), m, 0, cache; kw...)
        ih = P.integrate_master(d, x, (z, z), m, 1, cache; kw...)
        n = length(x)
        Ci, Co = P.wkb_boundary_coeffs(d, x, 1, n, m)
        combine(ih) = ((a, b) = P.homogeneous_coefficients(h1, h2, ih, Ci, Co, 1, n); a .* h1 .+ b .* h2 .+ ih)
        f = combine(ih)
        for _ in 1:2  # the two forced re-seeds at corotation
            ico = P.corotation_index(d, x, m)
            ih = P.integrate_master(d, x, (f[ico, 1], f[ico, 2]), m, 1, cache; kw...)
            f = combine(ih)
        end
        (d, x, m, h1, h2, ih, Ci, Co, 1, n), f
    end
    mr = NewtonianDisc(beta=80.0, lap_nphi=800)
    # FR22 setup (p = q = 0, h_p = 0.1, epsilon*h_p, beta = 0, N_phi = 1024, 1000 radii on [0.1, 5]).
    fr = NewtonianDisc(p=0.0, q=0.0, h0=0.1, beta=0.0, soft=0.3*0.1, indirect=0, lap_nphi=1024)
    cases = ((mr, 1, .05, 5., 2000, :backtracked),  # the unguarded step raises both X (obligation mr19-o-secant-m1-descent)
             (mr, 10, .05, 5., 2000, :accepted),
             (fr, 1, .1, 5., 1000, :rejected))       # full step: inner X falls, outer X rises
    for (d, m, rin, rout, n, expected) in cases
        x = loggrid(d, rin, rout, n)
        cache = PotentialCache(d, m, rin, rout; n)
        s = solve_m_mode(d, m, x, cache; phgrad_steps=1, diagnostics=true)
        @test s.indices == (1, n)
        args, baseline = shooting_basis(d, m, x, cache)
        @test isequal(baseline, solve_m_mode(d, m, x, cache))  # the unrefined default solution
        legacy = phase_gradient_refine(args...; baseline)
        @test isequal(legacy.mode, s.mode)
        g = phase_gradient_refine(args...; baseline, safeguard=true)
        a = g.safeguard
        Ci, Co = args[7], args[8]
        X0 = phase_error_diagnostic(d, baseline, x, m, 1, n)
        X_full = phase_error_diagnostic(d, legacy.mode, x, m, 1, n)
        @test a.status == expected
        @test a.X_before == X0
        @test a.X_after == phase_error_diagnostic(d, g.mode, x, m, 1, n)
        m == 1 && (@test g.C[1] === Ci)
        if expected == :rejected
            @test g.C === (Ci, Co)
            @test isequal(g.mode, baseline)
            @test X_full[1] < X0[1] && X_full[2] > X0[2]  # one edge improves, the other worsens
        else
            @test all(a.X_after .< X0)
            @test a.fraction == 2.0^-a.backtracks
            if expected == :accepted
                @test a.backtracks == 0
                @test g.C == legacy.C
                @test isequal(g.mode, legacy.mode)
            else
                @test a.backtracks >= 1
                @test g.C[2] == Co + a.fraction * (legacy.C[2] - Co)
                @test all(X_full .> X0)  # the full MR step ascends at both edges
                reject = phase_gradient_refine(args...; baseline, safeguard=true, max_backtracks=0)
                @test reject.safeguard.status == :rejected
                @test reject.C === (Ci, Co)
                @test isequal(reject.mode, baseline)
            end
        end
    end
end
