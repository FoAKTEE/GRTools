# Grounding for the :HeunC confluent-Heun radial backend (heunc.jl). CHARACTERIZED vs the LIGHT live
# oracle `Teukolsky` (closure KerrGeodesics/SWSH/SpecialFunctions/HypGeo — all light):
#   (1) scalar _heunc == Teukolsky._heunc at the legacy 1e-10 tolerance (series + continuation |z|≥1);
#   (2) public TeukolskyRadial(...;method=:HeunC) In/Up VALUE == legacy at the legacy HeunC rtols;
#   (3) HeunC DERIVATIVE (genuine confluent-Heun, analytic) == legacy teukolsky_radial_deriv on a
#       :HeunC carrier (which physically equals the MST derivative) + a central-FD self-check of the
#       green :HeunC value;
#   (4) STRONG cross-check :HeunC R(r) ≈ already-green :MST R(r) (same physical solution) + an
#       r-independent Wronskian W(In,Up).
# Replaces the :HeunC @test_throws that used to live in test_radial_carrier.jl.

include("fixtures/leg_radial_heunc.jl")  # REF_HEUNC_SCALAR / REF_HEUNC_LEG (baked Teukolsky HeunC oracle + inputs, Phase 7-B)
using Perturbation: TeukolskyRadial, teukolsky_radial_deriv

# HeunC(q,α,γ,δ,ε,z) input samples: series (|z|≤0.5), continuation (|z|>1), complex params + z, and a
# far real-negative z (the Teukolsky argument runs the negative real axis). Oracle = live legacy.
const _HEUNC_PTS = [
    (0.1, 0.2, 1.1, 1.2, 0.3, 0.05 + 0.0im),
    (0.5, 1.0, 2.0, 1.5, 0.7, 0.4 + 0.0im),
    (-0.3, 0.8, 1.3, -0.6, 1.1, 0.6 + 0.0im),
    (2.0, -1.0, 0.5, 0.9, -0.4, 0.9 + 0.0im),
    (0.5, 1.0, 2.0, 1.5, 0.7, -8.0 + 0.0im),        # analytic continuation, real-negative
    (0.3 + 0.2im, 0.7 - 0.1im, 1.4, 1.1, 0.5 + 0.3im, 0.7 + 0.5im),  # complex params + complex z
    (1.2, -0.5, 0.9, 1.3, 0.6, -25.0 + 0.0im),      # far continuation
]

@testset "Radial: scalar _heunc == legacy Teukolsky._heunc (series + continuation, baked, captured once)" begin
    worst = 0.0
    for (i, (q, al, ga, de, ep, z)) in enumerate(_HEUNC_PTS)
        got = Perturbation._heunc(q, al, ga, de, ep, z)
        ref = REF_HEUNC_SCALAR[i]                                 # baked LIGHT oracle
        worst = max(worst, abs(got - ref) / abs(ref))
        @test isapprox(got, ref; rtol=1e-10)
    end
    @test worst ≤ 1e-10
end

# Well-conditioned real-ν modes (same shape as the legacy REF_HEUNC_RAD oracle). λ,ν pinned so the
# GrTools + legacy backends are apples-to-apples (identical inputs → the HeunC math is the same port).
const _HEUNC_MODES = [
    (s=-2, l=2, m=2, a=0.5, w=0.3),
    (s=-2, l=2, m=1, a=0.9, w=0.2),
    (s=0,  l=2, m=1, a=0.5, w=0.1),
    (s=-1, l=2, m=1, a=0.0, w=0.2),
]

@testset "Radial: public :HeunC In/Up value == legacy Teukolsky :HeunC (baked, captured once)" begin
    for md in _HEUNC_MODES
        leg = REF_HEUNC_LEG[(md.s, md.l, md.m, md.a, md.w)]
        lam = leg.lam
        nu  = leg.nu
        Rg = TeukolskyRadial(md.s, md.l, md.m, md.a, md.w; method=:HeunC, lambda=lam, nu=nu)
        @test Rg.In.method === :HeunC && Rg.Up.method === :HeunC
        @test Rg.In.amplitudes.Transmission ≈ 1 && Rg.Up.amplitudes.Transmission ≈ 1
        for (j, rr) in enumerate((2.1, 3.0, 5.0, 10.0))
            # same confluent-Heun port on identical (λ,ν) → In agrees to ~machine; Up = Refl·In +
            # Inc·conj is a cancellation-limited difference at large r (identical conditioning in
            # both backends), so its GrTools-vs-legacy tolerance is relaxed with r.
            @test isapprox(Rg.In(rr), leg.vin[j]; rtol=1e-9)
            rtolUp = rr ≤ 3.0 ? 1e-8 : rr ≤ 5.0 ? 1e-7 : 1e-5
            @test isapprox(Rg.Up(rr), leg.vup[j]; rtol=rtolUp)
        end
    end
end

@testset "Radial: :HeunC derivative — analytic vs legacy (MST-equivalent, baked) + FD self-check" begin
    h = 1e-6
    for md in _HEUNC_MODES
        leg = REF_HEUNC_LEG[(md.s, md.l, md.m, md.a, md.w)]
        lam = leg.lam
        nu  = leg.nu
        Rg = TeukolskyRadial(md.s, md.l, md.m, md.a, md.w; method=:HeunC, lambda=lam, nu=nu)
        for (j, rr) in enumerate((3.0, 5.0))
            # legacy teukolsky_radial_deriv on a :HeunC carrier physically returns the MST derivative;
            # the GrTools analytic confluent-Heun derivative is the SAME physical solution → agrees to
            # the shared well-conditioned tolerance.
            @test isapprox(teukolsky_radial_deriv(Rg.In, rr), leg.din[j]; rtol=1e-7)
            rtolUp = rr ≤ 3.0 ? 1e-7 : 1e-5
            @test isapprox(teukolsky_radial_deriv(Rg.Up, rr), leg.dup[j]; rtol=rtolUp)
            # central-FD self-check of the green :HeunC value (oracle-free)
            dInFD = (Rg.In(rr + h) - Rg.In(rr - h)) / (2h)
            @test isapprox(teukolsky_radial_deriv(Rg.In, rr), dInFD; rtol=1e-6, atol=1e-7)
        end
    end
end

@testset "Radial: STRONG cross-check :HeunC ≡ :MST + r-independent Wronskian" begin
    for md in _HEUNC_MODES
        leg = REF_HEUNC_LEG[(md.s, md.l, md.m, md.a, md.w)]
        lam = leg.lam
        nu  = leg.nu
        Rh = TeukolskyRadial(md.s, md.l, md.m, md.a, md.w; method=:HeunC, lambda=lam, nu=nu)
        Rm = TeukolskyRadial(md.s, md.l, md.m, md.a, md.w; method=:MST,   lambda=lam, nu=nu)
        for rr in (3.0, 5.0, 10.0)
            @test isapprox(Rh.In(rr), Rm.In(rr); rtol=1e-7)
            rtolUp = rr ≤ 5.0 ? 1e-6 : 1e-5
            @test isapprox(Rh.Up(rr), Rm.Up(rr); rtol=rtolUp)
        end
        # Wronskian W = In·Up' - In'·Up of the homogeneous pair: the Teukolsky operator
        # Δ R'' + (2r-2)(1+s) R' + … = 0 is self-adjoint with weight Δ^{s+1}, so Δ^{s+1}·W is
        # r-independent. Check it agrees at two radii (a genuine homogeneous-solution invariant).
        wron(R, r) = begin
            Δ = r^2 - 2r + md.a^2
            Δ^(md.s + 1) * (R.In(r) * teukolsky_radial_deriv(R.Up, r) -
                            teukolsky_radial_deriv(R.In, r) * R.Up(r))
        end
        w1 = wron(Rh, 4.0); w2 = wron(Rh, 8.0)
        @test isapprox(w1, w2; rtol=1e-5)
    end
end
