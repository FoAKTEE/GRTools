# test_image_synthesis.jl — EGD::sw-image-synthesis: the mode-sum reconstruction of the spiral wake
# (Figs. 1-2 of arXiv:2608.17003, src/SpiralWaves/image_synthesis.jl). Gates:
#   (a) `mode_interp` is exact (interpolatory) at the mode's own nodes and zero outside its domain;
#   (b) ROUND TRIP: on synthetic modes with known amplitudes, `synthesize` reproduces the closed-form
#       sum f(r,Ψ) = Σ_m w_m Re[f̂_m e^{imΨ}] (w_0 = 1, w_{m≥1} = 2) to 1e-13;
#   (c) a single mode gives 2|f̂| cos(mΨ + arg f̂) — the real-field ×2 convention of (45a)/(73);
#   (d) the polar → Cartesian map and the fiducial image grid of (C1);
#   (e) a REAL solved mode synthesises a finite, m-fold symmetric wake; BigFloat inherits;
#   (f) `synthesize_cartesian` equals `synthesize` at the SAME (r, Ψ) to roundoff and blanks the
#       pixels outside the disc (the white region of Figs. 1-2);
#   (g) `arm_phase` recovers the known winding of an analytic spiral f̂_m = e^{-ikr} (crest at
#       Ψ = kr/m, so dΨ/dr = k/m) and unwraps across many turns;
#   (h) `wake_modes` defaults to the fiducial image domain (C1) and is bitwise thread-independent.
using Test
using AccretionDisc
using AccretionDisc.SpiralWaves

@testset "SpiralWaves image synthesis: δU(r, Ψ) = Σ_m w_m Re[δÛ_m e^{imΨ}] (45a)" begin

    # synthetic "mode-like" objects: anything with m, r and a 2×N amplitude matrix
    _mk(m, r, f1, f2) = (m = m, r = collect(r), y = permutedims(hcat(ComplexF64[f1(x) for x in r],
                                                                    ComplexF64[f2(x) for x in r])))
    rr = collect(range(6.0, 40.0; length = 257))
    mods = [_mk(1, rr, x -> cis(0.3x) / (1 + x), x -> ComplexF64(x, -x^2 / 100)),
            _mk(2, rr, x -> ComplexF64(cos(x / 3), sin(x / 5)), x -> ComplexF64(1 / x, x)),
            _mk(5, rr, x -> exp(-x / 20) * cis(-0.7x), x -> ComplexF64(0.0, x))]

    @testset "(a) mode_interp: exact at the nodes, zero outside the domain" begin
        mi = mode_interp(mods[2])
        @test mi.m == 2
        @test mi.r_lo == rr[1] && mi.r_hi == rr[end]
        worst = maximum(i -> abs(mi(rr[i]) - mods[2].y[1, i]), eachindex(rr))
        @test worst < 1e-13
        @test mi(rr[1] - 1e-9) == 0
        @test mi(rr[end] + 1e-9) == 0
        mi2 = mode_interp(mods[2]; component = 2)
        @test maximum(i -> abs(mi2(rr[i]) - mods[2].y[2, i]), eachindex(rr)) < 1e-13
        @test_throws ArgumentError mode_interp(mods[1]; component = 3)
    end

    @testset "(b) round trip against the closed-form mode sum" begin
        Ψq = collect(range(0, 2π; length = 37)[1:36])
        rq = rr[3:7:end]                                  # the modes' OWN nodes ⇒ interpolation exact
        for comp in (1, 2)
            F = synthesize(mods, rq, Ψq; component = comp)
            @test size(F) == (length(rq), length(Ψq))
            worst = 0.0
            for (i, r) in enumerate(rq), (j, Ψ) in enumerate(Ψq)
                k = findfirst(==(r), rr)
                ref = sum(2 * real(md.y[comp, k] * cis(md.m * Ψ)) for md in mods)
                worst = max(worst, abs(F[i, j] - ref) / max(1.0, abs(ref)))
            end
            @info "image synthesis round trip (component $comp): worst deviation = $worst"
            @test worst < 1e-13
        end
        @test_throws ArgumentError synthesize((), rq, Ψq)
    end

    @testset "(c) the ×2 real-field convention and m = 0" begin
        m = 3
        one_mode = [_mk(m, rr, x -> ComplexF64(0.4, -0.9), x -> ComplexF64(0, 0))]
        Ψq = collect(range(0, 2π; length = 61)[1:60])
        F = synthesize(one_mode, [rr[10]], Ψq)
        f = ComplexF64(0.4, -0.9)
        for (j, Ψ) in enumerate(Ψq)
            @test F[1, j] ≈ 2 * abs(f) * cos(m * Ψ + angle(f)) atol = 1e-13
        end
        # m-fold symmetry: Ψ → Ψ + 2π/m leaves the field invariant
        @test synthesize(one_mode, [rr[10]], [0.37])[1] ≈ synthesize(one_mode, [rr[10]], [0.37 + 2π / m])[1] atol = 1e-12
        # m = 0 is NOT doubled
        zero_mode = [_mk(0, rr, x -> ComplexF64(0.4, -0.9), x -> ComplexF64(0, 0))]
        @test synthesize(zero_mode, [rr[10]], [0.0])[1] ≈ 0.4 atol = 1e-13
    end

    @testset "(d) polar → Cartesian and the fiducial image grid" begin
        rq = [6.0, 10.0]; Ψq = [0.0, π / 2, π]
        X, Y = polar_to_cartesian(rq, Ψq)
        @test size(X) == size(Y) == (2, 3)
        @test X[1, 1] ≈ 6.0 && Y[1, 1] ≈ 0.0
        @test X[2, 2] ≈ 0.0 atol = 1e-14
        @test Y[2, 2] ≈ 10.0
        @test X[1, 3] ≈ -6.0
        bg = DiscBackground(; χ = 0.0, aEMRI = 15.0, h0 = 0.05)
        rg, Ψg = image_grid(bg; nr = 64, nΨ = 128)
        @test length(rg) == 64 && length(Ψg) == 128
        @test rg[1] ≈ r_isco(bg) && rg[end] ≈ 5 * bg.aEMRI
        @test Ψg[1] == 0 && Ψg[end] < 2π                   # periodic, endpoint excluded
        @test_throws ArgumentError image_grid(bg; nr = 1)
    end

    @testset "(e) a real solved mode, and the BigFloat path" begin
        bg = DiscBackground(; χ = 0.0, aEMRI = 15.0, h0 = 0.05, 𝔭 = 0.0, 𝔮 = 0.0)
        sols = [solve_mode(bg, m; grid = build_grid(bg, m; nL = 8, nR = 8, nC = 100, n_anc = 800,
                                                    extent = :fiducial)) for m in 2:4]
        rq, Ψq = image_grid(bg; nr = 48, nΨ = 96)
        F = synthesize(sols, rq, Ψq)
        @test size(F) == (48, 96)
        @test all(isfinite, F)
        @test maximum(abs, F) > 0
        # the wake is a genuine spiral: the azimuth of the maximum moves with radius
        jmax = [argmax(view(F, i, :)) for i in 1:size(F, 1)]
        @test length(unique(jmax)) > 3
        @info "wake image: max |δU| = $(maximum(abs, F)), distinct peak azimuths over $(size(F,1)) radii = $(length(unique(jmax)))"
        setprecision(BigFloat, 256) do
            bgB = DiscBackground{BigFloat}(; χ = big(0), aEMRI = big(15), h0 = big(0.05))
            solB = solve_mode(bgB, 3; grid = build_grid(bgB, 3; nL = 6, nR = 6, nC = 40, n_anc = 300,
                                                        extent = :fiducial), hmodes_N = 512)
            rqB, ΨqB = image_grid(bgB; nr = 8, nΨ = 16)
            FB = synthesize([solB], rqB, ΨqB)
            @test eltype(FB) === BigFloat
            @test all(isfinite, FB)
        end
    end

    @testset "(f) synthesize_cartesian: exact at the pixel, blanked outside the disc" begin
        xs = collect(range(-42.0, 42.0; length = 57))
        ys = collect(range(-42.0, 42.0; length = 53))
        C = synthesize_cartesian(mods, xs, ys)
        @test size(C) == (length(xs), length(ys))
        rlo, rhi = rr[1], rr[end]
        worst = 0.0
        nin = 0
        for (i, x) in enumerate(xs), (j, y) in enumerate(ys)
            r = hypot(x, y)
            if r < rlo || r > rhi
                @test isnan(C[i, j])
            else
                nin += 1
                ref = synthesize(mods, [r], [atan(y, x)])[1]
                worst = max(worst, abs(C[i, j] - ref) / max(1.0, abs(ref)))
            end
        end
        @info "synthesize_cartesian vs synthesize: $(nin) interior pixels, worst rel = $worst"
        @test nin > 100
        @test worst < 1e-14
        # explicit radial window + a finite blank value
        C2 = synthesize_cartesian(mods, xs, ys; r_lo = 20.0, r_hi = 30.0, blank = 0.0)
        for (i, x) in enumerate(xs), (j, y) in enumerate(ys)
            (hypot(x, y) < 20.0 || hypot(x, y) > 30.0) && @test C2[i, j] == 0.0
        end
        @test any(!iszero, C2)
        @test_throws ArgumentError synthesize_cartesian((), xs, ys)
        @test_throws ArgumentError synthesize_cartesian(mods, xs, ys; r_lo = 30.0, r_hi = 20.0)
        # component 2 is honoured
        @test synthesize_cartesian(mods, [25.0], [0.0]; component = 2)[1] ≈
              synthesize(mods, [25.0], [0.0]; component = 2)[1]
    end

    @testset "(g) arm_phase: the winding curve of an analytic spiral" begin
        Ψq = collect(range(0, 2π; length = 721)[1:720])
        # one arm (m = 1): 2 cos(Ψ − k r) ⇒ crest at Ψ = k r, six full turns over the range
        k = 0.8
        rq = collect(range(10.0, 60.0; length = 101))
        one_arm = [_mk(1, rq, x -> cis(-k * x), x -> ComplexF64(0, 0))]
        φ = arm_phase(synthesize(one_arm, rq, Ψq), Ψq)
        @test length(φ) == length(rq)
        slope = (φ[end] - φ[1]) / (rq[end] - rq[1])
        @info "arm_phase slope = $slope (exact k = $k); turns = $((φ[end] - φ[1]) / 2π)"
        @test slope ≈ k rtol = 1e-3
        @test maximum(i -> abs((φ[i] - φ[1]) - k * (rq[i] - rq[1])), eachindex(rq)) < 5e-3
        @test abs((φ[end] - φ[1]) / 2π) > 6       # genuinely unwrapped across several turns
        # a pure m = 3 pattern has m indistinguishable crests: the period is 2π/m
        m2, k2 = 3, 0.8
        rq2 = collect(range(10.0, 30.0; length = 41))
        three = [_mk(m2, rq2, x -> cis(-k2 * x), x -> ComplexF64(0, 0))]
        F3 = synthesize(three, rq2, Ψq)
        φ3 = arm_phase(F3, Ψq; period = 2π / m2)
        @test (φ3[end] - φ3[1]) / (rq2[end] - rq2[1]) ≈ k2 / m2 rtol = 1e-3
        @test_throws DimensionMismatch arm_phase(F3, Ψq[1:end-1])
        @test_throws ArgumentError arm_phase(F3, Ψq; period = 0)
    end

    @testset "(h) wake_modes: fiducial domain by default, thread-independent bits" begin
        bg = DiscBackground(; χ = 0.3, aEMRI = 15.0, h0 = 0.08)
        gk = (nL = 6, nR = 6, nC = 60, n_anc = 500)
        ms = 1:5
        par = wake_modes(bg, ms; grid_kwargs = gk, threads = true)
        ser = wake_modes(bg, ms; grid_kwargs = gk, threads = false)
        @test length(par) == length(ms)
        @test [s.m for s in par] == collect(ms)
        for (a, b) in zip(par, ser)
            @test a.r == b.r && a.y == b.y        # BITWISE
        end
        # the default extent is the image domain (C1): every mode shares [r_ISCO, 5a•]
        for s in par
            @test first(s.r) ≈ r_isco(bg)
            @test last(s.r) ≈ 5 * bg.aEMRI
        end
        # …and it can be overridden
        ad = wake_modes(bg, 3:3; grid_kwargs = merge(gk, (extent = :adaptive,)), threads = false)
        @test first(ad[1].r) > r_isco(bg) || last(ad[1].r) < 5 * bg.aEMRI
        @test_throws ArgumentError wake_modes(bg, Int[])
        @test_throws ArgumentError wake_modes(bg, [0, 1])
    end
end
