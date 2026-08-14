# Characterization of the interpolation library — ANALYTIC polynomial-exactness oracle (NO legacy dep).
#
# Each primitive is exact on the polynomial degree it interpolates; the expected value is the analytic
# polynomial / its analytic derivative, evaluated at OFF-NODE points (never the interpolant's output):
#   - not-a-knot cubic spline   exact on cubics  → reproduces p3 + p3'/p3''/p3''' to machine eps
#   - bicubic / tricubic        exact on bi/tri-cubic surfaces
#   - GslBicubic                exact on bilinear data / at nodes; natural derivs exact on lines
#   - cubic_hermite             reproduces the cubic through (f,f') endpoints
#   - fornberg_weights          weights·samples == analytic derivative; Richardson FD order
#   - barycentric / local       exact to its node degree
#   - cheb_diffmatrix           spectral derivative exact on polynomials of degree ≤ N
#   - fd1c/fd2c                 exact on cubics
#   - not-a-knot boundary identity: 3rd-derivative continuity across the 2nd/2nd-last knot

using NumericsCore
using NumericsCore: CubicSplineInterpolant, BicubicNotAKnot, TricubicNotAKnot, GslBicubic,
    natural_spline_node_derivs, gsl_bsearch, cubic_hermite, fornberg_weights, highorder_deriv,
    bary_weights, bary_eval, bary_interp, local_interp, cheb_nodes, cheb_diffmatrix, cheb_deriv,
    cheb_smooth_deriv, fd1c, fd2c

@testset "interp" begin
    @testset "CubicSplineInterpolant exact on cubics (value + derivatives)" begin
        p(x)  = 0.7 + 1.3x - 0.4x^2 + 0.25x^3
        dp(x) = 1.3 - 0.8x + 0.75x^2
        d2p(x) = -0.8 + 1.5x
        d3p(x) = 1.5
        t = [0.0, 0.6, 1.1, 1.9, 2.7, 3.5, 4.2]   # non-uniform
        s = CubicSplineInterpolant(t, p.(t))
        for xq in (0.3, 0.85, 1.5, 2.3, 3.9, 4.0)
            @test s(xq, 1) ≈ p(xq) atol=1e-11
            @test s(xq, 1; deriv_order=1) ≈ dp(xq) atol=1e-10
            @test s(xq, 1; deriv_order=2) ≈ d2p(xq) atol=1e-9
            @test s(xq, 1; deriv_order=3) ≈ d3p(xq) atol=1e-8
        end
    end

    @testset "three-node not-a-knot uses the unique parabola" begin
        q(x) = 0.7 - 1.2x + 0.45x^2
        dq(x) = -1.2 + 0.9x
        t = [0.0, 0.8, 2.3]
        s = CubicSplineInterpolant(t, q.(t))
        for xq in (0.2, 1.1, 2.0)
            @test s(xq, 1) ≈ q(xq) atol = 1e-13
            @test s(xq, 1; deriv_order = 1) ≈ dq(xq) atol = 1e-13
            @test s(xq, 1; deriv_order = 2) ≈ 0.9 atol = 1e-13
            @test iszero(s(xq, 1; deriv_order = 3))
        end

        w = [0.0, 0.7, 1.9]
        u = [0.0, 1.1, 2.4]
        f(x, y) = 1 + x - 2y + 0.3x^2 + 0.2x*y - 0.4y^2
        bs = BicubicNotAKnot(w, u, [f(x, y) for x in w, y in u])
        @test bs(0.4, 1.7) ≈ f(0.4, 1.7) atol = 1e-12

        tmix = [1.0, 1.0 + 1e-8, 1.0 + 2e-8]
        smix = CubicSplineInterpolant(tmix, Float32[0, 1, 4])
        @test eltype(smix.t) === Float64
        @test smix.t == tmix
        @test [smix(x, 1) for x in tmix] ≈ [0.0, 1.0, 4.0] atol = 1e-12
    end

    @testset "not-a-knot boundary identity (3rd-deriv continuity at 2nd / 2nd-last knot)" begin
        t = collect(0.0:0.5:4.0)
        s = CubicSplineInterpolant(t, sin.(t))            # non-cubic data
        # not-a-knot ⇒ the first two pieces are ONE cubic ⇒ equal 3rd derivative 6·c3
        @test s.c3[1, 1] ≈ s.c3[1, 2] atol=1e-12
        n = length(t)
        @test s.c3[1, n-2] ≈ s.c3[1, n-1] atol=1e-12
    end

    @testset "CubicSplineInterpolant multi-curve + BigFloat genericity" begin
        t = [0.0, 1.0, 2.0, 3.0, 4.0]
        Y = vcat((1.0 .+ t')[1, :]', (t .^ 2)')        # two curves
        s = CubicSplineInterpolant(t, Y)
        @test s(1.5, 1) ≈ 1 + 1.5 atol=1e-11
        @test s(1.5, 2) ≈ 1.5^2 atol=1e-11
        setprecision(BigFloat, 200) do
            tb = big.([0.0, 1.0, 2.0, 3.0, 4.0])
            pb(x) = big(1) + 2x - x^2 + x^3 / 3
            sb = CubicSplineInterpolant(tb, pb.(tb))
            @test sb(big(1.5)) isa Vector{BigFloat}
            @test abs(sb(big(1.5))[1] - pb(big(1.5))) < big(2.0)^-150
        end
    end

    @testset "BicubicNotAKnot exact on bicubic surfaces" begin
        f(x, y) = 1.0 + 2x - 0.5y + x*y - 0.3x^2*y + 0.2x^3 - 0.1y^3 + 0.05x^2*y^2 + 0.07x^3*y^3
        w = [0.0, 0.7, 1.3, 2.1, 2.9, 3.6]
        u = [0.0, 0.5, 1.2, 1.8, 2.6, 3.3]
        vals = [f(wi, uj) for wi in w, uj in u]
        bs = BicubicNotAKnot(w, u, vals)
        for (wq, uq) in ((0.4, 0.9), (1.0, 2.0), (2.5, 1.5), (3.1, 3.0))
            @test bs(wq, uq) ≈ f(wq, uq) atol=1e-9
        end
    end

    @testset "TricubicNotAKnot exact on tricubic surfaces" begin
        g(x, y, z) = 0.5 + x - 0.3y + 0.2z + x*y*z - 0.1x^3 + 0.05y^3 - 0.02z^3 + 0.03x^2*y*z
        u = [0.0, 0.6, 1.4, 2.2, 3.0]
        w = [0.0, 0.7, 1.5, 2.3, 3.1]
        z = [0.0, 0.8, 1.6, 2.4, 3.2]
        V = [g(ui, wj, zk) for ui in u, wj in w, zk in z]
        s = TricubicNotAKnot(u, w, z, V)
        for (uq, wq, zq) in ((0.3, 0.9, 1.1), (1.7, 2.0, 2.5), (2.8, 1.2, 0.4))
            @test s(uq, wq, zq) ≈ g(uq, wq, zq) atol=1e-9
        end
    end

    @testset "GslBicubic exact on bilinear data + at nodes; natural derivs exact on lines" begin
        # natural cubic spline derivs are the exact slope for linear data
        t = [0.0, 1.0, 2.5, 4.0, 5.0]
        line(x) = 3.0 - 1.5x
        d = natural_spline_node_derivs(t, line.(t))
        @test all(abs.(d .- (-1.5)) .< 1e-12)
        # bilinear surface exactly reproduced off-node
        bil(x, y) = 2.0 + 3x - y + 0.5x*y
        x = [0.0, 1.0, 2.0, 3.0]; y = [0.0, 1.0, 2.0, 3.0]
        Z = [bil(xi, yj) for xi in x, yj in y]
        gb = GslBicubic(x, y, Z)
        for (xq, yq) in ((0.5, 0.5), (1.5, 2.5), (2.7, 0.3))
            @test gb(xq, yq) ≈ bil(xq, yq) atol=1e-12
        end
        # at-node reproduction for a generic surface
        surf(x, y) = sin(x) + cos(y) + x*y
        Z2 = [surf(xi, yj) for xi in x, yj in y]
        gb2 = GslBicubic(x, y, Z2)
        @test gb2(2.0, 1.0) ≈ surf(2.0, 1.0) atol=1e-12
        # gsl_bsearch bracket index
        arr = [0.0, 1.0, 2.0, 3.0]
        @test gsl_bsearch(arr, 1.5) == 2
        @test gsl_bsearch(arr, -1.0) == 1
        @test gsl_bsearch(arr, 9.0) == 3
    end

    @testset "cubic_hermite reproduces the cubic through (f,f') endpoints" begin
        p(x)  = 1.0 - 0.5x + 2x^2 - 0.3x^3
        dp(x) = -0.5 + 4x - 0.9x^2
        @test cubic_hermite(1.0, 2.0, p(1.0), p(2.0), dp(1.0), dp(2.0), 1.37) ≈ p(1.37) atol=1e-13
        x = [0.0, 1.0, 2.0, 3.0]
        @test cubic_hermite(x, p.(x), dp.(x), 2.4) ≈ p(2.4) atol=1e-13
    end

    @testset "fornberg_weights == analytic derivative; Richardson FD order" begin
        nodes = [-2.0, -1.0, 0.0, 1.0, 2.0]
        p(x)  = 2.0 + x - 0.5x^2 + 0.25x^3 + 0.1x^4
        dp(x) = 1.0 - x + 0.75x^2 + 0.4x^3
        d2p(x) = -1.0 + 1.5x + 1.2x^2
        z = 0.3
        C = fornberg_weights(z, nodes, 2)
        d1 = sum(C[i, 2] * p(nodes[i]) for i in 1:5)
        d2 = sum(C[i, 3] * p(nodes[i]) for i in 1:5)
        @test d1 ≈ dp(z) atol=1e-10       # exact for deg ≤ 4
        @test d2 ≈ d2p(z) atol=1e-10
        # Richardson: a centred 2nd-order stencil error halves^2 (×1/4) when h halves on a non-poly fn
        fn(x) = exp(0.7x)
        deriv_h(h) = (fn(0.4 + h) - fn(0.4 - h)) / (2h)
        e1 = abs(deriv_h(0.1) - 0.7 * exp(0.7 * 0.4))
        e2 = abs(deriv_h(0.05) - 0.7 * exp(0.7 * 0.4))
        @test e1 / e2 ≈ 4 rtol=0.05       # 2nd-order convergence
    end

    @testset "highorder_deriv exact on polynomials over a grid" begin
        rgrid = collect(0.0:0.1:3.0)
        p(x)  = 1.0 + 0.5x - 0.2x^2 + 0.05x^3
        dp(x) = 0.5 - 0.4x + 0.15x^2
        f = p.(rgrid)
        for k in (5, 12, 20, 28)
            @test highorder_deriv(f, rgrid, k; npts=7, order=1) ≈ dp(rgrid[k]) atol=1e-9
        end
    end

    @testset "barycentric / local interpolation exact to node degree" begin
        nodes = [0.0, 0.5, 1.3, 2.1, 2.8, 3.5]          # 6 nodes ⇒ exact to degree 5
        p(x) = 1.0 - x + 0.5x^2 - 0.2x^3 + 0.1x^4 - 0.03x^5
        f = p.(nodes)
        w = bary_weights(nodes)
        for xq in (0.2, 1.0, 1.7, 3.0)
            @test bary_eval(nodes, w, f, xq) ≈ p(xq) atol=1e-10
            @test bary_interp(xq, nodes, f) ≈ p(xq) atol=1e-10
        end
        # local stencil exact on a low-degree poly within the window
        eqgrid = collect(0.0:0.25:5.0)
        q(x) = 2.0 + x - 0.3x^2 + 0.04x^3
        fe = q.(eqgrid)
        for xq in (0.6, 2.4, 4.1)
            @test local_interp(xq, eqgrid, fe; npts=8) ≈ q(xq) atol=1e-9
        end
        # complex-valued barycentric
        fc = (p.(nodes)) .+ im .* (q.(nodes))
        @test bary_eval(nodes, w, fc, 1.7) ≈ p(1.7) + im * q(1.7) atol=1e-9
    end

    @testset "cheb_diffmatrix spectral-exact on polynomials of degree ≤ N" begin
        N = 12; a = 0.0; b = 2.0
        D, x = cheb_diffmatrix(Float64, N, a, b)
        p(t)  = 1.0 - 0.5t + 0.3t^2 - 0.1t^3 + 0.05t^4
        dp(t) = -0.5 + 0.6t - 0.3t^2 + 0.2t^3
        df = D * p.(x)
        @test maximum(abs.(df .- dp.(x))) < 1e-9
        # cheb_nodes count + endpoints
        @test length(x) == N + 1
        @test x[1] ≈ b && x[end] ≈ a
    end

    @testset "cheb_deriv / cheb_smooth_deriv on smooth fields" begin
        rgrid = collect(0.0:0.05:3.0)
        f = sin.(rgrid)
        dnum = cheb_deriv(f, rgrid; N=40)
        @test maximum(abs.(dnum .- cos.(rgrid))[5:end-5]) < 1e-8
        # least-squares smoothing deriv is exact for a low-degree polynomial
        q(x)  = 1.0 + 0.5x - 0.2x^2 + 0.05x^3
        dq(x) = 0.5 - 0.4x + 0.15x^2
        ds = cheb_smooth_deriv(q.(rgrid), rgrid; deg=6)
        @test maximum(abs.(ds .- dq.(rgrid))) < 1e-9
    end

    @testset "fd1c / fd2c exact on cubics" begin
        p(x)   = 0.4 + 1.1x - 0.6x^2 + 0.3x^3
        dp(x)  = 1.1 - 1.2x + 0.9x^2
        d2p(x) = -1.2 + 1.8x
        for (r, h) in ((1.3, 0.1), (2.7, 0.05), (0.6, 0.2))
            @test fd1c(p, r, h) ≈ dp(r) atol=1e-10
            @test fd2c(p, r, h) ≈ d2p(r) atol=1e-9
        end
    end
end

# =====================================================================================================
# WL-faithful moving-window / tensor-Lagrange + log-Cheb interpolants
# GROUNDING: analytic — Lagrange/moving-window is EXACT on polynomials of its degree (independent oracle
# = the sampled polynomial + its analytic derivatives); log-Cheb is EXACT when log-mag is a low-degree
# polynomial in the Cheb variable and phase is affine (independent oracle = the closed-form R, R', R'').
# (WaSABI + KerrLorenzCirc are heavy READ-literal sources, NOT live-dep'd.)
# =====================================================================================================
using NumericsCore: Interp1WL, interp1_movingwindow, Interp2Hermite, interp2_hermite, partial
using NumericsCore: Interp3Lagrange, interp3_lagrange
using NumericsCore: LogChebInterp, log_cheb_interp, BigLogRInterp, biglogr_interp

@testset "interp_wl_logcheb" begin

    @testset "Interp1WL order-8 exact on deg≤8 (value + 1st + 2nd deriv)" begin
        c = [0.4, -1.1, 0.7, 0.05, -0.02, 0.011, -0.003, 0.0007, -0.00013]  # deg-8 coeffs
        P(x)   = sum(c[k] * x^(k-1) for k in 1:9)
        dP(x)  = sum((k-1) * c[k] * x^(k-2) for k in 2:9)
        d2P(x) = sum((k-1)*(k-2) * c[k] * x^(k-3) for k in 3:9)
        # uniform grid
        xu = collect(0.0:0.5:8.0)
        itp = interp1_movingwindow(xu, P.(xu); order=8)
        for xq in (0.3, 2.7, 4.1, 6.9, 7.8)
            @test itp(xq) ≈ P(xq) rtol=1e-10
            @test itp(xq; deriv_order=1) ≈ dP(xq) rtol=1e-9
            @test itp(xq; deriv_order=2) ≈ d2P(xq) rtol=1e-8
        end
        # steep log-spaced (non-uniform) grid — WL stencil boundary invariants
        xl = exp.(range(log(1.0), log(400.0); length=20))
        itl = interp1_movingwindow(xl, P.(xl); order=8)
        for xq in (1.5, 12.0, 90.0, 380.0)
            @test itl(xq) ≈ P(xq) rtol=1e-8
            @test itl(xq; deriv_order=1) ≈ dP(xq) rtol=1e-7
        end
        # complex samples
        g(x) = P(x) + im * (0.2 - 0.05x)
        itc = interp1_movingwindow(xu, g.(xu); order=8)
        @test itc(3.3) ≈ g(3.3) rtol=1e-10
    end

    @testset "Interp1WL window stays in-range / interval-centered" begin
        x = collect(0.0:1.0:12.0)
        itp = interp1_movingwindow(x, x.^2; order=6)   # deg-2 ⊂ deg-6 ⇒ exact everywhere
        for xq in (0.1, 0.9, 5.5, 11.4, 11.95)          # incl. both boundaries
            @test itp(xq) ≈ xq^2 rtol=1e-12
            @test itp(xq; deriv_order=1) ≈ 2xq rtol=1e-11
        end
    end

    @testset "Interp2Hermite exact on tensor bidegree + partials" begin
        xs = collect(range(-1.0, 1.0; length=5))   # per-axis degree ≤ 4
        ys = collect(range(0.0, 2.0; length=4))     # per-axis degree ≤ 3
        F(x, y)   = (1 + 2x - 0.5x^2 + 0.1x^3) * (0.3 - y + 0.4y^2)
        Fx(x, y)  = (2 - x + 0.3x^2) * (0.3 - y + 0.4y^2)
        Fy(x, y)  = (1 + 2x - 0.5x^2 + 0.1x^3) * (-1 + 0.8y)
        Fxy(x, y) = (2 - x + 0.3x^2) * (-1 + 0.8y)
        grid = [([x, y], F(x, y)) for x in xs for y in ys]
        itp = interp2_hermite(grid)
        for xq in (-0.7, 0.2, 0.85), yq in (0.4, 1.1, 1.7)
            @test itp(xq, yq) ≈ F(xq, yq) rtol=1e-11
            @test partial(itp, xq, yq, 1, 0) ≈ Fx(xq, yq) rtol=1e-10
            @test partial(itp, xq, yq, 0, 1) ≈ Fy(xq, yq) rtol=1e-10
            @test partial(itp, xq, yq, 1, 1) ≈ Fxy(xq, yq) rtol=1e-9
        end
    end

    @testset "Interp3Lagrange exact on tensor tridegree" begin
        xs = collect(range(0.0, 1.0; length=4))
        ys = collect(range(-1.0, 1.0; length=3))
        zs = collect(range(0.0, 2.0; length=3))
        G(x, y, z) = (1 - x + 0.5x^2) * (0.2 + y) * (1 + 0.3z)
        grid = [([x, y, z], G(x, y, z)) for x in xs for y in ys for z in zs]
        itp = interp3_lagrange(grid)
        for xq in (0.15, 0.8), yq in (-0.6, 0.7), zq in (0.5, 1.6)
            @test itp(xq, yq, zq) ≈ G(xq, yq, zq) rtol=1e-11
        end
    end

    @testset "LogChebInterp exact: log-mag polynomial + affine phase" begin
        # R(r) = exp(P(r)) · exp(i(α r + β)),  P a degree-4 polynomial ⇒ log|R| poly, phase affine
        α, β = 0.35, 0.7
        P(r)   = -0.3 + 0.4r - 0.15r^2 + 0.03r^3 - 0.002r^4
        dP(r)  = 0.4 - 0.3r + 0.09r^2 - 0.008r^3
        d2P(r) = -0.3 + 0.18r - 0.024r^2
        R(r)   = exp(P(r)) * exp(im * (α * r + β))
        Lp(r)  = dP(r) + im * α                       # (log R)'
        Lpp(r) = d2P(r)                               # (log R)''
        Rp(r)  = R(r) * Lp(r)
        Rpp(r) = R(r) * (Lpp(r) + Lp(r)^2)
        ci = log_cheb_interp(R, 1.0, 6.0; N=24)
        for r in (1.3, 2.6, 4.0, 5.7)
            @test ci(r) ≈ R(r) rtol=1e-10
            @test ci(r; deriv_order=1) ≈ Rp(r) rtol=1e-9
            @test ci(r; deriv_order=2) ≈ Rpp(r) rtol=1e-8
        end
    end

    @testset "BigLogRInterp exact in log-r + phase-unwrap round-trip" begin
        # R(r) = exp(P(u))·exp(i(α u + β)),  u=log r,  P degree-4 ⇒ exact in the log-r Cheb variable
        α, β = 1.5, 0.3
        P(u)   = 0.5 + 0.8u - 0.2u^2 + 0.05u^3 - 0.004u^4
        dP(u)  = 0.8 - 0.4u + 0.15u^2 - 0.016u^3
        d2P(u) = -0.4 + 0.3u - 0.048u^2
        R(r)   = (u = log(r); exp(P(u)) * exp(im * (α * u + β)))
        Lu(u)  = dP(u) + im * α
        Luu(u) = d2P(u)
        # closed-form r-derivatives via chain rule
        Rp(r)  = (u = log(r); R(r) * Lu(u) / r)
        Rpp(r) = (u = log(r); R(r) * (Luu(u) + Lu(u)^2 - Lu(u)) / r^2)
        ci = biglogr_interp(R, 1.0, 500.0; N=28)
        for r in (2.0, 15.0, 120.0, 480.0)
            @test ci(r) ≈ R(r) rtol=1e-9
            @test ci(r; deriv_order=1) ≈ Rp(r) rtol=1e-8
            @test ci(r; deriv_order=2) ≈ Rpp(r) rtol=1e-7
        end
        # phase-unwrap + exp-reconstruction round-trip: |R| and arg recovered
        for r in (3.0, 50.0, 300.0)
            @test abs(ci(r)) ≈ abs(R(r)) rtol=1e-9
        end
    end

    @testset "log-Cheb BigFloat genericity (setprecision 256)" begin
        setprecision(BigFloat, 256) do
            α = BigFloat(2) / 5
            P(u)  = BigFloat(1)/3 + BigFloat(3)/7 * u - BigFloat(1)/9 * u^2
            R(r)  = (u = log(r); exp(P(u)) * exp(im * α * u))
            Lu(u) = (BigFloat(3)/7 - BigFloat(2)/9 * u) + im * α
            Rp(r) = (u = log(r); R(r) * Lu(u) / r)
            ci = biglogr_interp(R, BigFloat(1), BigFloat(100); N=20)
            r = BigFloat(37)
            @test ci(r) isa Complex{BigFloat}
            @test abs(ci(r) - R(r)) < big(2.0)^-120 * abs(R(r))
            @test abs(ci(r; deriv_order=1) - Rp(r)) < big(2.0)^-110 * abs(Rp(r))
        end
    end
end
