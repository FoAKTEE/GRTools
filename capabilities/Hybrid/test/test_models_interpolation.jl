# Characterization of the WL-faithful interpolation leaf (Hybrid.Models, design §4.9).
#
# GROUNDING (no WaSABI live dep):
#   * Lagrange exactness identity (ANALYTIC closed form): an order-`k` moving-window
#     interpolant built from samples of a degree-≤k polynomial reproduces the polynomial AND
#     its 1st/2nd derivatives EXACTLY. Likewise the tensor-product 2D/3D global-Lagrange
#     interpolants reproduce a separable polynomial product and its analytic mixed partials.
#   * WL interval-centered stencil (CAPTURED FIXTURE): for a *non-polynomial* function the
#     window choice is load-bearing. Expected values are `Interpolation[data,
#     InterpolationOrder->8][q]` / `Derivative[n][..][q]` captured via a one-shot wolframscript
#     run on a uniform grid (exp(-x/7) sin x, x=1..20) and a log-spaced grid (1/(1+e^t),
#     x=e^t). A wrong stencil differs at O(1), so rtol 1e-9 robustly certifies the rule.
#   * interp_logexp closed form: itp interpolating the linear log-table (log r0, 2 log r0)
#     gives real(exp(itp(log r0))) == r0^2 exactly.
#   * kerr_grid_coords vs LIVE Spacetime.kerr_isco (light dep): (log(1-a), sqrt(isco(a)/r)).

using Test
using Hybrid.Models
using Spacetime: kerr_isco

@testset "Models.interpolation" begin

    # --- 1D order-8 Lagrange exactness (closed form) ---------------------------------------
    @testset "Interp1 polynomial exactness" begin
        p(x)  = 0.7 - 1.3x + 0.4x^2 - 0.05x^3 + 2e-3x^4 - 3e-4x^5 +
                1e-5x^6 - 4e-7x^7 + 9e-9x^8            # degree 8
        dp(x) = -1.3 + 0.8x - 0.15x^2 + 8e-3x^3 - 1.5e-3x^4 +
                6e-5x^5 - 2.8e-6x^6 + 7.2e-8x^7
        d2p(x)= 0.8 - 0.3x + 0.024x^2 - 6e-3x^3 + 3e-4x^4 -
                1.68e-5x^5 + 5.04e-7x^6
        xs = collect(1.0:12.0)
        itp = interp1_order8(xs, p.(xs))
        for q in (2.4, 6.5, 9.8, 11.2)
            @test itp(q)              ≈ p(q)   rtol=1e-8
            @test derivative(itp, q)  ≈ dp(q)  rtol=1e-7
            @test derivative2(itp, q) ≈ d2p(q) rtol=1e-6
        end
        # pairs-method builds the identical interpolant
        itp2 = interp1_order8([(xs[i], p(xs[i])) for i in eachindex(xs)])
        @test itp2(7.3) ≈ itp(7.3) rtol=0
    end

    # --- WL stencil characterization (captured wolframscript fixture) ----------------------
    @testset "Interp1 WL stencil (uniform grid)" begin
        xu = collect(1.0:20.0)
        yu = @. exp(-xu / 7) * sin(xu)
        itp = interp1_order8(xu, yu)
        # {q => (value, d1, d2)} from Interpolation[...,InterpolationOrder->8]
        fix = Dict(
            2.3  => (0.5377582784283734,  -0.5555332976230904,  -0.40083531951751405),
            5.5  => (-0.321456721124534,    0.368903116837352,    0.22157072272936812),
            11.7 => (-0.14318426851557914,  0.1420660129321463,   0.1050940423538892),
            18.4 => (-0.03139114393229938,  0.06942580163191488,  0.012226355346127135),
        )
        for (q, (v, d1, d2)) in fix
            @test itp(q)              ≈ v  rtol=1e-9
            @test derivative(itp, q)  ≈ d1 rtol=1e-9
            @test derivative2(itp, q) ≈ d2 rtol=1e-8
        end
    end

    @testset "Interp1 WL stencil (log-spaced grid)" begin
        t  = collect(range(0.0, 3.0; length = 16))
        xl = @. exp(t)
        yl = @. 1 / (1 + exp(t))
        itp = interp1_order8(xl, yl)
        fix = Dict(
            1.4  => (0.4166661997690849,   -0.17360703717076967,   0.144717030748148),
            3.9  => (0.2040752110240608,   -0.041626037895598364,  0.017218350638115262),
            9.2  => (0.09803618722811284,  -0.009630358719288705,  0.0018711441468946635),
            17.5 => (0.05349484742672044,  -0.0036863819640954934, -0.00007347403490494961),
        )
        for (q, (v, d1, d2)) in fix
            @test itp(q)              ≈ v  rtol=1e-9
            @test derivative(itp, q)  ≈ d1 rtol=1e-8
            @test derivative2(itp, q) ≈ d2 rtol=1e-7
        end
    end

    # --- log-exp composite (closed form) ---------------------------------------------------
    @testset "interp_logexp" begin
        r0s = collect(2.0:0.5:12.0)
        # table (log r0, log(r0^2)) = (log r0, 2 log r0): itp linear in log-space ⇒ exp = r0^2
        tbl = [(log(r), complex(2 * log(r))) for r in r0s]
        itp = interp1_order8(tbl)
        f = interp_logexp(itp)
        for r in (3.1, 5.7, 9.3)
            @test f(r) ≈ r^2 rtol=1e-10
            @test f(r) ≈ real(exp(itp(log(r))))  # definitional identity
        end
    end

    # --- 2D tensor-product Hermite (global Lagrange) + mixed partials -----------------------
    @testset "Interp2Hermite separable polynomial + partials" begin
        px(x)  = 1.0 + 0.5x - 0.2x^2 + 0.03x^3
        dpx(x) = 0.5 - 0.4x + 0.09x^2
        d2px(x)= -0.4 + 0.18x
        qy(y)  = 2.0 - 0.7y + 0.1y^2
        dqy(y) = -0.7 + 0.2y
        xs = collect(0.0:4.0)      # 5 nodes ⇒ degree-4 global Lagrange (exact for deg ≤4)
        ys = collect(0.0:3.0)      # 4 nodes ⇒ degree-3 global Lagrange (exact for deg ≤3)
        grid = [((x, y), px(x) * qy(y)) for x in xs for y in ys]
        itp = interp2_hermite(grid)
        for (xq, yq) in ((1.3, 2.1), (3.4, 0.6))
            @test itp(xq, yq)            ≈ px(xq) * qy(yq)   rtol=1e-10
            @test partial(itp, xq, yq, 1, 0) ≈ dpx(xq) * qy(yq)  rtol=1e-10
            @test partial(itp, xq, yq, 0, 1) ≈ px(xq) * dqy(yq)  rtol=1e-10
            @test partial(itp, xq, yq, 1, 1) ≈ dpx(xq) * dqy(yq) rtol=1e-10
            @test partial(itp, xq, yq, 2, 0) ≈ d2px(xq) * qy(yq) rtol=1e-9
        end
    end

    # --- 3D tensor-product global Lagrange -------------------------------------------------
    @testset "Interp3Lagrange separable polynomial" begin
        px(x) = 1.0 + 0.5x - 0.2x^2
        qy(y) = 2.0 - 0.3y
        rz(z) = 0.4 + 0.6z + 0.1z^2
        xs = collect(0.0:2.0); ys = collect(0.0:1.0); zs = collect(0.0:2.0)
        grid = [((x, y, z), px(x) * qy(y) * rz(z)) for x in xs for y in ys for z in zs]
        itp = interp3_lagrange(grid)
        for (xq, yq, zq) in ((1.4, 0.7, 1.1), (0.3, 0.2, 1.8))
            @test itp(xq, yq, zq) ≈ px(xq) * qy(yq) * rz(zq) rtol=1e-10
        end
    end

    # --- Kerr abscissa remap vs LIVE Spacetime.kerr_isco -----------------------------------
    @testset "kerr_grid_coords" begin
        @test kerr_grid_coords(0.0, 12.0) == (0.0, sqrt(6.0 / 12.0))   # isco(0)=6
        for (a, r) in ((0.0, 8.0), (0.5, 10.0), (0.9, 7.0), (-0.3, 9.0))
            c = kerr_grid_coords(a, r)
            @test c[1] ≈ log(1 - a)
            @test c[2] ≈ sqrt(kerr_isco(a) / r)
        end
    end

end
