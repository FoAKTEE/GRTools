# Characterization of the vendored adaptive-RK integrators — ANALYTIC ODE oracles (NO heavy live dep;
# H1Lorenz and FEW are both FORBIDDEN as live test-deps). Expected values are closed-form ODE solutions
# / conserved quantities, plus a captured-literal DP5 first-step from the published Butcher tableau:
#   - y'=y → exp(x); y'=−y → exp(−x) decay
#   - harmonic oscillator x''=−x: energy x²+v² conserved over many periods
#   - Kepler two-body: energy ½v²−1/r and angular momentum x·vy−y·vx conserved
#   - step-doubling order check: halving h shrinks the DP5 single-step error by ≈2^6 (local order 6)
#   - DOPR853 dense output matches exp at interior points; integrate accuracy within abstol
#   - captured literal: DP5 step (y'=y, y0=1, h=0.1) → y5 = 1.1051709183333333 (computed independently
#     from the H1Lorenz integrator.jl:9-19 tableau; err vs exp(0.1) = 2.6e-10 ⇒ order-6 local error)

using NumericsCore
using NumericsCore: dp5_step!, integrate_to_targets!, gbs_integrate_to_targets!, gbs_kmax, rk_tol
using NumericsCore: DOPR853, take_step, integrate, dense_eval

# y' = y  ⇒  f!(g,x,y) = (g .= y)
expand!(g, x, y) = (g .= y)
decay!(g, x, y) = (g .= .-y)
# harmonic oscillator state [x, v], x'=v, v'=-x
osc!(g, x, y) = (g[1] = y[2]; g[2] = -y[1])

@testset "ode" begin
    @testset "DP5 captured-literal first step (published tableau, y'=y)" begin
        y0 = [1.0]
        g = similar(y0); expand!(g, 0.0, y0); k1 = copy(g)
        err, y5, k7 = dp5_step!(expand!, 0.0, y0, k1, 0.1, 1e-12)
        @test y5[1] ≈ 1.1051709183333333 atol=1e-13      # independent tableau evaluation
        @test abs(y5[1] - exp(0.1)) < 3e-10              # order-6 local truncation error
        @test k7[1] ≈ y5[1] atol=1e-13                   # FSAL: last stage derivative = y5 for y'=y
    end

    @testset "DP5 step-doubling order (local error ∝ h^6)" begin
        function step_err(h)
            y0 = [1.0]; g = similar(y0); expand!(g, 0.0, y0); k1 = copy(g)
            _, y5, _ = dp5_step!(expand!, 0.0, y0, k1, h, 1e-12)
            abs(y5[1] - exp(h))
        end
        r = step_err(0.1) / step_err(0.05)
        @test r ≈ 2.0^6 rtol=0.1                          # halving h ⇒ ×1/64 error
    end

    @testset "DP5 integrate_to_targets! — y'=y → exp, y'=−y → decay" begin
        tol = rk_tol(Float64)
        Yp = integrate_to_targets!(expand!, [1.0], 0.0, [0.5, 1.0, 2.0], +1; tol=tol)
        @test Yp[1, 1] ≈ exp(0.5) rtol=1e-10
        @test Yp[2, 1] ≈ exp(1.0) rtol=1e-10
        @test Yp[3, 1] ≈ exp(2.0) rtol=1e-10
        Yd = integrate_to_targets!(decay!, [1.0], 0.0, [1.0, 3.0], +1; tol=tol)
        @test Yd[1, 1] ≈ exp(-1.0) rtol=1e-10
        @test Yd[2, 1] ≈ exp(-3.0) rtol=1e-9
    end

    @testset "DP5 harmonic-oscillator energy conservation (many periods)" begin
        tol = rk_tol(Float64)
        tf = 20 * 2π
        Y = integrate_to_targets!(osc!, [1.0, 0.0], 0.0, [tf], +1; tol=tol)
        x, v = Y[1, 1], Y[1, 2]
        @test x^2 + v^2 ≈ 1.0 rtol=1e-8                   # energy conserved
        @test x ≈ cos(tf) atol=1e-5                       # phase tracked
    end

    @testset "DP5 Kepler two-body energy + angular-momentum conservation" begin
        # circular orbit r=1, v=1: state [x,y,vx,vy] = [1,0,0,1]; accel = -r/|r|^3
        function kepler!(g, t, y)
            x, yy, vx, vy = y
            r3 = (x^2 + yy^2)^1.5
            g[1] = vx; g[2] = vy; g[3] = -x / r3; g[4] = -yy / r3
        end
        E(y) = 0.5 * (y[3]^2 + y[4]^2) - 1 / sqrt(y[1]^2 + y[2]^2)
        L(y) = y[1] * y[4] - y[2] * y[3]
        y0 = [1.0, 0.0, 0.0, 1.0]
        Y = integrate_to_targets!(kepler!, y0, 0.0, [10 * 2π], +1; tol=rk_tol(Float64))
        yf = Y[1, :]
        @test E(yf) ≈ E(y0) rtol=1e-8
        @test L(yf) ≈ L(y0) rtol=1e-9
    end

    @testset "GBS integrate_to_targets! — y'=y → exp (explicit kmax)" begin
        Y = gbs_integrate_to_targets!(expand!, [1.0], 0.0, [1.0, 2.0], +1; tol=1e-11, kmax=8)
        @test Y[1, 1] ≈ exp(1.0) rtol=1e-6                 # GBS Float64 path reaches ~1e-8
        @test Y[2, 1] ≈ exp(2.0) rtol=1e-6
        @test gbs_kmax(Float64) == 0                       # Float64 ⇒ use DP5
        # GBS is the BigFloat-feasible high-order path: its order ~2·kmax keeps the step count tractable
        setprecision(BigFloat, 128) do
            @test gbs_kmax(BigFloat) >= 6
            Yb = gbs_integrate_to_targets!(expand!, [big(1.0)], big(0.0), [big(1.0)], +1;
                                           tol=big(2.0)^-60, kmax=gbs_kmax(BigFloat))
            @test abs(Yb[1, 1] - exp(big(1.0))) < big(1.0e-11)
        end
    end

    @testset "DP5 BigFloat precision inheritance — type-generic, modest tol" begin
        # NOTE: DP5 is order 5 ⇒ step count ~tol^(-1/5); the full rk_tol(BigFloat)=2^(20-bits) is the
        # documented DP5 BigFloat throughput ceiling (GBS is the feasible high-precision path). Here we
        # verify the GENERIC element-type threading at a feasible tol, not the extreme-precision limit.
        setprecision(BigFloat, 120) do
            Y = integrate_to_targets!(expand!, [big(1.0)], big(0.0), [big(1.0)], +1; tol=big(1.0e-18))
            @test Y[1, 1] isa BigFloat
            @test abs(Y[1, 1] - exp(big(1.0))) < big(1.0e-15)
        end
    end

    @testset "DOPR853 integrate + dense output (y'=y)" begin
        ode!(dydx, x, y, args) = (dydx .= y)
        solver = DOPR853(ode!; abstol=1e-12)
        ts, ys, coeffs = integrate(solver, [1.0], 3.0; h0=0.05)
        @test ys[end][1] ≈ exp(ts[end]) rtol=1e-10
        # dense output at interior points matches exp
        for tq in (0.37, 1.25, 2.6)
            @test dense_eval(ts, coeffs, tq)[1] ≈ exp(tq) rtol=1e-10
        end
        @test_throws ArgumentError dense_eval(ts, coeffs, ts[end] + 1.0)
    end

    @testset "DOPR853 clips the final step to the requested endpoint" begin
        constant!(dydx, x, y, args) = (dydx .= 1.0)
        solver = DOPR853(constant!; abstol = 1e6)
        ts, ys, coeffs = integrate(solver, [0.0], 0.05; h0 = 1.0, max_step = 5)
        @test ts[end] === 0.05
        @test ys[end][1] ≈ 0.05 atol = 1e-14
        @test all(diff(ts) .> 0)
        @test length(coeffs) == length(ts) - 1
        @test dense_eval(ts, coeffs, 0.05)[1] ≈ 0.05 atol = 1e-14
        @test_throws ArgumentError integrate(solver, [0.0], -1.0)
        @test_throws ArgumentError integrate(solver, [0.0], 1.0; h0 = 0.0)
        @test_throws ArgumentError integrate(solver, [0.0], 1.0; max_step = 0)
    end

    @testset "DOPR853 take_step order (local error ∝ h^p, p≈8)" begin
        ode!(dydx, x, y, args) = (dydx .= y)
        function lerr(h)
            s = DOPR853(ode!; abstol=1e-13)
            ok, x1, y1, hnew, cont = take_step(s, 0.0, [1.0], h)
            abs(y1[1] - exp(h))
        end
        r = lerr(0.2) / lerr(0.1)
        @test r > 2.0^7                                    # high-order local error collapse
    end
end

# ---- integrate_grid: grid-emitting DP5 (dense + capped) ------------------------------------------
# Resolves the 4th vendored ODE copy (legacy `DiscModes.integrator`) into this family. Gates: exact
# solutions in both directions, dense == capped within tolerance, the two RHS calling conventions
# agree bitwise, and BigFloat genericity.
@testset "integrate_grid — grid-emitting DP5" begin
    # dy/dx = [y2, -y1] (harmonic oscillator) → y = [sin x, cos x] from y(0) = [0, 1]
    frhs(x, y) = [y[2], -y[1]]
    x = collect(range(0.0, 6.0; length = 61))
    Y = integrate_grid(frhs, [0.0, 1.0], x, 1, +1; rtol = 1e-12, atol = 1e-14)
    @test size(Y) == (61, 2)
    @test maximum(abs.(Y[:, 1] .- sin.(x))) < 1e-10
    @test maximum(abs.(Y[:, 2] .- cos.(x))) < 1e-10

    # inward integration from the last grid point reproduces the same solution
    Yin = integrate_grid(frhs, [sin(6.0), cos(6.0)], x, 61, -1; rtol = 1e-12, atol = 1e-14)
    @test maximum(abs.(Yin[:, 1] .- sin.(x))) < 1e-9

    # starting in the middle fills only the requested side (the other side stays zero)
    Ymid = integrate_grid(frhs, [sin(x[31]), cos(x[31])], x, 31, +1; rtol = 1e-12, atol = 1e-14)
    @test all(iszero, Ymid[1:30, :])
    @test maximum(abs.(Ymid[31:end, 1] .- sin.(x[31:end]))) < 1e-10

    # dense (cubic-Hermite) and capped (land-on-grid) paths agree
    Yc = integrate_grid(frhs, [0.0, 1.0], x, 1, +1; rtol = 1e-12, atol = 1e-14, dense = false)
    @test maximum(abs.(Yc[:, 1] .- sin.(x))) < 1e-10
    @test maximum(abs.(Yc .- Y)) < 1e-9

    # in-place RHS (3-arg) is detected and gives the SAME bits as the allocating 2-arg form
    struct _OscRHS end
    (::_OscRHS)(out::AbstractVector, x, y) = (out[1] = y[2]; out[2] = -y[1]; out)
    (::_OscRHS)(x, y) = [y[2], -y[1]]
    Yip = integrate_grid(_OscRHS(), [0.0, 1.0], x, 1, +1; rtol = 1e-12, atol = 1e-14)
    @test Yip == Y

    # exponential decay, BigFloat: y' = -y, y(0)=1
    setprecision(BigFloat, 128) do
        xb = collect(range(big(0.0), big(3.0); length = 31))
        Yb = integrate_grid((x, y) -> [-y[1]], [big(1.0)], xb, 1, +1;
                            rtol = big(1e-24), atol = big(1e-30))
        @test eltype(Yb) === BigFloat
        # the DENSE path interpolates with a cubic Hermite, so grid values carry its O(h⁴) error
        # (measured 2.5e-20 here) even when the STEPS are integrated to 1e-24 …
        @test maximum(abs.(Yb[:, 1] .- exp.(-xb))) < 1e-19
        # … while the capped path lands on the grid points and keeps the step tolerance
        Ybc = integrate_grid((x, y) -> [-y[1]], [big(1.0)], xb, 1, +1;
                             rtol = big(1e-24), atol = big(1e-30), dense = false)
        @test maximum(abs.(Ybc[:, 1] .- exp.(-xb))) < 1e-23
    end
end

# ---- dp5_step_buffered! == dp5_step! at rtol = atol = tol ----------------------------------------
@testset "dp5_step_buffered! matches dp5_step!" begin
    f!(g, x, y) = (g[1] = y[2]; g[2] = -y[1]; g)
    y = [0.3, 0.9]; k1 = similar(y); f!(k1, 0.0, y); h = 0.05; tol = 1e-10
    err_a, y5_a, k7_a = dp5_step!(f!, 0.0, copy(y), copy(k1), h, tol)
    kbuf = [Vector{Float64}(undef, 2) for _ in 1:7]; copyto!(kbuf[1], k1)
    ys = zeros(2); y5 = zeros(2); y4 = zeros(2)
    err_b = dp5_step_buffered!(f!, 0.0, copy(y), kbuf, ys, y5, y4, h; rtol = tol, atol = tol)
    @test err_b == err_a
    @test y5 == y5_a
    @test kbuf[7] == k7_a
end
