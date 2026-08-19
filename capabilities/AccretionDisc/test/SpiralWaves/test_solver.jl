# test_solver.jl — EGD::sw-mode-solver: the implicit-midpoint global solve of arXiv:2608.17003 Sec. II.E
# (src/SpiralWaves/solver.jl). Gates, in order:
#   (a) the generic banded LU (the BigFloat linear-algebra path) against a dense solve;
#   (b) `eigen2` against the defining equation C·Y = kY;
#   (c) MANUFACTURED EXACT solution: constant C and affine y_ex make the midpoint rule EXACT, so the
#       solver must return y_ex to roundoff on a NON-uniform grid;
#   (d) MANUFACTURED SMOOTH solution: second-order convergence (error ratio → 4 when Δr halves);
#   (e) the physical mode (𝔭 = 𝔮 = 0, χ = 0, h0 = 0.05, a• = 15, m = 2…5): finite everywhere, the
#       Landau-regularised corotation resolved, residual of (64) at roundoff, boundary rows satisfied;
#   (f) the boundary branch: the selected wave carries its group velocity AWAY from the orbit and is a
#       TRAILING spiral (Im k > 0 at BOTH edges with the paper's own Fourier convention) — ERRATUM E-M4;
#       `bc = :paper` returns the other root, as printed;
#   (g) determinism: threaded == serial and BLAS-thread-independent, BITWISE; :sparse == :banded;
#   (h) the BigFloat path end-to-end on a 200-point grid.
using Test
using LinearAlgebra
using SparseArrays
using AccretionDisc
using AccretionDisc.SpiralWaves
using AccretionDisc.SpiralWaves: _band_from_sparse, _band_lu!, _band_solve!, _solve_linear
using Spacetime: kerr_equatorial_potentials

# reference disc of the paper's nominal configuration (M:876-884)
_bg_ref(::Type{T} = Float64; χ = 0, aEMRI = 15, h0 = 0.05, 𝔭 = 0, 𝔮 = 0) where {T} =
    DiscBackground{T}(; χ = T(χ), aEMRI = T(aEMRI), h0 = T(h0), 𝔭 = T(𝔭), 𝔮 = T(𝔮))
# small but genuinely adaptive grid (the reference triple is (240, 240, 24000), A:140)
_grid_small(bg, m; nL = 30, nR = 30, nC = 600, n_anc = 2000) =
    build_grid(bg, m; nL = nL, nR = nR, nC = nC, n_anc = n_anc)

@testset "SpiralWaves solver: implicit midpoint (64) + outgoing-wave BCs (66)-(68)" begin

    @testset "(a) generic banded LU == dense solve" begin
        for (T, tol) in ((Float64, 1e-10), (BigFloat, 1e-40))
            setprecision(BigFloat, 256) do
                n, kl, ku = 24, 2, 2
                A = spzeros(Complex{T}, n, n)
                for i in 1:n, j in max(1, i - kl):min(n, i + ku)
                    A[i, j] = Complex{T}(T(sin(3.1 * i + 1.7 * j)), T(cos(0.9 * i - 2.3 * j)))
                end
                for i in 1:n                       # make it comfortably nonsingular
                    A[i, i] += Complex{T}(T(6), T(0))
                end
                b = Complex{T}[Complex{T}(T(cos(i)), T(sin(2i))) for i in 1:n]
                xref = Matrix(A) \ b
                AB = _band_from_sparse(A, kl, ku)
                ipiv = _band_lu!(AB, kl, ku)
                x = _band_solve!(AB, ipiv, kl, ku, copy(b))
                @test maximum(abs, x - xref) / maximum(abs, xref) < tol
                @test maximum(abs, A * x - b) / maximum(abs, b) < tol
                # a genuinely pivot-needing case: zero leading diagonal entry
                A2 = copy(A); A2[1, 1] = zero(Complex{T})
                AB2 = _band_from_sparse(A2, kl, ku)
                ip2 = _band_lu!(AB2, kl, ku)
                @test ip2[1] != 1                                   # the first column had to pivot
                x2 = _band_solve!(AB2, ip2, kl, ku, copy(b))
                @test maximum(abs, A2 * x2 - b) / maximum(abs, b) < tol
                # an entry outside the declared band must be refused, never silently dropped
                A3 = copy(A); A3[1, 1 + kl + ku + 1] = one(Complex{T})
                @test_throws ArgumentError _band_from_sparse(A3, kl, ku)
            end
        end
    end

    @testset "(b) eigen2 solves C·Y = kY" begin
        for T in (Float64, BigFloat)
            setprecision(BigFloat, 256) do
                for (a, b, c, d) in ((1.0, 2.0, -0.5, 0.3), (0.0, 1.0, -1.0, 0.0), (2.0, 0.0, 0.0, -3.0),
                                     (0.1, 1e-8, 1e8, 0.1))
                    C = Complex{T}[Complex{T}(T(a), T(0.2)) Complex{T}(T(b), T(-0.1));
                                   Complex{T}(T(c), T(0.05)) Complex{T}(T(d), T(0.3))]
                    k1, k2, Y1, Y2 = eigen2(C)
                    for (k, Y) in ((k1, Y1), (k2, Y2))
                        res = (C[1, 1] * Y[1] + C[1, 2] * Y[2] - k * Y[1],
                               C[2, 1] * Y[1] + C[2, 2] * Y[2] - k * Y[2])
                        scale = max(maximum(abs, C), abs(k))
                        @test max(abs(res[1]), abs(res[2])) / scale < (T === Float64 ? 1e-12 : 1e-40)
                    end
                    @test abs(k1 + k2 - (C[1, 1] + C[2, 2])) / max(abs(k1), abs(k2)) < (T === Float64 ? 1e-13 : 1e-40)
                end
                # Y⊥ is orthogonal to Y in the UNCONJUGATED bilinear form (the (68) row is "y ∥ Y")
                C = Complex{T}[Complex{T}(T(0.3), T(0)) Complex{T}(T(0), T(-1.2));
                               Complex{T}(T(0.7), T(0)) Complex{T}(T(-0.2), T(0))]
                k, Y, W = outgoing_branch(C, T(1); outward = true)
                @test abs(W[1] * Y[1] + W[2] * Y[2]) < (T === Float64 ? 1e-15 : 1e-45)
            end
        end
    end

    # manufactured problem: y' = C(r) y + S(r) with y_ex known, S := y_ex' − C y_ex
    _mansolve(r, Cf, Sf, yex; linsolve = :auto) = begin
        rm = [(r[i] + r[i+1]) / 2 for i in 1:length(r)-1]
        Cmid = [Cf(x) for x in rm]
        Smid = [Sf(x) for x in rm]
        solve_midpoint(r, Cmid, Smid, (1, 0), yex(r[1])[1], (0, 1), yex(r[end])[2]; linsolve = linsolve)
    end

    @testset "(c) manufactured EXACT: constant C, affine y_ex, non-uniform grid" begin
        C0 = ComplexF64[0.1+0.0im 0.2+0.3im; -0.4+0.1im 0.05-0.2im]
        yex(r) = (ComplexF64(1 + 2r, 0.5 - r), ComplexF64(3 - r, 2 + 0.25r))
        dy = (ComplexF64(2, -1), ComplexF64(-1, 0.25))
        Cf(r) = C0
        Sf(r) = begin
            y = yex(r)
            [dy[1] - (C0[1, 1] * y[1] + C0[1, 2] * y[2]), dy[2] - (C0[2, 1] * y[1] + C0[2, 2] * y[2])]
        end
        r = 1.0 .+ 4.0 .* (range(0, 1; length = 51) .^ 1.7)          # strictly increasing, non-uniform
        for ls in (:auto, :banded)
            y = _mansolve(r, Cf, Sf, yex; linsolve = ls)
            err = maximum(i -> max(abs(y[1, i] - yex(r[i])[1]), abs(y[2, i] - yex(r[i])[2])), eachindex(r))
            @info "manufactured exact (linsolve = $ls): max |y − y_ex| = $err"
            @test err < 1e-12
        end
    end

    @testset "(d) manufactured SMOOTH: second-order convergence of the midpoint rule" begin
        Cf(r) = ComplexF64[0.3+0.1r  (-0.5*sin(r))+0.2im; 1/(1+r)  0.2*cos(r)-0.1im]
        yex(r) = (ComplexF64(cos(2r), sin(3r)), ComplexF64(exp(-r/3), r^2 / 10))
        dyex(r) = (ComplexF64(-2sin(2r), 3cos(3r)), ComplexF64(-exp(-r/3)/3, r/5))
        Sf(r) = begin
            y = yex(r); C = Cf(r); d = dyex(r)
            [d[1] - (C[1, 1] * y[1] + C[1, 2] * y[2]), d[2] - (C[2, 1] * y[1] + C[2, 2] * y[2])]
        end
        errs = Float64[]
        for N in (51, 101, 201, 401, 801)
            r = collect(range(1.0, 5.0; length = N))
            y = _mansolve(r, Cf, Sf, yex)
            push!(errs, maximum(i -> max(abs(y[1, i] - yex(r[i])[1]), abs(y[2, i] - yex(r[i])[2])), eachindex(r)))
        end
        ratios = errs[1:end-1] ./ errs[2:end]
        @info "manufactured smooth: errors $(errs), halving ratios $(ratios) (2nd order ⇒ 4)"
        @test all(3.9 .< ratios .< 4.1)
        @test errs[end] < 5e-5
        @test errs[end] * 800^2 < 30                       # the O(Δr²) constant is O(1) on this problem
    end

    @testset "(e) physical mode 𝔭=𝔮=0, χ=0, h0=0.05, a•=15, m = 2…5" begin
        bg = _bg_ref()
        for m in 2:5
            g = _grid_small(bg, m)
            sol = solve_mode(bg, m; grid = g)
            @test sol isa ModeSolution{Float64}
            @test length(sol) == length(g.r)
            @test all(isfinite, sol.y)
            @test sol.residual < 1e-10
            @test sol.bc_residual < 1e-10 * maximum(abs, sol.y)
            # corotation is REGULARISED and RESOLVED: |y| stays finite through the central block and the
            # node-to-node variation there is a small fraction of the block amplitude (the sinh block of
            # (C4)-(C5) resolves the Landau-shifted pole). The GLOBAL max sits at the inner edge, not at
            # r_c — the short-wavelength wave amplifies inward; that is physics, not a pole.
            ic = argmax(vec(maximum(abs, sol.y; dims = 1)))
            blk = g.blocks.C
            ampC = maximum(abs, view(sol.y, :, blk))
            jump = maximum(abs, diff(view(sol.y, 1, blk))) / ampC
            @info "m = $m: N = $(length(g.r)), max|y| = $(maximum(abs, sol.y)) at r = $(sol.r[ic]) (r_c = $(g.r_c)), corotation-block max = $ampC, max node-to-node jump there = $jump, residual(64) = $(sol.residual), bc = $(sol.bc_residual), $(round(sol.seconds; digits = 3)) s"
            @test isfinite(ampC) && ampC > 0
            @test jump < 0.05                              # the corotation block resolves the structure
            @test all(isfinite, view(sol.y, :, blk))
            # δû_φ at the nodes (route R2) is finite as well
            duphi = delta_u_phi(bg, sol)
            @test length(duphi) == length(sol.r)
            @test all(isfinite, duphi)
        end
    end

    @testset "(f) boundary branch: group velocity away from the orbit (ERRATUM E-M4)" begin
        bg = _bg_ref()
        m = 3
        g = _grid_small(bg, m)
        sol = solve_mode(bg, m; grid = g)
        ωd_in = omega_d(bg, m, g.r[1]); ωd_out = omega_d(bg, m, g.r[end])
        @test ωd_in < 0 && ωd_out > 0                       # inside / outside corotation
        # the WKB radial wavenumber of the selected branch and the group velocity c_s² k_r/ω_d
        kr_in = imag(sol.k_in); kr_out = imag(sol.k_out)
        vg_in = csound(bg, g.r[1])^2 * kr_in / ωd_in
        vg_out = csound(bg, g.r[end])^2 * kr_out / ωd_out
        @info "BC branch: k_in = $(sol.k_in), k_out = $(sol.k_out); v_g(in) = $vg_in (want < 0), v_g(out) = $vg_out (want > 0)"
        @test vg_in < 0                                      # away from the orbit ⇒ inward at r_in
        @test vg_out > 0                                     # away from the orbit ⇒ outward at r_out
        @test kr_in > 0 && kr_out > 0                        # TRAILING spiral dφ/dr = −k_r/m < 0
        @test imag(sol.k_in) > 0 && imag(sol.k_out) > 0      # ⇒ Im[k] > 0, NOT the printed Im[k] < 0
        @test abs(imag(sol.k_out)) > 10 * abs(real(sol.k_out))   # the edge is in a wave zone
        # the RELATIVISTIC dispersion relation k_r = e^{μ}u^t√(ω_d²−κ²)/c_s reproduces the eigenvalue;
        # the NEWTONIAN one of grid.jl (App. C) is too small by exactly the factor e^{μ}u^t [FINDING]
        for (r, kr) in ((g.r[1], kr_in), (g.r[end], kr_out))
            rel = abs(kr - wkb_wavenumber_gr(bg, m, r)) / kr
            fac = kr / wkb_wavenumber(bg, m, r)
            emu_ut = exp(kerr_equatorial_potentials(bg.χ, r).μ) * rotation_profile(bg, r).ut
            @info "edge r = $r: Im k = $kr, relativistic WKB = $(wkb_wavenumber_gr(bg, m, r)) (rel $rel), Newtonian WKB = $(wkb_wavenumber(bg, m, r)), ratio = $fac vs e^μ u^t = $emu_ut"
            @test rel < 2e-3
            @test isapprox(fac, emu_ut; rtol = 5e-3)
        end
        @test wkb_wavenumber_gr(bg, m, g.r[1]) > 1.7 * wkb_wavenumber(bg, m, g.r[1])   # ISCO: factor √3
        # the literal reading of M:776 selects the OTHER root (the ingoing/leading branch)
        solp = solve_mode(bg, m; grid = g, bc = :paper)
        @test imag(solp.k_in) < 0 && imag(solp.k_out) < 0
        @test isapprox(solp.k_in, conj(sol.k_in); rtol = 1e-6)
        @test !isapprox(solp.y[1, end], sol.y[1, end]; rtol = 1e-3)
        @test_throws ArgumentError solve_mode(bg, m; grid = g, bc = :nonsense)
    end

    @testset "(g) determinism: threads, BLAS threads, and the two linear-algebra backends" begin
        bg = _bg_ref()
        m = 3
        g = _grid_small(bg, m; nL = 20, nR = 20, nC = 200, n_anc = 1000)
        s_serial = solve_mode(bg, m; grid = g, threads = false)
        s_thread = solve_mode(bg, m; grid = g, threads = true)
        @test s_serial.y == s_thread.y                                  # BITWISE
        nb = BLAS.get_num_threads()
        try
            BLAS.set_num_threads(1); s1 = solve_mode(bg, m; grid = g)
            BLAS.set_num_threads(4); s4 = solve_mode(bg, m; grid = g)
            @test s1.y == s4.y                                          # BITWISE
        finally
            BLAS.set_num_threads(nb)
        end
        s_band = solve_mode(bg, m; grid = g, linsolve = :banded)
        rel = maximum(abs, s_band.y - s_serial.y) / maximum(abs, s_serial.y)
        @info "linsolve :sparse (UMFPACK) vs :banded (generic LU): max rel difference = $rel; threads = $(Threads.nthreads())"
        @test rel < 1e-8
        @test s_band.residual < 1e-10
    end

    @testset "(h) BigFloat path on a 200-point grid" begin
        setprecision(BigFloat, 256) do
            bgB = _bg_ref(BigFloat)
            m = 3
            gB = build_grid(bgB, m; nL = 8, nR = 8, nC = 60, n_anc = 400)
            @test length(gB.r) ≥ 60
            solB = solve_mode(bgB, m; grid = gB, hmodes_N = 512)
            @test solB isa ModeSolution{BigFloat}
            @test eltype(solB.y) === Complex{BigFloat}
            @test all(isfinite, solB.y)
            @info "BigFloat(256): N = $(length(gB.r)), residual(64) = $(Float64(solB.residual)), bc = $(Float64(solB.bc_residual)), $(round(solB.seconds; digits = 2)) s"
            @test solB.residual < 1e-25
            # the same physical branch selection as Float64
            @test imag(solB.k_in) > 0 && imag(solB.k_out) > 0
            # Float64 agreement on the SAME grid
            bg64 = _bg_ref()
            g64 = build_grid(bg64, m; nL = 8, nR = 8, nC = 60, n_anc = 400)
            if length(g64.r) == length(gB.r)
                sol64 = solve_mode(bg64, m; grid = g64, hmodes_N = 512)
                rel = maximum(abs, ComplexF64.(solB.y) - sol64.y) / maximum(abs, sol64.y)
                @info "BigFloat vs Float64 on the same grid: max rel difference = $rel"
                @test rel < 1e-6
            end
        end
    end
end
