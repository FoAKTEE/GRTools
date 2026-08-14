# Characterization of the time-domain ±m mode-sum. NO live heavy dep (FEW/KerrBackwards forbidden).
# Grounding:
#   - direct_mode_sum: HAND-COMPUTED reference sums from h_i = Σ_j y_j A_ij cis(-ψ) + [m_j>0] y'_j conj(A_ij) cis(+ψ),
#     ψ = m_j Φ_φ + n_j Φ_r. Phases chosen as multiples of π/2 so every value is an exact closed form:
#       * K=1, m=0  (no partner)           → h = 2+3im
#       * K=1, m>0  (partner branch)        → h = 9+2im
#       * K=2, mixed (m=0 + m=2 partner)    → h = 2-1im
#   - interpolated_mode_sum: spline-resample == direct_mode_sum identity at the original nodes. The
#     trajectory is a global cubic in t and the dense grid coincides with the nodes, so the not-a-knot
#     CubicSplineInterpolant reproduces the data exactly ⇒ equals direct_mode_sum on the raw modes.
#   - PARALLEL == SERIAL bit-identity (===) for both kernels.
#   - FD-SPA (fd_interpolated_mode_sum/spa_func/fd_flip_average, S4-a hoist from EMRIWaveforms.ModeSum.FEW):
#       * spa_func / _fd_seg_freq_range: legacy-KB captured literals (pre-change grounded) + the |x|=7
#         series/asymptotic continuity identity; _cube_roots: analytic roots of (x−1)(x−2)(x−3).
#       * fd_flip_average: ±f Hermitian-flip symmetry (analytic self-consistency identity).
#       * fd_interpolated_mode_sum: EXACT (bit-level ==) match of every output bin against the BAKED
#         pre-hoist fixture fixtures/fdspa_ref.jl (captured from the committed EMRIWaveforms.ModeSum.FEW
#         local kernel @ f409fc7 on the shared fixtures/fdspa_inputs.jl generators) — hp/hc, the
#         include_minus_m=false and separate_modes=true variants, and a zero-crossing mode set whose
#         single-(mode,segment) bin ranges straddle the zero bin (exercising the direct-before-partner,
#         partner-before-direct AND i==zidx interleaved ±m accumulation orders of the legacy scatter).
#       * PARALLEL == SERIAL bit-identity (isequal — NaN-aware) for combined + separate-modes paths.

using Test
using Kernels

@testset "modesum" begin
    @testset "direct_mode_sum: hand-computed closed forms" begin
        # K=1, m=0: ψ=0, cis(0)=1, no partner. h = y*A = 1*(2+3im).
        h = direct_mode_sum(reshape(ComplexF64[2+3im], 1, 1), ComplexF64[1, 999],
                            [0.0], [0.0], [0], [0])
        @test h == ComplexF64[2+3im]

        # K=1, m>0: Φ_φ=π/2 ⇒ ψ=π/2, cis(-π/2)=-i, cis(π/2)=i. y=1, y'=2, A=2+3im.
        # h = (2+3im)(-i) + 2 conj(2+3im)(i) = (3-2i) + 2(3+2i) = 9+2i.
        h = direct_mode_sum(reshape(ComplexF64[2+3im], 1, 1), ComplexF64[1, 2],
                            [π/2], [0.0], [1], [0])
        @test h ≈ ComplexF64[9+2im] atol=1e-13

        # K=2: mode1 (m=0,n=1), mode2 (m=2,n=0). Φ_φ=π/4, Φ_r=π/2.
        #   ψ1 = π/2 ⇒ contrib1 = 1*(1+0im)*(-i) = -i  (no partner).
        #   ψ2 = π/2 ⇒ contrib2 = (1+1im)(-i) + conj(1+1im)(i) = (1-i)+(1+i) = 2.
        #   h = 2 - i.
        teuk = ComplexF64[1+0im 1+1im]              # 1 x 2
        ylms = ComplexF64[1, 1, 0, 1]               # base y1,y2 ; partner y1',y2'
        h = direct_mode_sum(teuk, ylms, [π/4], [π/2], [0, 2], [1, 0])
        @test h ≈ ComplexF64[2-1im] atol=1e-13
    end

    @testset "direct_mode_sum: parallel === serial (bit-identity)" begin
        nt, K = 37, 4
        teuk = ComplexF64[cos(0.3i + 0.7j) + im*sin(0.2i - 0.5j) for i in 1:nt, j in 1:K]
        ylms = ComplexF64[cis(0.11k) for k in 1:2K]
        Phi_phi = Float64[0.05i^2 - 0.3i for i in 1:nt]
        Phi_r = Float64[0.02i^2 + 0.1i for i in 1:nt]
        m_arr = [0, 1, 2, 3]; n_arr = [0, -1, 1, 2]
        hs = direct_mode_sum(teuk, ylms, Phi_phi, Phi_r, m_arr, n_arr; threaded = false)
        ht = direct_mode_sum(teuk, ylms, Phi_phi, Phi_r, m_arr, n_arr; threaded = true)
        @test ht == hs
        @test all(ht[i] === hs[i] for i in 1:nt)
    end

    @testset "interpolated_mode_sum: spline-resample == direct at nodes" begin
        # cubic trajectory in t; nodes = dense grid (dt=10) so the spline is exact.
        N = 8
        t = Float64[10.0 * (i - 1) for i in 1:N]
        K = 2
        cubic(c0, c1, c2, c3, x) = c0 + c1*x + c2*x^2 + c3*x^3
        teuk = Matrix{ComplexF64}(undef, N, K)
        for i in 1:N, j in 1:K
            x = t[i] / 100
            re = cubic(0.5j, 0.3, -0.2, 0.05, x)
            im_ = cubic(-0.1, 0.4j, 0.1, -0.03, x)
            teuk[i, j] = complex(re, im_)
        end
        Phi_phi = Float64[cubic(0.0, 1.0, -0.1, 0.02, t[i] / 100) for i in 1:N]
        Phi_r = Float64[cubic(0.0, 0.7, 0.05, -0.01, t[i] / 100) for i in 1:N]
        ylms = ComplexF64[cis(0.2k) for k in 1:2K]
        m_arr = [0, 1]; n_arr = [1, 0]

        hi = interpolated_mode_sum(t, teuk, ylms, Phi_phi, Phi_r, m_arr, n_arr;
                                   dt = 10.0, out_len = N)
        hd = direct_mode_sum(teuk, ylms, Phi_phi, Phi_r, m_arr, n_arr)
        @test hi ≈ hd rtol=1e-8

        # parallel === serial
        hs = interpolated_mode_sum(t, teuk, ylms, Phi_phi, Phi_r, m_arr, n_arr;
                                   dt = 10.0, out_len = N, threaded = false)
        ht = interpolated_mode_sum(t, teuk, ylms, Phi_phi, Phi_r, m_arr, n_arr;
                                   dt = 10.0, out_len = N, threaded = true)
        @test all(ht[i] === hs[i] for i in eachindex(hs))
    end
end

include(joinpath(@__DIR__, "fixtures", "fdspa_inputs.jl"))
include(joinpath(@__DIR__, "fixtures", "fdspa_ref.jl"))

# Exact bin-by-bin comparison against a baked (bin, re, im) fixture: nonzero finite bins must match
# bit-exactly (== on repr-round-tripped Float64), NaN bins must be NaN, every other bin exactly zero.
function fdspa_check(v::AbstractVector, ref, nanbins)
    nz = Dict(b => complex(re, im_) for (b, re, im_) in ref)
    nanset = Set(nanbins)
    ok_nz = all(v[b] == z for (b, z) in nz)
    ok_nan = all(isnan(v[b]) for b in nanset)
    ok_zero = all(iszero(v[i]) for i in eachindex(v) if !haskey(nz, i) && !(i in nanset))
    return (ok_nz, ok_nan, ok_zero)
end
function fdspa_check(M::AbstractMatrix, ref, nanbins)
    nz = Dict((mi, b) => complex(re, im_) for (mi, b, re, im_) in ref)
    nanset = Set(nanbins)
    ok_nz = all(M[k[1], k[2]] == z for (k, z) in nz)
    ok_nan = all(isnan(M[k[1], k[2]]) for k in nanset)
    ok_zero = all(iszero(M[I]) for I in CartesianIndices(M)
                  if !haskey(nz, (I[1], I[2])) && !((I[1], I[2]) in nanset))
    return (ok_nz, ok_nan, ok_zero)
end
fdspa_run(c; kw...) = fd_interpolated_mode_sum(c.t, c.teuk, c.ylms, c.Phi_phi, c.Phi_r,
                                               c.f_phi, c.f_r, c.m_arr, c.n_arr;
                                               f_arr = c.f_arr, kw...)

@testset "modesum FD-SPA" begin
    @testset "spa_func: legacy-KB captured + series/asymptotic continuity" begin
        # captured legacy KerrBackwards spa_func values (pre-change grounded literals)
        @test spa_func(7.0)  ≈ 0.8768753314049489  - 0.8942724953516838im rtol = 1e-13
        @test spa_func(0.5)  ≈ 0.7691512654746617  - 0.9158084090089567im rtol = 1e-13
        @test spa_func(-3.0) ≈ 0.9025394332751743  - 0.8634854343984399im rtol = 1e-13
        @test spa_func(20.0) ≈ 0.8830719695332906  - 0.8892179686552646im rtol = 1e-13
        # ascending series (|x|≤7) and asymptotic series (|x|>7) agree across the branch switch
        @test spa_func(7.0) ≈ spa_func(7.0001) rtol = 1e-5
    end

    @testset "_cube_roots / _fd_seg_freq_range" begin
        r = Kernels._cube_roots(1.0, -6.0, 11.0, -6.0)   # (x−1)(x−2)(x−3)
        @test sort(real.(collect(r))) ≈ [1.0, 2.0, 3.0] rtol = 1e-12
        @test all(abs.(imag.(collect(r))) .< 1e-12)
        # Exact spline degeneracies dispatch by polynomial degree instead of dividing by a=0.
        rq = Kernels._cube_roots(0.0, 1.0, -3.0, 2.0) # (x-1)(x-2)
        @test sort(real.(collect(rq))) ≈ [1.0, 2.0] rtol = 1e-14
        rscale = Kernels._cube_roots(0.0, 1.0, -1e10, 1.0)
        scaled_residual(r) = abs(evalpoly(r, (1.0, -1e10, 1.0))) /
                             (1 + 1e10 * abs(r) + abs2(r))
        @test all(r -> scaled_residual(r) < 10eps(Float64), rscale)
        @test minimum(abs, rscale) ≈ 1e-10 rtol = 1e-14
        @test length(Kernels._cube_roots(0.0, 0.0, 2.0, -4.0)) == 1
        @test only(Kernels._cube_roots(0.0, 0.0, 2.0, -4.0)) == 2.0 + 0.0im
        @test isempty(Kernels._cube_roots(0.0, 0.0, 0.0, 1.0))
        @test isempty(Kernels._cube_roots(0.0, 0.0, 0.0, 0.0))
        @test all(iszero, Kernels._cube_roots(1.0, 0.0, 0.0, 0.0))
        # Cardano's +sqrt branch cancels for x^3-1; the alternate branch must
        # recover all three roots rather than silently returning a triple zero.
        rpure = Kernels._cube_roots(1.0, 0.0, 0.0, -1.0)
        @test any(r -> isapprox(r, 1.0; rtol = 1e-14, atol = 1e-14), rpure)
        @test all(r -> abs(r^3 - 1) < 1e-13, rpure)
        # captured turnover-aware segment frequency range (pre-change grounded literal)
        @test all(Kernels._fd_seg_freq_range(1.0, -2.0, 0.5, 0.1, 2.0) .≈
                  (-0.4520627463118869, 1.0999999999999996))
        # Quadratic F(x)=(x-1)^2+2 on [0,2] has an interior minimum missed by cubic-only logic.
        @test Kernels._fd_seg_freq_range(0.0, 1.0, -2.0, 3.0, 2.0) == (2.0, 3.0)
    end

    @testset "fd_flip_average: ±f Hermitian symmetry" begin
        nf = 9
        freq = [1.0 * (i - 5) for i in 1:nf]
        wf = ComplexF64[complex(0.3 * i, -0.2 * i + 0.1) for i in 1:nf]
        hp, hc = fd_flip_average(wf, freq)
        for i in 1:nf
            j = nf + 1 - i
            @test hp[i] ≈ conj(hp[j]) rtol = 1e-13
            @test hc[i] ≈ conj(hc[j]) rtol = 1e-13
        end
    end

    @testset "fd_interpolated_mode_sum == baked pre-hoist fixture (case 1)" begin
        c1 = fdspa_case1_inputs()
        hp, hc = fdspa_run(c1)
        @test fdspa_check(hp, FDSPA_C1_HP, FDSPA_C1_HP_NAN) == (true, true, true)
        @test fdspa_check(hc, FDSPA_C1_HC, FDSPA_C1_HC_NAN) == (true, true, true)
        hpn, hcn = fdspa_run(c1; include_minus_m = false)
        @test fdspa_check(hpn, FDSPA_C1_NOMM_HP, FDSPA_C1_NOMM_HP_NAN) == (true, true, true)
        @test fdspa_check(hcn, FDSPA_C1_NOMM_HC, FDSPA_C1_NOMM_HC_NAN) == (true, true, true)
        sep = fdspa_run(c1; separate_modes = true)
        @test size(sep) == (size(c1.teuk, 2), length(c1.f_arr))
        @test fdspa_check(sep, FDSPA_C1_SEP, FDSPA_C1_SEP_NAN) == (true, true, true)
    end

    @testset "fd_interpolated_mode_sum == baked pre-hoist fixture (case 2, zero-crossing ±m orders)" begin
        c2 = fdspa_case2_inputs()
        hp, hc = fdspa_run(c2)
        @test fdspa_check(hp, FDSPA_C2_HP, FDSPA_C2_HP_NAN) == (true, true, true)
        @test fdspa_check(hc, FDSPA_C2_HC, FDSPA_C2_HC_NAN) == (true, true, true)
    end

    @testset "fd_interpolated_mode_sum: parallel == serial (bit-identity, NaN-aware)" begin
        for c in (fdspa_case1_inputs(), fdspa_case2_inputs())
            hps, hcs = fdspa_run(c; threaded = false)
            hpt, hct = fdspa_run(c; threaded = true)
            @test isequal(hpt, hps)   # isequal: elementwise bit-level incl. NaN bins
            @test isequal(hct, hcs)
        end
        c1 = fdspa_case1_inputs()
        ss = fdspa_run(c1; separate_modes = true, threaded = false)
        st = fdspa_run(c1; separate_modes = true, threaded = true)
        @test isequal(st, ss)
    end

    @testset "fd_interpolated_mode_sum: argument validation" begin
        c1 = fdspa_case1_inputs()
        # no exact-0.0 frequency bin
        @test_throws ArgumentError fd_interpolated_mode_sum(c1.t, c1.teuk, c1.ylms, c1.Phi_phi,
                                                            c1.Phi_r, c1.f_phi, c1.f_r, c1.m_arr,
                                                            c1.n_arr; f_arr = c1.f_arr .+ 1e-9)
        # ylms must have length 2K
        @test_throws ArgumentError fd_interpolated_mode_sum(c1.t, c1.teuk, c1.ylms[1:3], c1.Phi_phi,
                                                            c1.Phi_r, c1.f_phi, c1.f_r, c1.m_arr,
                                                            c1.n_arr; f_arr = c1.f_arr)
    end
end
