# ModeSum characterization.
#
# Grounding:
#  (A) ANALYTIC hand-computed h = Σ A·sYlm·e^{-iΦ} — a single-sample, K=2 toy whose strain is written out
#      term by term (closed-form complex sum, independent of the kernel's loop) and matched to full
#      precision, for both the KB fork (no (−1)^l) and the FEW fork (explicit (−1)^l partner sign).
#  (B) parallel == serial bit-identical (threaded vs serial Kernels reduction).
#  (C) interpolated ≡ direct at coincident dense nodes (spline reproduces the trajectory).
#  (D) FD-SPA forks error as documented [HOLE]s.
# (The underlying Kernels.direct/interpolated_mode_sum are themselves characterized vs legacy in Kernels.)

using Test
using Random
import EMRIWaveforms.ModeSum
const KBsum = ModeSum.KerrBackwards
const FEWsum = ModeSum.FEW

@testset "ModeSum" begin
    # ---- (A) hand-computed strain, KB fork ----
    @testset "hand-computed direct sum (KB fork)" begin
        a1 = 0.7 + 0.2im; a2 = -0.3 + 0.9im
        y1 = 1.1 + 0.0im; y2 = 0.5 - 0.4im; y2p = 0.8 + 0.1im
        φ = 0.37; ρ = 1.21
        teuk = reshape([a1, a2], 1, 2)
        ylms = ComplexF64[y1, y2, 0.0, y2p]   # [base(2); partners(2)] — partner of m=0 mode unused
        m = [0, 2]; n = [0, 1]
        h = KBsum.direct_mode_sum(teuk, ylms, [φ], [ρ], m, n)
        href = y1 * a1 + y2 * a2 * cis(-(2φ + ρ)) + y2p * conj(a2) * cis(2φ + ρ)
        @test h[1] ≈ href rtol=1e-14
    end

    # ---- (A) hand-computed strain, FEW fork (explicit (−1)^l) ----
    @testset "hand-computed direct sum (FEW fork, (−1)^l)" begin
        a1 = 0.7 + 0.2im; a2 = -0.3 + 0.9im
        y1 = 1.1 + 0.0im; y2 = 0.5 - 0.4im; y2p = 0.8 + 0.1im
        φ = 0.37; ρ = 1.21
        teuk = reshape([a1, a2], 1, 2)
        ylms = ComplexF64[y1, y2, 0.0, y2p]
        l = [2, 3]; m = [0, 2]; k = [0, 0]; n = [0, 1]   # l2=3 ⇒ (−1)^l = −1
        h = FEWsum.direct_mode_sum(teuk, ylms, [φ], [0.0], [ρ], l, m, k, n)
        href = y1 * a1 + y2 * a2 * cis(-(2φ + ρ)) + (-1)^3 * y2p * conj(a2) * cis(2φ + ρ)
        @test h[1] ≈ href rtol=1e-14
        # nonzero k (generic-Kerr polar term) rejected by the equatorial consolidation.
        @test_throws ArgumentError FEWsum.direct_mode_sum(teuk, ylms, [φ], [0.0], [ρ], l, m, [1, 0], n)
    end

    # ---- (B) parallel == serial bit-identical ----
    @testset "parallel == serial" begin
        Random.seed!(20240630)
        nt, K = 64, 5
        teuk = randn(ComplexF64, nt, K)
        ylms = randn(ComplexF64, 2K)
        pp = cumsum(rand(nt)); pr = cumsum(rand(nt))
        m = [0, 1, 2, 3, 4]; n = [0, -1, 1, -2, 2]
        hs = KBsum.direct_mode_sum(teuk, ylms, pp, pr, m, n; threaded = false)
        hp = KBsum.direct_mode_sum(teuk, ylms, pp, pr, m, n; threaded = true)
        @test hs == hp   # bit-identical
        l = [2, 3, 4, 5, 6]; k = zeros(Int, K)
        fs = FEWsum.direct_mode_sum(teuk, ylms, pp, [0.0 for _ in 1:nt], pr, l, m, k, n; threaded = false)
        fp = FEWsum.direct_mode_sum(teuk, ylms, pp, [0.0 for _ in 1:nt], pr, l, m, k, n; threaded = true)
        @test fs == fp
    end

    # ---- (C) interpolated ≡ direct at coincident dense nodes ----
    @testset "interpolated at coincident nodes ≡ direct" begin
        dt = 10.0
        nt = 6
        t = [dt * (i - 1) for i in 1:nt]
        teuk = ComplexF64[cis(0.1 * i) * (1 + 0.0im) + 0.0im for i in 1:nt, _ in 1:2]
        # smooth trajectories so the spline reproduces near-exactly at nodes
        teuk = ComplexF64[complex(cos(0.05 * i * j), sin(0.03 * i * j)) for i in 1:nt, j in 1:2]
        ylms = ComplexF64[1.0, 0.5 - 0.2im, 0.0, 0.7 + 0.1im]
        m = [0, 2]; n = [0, 1]
        pp = [0.2 * i for i in 1:nt]; pr = [0.15 * i for i in 1:nt]
        hd = KBsum.direct_mode_sum(teuk, ylms, pp, pr, m, n)
        hi = KBsum.interpolated_mode_sum(t, teuk, ylms, pp, pr, m, n; dt = dt, out_len = nt)
        @test hi ≈ hd rtol=1e-10
    end

    # ---- (D) FD-SPA stationary-phase mode-sum ----
    # Grounding:
    #  · _spa_func: ASCENDING-series vs ASYMPTOTIC-series agree across the |x|=7 switch (continuity),
    #    plus captured legacy KB values.
    #  · _cube_roots: roots of (x−1)(x−2)(x−3) = x³−6x²+11x−6 are {1,2,3} (analytic).
    #  · fd_flip_average: ±f Hermitian-flip symmetry (self-consistency identity).
    #  · fd_interpolated_mode_sum: CAPTURED legacy KB output (one-shot KerrBackwards run, baked literals).
    @testset "FD-SPA _spa_func" begin
        # series/asymptotic continuity across the |x|=7 branch switch
        @test FEWsum._spa_func(7.0) ≈ FEWsum._spa_func(7.0001) rtol = 1e-5
        # captured legacy KB _spa_func values
        @test FEWsum._spa_func(7.0)  ≈ 0.8768753314049489  - 0.8942724953516838im rtol = 1e-13
        @test FEWsum._spa_func(0.5)  ≈ 0.7691512654746617  - 0.9158084090089567im rtol = 1e-13
        @test FEWsum._spa_func(-3.0) ≈ 0.9025394332751743  - 0.8634854343984399im rtol = 1e-13
        @test FEWsum._spa_func(20.0) ≈ 0.8830719695332906  - 0.8892179686552646im rtol = 1e-13
    end

    @testset "FD-SPA _cube_roots" begin
        r = FEWsum._cube_roots(1.0, -6.0, 11.0, -6.0)   # (x−1)(x−2)(x−3)
        rr = sort(real.(collect(r)))
        @test rr ≈ [1.0, 2.0, 3.0] rtol = 1e-12
        @test all(abs.(imag.(collect(r))) .< 1e-12)
        # captured turnover-aware seg freq range
        @test all(FEWsum._fd_seg_freq_range(1.0, -2.0, 0.5, 0.1, 2.0) .≈ (-0.4520627463118869, 1.0999999999999996))
    end

    @testset "FD-SPA fd_flip_average ±f symmetry" begin
        # h₊,h₍×₎ from flip-average have Hermitian ±f structure: hp(−f)=conj(hp(f)), hc(−f)=conj(hc(f)).
        nf = 9
        freq = [1.0 * (i - 5) for i in 1:nf]   # symmetric, zero-centred
        wf = ComplexF64[complex(0.3 * i, -0.2 * i + 0.1) for i in 1:nf]
        hp, hc = FEWsum.fd_flip_average(wf, freq)
        for i in 1:nf
            j = nf + 1 - i            # mirror bin (−f)
            @test hp[i] ≈ conj(hp[j]) rtol = 1e-13
            @test hc[i] ≈ conj(hc[j]) rtol = 1e-13
        end
    end

    @testset "FD-SPA fd_interpolated_mode_sum vs captured KB" begin
        nt = 5; K = 2
        t = collect(0.0:100.0:400.0)
        teuk = ComplexF64[complex(0.5*cos(0.01*i*j), 0.3*sin(0.013*i*j)) for i in 1:nt, j in 1:K]
        ylms = ComplexF64[1.0+0.0im, 0.6-0.2im, 0.8+0.1im, 0.4+0.3im]
        Phi_phi = [0.02*(i-1)*100 for i in 1:nt]
        Phi_r   = [0.013*(i-1)*100 for i in 1:nt]
        f_phi = fill(2.0e-3, nt) .+ [1e-5*(i-1) for i in 1:nt]
        f_r   = fill(1.3e-3, nt) .+ [0.7e-5*(i-1) for i in 1:nt]
        m_arr = [2, 2]; n_arr = [0, 1]
        df = 2.0e-5; nf = 601
        f_arr = [df*(i - (nf+1)÷2) for i in 1:nf]
        hp, hc = FEWsum.fd_interpolated_mode_sum(t, teuk, ylms, Phi_phi, Phi_r, f_phi, f_r,
                                                 m_arr, n_arr; f_arr = f_arr)
        # captured legacy KB nonzero bins: (idx, Re hp, Im hp, Re hc, Im hc)
        cap = [
            (31, -2.01635917012076618e+02, -4.97686464926241911e+02,  9.19841225130812745e+01,  5.54335243467494934e+01),
            (32,  4.45864065321739133e+02, -2.99741828898832580e+02, -3.53692309836183654e+01,  1.01462380462421891e+02),
            (33,  3.78069974186625132e+02,  3.82026586753001652e+02, -1.06334386280220173e+02, -1.57558569780852338e+01),
            (34, -1.19272241716505823e+03,  1.05549734672761178e+03,  6.41759051390957893e+01, -3.12006265536224760e+02),
            (35, -4.79517980291258027e+02, -2.43468715867075616e+02,  1.05939122750650313e+02, -1.85871630962188874e+01),
            (98,  1.45193900453508786e+03,  1.39696080341560310e+03, -2.39792025262301081e+02,  7.03958433551262033e+01),
            (504, 1.45193900453508786e+03, -1.39696080341560310e+03, -2.39792025262301081e+02, -7.03958433551262033e+01),
            (567, -4.79517980291258027e+02, 2.43468715867075616e+02,  1.05939122750650313e+02,  1.85871630962188874e+01),
            (568, -1.19272241716505823e+03, -1.05549734672761178e+03, 6.41759051390957893e+01,  3.12006265536224760e+02),
            (569, 3.78069974186625132e+02, -3.82026586753001652e+02, -1.06334386280220173e+02, 1.57558569780852338e+01),
            (570, 4.45864065321739133e+02,  2.99741828898832580e+02, -3.53692309836183654e+01, -1.01462380462421891e+02),
            (571, -2.01635917012076618e+02, 4.97686464926241911e+02,  9.19841225130812745e+01, -5.54335243467494934e+01),
        ]
        nzset = Set(c[1] for c in cap)
        # legacy KB also emits NaN at the turnover bins (fdot→0 SPA singularity); GrTool reproduces them.
        nanset = Set([100, 101, 501, 502])
        for (idx, hpr, hpi, hcr, hci) in cap
            @test hp[idx] ≈ complex(hpr, hpi) rtol = 1e-10
            @test hc[idx] ≈ complex(hcr, hci) rtol = 1e-10
        end
        for i in nanset
            @test isnan(hp[i]) || isnan(hc[i])
        end
        # all other bins are exactly zero
        for i in 1:nf
            if !(i in nzset) && !(i in nanset)
                @test hp[i] == 0.0 && hc[i] == 0.0
            end
        end
    end
end
