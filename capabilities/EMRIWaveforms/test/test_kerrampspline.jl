# KerrAmpSpline characterization — clamped tensor cubic B-spline evaluator (FITPACK fpbspl/splev /
# scipy fpbisp-bisplev equivalent) + the ZNAmps Kerr-grid HDF5 loader.
#
# Grounding (NO Python — pure analytic B-spline identities + structural loader checks):
#  (1) partition of unity   — coeffs ≡ 1 ⇒ spline ≡ 1 everywhere (clamped knots sum-to-one), rtol~eps.
#  (2) linear reproduction  — coeffs = Greville abscissae τ_i = mean(t[i+1..i+k]) ⇒ spline(x) = x exactly.
#  (3) cubic exact repro    — coeffs = polar form (blossom) of a degree-≤3 polynomial P ⇒ spline = P
#      (blossom of x^n on knots (t_{i+1},t_{i+2},t_{i+3}): 1, mean, e2/3, product), and the tensor
#      product P(w)Q(u) for bisplev_clamped.
#  (4) C2 / derivative continuity across an interior knot (cubic B-spline is C^2): left/right finite-
#      difference 1st & 2nd derivatives agree.
#  Loader (read_kerrampgrid) grounded STRUCTURALLY: n_coeff = n_knots − 4 (37→33, 21→17), 33²=1089,
#  17²=289, root attrs lmax=10/mmax=10/nmax=55 ⇒ 6993 (l,m,n) modes; HDF5.jl reverses dims so the
#  CoeffsRegion axis layout is (coefflat, ReIm, mode, z); a read slab → cmat → bisplev is finite.

using Test
using EMRIWaveforms.AmplitudeModels.KerrAmpSpline
using EMRIWaveforms.AmplitudeModels.KerrAmpSpline: _find_interval, _flat_to_cmat, kerr_mode_table
using TestKit: oracle_root, oracle_or_skip

const _FEWDATA_KAS = oracle_root("GRTOOLS_FEWDATA_DIR",
    "/data/haiyangw/claude/emri-disk/ref-code/FastEMRIWaveforms/.fewdata")
const _ZNAMPS = joinpath(_FEWDATA_KAS, "ZNAmps_l10_m10_n55_DS2Outer.h5")

# --- analytic helpers (oracle: closed-form B-spline polar form) -----------------------------
# Greville abscissa of coeff i (1-based) for degree k: mean of k interior knots t[i+1..i+k].
_greville(t, i; k = 3) = sum(@view t[(i + 1):(i + k)]) / k

# polar form (blossom) of a cubic P(x)=a3 x^3 + a2 x^2 + a1 x + a0 evaluated at the three knot-averages
# of coeff i — the exact B-spline coefficient that reproduces P (degree ≤ 3, clamped cubic).
function _cubic_blossom_coeff(t, i, a0, a1, a2, a3)
    ta = t[i + 1]; tb = t[i + 2]; tc = t[i + 3]
    e1 = ta + tb + tc
    e2 = ta * tb + ta * tc + tb * tc
    e3 = ta * tb * tc
    return a3 * e3 + a2 * (e2 / 3) + a1 * (e1 / 3) + a0
end

# a representative clamped cubic knot vector: 4×0, 8 interior, 4×1 ⇒ 16 knots ⇒ 12 coeffs
const _TK = vcat(fill(0.0, 4), [0.12, 0.24, 0.36, 0.48, 0.6, 0.72, 0.84, 0.92], fill(1.0, 4))

@testset "KerrAmpSpline" begin

    ncoef = length(_TK) - 4
    @test ncoef == 12

    # =========================================================================
    # (1) partition of unity: coeffs ≡ 1 ⇒ spline ≡ 1 on [t_{k+1}, t_{n-k}] = [0,1]
    # =========================================================================
    @testset "deboor partition-of-unity" begin
        c1 = ones(ncoef)
        for x in 0.0:0.037:1.0
            @test deboor_cubic(_TK, c1, x) ≈ 1.0 atol = 1e-13
        end
        @test deboor_cubic(_TK, c1, 1.0) ≈ 1.0 atol = 1e-13
    end

    # =========================================================================
    # (2) linear reproduction: coeffs = Greville abscissae ⇒ spline(x) = x
    # =========================================================================
    @testset "deboor linear reproduction" begin
        cg = [_greville(_TK, i) for i in 1:ncoef]
        for x in 0.0:0.041:1.0
            @test deboor_cubic(_TK, cg, x) ≈ x atol = 1e-12
        end
    end

    # =========================================================================
    # (3) cubic-polynomial exact reproduction (degree ≤ 3 closed form via blossom)
    # =========================================================================
    @testset "deboor cubic exact reproduction" begin
        a0, a1, a2, a3 = -0.7, 1.3, -2.1, 0.9
        P(x) = a0 + a1 * x + a2 * x^2 + a3 * x^3
        cc = [_cubic_blossom_coeff(_TK, i, a0, a1, a2, a3) for i in 1:ncoef]
        for x in 0.0:0.033:1.0
            @test deboor_cubic(_TK, cc, x) ≈ P(x) atol = 1e-12
        end
    end

    # =========================================================================
    # (4) C2 / derivative continuity across an interior knot
    # =========================================================================
    @testset "deboor C2 continuity" begin
        cc = [sin(2.3 * i) + 0.4 * i for i in 1:ncoef]   # arbitrary smooth-ish coeffs
        xk = _TK[8]                                       # an interior knot
        f(x) = deboor_cubic(_TK, cc, x)
        h = 1e-4
        # 2nd-order one-sided derivative limits AT the knot (both estimate f'(xk)/f''(xk));
        # for a C^2 spline the left and right limits coincide.
        d1L = (3f(xk) - 4f(xk - h) + f(xk - 2h)) / (2h)
        d1R = (-3f(xk) + 4f(xk + h) - f(xk + 2h)) / (2h)
        # continuous f' ⇒ limits agree to truncation O(h²·Δf'''); a jump would be O(1).
        @test d1L ≈ d1R atol = 1e-3
        f2L = (f(xk) - 2f(xk - h) + f(xk - 2h)) / h^2
        f2R = (f(xk) - 2f(xk + h) + f(xk + 2h)) / h^2
        @test f2L ≈ f2R atol = 1e-2 * max(1, abs(f2L))
    end

    # =========================================================================
    # bisplev_clamped: tensor identities
    # =========================================================================
    @testset "bisplev tensor identities" begin
        tw = _TK; tu = vcat(fill(0.0, 4), [0.2, 0.4, 0.6, 0.8], fill(1.0, 4))  # 12 knots ⇒ 8 u-coeffs
        nwc = length(tw) - 4; nuc = length(tu) - 4
        @test (nwc, nuc) == (12, 8)

        # partition of unity (tensor): coeffs ≡ 1 ⇒ surface ≡ 1
        C1 = ones(nwc, nuc)
        for w in (0.05, 0.37, 0.73, 0.99), u in (0.1, 0.5, 0.9)
            @test bisplev_clamped(tw, tu, C1, w, u) ≈ 1.0 atol = 1e-13
        end

        # separable reproduction: c[i,j] = greville_w_i ⇒ surface = w ; c[i,j] = greville_u_j ⇒ surface = u
        Cw = [_greville(tw, i) for i in 1:nwc, j in 1:nuc]
        Cu = [_greville(tu, j) for i in 1:nwc, j in 1:nuc]
        for w in (0.12, 0.55, 0.88), u in (0.15, 0.65, 0.95)
            @test bisplev_clamped(tw, tu, Cw, w, u) ≈ w atol = 1e-12
            @test bisplev_clamped(tw, tu, Cu, w, u) ≈ u atol = 1e-12
        end

        # tensor cubic product P(w)*Q(u): c[i,j] = blossomP_i * blossomQ_j ⇒ surface = P(w)Q(u)
        pa = (0.3, -1.1, 0.7, 0.5); qa = (-0.2, 0.9, -1.4, 0.6)
        P(x) = pa[1] + pa[2]x + pa[3]x^2 + pa[4]x^3
        Q(x) = qa[1] + qa[2]x + qa[3]x^2 + qa[4]x^3
        cP = [_cubic_blossom_coeff(tw, i, pa...) for i in 1:nwc]
        cQ = [_cubic_blossom_coeff(tu, j, qa...) for j in 1:nuc]
        Cpq = [cP[i] * cQ[j] for i in 1:nwc, j in 1:nuc]
        for w in (0.18, 0.5, 0.82), u in (0.22, 0.6, 0.91)
            @test bisplev_clamped(tw, tu, Cpq, w, u) ≈ P(w) * Q(u) atol = 1e-11
        end
    end

    # =========================================================================
    # _find_interval / _flat_to_cmat unit checks
    # =========================================================================
    @testset "interval + reshape helpers" begin
        # clamped: valid l ∈ [k+1, n-k-1]
        @test _find_interval(_TK, 3, 0.0) == 4
        @test _find_interval(_TK, 3, 1.0) == length(_TK) - 4
        # C-order reshape: cmat[i,j] = flat[(i-1)*nu + j]
        nw, nu = 3, 4
        flat = collect(1.0:(nw * nu))
        cm = _flat_to_cmat(flat, nw, nu)
        @test size(cm) == (nw, nu)
        for i in 1:nw, j in 1:nu
            @test cm[i, j] == flat[(i - 1) * nu + j]
        end
    end

    # =========================================================================
    # type-genericity (BigFloat inheritance)
    # =========================================================================
    @testset "BigFloat genericity" begin
        tb = BigFloat.(_TK)
        cb = ones(BigFloat, length(_TK) - 4)
        v = deboor_cubic(tb, cb, BigFloat("0.371"))
        @test v isa BigFloat
        @test v ≈ one(BigFloat) atol = BigFloat("1e-40")
    end

    # =========================================================================
    # read_kerrampgrid — structural grounding vs the 4.8GB ZNAmps grid
    # =========================================================================
    @testset "read_kerrampgrid structural" begin
        if oracle_or_skip(isfile(_ZNAMPS), "read_kerrampgrid structural vs the ZNAmps grid"; env = "GRTOOLS_FEWDATA_DIR")
            # mode table replicates FEW KerrEccentricEquatorial ordering (m=0 first, then m>0)
            modes_all, index = kerr_mode_table(; lmax = 10, nmax = 55)
            @test length(modes_all) == 6993
            @test modes_all[1] == (2, 0, -55)          # first m=0 block, l=2, n=-55
            @test index[(2, 2, 0)] isa Int
            @test 1 <= index[(2, 2, 0)] <= 6993

            req = [(2, 2, 0), (3, 2, 1), (4, 0, -3)]
            g = read_kerrampgrid(_ZNAMPS; modes = req)

            @test g.lmax == 10 && g.mmax == 10 && g.nmax == 55
            @test g.num_total_modes == 6993
            @test g.modes == req
            @test length(g.z_knots) == 33
            @test g.z_knots[1] == 0.0 && g.z_knots[end] == 1.0
            # n_coeff = n_knots − 4
            @test length(g.twA) == 37 && length(g.tuA) == 37
            @test length(g.twB) == 21 && length(g.tuB) == 21
            # raw coeff slabs: (ncoefflat, 2, nmodes, nz)
            @test size(g.coeffsA) == (1089, 2, length(req), 33)
            @test size(g.coeffsB) == (289, 2, length(req), 33)

            # a read slab → cmat → bisplev is finite (mechanics sanity at an interior (w,u))
            nwcA = length(g.twA) - 4; nucA = length(g.tuA) - 4
            rec = _flat_to_cmat(g.coeffsA[:, 1, 1, 17], nwcA, nucA)
            imc = _flat_to_cmat(g.coeffsA[:, 2, 1, 17], nwcA, nucA)
            vr = bisplev_clamped(g.twA, g.tuA, rec, 0.3, 0.4)
            vi = bisplev_clamped(g.twA, g.tuA, imc, 0.3, 0.4)
            @test isfinite(vr) && isfinite(vi)
        end
    end
end
