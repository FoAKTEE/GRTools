# test_hyperboloidal.jl — ESSF node `hyperboloidal-transformation` (schedule W5-C), written
# BEFORE the implementation (TDD). The node delivers the hyperboloidal slicing KINEMATICS
# (arXiv:2606.04998 Sec. II.B.2, main.tex:255-262 tortoise + 267-314 height/boost/conjugation):
# the Kerr tortoise chart r_*(r) as a thin wrap of Spacetime.tortoise_kerr, and the
# height-function carrier (h, H = dh/dr_*, H′ = dH/dr_*) evaluated at Boyer–Lindquist r.
# Predecessor scope (NOT re-delivered here, only closure-tested against): the phase map
# Ψ_m = Ψ̃_m e^{−imΩh} + inverse landed in ansatz.jl (`psi_from_tilde`/`tilde_from_psi`,
# node mode-ansatz-decomposition) and the (H, H′)-parametrized flux-form coefficients landed
# in tsliced.jl (node elliptic-operator-tsliced).
#
# Gates (node contract, GLOBAL_DAG row 27):
#  (a) r_*(r) literals — `rstar` equals an independent test-side transcription of the
#      main.tex:257 log form r + (2M/(r₊−r₋))[r₊ln((r−r₊)/2M) − r₋ln((r−r₋)/2M)] (M = 1,
#      horizons recomputed here from a), 1e-13, incl. the a = 0 degeneracy (r₋ = 0 term
#      drops; Spacetime delegates to tortoise_schw) — plus the −∞ log approach at r₊.
#  (b) dr_*/dr = (r²+a²)/Δ (main.tex:260-262) — BigFloat(256) central FD of `rstar` vs the
#      direct transcription AND vs 1/flux_coeff_r (A^rr = Δ/(r²+a²) IS dr/dr_*, binding the
#      kinematics to the coefficient layer), 1e-13.
#  (c) H = 0 ⇒ transformed operator ≡ t-sliced operator — through the t_slicing()
#      composition path the coefficients and the phase map are BITWISE the t-sliced
#      defaults (zero tolerance).
#  (d) boost/height relation — the carrier invariant H = dh/dr_*, H′ = dH/dr_* holds for
#      the smooth tanh test profile by BigFloat(256) FD (1e-30; truncation ~h² = 1e-40),
#      H(±∞) = ±1 (main.tex:275-279), and t_slicing() is exact zeros.
#  (e) CONJUGATION IDENTITY (the row's oracle) — △_m[Ψ̃ e^{−imΩh}] = e^{−imΩh} △̃_m[Ψ̃]
#      at 1e-10 for the smooth tanh boost, m ∈ {0, 1, 2, 5}, a ∈ {0, 0.5, −0.7}, scattered
#      fixed (r, z) draws; both sides assembled from the EXPORTED coefficient layer with
#      BigFloat(256) central-FD divergences (h = 1e-20). The transformed field's analytic
#      ∂_r uses H·dr_*/dr for the phase derivative, so the identity genuinely binds the
#      h ↔ H ↔ r_* consistency of the kinematics to the (H, H′) slots of tsliced.jl.
#  (f) BigFloat(256) smoke — precision-256 rstar/kinematics through a full-BigFloat orbit,
#      identity (e) at 1e-40 (h = 1e-23: truncation ~h² = 1e-46 dominates roundoff).
#
# All local names are hyp_-prefixed: sibling test files share the runtests.jl include scope
# (test_operator_tsliced.jl owns cdiff/OPFIELDS/apply_flux_form), and this file must also
# run standalone as the node's verification command.
using Test
using ScalarSelfForce

# ---- BigFloat central FD ----------------------------------------------------------------
hyp_cdiff(f, x, h) = (f(x + h) - f(x - h)) / (2h)

# ---- smooth test boost: tanh profile with EXACT height antiderivative --------------------
# H(r_*) = tanh(k(r_*−c)) satisfies the hyperboloidal condition H(±∞) = ±1 (main.tex:275-279);
# h(r_*) = log(cosh(k(r_*−c)))/k is its exact antiderivative and H′ = k/cosh²(k(r_*−c)) its
# exact derivative — the carrier invariant is FD-gated in (d) before the oracle (e) uses it.
hyp_tanh_height(k, c) = HeightFunction(
    rs -> log(cosh(k * (rs - c))) / k,
    rs -> tanh(k * (rs - c)),
    rs -> k / cosh(k * (rs - c))^2,
)

# ---- smooth complex t-sliced test fields with EXACT first derivatives --------------------
# Ψ̃(r, z) = P(z)·g(r) + Q(z)·k(r) (complex z-polynomials × Gaussian + decaying exponential);
# only FIRST derivatives are analytic — second derivatives enter via hyp_cdiff of products.
# DIVISION-CONSISTENT dg: the conjugation identity (e)/(f) is sensitive to fr − ∂_r f (the
# residual is iκA·(fr − ∂_r f), κ = mΩH·dr_*/dr), so dg must be the EXACT derivative of the
# coded g at working precision — `sr^2` would be rounded in Float64 (≠ the sr-division chain
# of g at 1e-17 rel, breaking the 1e-40 BigFloat gate); divide by the same `sr` twice instead.
function hyp_field(pc, qc, rc, sr)
    dpc = ntuple(i -> i * pc[i + 1], length(pc) - 1)
    dqc = ntuple(i -> i * qc[i + 1], length(qc) - 1)
    g(r) = exp(-((r - rc) / sr)^2)
    dg(r) = -2 * ((r - rc) / sr) / sr * g(r)
    k(r) = exp(-r / (2rc))
    dk(r) = -k(r) / (2rc)
    (f  = (r, z) -> evalpoly(z, pc) * g(r) + evalpoly(z, qc) * k(r),
     dr = (r, z) -> evalpoly(z, pc) * dg(r) + evalpoly(z, qc) * dk(r),
     dz = (r, z) -> evalpoly(z, dpc) * g(r) + evalpoly(z, dqc) * k(r))
end

# deterministic "random-ish" complex coefficients (no Random dep in the test target)
const HYPFIELDS = (
    hyp_field((complex(0.71, -0.33), complex(0.28, 0.84), complex(-0.41, 0.17),
               complex(0.209, -0.114)),
              (complex(0.093, -0.502), complex(0.77, 0.36), complex(-0.25, 0.131)),
              8.0, 2.8),
    hyp_field((complex(-0.62, 0.45), complex(0.19, -0.57), complex(0.33, 0.08),
               complex(-0.147, 0.222)),
              (complex(0.54, 0.21), complex(-0.31, -0.66), complex(0.118, 0.29)),
              6.0, 3.2),
)

# ---- hyperboloidally-transformed field Ψ = Ψ̃·e^{−imΩh(r_*(r))}, analytic ∂_r/∂_z --------
# ∂_rΨ = (∂_rΨ̃ − imΩ·H(r_*(x))·(dr_*/dr)·Ψ̃)·e^{−imΩh}, dr_*/dr = (x²+a²)/Δ in the SAME
# factored-horizon form as the coefficient layer. Uses H = dh/dr_* — the boost/height
# relation gated in (d) — so identity (e) fails if the carrier's h and H are inconsistent.
function hyp_transformed_field(orb, m, hf, tld)
    Om = orb.Omega
    drs(x) = (x^2 + orb.a^2) / ((x - orb.rp) * (x - orb.rm))
    ph(x) = cis(-m * Om * hf.h(rstar(orb, x)))
    (f  = (x, y) -> tld.f(x, y) * ph(x),
     dr = (x, y) -> (tld.dr(x, y) -
                     im * m * Om * hf.H(rstar(orb, x)) * drs(x) * tld.f(x, y)) * ph(x),
     dz = (x, y) -> tld.dz(x, y) * ph(x))
end

# ---- flux-form application (uses ONLY the exported coefficient layer) ---------------------
function hyp_apply(orb, m, fld, r, z, H, Hp, h)
    -hyp_cdiff(x -> flux_coeff_r(orb, x) * fld.dr(x, z), r, h) -
        hyp_cdiff(y -> flux_coeff_z(orb, r, y) * fld.dz(r, y), z, h) +
        beta_m(orb, m, r, z, H, Hp) * fld.f(r, z) +
        gamma_r_m(orb, m, r, H) * fld.dr(r, z) +
        gamma_z_m(orb, m, r, z) * fld.dz(r, z)
end

@testset "MModeOperator slicing kinematics (hyperboloidal-transformation)" begin

    @testset "exports" begin
        for s in (:rstar, :HeightFunction, :t_slicing, :slicing_kinematics)
            @test Base.isexported(ScalarSelfForce, s)
        end
    end

    @testset "(a) r_*(r) log form vs main.tex:257 literals (1e-13)" begin
        for a in (0.0, 0.5, 0.998, -0.7)
            orb = CircularOrbit(a, a == 0.998 ? 6.0 : 10.0)
            # test-side transcription, M = 1 horizons recomputed from a (NOT orbit fields);
            # at a = 0 the r₋ term drops exactly (r₋ = 0), covering the tortoise_schw
            # delegation branch of the Spacetime wrap
            rp = 1 + sqrt(1 - a^2)
            rm = 1 - sqrt(1 - a^2)
            for r in (rp + 0.3, 4.7, 11.3, 28.0)
                ref = r + (2 / (rp - rm)) * (rp * log((r - rp) / 2) - rm * log((r - rm) / 2))
                @test abs(rstar(orb, r) - ref) <= 1e-13 * max(1.0, abs(ref))
            end
            # structural pins: log approach to −∞ at r₊ (this is why hyperboloidal slicing
            # is formulated in r_*, main.tex:265) and monotonicity
            @test rstar(orb, rp + 1e-9) < -30
            @test rstar(orb, 5.0) < rstar(orb, 5.1)
        end
    end

    @testset "(b) dr_*/dr = (r²+a²)/Δ (main.tex:260-262, 1e-13)" begin
        worst = 0.0
        setprecision(BigFloat, 256) do
            hstep = big"1e-20"
            for (a, r0) in ((0.0, 10.0), (0.5, 10.0), (0.998, 6.0), (-0.7, 10.0))
                orb = CircularOrbit(a, r0)
                for rr in (orb.rp + 0.4, orb.rp + 2.2, 9.4, 31.0)
                    r = BigFloat(rr)
                    fd = hyp_cdiff(x -> rstar(orb, x), r, hstep)
                    ref = (r^2 + a^2) / (r^2 - 2r + a^2)     # direct transcription, M = 1
                    dev1 = Float64(abs(fd - ref) / abs(ref))
                    # binds the kinematics to the coefficient layer: A^rr = Δ/(r²+a²)
                    # (flux_coeff_r) IS dr/dr_*, the inverse of this identity
                    dev2 = Float64(abs(fd * flux_coeff_r(orb, r) - 1))
                    worst = max(worst, dev1, dev2)
                    @test dev1 <= 1e-13
                    @test dev2 <= 1e-13
                end
            end
        end
        println("hyperboloidal gate (b) worst dr*/dr rel dev: ", worst)
    end

    @testset "(c) H = 0 ⇒ transformed ≡ t-sliced (bitwise)" begin
        ts = t_slicing()
        for (a, r0) in ((0.0, 10.0), (0.5, 10.0), (-0.7, 10.0))
            orb = CircularOrbit(a, r0)
            for m in (0, 1, 2, 5), (r, z) in ((orb.rp + 0.9, 0.41), (13.7, -0.28))
                kin = slicing_kinematics(ts, orb, r)
                @test kin.h == 0 && kin.H == 0 && kin.Hprime == 0
                # coefficient layer through the composition path == t-sliced defaults
                @test beta_m(orb, m, r, z, kin.H, kin.Hprime) == beta_m(orb, m, r, z)
                @test gamma_r_m(orb, m, r, kin.H) == gamma_r_m(orb, m, r)
                # phase map (ansatz.jl predecessor): e^{−imΩ·0} = 1 ⇒ Ψ ≡ Ψ̃ exactly
                psi = complex(0.83, -0.44)
                @test psi_from_tilde(m, psi, orb.Omega, kin.h) == psi
                @test tilde_from_psi(m, psi, orb.Omega, kin.h) == psi
            end
        end
    end

    @testset "(d) boost/height relation of the carrier (FD, 1e-30)" begin
        setprecision(BigFloat, 256) do
            hstep = big"1e-20"
            hf = hyp_tanh_height(0.35, 9.0)
            for rs in (big"-4.2", big"7.9", big"9.0", big"23.6")
                @test abs(hyp_cdiff(hf.h, rs, hstep) - hf.H(rs)) <= big"1e-30"
                @test abs(hyp_cdiff(hf.H, rs, hstep) - hf.Hprime(rs)) <= big"1e-30"
            end
            # hyperboloidal condition H(±∞) = ±1 (main.tex:275-279) for the test profile
            @test abs(hf.H(big"1e3") - 1) <= big"1e-70"
            @test abs(hf.H(big"-1e3") + 1) <= big"1e-70"
        end
        # t_slicing: exact zeros of the argument's type (h ≡ 0 ⇒ H ≡ 0 ⇒ H′ ≡ 0)
        ts = t_slicing()
        @test ts.h(3.7) === 0.0 && ts.H(3.7) === 0.0 && ts.Hprime(3.7) === 0.0
        @test ts.h(big"3.7") isa BigFloat
    end

    @testset "(e) conjugation identity △_m[Ψ̃e^{−imΩh}] = e^{−imΩh}△̃_m[Ψ̃] (1e-10)" begin
        worst = 0.0
        setprecision(BigFloat, 256) do
            hstep = big"1e-20"
            hf = hyp_tanh_height(0.35, 9.0)
            for (a, r0) in ((0.0, 10.0), (0.5, 10.0), (-0.7, 10.0))
                orb = CircularOrbit(a, r0)
                # scattered fixed (r, z) draws (deterministic; near-horizon, mid, far)
                pts = ((orb.rp + 0.8, 0.37), (orb.rp + 4.6, -0.63), (17.3, 0.21))
                for m in (0, 1, 2, 5), tld in HYPFIELDS, (rr, zz) in pts
                    r, z = BigFloat(rr), BigFloat(zz)
                    kin = slicing_kinematics(hf, orb, r)
                    psi = hyp_transformed_field(orb, m, hf, tld)
                    lhs = hyp_apply(orb, m, psi, r, z, kin.H, kin.Hprime, hstep)
                    rhs = cis(-m * orb.Omega * kin.h) *
                          hyp_apply(orb, m, tld, r, z, zero(r), zero(r), hstep)
                    dev = Float64(abs(lhs - rhs) / max(one(BigFloat), abs(rhs)))
                    worst = max(worst, dev)
                    @test dev <= 1e-10
                end
            end
        end
        println("hyperboloidal gate (e) worst conjugation-identity rel dev: ", worst)
    end

    @testset "(f) BigFloat(256) smoke: precision-256 kinematics + tight identity" begin
        setprecision(BigFloat, 256) do
            orb = CircularOrbit(big"0.5", big"10.0")
            r, z = big"7.3", big"0.44"
            rs = rstar(orb, r)
            @test rs isa BigFloat && precision(rs) == 256
            hf = hyp_tanh_height(0.35, 9.0)
            kin = slicing_kinematics(hf, orb, r)
            @test kin.h isa BigFloat && precision(kin.h) == 256
            @test kin.H isa BigFloat && precision(kin.H) == 256
            @test kin.Hprime isa BigFloat && precision(kin.Hprime) == 256
            @test kin.h == hf.h(rs) && kin.H == hf.H(rs)   # composition = h∘r_*, H∘r_*
            # identity (e) with the fully consistent BigFloat orbit: truncation ~h² = 1e-46
            # dominates roundoff ~eps₂₅₆/h; gate 1e-40 (30 orders below the Float64 gate)
            hstep = big"1e-23"
            for m in (0, 5)
                tld = HYPFIELDS[1]
                psi = hyp_transformed_field(orb, m, hf, tld)
                lhs = hyp_apply(orb, m, psi, r, z, kin.H, kin.Hprime, hstep)
                rhs = cis(-m * orb.Omega * kin.h) *
                      hyp_apply(orb, m, tld, r, z, zero(r), zero(r), hstep)
                @test abs(lhs - rhs) <= big"1e-40" * max(one(BigFloat), abs(rhs))
            end
        end
    end
end
