# test_operator_tsliced.jl — ESSF node `elliptic-operator-tsliced` (schedule W4-C), written
# BEFORE the implementation (TDD). The node delivers the m-mode elliptic operator COEFFICIENT
# layer in first-order flux form (eq:flux main.tex:366-368, eq:first-order main.tex:370-373,
# β_m/γ^i_m main.tex:375-383) with a CALLER-supplied boost (H, H′ = dH/dr_*) so one code path
# serves t/v/u slicing — t-slicing is the (H, H′) = (0, 0) default; the vtu piecewise assembly
# + junction conditions are a LATER node — plus the z² equatorial-symmetry variant
# 𝓕^{z²} = 4z²𝓕^z, γ^{z²} = 2zγ^z + 2(1−z²)/(r²+a²) (main.tex:389-394, factors READ from tex).
#
# Gates (node contract):
#  (a) INDEPENDENT cross-check — the flux-form operator −∂_i𝓕^i + β_mΨ + γ^i_m∂_iΨ (divergences
#      by BigFloat(256) central FD of coefficient×derivative products, h = 1e-20: truncation
#      ~h² = 1e-40, roundoff ~1e-57) equals a DIRECT test-side transcription of eq:triangle
#      (main.tex:247-251) + the hyperboloidal replacement rules (main.tex:306-308), with
#      ∂_{r_*} = Δ/(r²+a²)∂_r (main.tex:356), on smooth complex fields — gate 1e-9 relative
#      (prompted node gate; measured ~1e-15). The replacement-rule identity is algebraic and
#      POINTWISE in (H, H′), so the constant samples H ∈ {0, −1, +1, 0.3} plus (H, H′) =
#      (0.3, 0.2) cover general boost functions. (a, r0) ∈ {(0,10), (0.5,10), (0.998,6),
#      (−0.7,10)} incl. a < 0; m ∈ {0, 1, 2, 5}.
#  (b) horizon regularity (main.tex:386; DAG claim #12) — with H = −1 the 1/Δ pieces of β_m
#      cancel: β_m(r₊ + 10^−k), k = 4..10, stays within 1e-3 of the CLOSED-FORM horizon limit
#      transcribed here; the t-sliced (H = 0) β_m DIVERGES ~1/Δ on the same approach
#      (discriminator showing the cancellation gate is not vacuous).
#  (c) z² equatorial variant — on an even-in-z field the operator assembled in the w = z²
#      coordinate (flux_coeff_z2 + gamma_z2_m, w-FD divergence) equals the z-form at 1e-10.
#  (d) m → −m conjugation — coefficient_{−m} = conj(coefficient_m) for β_m/γ^r_m; γ^z_m and
#      γ^{z²}_m real and even in m. This is the |m| reading of the angular-origin terms (the
#      printed m(m+1), 2mz are the m ≥ 0 branch of |m|(|m|+1), 2|m|z coming from the
#      sin^{|m|}θ reduction), consistent with the ansatz layer's Ψ_{−m} = conj(Ψ_m).
#  (e) BigFloat(256) smoke — full-BigFloat orbit + coefficients, identity (a) at 1e-45
#      (h = 1e-23: FD truncation ~h² = 1e-46 dominates, roundoff ~eps₂₅₆/h = 1.7e-54),
#      Complex{BigFloat}/BigFloat outputs with precision 256.
#
# The eq:triangle reference transcription lives ONLY in this file (independent of the
# implementation, which transcribes the eq:flux / β / γ printed forms instead).
using Test
using ScalarSelfForce

# ---- BigFloat central FD ----------------------------------------------------------------
cdiff(f, x, h) = (f(x + h) - f(x - h)) / (2h)

# ---- smooth complex test fields with EXACT first derivatives -----------------------------
# Ψ(r, z) = P(z)·g(r) + Q(z)·k(r): complex z-polynomials times a Gaussian and a decaying
# exponential radial profile. Only FIRST derivatives are coded analytically; every second
# derivative in the operators below enters through cdiff of a first-derivative product.
function smoothfield(pc, qc, rc, sr)
    dpc = ntuple(i -> i * pc[i + 1], length(pc) - 1)
    dqc = ntuple(i -> i * qc[i + 1], length(qc) - 1)
    g(r) = exp(-((r - rc) / sr)^2)
    dg(r) = -2 * (r - rc) / sr^2 * g(r)
    k(r) = exp(-r / (2rc))
    dk(r) = -k(r) / (2rc)
    (f  = (r, z) -> evalpoly(z, pc) * g(r) + evalpoly(z, qc) * k(r),
     dr = (r, z) -> evalpoly(z, pc) * dg(r) + evalpoly(z, qc) * dk(r),
     dz = (r, z) -> evalpoly(z, dpc) * g(r) + evalpoly(z, dqc) * k(r))
end

# deterministic "random-ish" complex coefficients (no Random dep in the test target)
const OPFIELDS = (
    smoothfield((complex(0.83, -0.21), complex(0.40, 0.90), complex(-0.37, 0.62),
                 complex(0.155, 0.201)),
                (complex(-0.058, -0.113), complex(0.91, -0.44), complex(0.21, 0.077)),
                9.0, 2.5),
    smoothfield((complex(-0.50, 0.70), complex(0.31, -0.11), complex(0.09, 0.44),
                 complex(-0.23, 0.05)),
                (complex(0.66, 0.13), complex(-0.40, 0.52), complex(0.12, -0.31)),
                5.0, 3.5),
)

# ---- flux-form application (uses ONLY the exported coefficient layer) ---------------------
function apply_flux_form(orb, m, fld, r, z, H, Hp, h)
    -cdiff(x -> flux_coeff_r(orb, x) * fld.dr(x, z), r, h) -
        cdiff(x -> flux_coeff_z(orb, r, x) * fld.dz(r, x), z, h) +
        beta_m(orb, m, r, z, H, Hp) * fld.f(r, z) +
        gamma_r_m(orb, m, r, H) * fld.dr(r, z) +
        gamma_z_m(orb, m, r, z) * fld.dz(r, z)
end

# ---- INDEPENDENT reference: eq:triangle (main.tex:247-251) + replacement rules ------------
# (main.tex:306-308: ∂_{r_*} → ∂_{r_*} + imΩH, ∂²_{r_*} → ∂²_{r_*} + 2imΩH∂_{r_*} + imΩH′
# − m²Ω²H²), with ∂_{r_*} = Δ/(r²+a²)∂_r (main.tex:356). Transcribed AS PRINTED, so valid for
# m ≥ 0 (the m < 0 branch is covered by conjugation gate (d)). Δ is deliberately computed as
# r² − 2Mr + a² here, NOT the implementation's factored (r − r₊)(r − r₋).
function apply_triangle_ref(orb, m, fld, r, z, H, Hp, h)
    a, M, Om = orb.a, orb.M, orb.Omega
    Dl(x) = x^2 - 2M * x + a^2
    w(x) = x^2 + a^2
    psi = fld.f(r, z)
    dpsi_r = fld.dr(r, z)
    dpsi_z = fld.dz(r, z)
    D1 = Dl(r) / w(r) * dpsi_r                                          # ∂Ψ/∂r_*
    D2 = Dl(r) / w(r) * cdiff(x -> Dl(x) / w(x) * fld.dr(x, z), r, h)   # ∂²Ψ/∂r_*²
    D1h = D1 + im * m * Om * H * psi
    D2h = D2 + 2im * m * Om * H * D1 + (im * m * Om * Hp - m^2 * Om^2 * H^2) * psi
    Sig = r^2 + a^2 * z^2
    Zang = cdiff(x -> (1 - x^2) * fld.dz(r, x), z, h) - 2m * z * dpsi_z - m * (m + 1) * psi
    -(w(r) / Dl(r)) * D2h + (2a^2 / (r * w(r)) - 2im * m * a / Dl(r)) * D1h -
        (Zang + m^2 * Om * (2M * r * (Om * w(r) - 2a) + Om * Dl(r) * Sig) / Dl(r) * psi +
         2 * (a^2 - im * m * a * r - M * r) / r^2 * psi) / w(r)
end

@testset "MModeOperator flux-form coefficients (elliptic-operator-tsliced)" begin

    @testset "exports" begin
        for s in (:flux_coeff_r, :flux_coeff_z, :flux_coeff_z2,
                  :beta_m, :gamma_r_m, :gamma_z_m, :gamma_z2_m)
            @test Base.isexported(ScalarSelfForce, s)
        end
    end

    @testset "(a) flux form == eq:triangle + replacement rules (1e-9)" begin
        worst = 0.0
        setprecision(BigFloat, 256) do
            h = big"1e-20"
            for (a, r0) in ((0.0, 10.0), (0.5, 10.0), (0.998, 6.0), (-0.7, 10.0))
                orb = CircularOrbit(a, r0)
                pts = ((orb.rp + 0.6, 0.33), (orb.rp + 3.1, -0.71), (18.7, 0.52))
                for m in (0, 1, 2, 5),
                    (H, Hp) in ((0.0, 0.0), (-1.0, 0.0), (1.0, 0.0), (0.3, 0.0), (0.3, 0.2)),
                    fld in OPFIELDS, (rr, zz) in pts

                    r, z = BigFloat(rr), BigFloat(zz)
                    lhs = apply_flux_form(orb, m, fld, r, z, H, Hp, h)
                    ref = apply_triangle_ref(orb, m, fld, r, z, H, Hp, h)
                    dev = Float64(abs(lhs - ref) / max(one(BigFloat), abs(ref)))
                    worst = max(worst, dev)
                    @test dev <= 1e-9
                end
            end
        end
        println("operator-tsliced gate (a) worst flux-vs-triangle rel dev: ", worst)
    end

    @testset "(b) horizon regularity of beta_m at H = -1 (main.tex:386)" begin
        worst = 0.0
        for (a, r0, m) in ((0.5, 10.0, 1), (0.998, 6.0, 5), (-0.7, 10.0, 2))
            orb = CircularOrbit(a, r0)
            z = 0.37
            # closed-form horizon limit of β_m at (H, H′) = (−1, 0), M = 1, transcribed here:
            # the (w/Δ)·m²Ω²(H²−1) piece is exactly 0; (w/Δ)·(2am²Ω/w)(2Mr/w − 1) → −2am²Ω/w₊
            # via 2Mr₊ − r₊² − a² = −Δ(r₊) cancellation; the a²Δ(1−z²)/w² piece → finite.
            wp = orb.rp^2 + a^2
            blim = (m * (m + 1) + 2 / orb.rp - 2a^2 / orb.rp^2 +
                    2im * a * m / orb.rp * (1 - a * orb.Omega)) / wp +
                   m^2 * orb.Omega^2 * a^2 * (1 - z^2) / wp - 2a * m^2 * orb.Omega / wp
            for k in 4:10
                r = orb.rp + 10.0^-k
                b = beta_m(orb, m, r, z, -1.0, 0.0)
                @test isfinite(b)
                dev = abs(b - blim) / (1 + abs(blim))
                worst = max(worst, dev)
                @test dev <= 1e-3
            end
            # discriminator: WITHOUT the boost (t-slicing H = 0) β_m blows up ~ 1/Δ there,
            # so the H = −1 boundedness above is a genuine cancellation, not a weak gate
            @test abs(beta_m(orb, m, orb.rp + 1e-8, z, 0.0, 0.0)) > 1e3
        end
        println("operator-tsliced gate (b) worst |beta - horizon limit| rel dev: ", worst)
    end

    @testset "(c) z² equatorial-symmetry variant (main.tex:389-394)" begin
        # even field u(r, w), w = z²: u = (c0 + c2 w + c4 w²)·g(r), exact ∂_r, ∂_w
        c = (complex(0.70, -0.20), complex(-0.45, 0.31), complex(0.18, 0.09))
        dc = (c[2], 2 * c[3])
        g(r) = exp(-((r - 8) / 3)^2)
        dg(r) = -2 * (r - 8) / 9 * g(r)
        uf(r, w) = evalpoly(w, c) * g(r)
        ur(r, w) = evalpoly(w, c) * dg(r)
        uw(r, w) = evalpoly(w, dc) * g(r)
        worst = 0.0
        setprecision(BigFloat, 256) do
            h = big"1e-20"
            for (a, r0) in ((0.5, 10.0), (0.998, 6.0)), m in (0, 2, 5),
                (H, Hp) in ((0.0, 0.0), (-1.0, 0.0)), zz in (0.3, 0.7)

                orb = CircularOrbit(a, r0)
                r, z = big"7.3", BigFloat(zz)
                wv = z^2
                # z-form of the SAME field: Ψ(r, z) = u(r, z²), ∂_zΨ = 2z ∂_w u
                fldz = (f = (x, y) -> uf(x, y^2), dr = (x, y) -> ur(x, y^2),
                        dz = (x, y) -> 2y * uw(x, y^2))
                Lz = apply_flux_form(orb, m, fldz, r, z, H, Hp, h)
                # w-form: divergence in w with the z²-variant coefficients (∂_w slot)
                Lw = -cdiff(x -> flux_coeff_r(orb, x) * ur(x, wv), r, h) -
                     cdiff(x -> flux_coeff_z2(orb, r, sqrt(x)) * uw(r, x), wv, h) +
                     beta_m(orb, m, r, z, H, Hp) * uf(r, wv) +
                     gamma_r_m(orb, m, r, H) * ur(r, wv) +
                     gamma_z2_m(orb, m, r, z) * uw(r, wv)
                dev = Float64(abs(Lz - Lw) / max(one(BigFloat), abs(Lz)))
                worst = max(worst, dev)
                @test dev <= 1e-10
            end
        end
        println("operator-tsliced gate (c) worst z2-vs-z rel dev: ", worst)
    end

    @testset "(d) m → −m conjugation of the coefficient layer" begin
        for (a, r0) in ((0.5, 10.0), (-0.7, 10.0)), m in (1, 2, 5)
            orb = CircularOrbit(a, r0)
            for (r, z, H, Hp) in ((7.7, 0.41, 0.0, 0.0), (3.2, -0.6, -1.0, 0.0),
                                  (12.0, 0.11, 0.3, 0.2))
                @test isapprox(beta_m(orb, -m, r, z, H, Hp),
                               conj(beta_m(orb, m, r, z, H, Hp)); rtol = 1e-15)
                @test isapprox(gamma_r_m(orb, -m, r, H),
                               conj(gamma_r_m(orb, m, r, H)); rtol = 1e-15)
                # angular-origin coefficients are REAL and even in m (|m| reading)
                @test gamma_z_m(orb, -m, r, z) == gamma_z_m(orb, m, r, z)
                @test gamma_z2_m(orb, -m, r, z) == gamma_z2_m(orb, m, r, z)
            end
        end
        # m = 0: no azimuthal phase, coefficients purely real (even with H′ ≠ 0)
        orb0 = CircularOrbit(0.9, 8.0)
        @test imag(beta_m(orb0, 0, 6.6, 0.5, 0.3, 0.2)) == 0
        @test imag(gamma_r_m(orb0, 0, 6.6, 0.3)) == 0
    end

    @testset "structural pins: boundary zeros, signs, t-slicing defaults" begin
        for (a, r0) in ((0.0, 10.0), (0.5, 10.0), (0.998, 6.0), (-0.7, 10.0))
            orb = CircularOrbit(a, r0)
            # 𝓕^r coefficient vanishes AT the horizon — exact zero via factored Δ off the
            # cached r₊ (main.tex:385 + 399-400: no explicit horizon BC needed)
            @test flux_coeff_r(orb, orb.rp) == 0
            # 𝓕^z coefficient vanishes at the poles exactly (sin^{|m|}θ regularity)
            @test flux_coeff_z(orb, 7.0, 1.0) == 0
            @test flux_coeff_z(orb, 7.0, -1.0) == 0
            # ellipticity: both principal coefficients positive in the open domain
            @test flux_coeff_r(orb, orb.rp + 1.3) > 0
            @test flux_coeff_z(orb, 7.0, 0.3) > 0
            # z²-variant: A^{z²} = 4z²A^zz (main.tex:391), zero at the equator w = 0 AND
            # the pole w = 1 (the even half-domain needs no explicit boundary data either)
            @test flux_coeff_z2(orb, 7.0, 0.5) == 4 * 0.25 * flux_coeff_z(orb, 7.0, 0.5)
            @test flux_coeff_z2(orb, 7.0, 0.0) == 0
            @test flux_coeff_z2(orb, 7.0, 1.0) == 0
            # γ^{z²} = 2zγ^z + 2(1−z²)/(r²+a²) (main.tex:392)
            z = 0.62
            @test isapprox(gamma_z2_m(orb, 3, 7.0, z),
                           2z * gamma_z_m(orb, 3, 7.0, z) + 2 * (1 - z^2) / (49 + orb.a^2);
                           rtol = 1e-15)
            # t-slicing = the (H, H′) = (0, 0) defaults, bitwise
            @test beta_m(orb, 2, 7.0, z) == beta_m(orb, 2, 7.0, z, 0.0, 0.0)
            @test gamma_r_m(orb, 2, 7.0) == gamma_r_m(orb, 2, 7.0, 0.0)
        end
    end

    @testset "(e) BigFloat(256) smoke: full-precision orbit + identity" begin
        setprecision(BigFloat, 256) do
            orb = CircularOrbit(big"0.5", big"10.0")
            @test orb isa CircularOrbit{BigFloat}
            r, z = big"7.1", big"0.44"
            H, Hp = big"0.3", big"0.2"
            b = beta_m(orb, 3, r, z, H, Hp)
            gr = gamma_r_m(orb, 3, r, H)
            gz = gamma_z_m(orb, 3, r, z)
            @test b isa Complex{BigFloat} && precision(real(b)) == 256
            @test gr isa Complex{BigFloat} && precision(real(gr)) == 256
            @test gz isa BigFloat && precision(gz) == 256
            @test flux_coeff_r(orb, r) isa BigFloat
            @test flux_coeff_z2(orb, r, z) isa BigFloat
            # h chosen for the 2nd-order-FD error floor at 256 bits: truncation ~h² = 1e-46
            # dominates roundoff ~eps₂₅₆/h = 1.7e-54; gate 1e-45 (30 orders below Float64)
            h = big"1e-23"
            for fld in OPFIELDS, m in (0, 3)
                lhs = apply_flux_form(orb, m, fld, r, z, H, Hp, h)
                ref = apply_triangle_ref(orb, m, fld, r, z, H, Hp, h)
                @test abs(lhs - ref) <= big"1e-45" * max(one(BigFloat), abs(ref))
            end
        end
    end
end
