# test_mode_ansatz.jl — ESSF node `mode-ansatz-decomposition` (schedule W3-C), written BEFORE
# the implementation (TDD). Gates (GLOBAL_DAG #25 oracle column + node contract):
#   * round-trip mode-decompose -> reconstruct at 1e-12: smooth axisymmetric field, single-mode
#     fields m in {1, 2, 5}, and a multi-mode random-coefficient REAL field (m = -5..5), all
#     decomposed with a TEST-SIDE trapezoidal/Fourier azimuthal quadrature (exact for trig
#     polynomials of degree < N/2 on the periodic grid) — the module only strips/applies the
#     ansatz factors (eq:mode_sum main.tex:228-231, eq:modeIntegral main.tex:233-236);
#   * conjugate symmetry on a random real field: Psi_{-m} = conj(Psi_m) in the locked
#     convention, both t-sliced and after the hyperboloidal phase e^{-i m Omega h}
#     (eq:toTilde main.tex:286-289);
#   * sin^{|m|}(theta) z-derivative closed forms vs BigFloat(256) central finite differences
#     (h = 1e-20: truncation ~h^2 = 1e-40, roundoff ~eps/h^2 = 1e-37 — gate 1e-30), plus the
#     exact edge pins (|m| = 2 second derivative == -2 including at the poles where the naive
#     s^{|m|-4} form is 0*Inf; |m| = 1 pole divergence signs; even-|m| finite pole values);
#   * a = 0 azimuthal-chart degeneracy: kerr_varphi(0, r, phi) == phi exactly, and for a != 0
#     the offset equals the Spacetime.kerr_varphi_offset closed form (main.tex:223-226);
#   * effsource fixture bridge (W2-D headers, fixtures/effsource/): the comoving Psi^P_m
#     fixture obeys the module's documented conjugate-symmetry relation and the
#     e^{-i m phi_p} phase structure, and the module's OWN forward integral produces the SAME
#     e^{-i m phi_p} structure on an analytic real field reflection-symmetric about the
#     particle azimuth (a = 0 config: the phi/varphi charts degenerate, no chart ambiguity).
#     The OPPOSITE rotation is asserted NOT to realign the phase — the spectre
#     CircularOrbit.hpp doxygen sign trap (audit-spectre-ssf §9: doc says e^{+im(phi-phi_p)}
#     for Psi^P_m but e^{-im(phi-phi_p)} for S^eff_m; the code applies ONE rotation to both);
#   * type-genericity: BigFloat(256) inheritance through all ansatz factors (gate 1e-70).
using Test
using ScalarSelfForce
using Spacetime: kerr_varphi_offset

# ---- test-side azimuthal quadrature (the module NEVER integrates) ----------------------
# Periodic trapezoid == discrete Fourier analysis: c_m = (1/2π)∫ Φ e^{-i m φ̄} dφ̄ sampled on
# N uniform nodes; exact (to roundoff) for trig polynomials of degree < N/2.
const NPHI = 64
azgrid(N = NPHI) = [-π + 2π * (j - 1) / N for j in 1:N]
fourier_cm(vals, m, phis = azgrid()) =
    sum(vals[j] * cis(-m * phis[j]) for j in eachindex(vals, phis)) / length(vals)

# deterministic "random-ish" complex coefficients (no Random dep in the test target)
const PSI_PLANTED = Dict(
    0 => complex(0.83, 0.0),          # m = 0 real for a real field
    1 => complex(0.91, -0.44),
    2 => complex(-0.37, 0.62),
    3 => complex(0.155, 0.201),
    4 => complex(-0.058, -0.113),
    5 => complex(0.021, 0.077),
)

@testset "MModeOperator ansatz (mode-ansatz-decomposition)" begin

    @testset "exports" begin
        for s in (:sintheta_of_z, :sinm_theta, :sinm_theta_dz, :sinm_theta_d2z,
                  :kerr_varphi, :mode_term, :mode_amplitude, :reconstruct_field,
                  :psi_from_tilde, :tilde_from_psi)
            @test Base.isexported(ScalarSelfForce, s)
        end
    end

    @testset "z-chart: sinθ branch + sin^{|m|} factors" begin
        # locked chart: z = cosθ, θ = acos(z) ∈ [0, π], sinθ = +√(1-z²) ≥ 0
        for th in (0.13, 0.7, π / 2, 2.2, 3.0)
            @test isapprox(sintheta_of_z(cos(th)), sin(th); rtol = 1e-14)
            for m in (0, 1, 2, 5, -3)
                @test isapprox(sinm_theta(m, cos(th)), sin(th)^abs(m); rtol = 1e-13)
            end
        end
        @test sintheta_of_z(1.0) == 0.0
        @test sintheta_of_z(-1.0) == 0.0
        @test sintheta_of_z(0.0) == 1.0
        # sin^0 θ ≡ 1 everywhere, poles included (the m = 0 ansatz has no angular factor)
        @test sinm_theta(0, 1.0) == 1.0
        @test sinm_theta(0, -1.0) == 1.0
        # |m| ≥ 1 factors vanish at the poles
        @test sinm_theta(1, 1.0) == 0.0
        @test sinm_theta(5, -1.0) == 0.0
        # BigFloat inheritance (values + precision)
        setprecision(BigFloat, 256) do
            zb = big"0.4"
            sb = sintheta_of_z(zb)
            @test sb isa BigFloat && precision(sb) == 256
            @test abs(sb - sin(acos(zb))) < big"1e-70"
        end
    end

    @testset "sin^{|m|} z-derivatives: closed forms vs BigFloat(256) FD" begin
        setprecision(BigFloat, 256) do
            h = big"1e-20"
            for m in (0, 1, 2, 3, 5, 8, -2, -5),
                z in BigFloat.((-95 // 100, -1 // 2, 0, 3 // 10, 9 // 10))
                f(x) = sinm_theta(m, x)
                d1_fd = (f(z + h) - f(z - h)) / (2h)
                d2_fd = (f(z + h) - 2f(z) + f(z - h)) / h^2
                d1 = sinm_theta_dz(m, z)
                d2 = sinm_theta_d2z(m, z)
                @test abs(d1 - d1_fd) <= big"1e-30" * max(one(BigFloat), abs(d1))
                @test abs(d2 - d2_fd) <= big"1e-30" * max(one(BigFloat), abs(d2))
            end
        end
        # exact structural pins (Float64)
        for z in (-0.9, 0.0, 0.37, 1.0)
            @test sinm_theta_dz(0, z) == 0.0
            @test sinm_theta_d2z(0, z) == 0.0
            @test sinm_theta_dz(2, z) == -2z          # d/dz (1-z²) = -2z exactly
            @test sinm_theta_d2z(2, z) == -2.0        # incl. poles: naive s^{|m|-4} is 0*Inf
        end
        # |m| = 1 pole divergences with the correct signs: d/dz √(1-z²) = -z/√(1-z²)
        @test sinm_theta_dz(1, 1.0) == -Inf
        @test sinm_theta_dz(1, -1.0) == Inf
        # even |m| ≥ 4 finite pole value: f = (1-z²)², f''(±1) = 8
        @test sinm_theta_d2z(4, 1.0) == 8.0
        @test sinm_theta_d2z(4, -1.0) == 8.0
    end

    @testset "azimuthal chart: φ̄ composed from Spacetime.kerr_varphi_offset" begin
        # a = 0 degeneracy: φ̄ ≡ φ EXACTLY (offset is a hard zero, not roundoff)
        for r in (2.1, 6.0, 50.0), phi in (0.0, 1.3, -2.0)
            @test kerr_varphi(0.0, r, phi) == phi
        end
        # a ≠ 0: composition equals offset, and offset equals the main.tex:225 closed form
        a, r, phi = 0.7, 8.0, 1.234
        rp = 1 + sqrt(1 - a^2)
        rm = 1 - sqrt(1 - a^2)
        lit = a / (rp - rm) * log((r - rp) / (r - rm))
        @test isapprox(kerr_varphi(a, r, phi) - phi, kerr_varphi_offset(a, r); rtol = 1e-15)
        @test isapprox(kerr_varphi(a, r, phi), phi + lit; rtol = 1e-14)
        # sign: prograde a > 0 ⇒ offset < 0 near the horizon (log term negative), odd in a
        @test kerr_varphi(a, rp + 1e-3, 0.0) < -3.0
        @test isapprox(kerr_varphi(-a, r, 0.0), -kerr_varphi(a, r, 0.0); rtol = 1e-14)
    end

    @testset "round-trip: decompose → mode_amplitude → reconstruct (1e-12)" begin
        orbit = CircularOrbit(0.7, 8.0)
        Omega = orbit.Omega
        t, r, z = 2.3, 8.5, 0.4
        s = sqrt(1 - z^2)
        phis = azgrid()

        @testset "smooth axisymmetric field" begin
            fval = exp(-(r - 8.0)^2 / 20) * (1 + z^2 / 3)     # Φ = f(r, z), no φ̄/t dependence
            vals = fill(complex(fval), NPHI)
            psi0 = mode_amplitude(0, fourier_cm(vals, 0); Omega, ts = t, r, z)
            @test isapprox(psi0, r * fval; rtol = 1e-13)
            for m in (1, 2, 5)
                @test abs(fourier_cm(vals, m)) <= 1e-13 * abs(fval)
            end
            recon = reconstruct_field([0], [psi0]; Omega, ts = t, r, z, varphi = 0.77)
            @test isapprox(recon, fval; rtol = 1e-13)
        end

        @testset "single-mode fields m ∈ {1, 2, 5}" begin
            for m0 in (1, 2, 5)
                psian = PSI_PLANTED[m0]
                # INDEPENDENT test-side construction of the ansatz (eq:mode_sum transcribed)
                vals = [psian * s^m0 / r * exp(im * m0 * (ph - Omega * t)) for ph in phis]
                got = mode_amplitude(m0, fourier_cm(vals, m0); Omega, ts = t, r, z)
                @test isapprox(got, psian; rtol = 1e-12)
                for m in (0, 1, 2, 5)
                    m == m0 && continue
                    @test abs(fourier_cm(vals, m)) <= 1e-14
                end
                for ph in (-2.9, -0.4, 0.0, 1.1, 2.8)
                    recon = reconstruct_field([m0], [got]; Omega, ts = t, r, z, varphi = ph)
                    direct = psian * s^m0 / r * exp(im * m0 * (ph - Omega * t))
                    @test isapprox(recon, direct; rtol = 1e-12)
                end
            end
        end

        @testset "multi-mode real field m = -5..5" begin
            ms = collect(-5:5)
            planted = [m >= 0 ? PSI_PLANTED[m] : conj(PSI_PLANTED[-m]) for m in ms]
            field(ph) = real(sum(planted[k] * s^abs(ms[k]) / r *
                                 exp(im * ms[k] * (ph - Omega * t)) for k in eachindex(ms)))
            vals = complex.(field.(phis))
            got = [mode_amplitude(m, fourier_cm(vals, m); Omega, ts = t, r, z) for m in ms]
            for k in eachindex(ms)
                @test isapprox(got[k], planted[k]; rtol = 1e-12, atol = 1e-14)
            end
            for ph in (-3.0, -1.2, 0.3, 0.9, 2.5)
                recon = reconstruct_field(ms, got; Omega, ts = t, r, z, varphi = ph)
                @test isapprox(recon, complex(field(ph)); rtol = 1e-12)
                @test abs(imag(recon)) <= 1e-12 * abs(recon)   # real field reconstructs real
            end
        end
    end

    @testset "conjugate symmetry: Ψ_{-m} = conj(Ψ_m), both slicings" begin
        orbit = CircularOrbit(0.7, 8.0)
        Omega = orbit.Omega
        t, r, z, h = 1.7, 6.3, -0.25, 3.9
        phis = azgrid()
        # random-ish REAL azimuthal trig polynomial (degree 5)
        coefs = [(0.7, 0.0), (1.1, -0.3), (-0.4, 0.8), (0.25, 0.15), (-0.09, -0.31), (0.05, 0.12)]
        field(ph) = sum(c[1] * cos(k * ph) + c[2] * sin(k * ph)
                        for (k, c) in zip(0:5, coefs))
        vals = complex.(field.(phis))
        tpsi = Dict(m => mode_amplitude(m, fourier_cm(vals, m); Omega, ts = t, r, z)
                    for m in -5:5)
        for m in 1:5
            # t-sliced: Ψ̃_{-m} = conj(Ψ̃_m)
            @test isapprox(tpsi[-m], conj(tpsi[m]); rtol = 1e-12)
            # hyperboloidal phase preserves it: Ψ_{-m} = conj(Ψ_m)
            @test isapprox(psi_from_tilde(-m, tpsi[-m], Omega, h),
                           conj(psi_from_tilde(m, tpsi[m], Omega, h)); rtol = 1e-12)
            # round trip tilde ↔ psi is the identity
            @test isapprox(tilde_from_psi(m, psi_from_tilde(m, tpsi[m], Omega, h), Omega, h),
                           tpsi[m]; rtol = 1e-15)
        end
        @test abs(imag(tpsi[0])) <= 1e-13 * abs(tpsi[0])   # m = 0 mode of a real field is real
    end

    @testset "hyperboloidal slicing consistency (eq:toTilde)" begin
        # Ψ_m e^{im(φ̄-Ωs)} == Ψ̃_m e^{im(φ̄-Ωt)} with s = t - h: the SAME physical Φ
        orbit = CircularOrbit(0.5, 10.0)
        Omega = orbit.Omega
        t, r, z, varphi = 4.1, 12.0, 0.6, -1.9
        for m in (-3, -1, 0, 2, 5), h in (0.0, -2.7, 5.3)
            tpsi = PSI_PLANTED[abs(m)]
            psi = psi_from_tilde(m, tpsi, Omega, h)
            @test isapprox(mode_term(m, psi; Omega, ts = t - h, r, z, varphi),
                           mode_term(m, tpsi; Omega, ts = t, r, z, varphi); rtol = 1e-13)
        end
    end

    @testset "effsource fixture bridge: Ψ^P_m phase convention (W2-D headers)" begin
        fxdir = joinpath(@__DIR__, "fixtures", "effsource")
        isdefined(Main, :EFFSOURCE_CONFIGS) || include(joinpath(fxdir, "effsource_configs.jl"))
        isdefined(Main, :EFFSOURCE_PHIS_M_COMOVING) ||
            include(joinpath(fxdir, "phis_m_comoving.jl"))
        M = EFFSOURCE_PHIS_M_COMOVING                     # 7 × N, records are columns
        lut = Dict{Tuple{Int, Int, Int}, ComplexF64}()
        for k in axes(M, 2)
            lut[(Int(M[1, k]), Int(M[2, k]), Int(M[3, k]))] = complex(M[6, k], M[7, k])
        end
        # (i) the fixture obeys the module's documented relation Ψ^P_{-m} = conj(Ψ^P_m)
        # (ii) Ψ^P_m e^{+im φ_p} is REAL: the e^{-im φ_p} structure of the locked forward
        #      integral for a real field reflection-symmetric about the particle azimuth
        for cfg in (1, 4), ipt in (1, 5, 20), m in (1, 2)
            psim = lut[(cfg, m, ipt)]
            @test isapprox(lut[(cfg, -m, ipt)], conj(psim); rtol = 1e-12)
            rot = psim * cis(m * EFFSOURCE_PHIP)
            @test abs(imag(rot)) <= 1e-12 * abs(rot)
        end
        # (iii) module tie-in at a = 0 (cfg 1: φ̄ ≡ φ, no chart ambiguity): OUR forward
        # integral on an analytic real field reflection-symmetric about φ_p = 0.31 yields
        # EXACTLY the fixture's e^{-im φ_p} phase structure. The OPPOSITE rotation (the
        # spectre doxygen trap sign) does NOT realign the phase.
        orbit = CircularOrbit(0.0, 10.0)
        Omega = orbit.Omega
        t, r, z = 0.0, 10.0, 0.3                          # fixture particle is at t = 0
        s = sqrt(1 - z^2)
        gam = [0.7, 1.1, -0.4, 0.25]                      # γ_0 + Σ γ_k cos(k(φ̄ - φ_p))
        field(ph) = gam[1] + sum(gam[k + 1] * cos(k * (ph - EFFSOURCE_PHIP)) for k in 1:3)
        vals = complex.(field.(azgrid()))
        for m in (-3, -2, -1, 1, 2, 3)
            got = mode_amplitude(m, fourier_cm(vals, m); Omega, ts = t, r, z)
            expected = gam[abs(m) + 1] / 2 * cis(-m * EFFSOURCE_PHIP) * r / s^abs(m)
            @test isapprox(got, expected; rtol = 1e-12)
            wrong = got * cis(-m * EFFSOURCE_PHIP)        # opposite-convention rotation
            @test abs(imag(wrong)) > 0.4 * abs(wrong)     # |sin(2m·0.31)| ≥ 0.58 for these m
        end
    end

    @testset "type-genericity: BigFloat(256) inheritance" begin
        setprecision(BigFloat, 256) do
            Omega = 1 / (big"0.7" + sqrt(big"8.0"^3))     # closed form, main.tex:196
            t, r, z, varphi = big"2.3", big"8.5", big"0.4", big"1.1"
            psian = complex(big"0.91", big"-0.44")
            for m in (0, 1, 5)
                term = mode_term(m, psian; Omega, ts = t, r, z, varphi)
                @test term isa Complex{BigFloat} && precision(real(term)) == 256
                # algebraic inverse: strip e^{imφ̄} → azimuthal coefficient → amplitude
                back = mode_amplitude(m, term * cis(-m * varphi); Omega, ts = t, r, z)
                @test abs(back - psian) <= big"1e-70" * abs(psian)
            end
            h = big"3.9"
            @test abs(tilde_from_psi(2, psi_from_tilde(2, psian, Omega, h), Omega, h) -
                      psian) <= big"1e-70" * abs(psian)
            @test sinm_theta(3, z) isa BigFloat && precision(sinm_theta(3, z)) == 256
        end
    end
end
