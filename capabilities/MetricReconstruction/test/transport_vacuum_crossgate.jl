# Grounding for MetricReconstruction.Transport — the vacuum cross-gate G-E-vac (dossier §4.3).
#
# The seven LY transport equations reduce IN VACUUM (h_mm̄ = Φ_ij = μ^{(1)} = 0) to the RLGP-F /
# D4 assembler algebra. The HONEST, chart-independent comparable (the charts differ — LY is BL-
# Kinnersley with explicit NP-γ, D4 is HPHC-cl with γ absorbed into the GHP þ′ boost weight, so a
# full POINTWISE cross-chart coefficient identity needs the GHP boost-map derivation = [FUTURE],
# same reasoning as D4 G-D4f) is:
#   (1) each Δ_{[a,b,c,d]} operator expands to the published RLGP-F vacuum operator
#       (a·μ+b·μ̄+c·γ+d·γ̄), verified via `weight_combo` on the NPGHP table;
#   (2) the generalization SOURCE terms (Φ_ij, R_g, R_h, −3Ψ2 μ^{(1)}) all VANISH in vacuum, so the
#       chain collapses to the RLGP-F 7-equation form (= the D4 recon_set_*_k! sources).

using Test
using MetricReconstruction.Transport
import Spacetime.NPGHP as NPGHP

@testset "Transport — vacuum cross-gate (G-E-vac, dossier §4.3)" begin

    # The RLGP-F vacuum operator coefficient combinations (a·μ+b·μ̄+c·γ+d·γ̄), quoted from the
    # Formalism.tex forms cited in the dossier §2 (F: line numbers). Keyed by LY equation label.
    rlgp_combo(np, a, b, c, d) = a * np.μ + b * conj(np.μ) + c * np.γ + d * conj(np.γ)
    expected = (
        mu_transport_final = (2, 0, 1, 1),    # Δ+2μ+γ+γ̄            (LY:230)
        trace_transport    = (-1, 1, 0, 0),   # Δ−μ+μ̄               (LY:235)
        lambda_transport   = (1, 1, 3, -1),   # Δ+μ+μ̄+3γ−γ̄  (F:1162, B1)
        hmm_transport      = (-1, 1, 2, -2),  # Δ−μ+μ̄+2γ−2γ̄ (F:1168, B2)
        psi3_transport     = (4, 0, 2, 0),    # Δ+4μ+2γ      (F:1182, B3)
        pi_transport       = (0, 0, 1, -1),   # Δ+γ−γ̄       (F:1204, B4)
        hlm_transport      = (0, 1, 0, -2),   # Δ+μ̄−2γ̄      (F:1213, B5)
        psi2_transport     = (3, 0, 0, 0),    # Δ+3μ        (F:1222, B6)
        hll_transport      = (-1, 1, -2, -2), # Δ−μ+μ̄−2γ−2γ̄ (F:2437, B7)
    )

    @testset "G-E-vac-1  Δ-operators expand to the published RLGP-F vacuum forms" begin
        # structural (integer) quartets match the dossier §2 documented labels
        for (name, quad) in pairs(expected)
            op = getproperty(CHAIN_DELTA_OPS, name)
            @test (op.a, op.b, op.c, op.d) == quad
        end
        # and the spin-coefficient EXPANSION matches (a·μ+b·μ̄+c·γ+d·γ̄) on the NPGHP Kerr table,
        # sweeping a≠0 and several (r,θ) — chart-independent operator-coefficient identity.
        for a in (0.0, 0.5, 0.9), r in (3.0, 6.0), θ in (0.6, 1.4)
            np = NPGHP.np_spincoeff(a, r, θ)
            for (name, quad) in pairs(expected)
                op = getproperty(CHAIN_DELTA_OPS, name)
                @test isapprox(weight_combo(op, np), rlgp_combo(np, quad...); atol = 1e-13)
            end
        end
        # genuine (non-tautological) check that weight_combo implements the FULL leg map of LY
        # eq:derivative-notation (LY:204–211) for ALL FOUR legs — the quartet-to-spin-coeff dict:
        #   D→{ε,ε̄,ρ,ρ̄}  Δ→{μ,μ̄,γ,γ̄}  δ→{ᾱ,β,π̄,τ}  δ̄→{α,β̄,π,τ̄}.
        np = NPGHP.np_spincoeff(0.7, 4.0, 1.1)
        (A, B, C, D_) = (2, -1, 3, -2)
        @test isapprox(weight_combo(Dw(A, B, C, D_), np),
                       A * np.ε + B * conj(np.ε) + C * np.ρ + D_ * conj(np.ρ); atol = 1e-13)
        @test isapprox(weight_combo(Δw(A, B, C, D_), np),
                       A * np.μ + B * conj(np.μ) + C * np.γ + D_ * conj(np.γ); atol = 1e-13)
        @test isapprox(weight_combo(δw(A, B, C, D_), np),
                       A * conj(np.α) + B * np.β + C * conj(np.π) + D_ * np.τ; atol = 1e-13)
        @test isapprox(weight_combo(δbarw(A, B, C, D_), np),
                       A * np.α + B * conj(np.β) + C * np.π + D_ * conj(np.τ); atol = 1e-13)
    end

    @testset "G-E-vac-2  generalization source terms vanish in vacuum" begin
        T = Float64
        R = 6.0
        zp = ShellProfile{T}(lzero(T), lzero(T), 0.0, 0.0, T(R))   # a zero (vacuum) NP-Ricci source
        # np_ricci_from_T of a zero stress-energy is zero (Φ_ij = Λ = 0)
        φ = np_ricci_from_T(0.0, 0.0, 0.0, 0.0)
        @test φ.Φ00 == 0 && φ.Φ22 == 0 && φ.Φ11 == 0 && φ.Λ == 0
        # R_g reduced from ZERO Φ sources vanishes (−D_{[2,2,0,-1]}Φ22 − 2ΔΛ − 2μΦ11 = 0)
        np = schw_np_laurent(T, 1.0)
        f = schw_f(T, 1.0)
        Rg0 = pscale(-1.0, pderiv(zp) - psmooth(np.ρ, zp)) -
              pscale(-1.0, psmooth(f, pderiv(zp))) - psmooth(lscale(2, np.μ), zp)
        @test isempty(Rg0.inside.c) && isempty(Rg0.outside.c) && Rg0.dc == 0 && Rg0.dpc == 0
        # R_h (monopole) is identically zero (Φ21=0, spherical) — the vacuum reduction, ANCHOR-5
        @test let rh = Rh_shell(T, 1.0, R)
            isempty(rh.inside.c) && isempty(rh.outside.c) && rh.dc == 0 && rh.dpc == 0
        end
        # the B6 generalization term −3Ψ2·μ^{(1)} vanishes when μ^{(1)}=0 (vacuum): trivially,
        # −3·ψ2_bg·(zero μ) = 0
        ψ2bg = np.ψ₂
        @test isempty((psmooth(lscale(-3, ψ2bg), zp)).inside.c)
    end

    # DOCUMENTED correspondence (dossier §4.3 table) — the vacuum Transport equation ↦ D4 assembler:
    #   B1 lambda_transport ↦ recon_set_lam_k!   (hphc_reconstruction.jl:236)
    #   B2 hmm_transport    ↦ recon_set_hmbmb_k! (:255)     B3 psi3 ↦ recon_set_psi3_k! (:227)
    #   B4 pi_transport     ↦ recon_set_pi_k!    (:267)     B5 hlm  ↦ recon_set_hlmb_k! (:278)
    #   B6 psi2_transport   ↦ recon_set_psi2_k!  (:245)     B7 hll  ↦ recon_set_muhll_k! (:293)
    # The μ-BLOCK structure matches (e.g. B1: LY (a,b)=(1,1) ↦ D4 −(μ₀+μ̄₀)); the full POINTWISE
    # cross-chart coefficient identity (γ-block absorbed into the GHP þ′ boost weight; the HPHC R
    # Jacobian) needs the GHP boost-map derivation and is [FUTURE] (same reasoning as D4 G-D4f).
    @testset "G-E-vac-3  μ-block structure matches D4 recon_set_lam_k! (honest spot-check)" begin
        # LY B1 λ: isolated coupling −(μ+μ̄+3γ−γ̄); its μ-block (coeff of μ, of μ̄) is (1,1),
        # matching D4 recon_set_lam_k!'s −R(μ₀+conj μ₀) [= −(μ₀+μ̄₀)].  [FUTURE]: γ-block boost-map.
        op = CHAIN_DELTA_OPS.lambda_transport
        @test (op.a, op.b) == (1, 1)        # μ, μ̄ coefficients == D4 (μ₀+μ̄₀) structure
    end
end
