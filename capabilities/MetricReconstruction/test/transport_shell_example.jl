# Grounding for MetricReconstruction.Transport — the DECISIVE Move-E gate: the Li–Yunes
# (arXiv:2605.11080v1) Schwarzschild static thin-shell worked example (Supplemental Material).
#
# EVERY expected value is a PAPER ANCHOR — a closed form quoted from the tex with its line number —
# NOT whatever the code emits. The chain is developed from the transport equations + BC/jump/
# continuity and CHECKED against the independent closed forms (and, at G-E6, against the linearized
# exact Frauendiener metric 𝔥, which is derived independently of the transport chain).
#
# Gates: G-E1 sources (T→Φ,Λ; R_g) · G-E2 vanishing · G-E3 trace branch · G-E4 Ψ2 + δ-cancellation
#        · G-E5 h_ll + continuity · G-E6 end-to-end metric match after ξ_t gauge.
# {T}-generic: the whole example re-runs at BigFloat-256 (anchors at rel 1e-40).

using Test
using MetricReconstruction.Transport
import Spacetime.NPGHP as NPGHP

# tolerance helper: symbolic-coefficient equality of two Laurents against an anchor Laurent
lmatch(a, b; atol) = Transport.lapprox(a, b; atol = atol)

# Build the anchor closed forms (from the tex) for a given (T, M, R). Returns a NamedTuple of the
# Laurent/scalar anchors. M=area mass; f(R)=1−2M/R.
function shell_anchors(::Type{T}, M, R) where {T}
    Mt, Rt = T(M), T(R)
    fR = 1 - 2Mt / Rt
    πc = T(pi)
    (
     # ANCHOR-1 EOS (LY:718–721)
     σ = Mt / (4πc * Rt^2 * sqrt(fR)),
     p = Mt^2 / (8πc * Rt^3 * fR^(T(3) / 2)),
     # ANCHOR-3 T_ab δ-coefficients (LY:749–756)
     Tll = Mt / (4πc * Rt^2 * fR), Tnn = Mt * fR / (16πc * Rt^2),
     Tln = Mt / (8πc * Rt^2), Tmm̄ = Mt^2 / (8πc * Rt^3 * fR),
     # ANCHOR-4 NP-Ricci δ-coefficients (LY:758–769)
     Φ00 = Mt / (Rt^2 * fR), Φ22 = Mt * fR / (4Rt^2),
     Φ11 = Mt * (Rt - Mt) / (4Rt^3 * fR), Λ = Mt * (Rt - 3Mt) / (12Rt^3 * fR),
     # R_g δ′/δ (LY:810–814)
     Rg_dpc = -Mt * (2Rt - 3Mt) / (12Rt^3), Rg_dc = Mt^2 / (12Rt^4 * fR),
     # ANCHOR-7 jump (LY:795) ; ANCHOR-8 h_mm̄ (LY:799–802) ; ANCHOR-9 μ (LY:816–818)
     hmm_jump = -2Mt / fR,
     hmm_in = laurent(T, 0 => 2Mt / (Rt * fR), -1 => -2Mt / fR),
     μ_in = laurent(T, -2 => -Mt / (2fR), -3 => Mt^2 / fR),
     # ANCHOR-10/psi2 (LY:820–828)
     Ψ2_in = laurent(T, -3 => 3Mt^2 / (Rt * fR), -4 => -3Mt^2 / fR),
     Ψ2_out = laurent(T, -3 => -Mt),
     Ψ2_dc = Mt * (2Rt - 3Mt) / (6Rt^3 * fR),
     # ANCHOR-11 δ-cancellation Θ residual (LY:850–851)
     cancel_in = laurent(T, -3 => Mt * (Rt - 3Mt) / (fR * Rt)),
     cancel_out = laurent(T, -3 => Mt),
     # ANCHOR-12 h_ll numerators (LY:863–869)  h_ll = numer/f²
     hll_numer_out = laurent(T, -1 => 8Mt),
     hll_numer_in = laurent(T, 0 => 8Mt / (Rt * fR), -1 => -24Mt^2 / (Rt * fR), -2 => 8Mt^2 / fR),
     hllR = 8Mt / (Rt * fR^2),
     fR = fR,
    )
end

# Run all closed-form (symbolic-coefficient) gates for one (T, M, R) at symbolic tolerance `atol`.
function run_shell_symbolic(::Type{T}, M, R; atol) where {T}
    a = shell_anchors(T, M, R)
    sol = transport_solve(ShellSource(M = M, R = R), T)

    @testset "G-E1 sources  (M=$M R=$R  $(T))" begin
        # G-E1a: shell_T_np == ANCHOR-3 (EOS → tetrad-contracted NP T, LY:749–756)
        Tn = shell_T_np(T, M, R)
        @test isapprox(Tn.Tll, a.Tll; atol = atol)
        @test isapprox(Tn.Tnn, a.Tnn; atol = atol)
        @test isapprox(Tn.Tln, a.Tln; atol = atol)
        @test isapprox(Tn.Tmm̄, a.Tmm̄; atol = atol)
        # G-E1c: np_ricci_from_T == ANCHOR-4 (on-shell projection, LY:757)
        φ = np_ricci_from_T(Tn.Tll, Tn.Tnn, Tn.Tln, Tn.Tmm̄)
        @test isapprox(φ.Φ00, a.Φ00; atol = atol)
        @test isapprox(φ.Φ22, a.Φ22; atol = atol)
        @test isapprox(φ.Φ11, a.Φ11; atol = atol)
        @test isapprox(φ.Λ, a.Λ; atol = atol)
        # R_g reduced from the sources == LY:810–814 (δ′ and δ)
        Rg = Rg_shell(T, M, R)
        @test isapprox(Rg.dpc, a.Rg_dpc; atol = atol)
        @test isapprox(Rg.dc, a.Rg_dc; atol = atol)
    end

    @testset "G-E2 vanishing (ANCHOR-5)  (M=$M R=$R  $(T))" begin
        van = sol.vanishing
        @test van.all_zero
        for fld in (van.λ, van.h_m̄m̄, van.Ψ3, van.π, van.h_lm̄)
            @test isempty(fld.inside.c) && isempty(fld.outside.c) && fld.dc == 0 && fld.dpc == 0
        end
    end

    @testset "G-E3 trace branch (ANCHOR-6..9)  (M=$M R=$R  $(T))" begin
        # ANCHOR-8 h_mm̄
        @test lmatch(sol.h_mm̄.inside, a.hmm_in; atol = atol)
        @test isempty(sol.h_mm̄.outside.c)
        # ANCHOR-7 jump  [r²∂_r h_mm̄] = −2M/f(R)
        r2 = laurent(T, 2 => one(T))
        r2dh = psmooth(r2, pderiv(sol.h_mm̄))
        @test isapprox(pjump(r2dh), a.hmm_jump; atol = atol)
        # ANCHOR-9 μ
        @test lmatch(sol.μ.inside, a.μ_in; atol = atol)
        @test isempty(sol.μ.outside.c)
        # A1 residual: mu_transport_final  Δ_{[2,0,1,1]}μ == −Φ22  (LY:229) — pure δ, zero interior
        resid = delta_static(Δw(2, 0, 1, 1), sol.μ, M)
        @test lmatch(resid.inside, lzero(T); atol = atol)
        @test lmatch(resid.outside, lzero(T); atol = atol)
        @test isapprox(resid.dc, -sol.Φ22.dc; atol = atol)
    end

    @testset "G-E4 Ψ2 + δ-cancellation (ANCHOR-10/11)  (M=$M R=$R  $(T))" begin
        @test lmatch(sol.Ψ2.inside, a.Ψ2_in; atol = atol)
        @test lmatch(sol.Ψ2.outside, a.Ψ2_out; atol = atol)
        @test isapprox(sol.Ψ2.dc, a.Ψ2_dc; atol = atol)
        # THE δ-cancellation: the δ coefficients of ∂_rμ − Ψ2 − 2Λ cancel EXACTLY (LY:853)
        @test isapprox(sol.cancel.dc, zero(T); atol = atol)
        @test lmatch(sol.cancel.inside, a.cancel_in; atol = atol)
        @test lmatch(sol.cancel.outside, a.cancel_out; atol = atol)
    end

    @testset "G-E5 h_ll + continuity (ANCHOR-12/13)  (M=$M R=$R  $(T))" begin
        @test lmatch(sol.h_ll.numer_out, a.hll_numer_out; atol = atol)
        @test lmatch(sol.h_ll.numer_in, a.hll_numer_in; atol = atol)
        # ANCHOR-13 continuity at r=R  h_ll⁺(R)=h_ll⁻(R)=8M/(Rf²(R))
        fR = a.fR
        hllR_out = leval(sol.h_ll.numer_out, T(R)) / fR^2
        hllR_in = leval(sol.h_ll.numer_in, T(R)) / fR^2
        @test isapprox(hllR_out, a.hllR; atol = atol)
        @test isapprox(hllR_in, a.hllR; atol = atol)
    end

    return sol
end

# G-E6 end-to-end: reconstructed h_μν == exact 𝔥_μν (tt,rr,φφ) directly; the residual h_tr is
# exactly removable by ξ_t (h_tr == gauge_ode_rhs), and ξ_t integrates that RHS.
function run_shell_endtoend(::Type{T}, M, R; rtol, θ = T(1)) where {T}
    sol = transport_solve(ShellSource(M = M, R = R), T)
    Rt = T(R)
    samples = (Rt - Rt / 4, Rt - Rt / 8, Rt + Rt / 8, Rt + Rt / 2)   # both sides, off-horizon
    @testset "G-E6 end-to-end (ANCHOR-14)  (M=$M R=$R  $(T))" begin
        for r in samples
            @test isapprox(recon_htt(sol.h_ll, r), exact_htt(M, R, r); rtol = rtol)
            @test isapprox(recon_hrr(sol.h_ll, r), exact_hrr(M, R, r); rtol = rtol)
            @test isapprox(recon_hφφ_over_sin2(sol.h_mm̄, r) * sin(θ)^2,
                           exact_hφφ_over_sin2(M, R, r) * sin(θ)^2; rtol = rtol, atol = eps(T))
            # h_tr is exactly the removable gauge form (⇒ after ξ_t transform h_tr→0=𝔥_tr, LY:968)
            @test isapprox(recon_htr(sol.h_ll, r), gauge_ode_rhs(M, R, r); rtol = rtol)
        end
        # independent check: ξ_t (closed form, LY:941–967) integrates its ODE:
        # (∂_r − f′/f)ξ_t == gauge_ode_rhs.  High-precision central difference.
        h = T(R) * (T <: BigFloat ? T(10)^(-20) : T(1e-6))
        for r in samples
            fr = 1 - 2T(M) / r
            fpr = 2T(M) / r^2
            dξ = (gauge_xi_t(M, R, r + h) - gauge_xi_t(M, R, r - h)) / (2h)
            lhs = dξ - (fpr / fr) * gauge_xi_t(M, R, r)
            @test isapprox(lhs, gauge_ode_rhs(M, R, r); rtol = (T <: BigFloat ? T(1e-30) : 1e-6))
        end
    end
    return sol
end

@testset "Transport — Schwarzschild static thin shell (LY 2605.11080 Supplement C)" begin
    # G-E1b: the a=0 Schwarzschild NP table is the NPGHP kernel, not a re-hardcode (LY:743–748).
    @testset "G-E1b  schw_np_laurent == Spacetime.NPGHP.np_spincoeff(0,r,θ)" begin
        np = schw_np_laurent(Float64, 1.0)
        for r in (2.5, 4.0, 7.0), θ in (0.7, 1.3, 2.4)
            ref = NPGHP.np_spincoeff(0.0, r, θ)
            @test isapprox(leval(np.ρ, r), real(ref.ρ); atol = 1e-13)
            @test isapprox(leval(np.μ, r), real(ref.μ); atol = 1e-13)
            @test isapprox(leval(np.γ, r), real(ref.γ); atol = 1e-13)
            @test isapprox(leval(np.ψ₂, r), real(ref.ψ₂); atol = 1e-13)
            @test isapprox(imag(ref.μ), 0.0; atol = 1e-13)   # real at a=0
        end
        # weight_combo of Δ_{[2,0,1,1]} on the NPGHP a=0 table == 2μ+2γ  (LY:207)
        ref = NPGHP.np_spincoeff(0.0, 5.0, 1.0)
        @test isapprox(weight_combo(Δw(2, 0, 1, 1), ref), 2ref.μ + 2ref.γ; atol = 1e-14)
    end

    # Float64 closed-form + end-to-end, for (M,R) in {(1,4),(1,10)}
    for (M, R) in ((1.0, 4.0), (1.0, 10.0))
        run_shell_symbolic(Float64, M, R; atol = 1e-11)
        run_shell_endtoend(Float64, M, R; rtol = 1e-13)
    end

    # {T}-generic: the whole example re-runs at BigFloat-256 (anchors at rel 1e-40)
    setprecision(BigFloat, 256) do
        for (M, R) in ((1, 4), (1, 10))
            run_shell_symbolic(BigFloat, big(M), big(R); atol = big(10.0)^(-40))
            run_shell_endtoend(BigFloat, big(M), big(R); rtol = big(10.0)^(-40))
        end
    end
end
