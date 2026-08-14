# Grounding for MetricReconstruction.Transport — the F6-B DECISIVE gates: the static ℓ=1 m=0 AXIAL
# angular-momentum (Lense–Thirring) completion (F6B_DERIVATION.md; referee F6B_referee.md).
#
# EVERY target is an EXACT-KERR referee anchor:
#   • the metric  d/da kerr_metric(a,r,θ;:BL)|₀  (computed IN-TEST by central difference, h=1e-6, as
#     the referee did) and its analytic O(a) form  −2δJ sin²θ/r  (g_tφ only);
#   • the frame-invariant Weyl scalar  Ψ₂⁽¹⁾ = −3i δJ cosθ/r⁴  (referee ITEM-3).
# The chain ODE residuals are checked against the IMPLEMENTED operators (CHAIN_DELTA_OPS + the a=0
# NPGHP table), never the note's restatement.  {T}-generic: everything re-runs at BigFloat-256.
#
# Gates: G-F6B-1 metric (DECISIVE) · G-F6B-2 Weyl (INVARIANT) · G-F6B-3 zero sectors ·
#        G-F6B-4 {T} BigFloat-256 · plus the verify-the-algebra chain residuals (B3/B4/B5/B6).

using Test
using MetricReconstruction.Transport
import Spacetime
import Spacetime.NPGHP as NPGHP

# central-difference d/da kerr_metric|₀ (the referee's own FD-oracle), scaled by δJ
function fd_metric_pert(::Type{T}, δJ, r, θ, h) where {T}
    gp = Spacetime.kerr_metric(T(h), T(r), T(θ); coords = :BL)
    gm = Spacetime.kerr_metric(-T(h), T(r), T(θ); coords = :BL)
    T(δJ) .* (gp .- gm) ./ (2 * T(h))
end

# analytic O(δJ) BL perturbation: only h_tφ = h_φt = −2δJ sin²θ/r
function analytic_metric_pert(::Type{T}, δJ, r, θ) where {T}
    M = zeros(Complex{T}, 4, 4)
    M[1, 4] = M[4, 1] = -2 * T(δJ) * sin(T(θ))^2 / T(r)
    M
end

# ================================================================================================
# G-F6B-1 (metric, DECISIVE) : assembled h_μν (ORG assembly + axial gauge) == d/da kerr_metric|₀
# ================================================================================================
function run_metric_gate(::Type{T}, δJ; atol_fd, atol_ana, use_fd) where {T}
    dc = DipoleCompletion(T(δJ))
    grid = [(T(r), T(θ)) for r in (3, 5, 10, 30) for θ in (T(6)/10, T(11)/10, T(pi)/2, T(23)/10)]
    @testset "G-F6B-1 metric (M=1 δJ=$δJ  $(T))" begin
        for (r, θ) in grid
            H = dipole_h(dc, r, θ)
            @test eltype(H) == Complex{T}
            # vs the analytic O(δJ) BL form (g_tφ only) — the tight, {T}-clean target
            tgt = analytic_metric_pert(T, δJ, r, θ)
            @test maximum(abs.(H .- tgt)) ≤ atol_ana
            # vs the central-difference kerr_metric FD-oracle (Float64 only; FD-truncation-limited)
            if use_fd
                ref = fd_metric_pert(T, δJ, r, θ, T(1e-6))
                @test maximum(abs.(H .- ref)) ≤ atol_fd
            end
            # the ONLY nonzero component is (t,φ); imaginary part is zero everywhere
            @test abs(imag(H[1, 4])) ≤ atol_ana
        end
        # structural: the ORG assembly carries h_tφ AND h_rφ = h_tφ/f; the gauge removes only h_rφ.
        r, θ = T(5), T(11)/10
        Horg = dipole_assemble_org(dc, r, θ)
        f = 1 - 2/r
        @test isapprox(real(Horg[1, 4]), -2*T(δJ)*sin(θ)^2/r; atol = atol_ana)      # h_tφ^ORG
        @test isapprox(real(Horg[2, 4]), -2*T(δJ)*sin(θ)^2/(r*f); atol = atol_ana)   # h_rφ^ORG = h_tφ/f
        @test isapprox(dipole_gauge_hrphi(dc, r, θ), +2*T(δJ)*sin(θ)^2/(r*f); atol = atol_ana)
        @test isapprox(dipole_h(dc, r, θ)[2, 4], zero(Complex{T}); atol = atol_ana)   # h_rφ removed
        # the gauge generator c(r) integrates its ODE: c′ − 2c/r == 2δJ/(rf)  (central diff)
        hh = T <: BigFloat ? T(10)^(-20) : T(1e-6)
        for rr in (T(4), T(7))
            dcp = (dipole_gauge_c(dc, rr + hh) - dipole_gauge_c(dc, rr - hh)) / (2hh)
            lhs = dcp - 2 * dipole_gauge_c(dc, rr) / rr
            @test isapprox(lhs, 2*T(δJ)/(rr*(1 - 2/rr)); rtol = (T <: BigFloat ? T(10)^(-25) : 1e-6))
        end
    end
    return dc
end

# ================================================================================================
# G-F6B-2 (Weyl, INVARIANT) : chain Ψ₂⁽¹⁾ == δJ·(−3i cosθ/r⁴).
#
# INVARIANCE ARGUMENT (restated + cited): the a=0 Schwarzschild background is type-D with
# Ψ₀=Ψ₁=Ψ₃=Ψ₄=0.  Under any O(a) tetrad rotation (class-I about l: Ψ₂→Ψ₂+2bΨ₁+b²Ψ₀; class-II about
# n: Ψ₂→Ψ₂+2aΨ₃+a²Ψ₄; spin/boost: Ψ₂ has weight 0) the first-order δΨ₂ is INVARIANT because every
# background scalar it would mix with vanishes.  Hence the chain (ORG-frame) Ψ₂⁽¹⁾ equals the
# Kinnersley-frame referee value δΨ₂ = −3i δJ cosθ/r⁴ (referee ITEM-3; = the ℓ=1 part of −1/ζ³,
# ζ=r−iaδJcosθ).  This is a NON-tautological cross-check: C₃ was fixed by the metric normalization.
# ================================================================================================
function run_weyl_gate(::Type{T}, δJ; atol) where {T}
    dc = DipoleCompletion(T(δJ))
    grid = [(T(r), T(θ)) for r in (3, 5, 10, 30) for θ in (T(6)/10, T(11)/10, T(23)/10)]
    @testset "G-F6B-2 Weyl invariant (δJ=$δJ  $(T))" begin
        for (r, θ) in grid
            Ψ2 = dipole_psi2(dc, r, θ)
            @test eltype(Ψ2) == Complex{T}
            @test isapprox(Ψ2, Complex{T}(0, -3*T(δJ)*cos(θ)/r^4); atol = atol)
        end
        # radial Laurent shape is exactly r⁻⁴ with coefficient 2√2 C₃ = −3i δJ
        L = dipole_psi2_radial(dc)
        @test isapprox(leval(L, T(1)), Complex{T}(0, -3*T(δJ)); atol = atol)   # coeff at r=1
    end
    return dc
end

# ================================================================================================
# G-F6B-3 (zero sectors) : spin-2 fields, trace branch, and h_ll identically zero.
# ================================================================================================
function run_zero_gate(::Type{T}, δJ; atol) where {T}
    dc = DipoleCompletion(T(δJ))
    @testset "G-F6B-3 zero sectors (δJ=$δJ  $(T))" begin
        for (r, θ) in [(T(4), T(1)), (T(9), T(2))]
            @test dipole_lambda(dc, r, θ) == 0            # spin-2 λ (ℓ=1<|s|=2)
            @test dipole_hmbarmbar(dc, r, θ) == 0          # spin-2 h_m̄m̄
            @test dipole_hmmbar(dc, r, θ) == 0             # trace/polar spin-0 h_mm̄
            @test dipole_hll(dc, r, θ) == 0                # polar spin-0 h_ll (AXIOM-F6B-3)
        end
        # DIRECT consistency: the assembled metric really has NO spin-0 (polar) content — project
        # the reconstructed BL h_μν onto l⊗l, n⊗n, l⊗n, m⊗m̄.  All must vanish (referee ITEM-2).
        for (r, θ) in [(T(5), T(11)/10), (T(12), T(7)/10)]
            H = dipole_h(dc, r, θ)
            lu, nu, mu, mbu = Spacetime.kinnersley_tetrad(zero(T), r, θ; coords = :BL)
            proj(ea, eb) = sum(H[μ, ν] * ea[μ] * eb[ν] for μ in 1:4, ν in 1:4)
            @test abs(proj(lu, lu)) ≤ atol     # h_ll = 0
            @test abs(proj(nu, nu)) ≤ atol     # h_nn = 0
            @test abs(proj(lu, nu)) ≤ atol     # h_ln = 0
            @test abs(proj(mu, mbu)) ≤ atol    # h_mm̄ = 0
        end
    end
    return dc
end

# ================================================================================================
# VERIFY-THE-ALGEBRA chain residuals — B3/B4/B5/B6 solved against the IMPLEMENTED operators
# (CHAIN_DELTA_OPS quartets + the a=0 NPGHP table via weight_combo), NOT the note's restatement.
# Static reduction:  Δ_{[a,b,c,d]} X = −(f/2)∂_r X + weight_combo(op, np_a0)·X.
# ================================================================================================
function run_chain_residuals(::Type{T}, δJ; atol) where {T}
    dc = DipoleCompletion(T(δJ))
    grid = [(T(r), T(θ)) for r in (3, 5, 10) for θ in (T(6)/10, T(11)/10, T(23)/10)]
    # a=0 NPGHP background table (real μ,γ,β,τ) — the implemented ground-truth coefficients
    np0(r, θ) = NPGHP.np_spincoeff(zero(T), r, θ)
    fdd(g, r) = begin hh = T(1e-6); (g(r + hh) - g(r - hh)) / (2hh) end   # central diff for non-Laurent
    @testset "verify-the-algebra: chain ODE residuals (δJ=$δJ  $(T))" begin
        for (r, θ) in grid
            f = 1 - 2/r
            np = np0(r, θ)
            # B3: Δ_{[4,0,2,0]}Ψ₃ = 0   (homogeneous).  Ψ₃, Ψ₃′ from the Laurent (exact).
            L3 = dipole_psi3_radial(dc)
            alg3 = weight_combo(CHAIN_DELTA_OPS.psi3_transport, np)          # = 4μ+2γ = −2/r+5/r²
            @test isapprox(real(alg3), -2/r + 5/r^2; atol = atol)
            res3 = -(f/2)*leval(lderiv(L3), r) + alg3*leval(L3, r)
            @test isapprox(res3, zero(Complex{T}); atol = atol)
            # B4: Δ_{[0,0,1,-1]}π = −Ψ₃   (alg 0).
            L4 = dipole_pi_radial(dc)
            alg4 = weight_combo(CHAIN_DELTA_OPS.pi_transport, np)            # = 0
            @test isapprox(real(alg4), 0; atol = atol)
            res4 = (-(f/2)*leval(lderiv(L4), r) + alg4*leval(L4, r)) - (-leval(L3, r))
            @test isapprox(res4, zero(Complex{T}); atol = atol)
            # B5: Δ_{[0,1,0,-2]}h_lm̄ = +2π   (SOURCE +2π — the referee-consistent sign, header).
            alg5 = weight_combo(CHAIN_DELTA_OPS.hlm_transport, np)           # = μ−2γ = −1/(2r)
            @test isapprox(real(alg5), -1/(2r); atol = atol)
            h_of(rr) = dipole_hlmbar_radial(dc, rr)
            res5 = (-(f/2)*fdd(h_of, r) + alg5*h_of(r)) - (2*leval(L4, r))
            @test isapprox(res5, zero(Complex{T}); atol = (T <: BigFloat ? T(1e-12) : 1e-6))  # FD-limited
            # B6: Δ_{[3,0,0,0]}Ψ₂ = (δ+2β−2τ)Ψ₃.   δ-source built from the Spacetime tetrad + NPGHP β,τ.
            L2 = dipole_psi2_radial(dc)
            alg6 = weight_combo(CHAIN_DELTA_OPS.psi2_transport, np)          # = 3μ = −3/(2r)+3/r²
            @test isapprox(real(alg6), -3/(2r) + 3/r^2; atol = atol)
            lu, nu, mu, mbu = Spacetime.kinnersley_tetrad(zero(T), r, θ; coords = :BL)
            # Ψ₃(θ) = leval(L3,r)·sinθ ; ∂_θ Ψ₃ = leval(L3,r)·cosθ ; δ = m^θ∂_θ (m=0 ⇒ ∂_φ term = 0)
            δΨ3 = mu[3] * (leval(L3, r) * cos(θ))
            src6 = δΨ3 + 2*np.β*(leval(L3, r)*sin(θ)) - 2*np.τ*(leval(L3, r)*sin(θ))
            resid6 = (-(f/2)*leval(lderiv(L2), r) + alg6*leval(L2, r)) * cos(θ) - src6
            @test isapprox(resid6, zero(Complex{T}); atol = atol)
        end
    end
    return dc
end

# ================================================================================================
@testset "Transport — F6-B dipole (Lense–Thirring) completion (F6B_DERIVATION)" begin
    # sanity: the derived amplitude C₃ = −3√2 i δJ/4
    let dc = DipoleCompletion(1.0)
        @test isapprox(dc.C₃, -3*sqrt(2)*im/4; atol = 1e-15)
        @test dc.δJ == 1.0 && dc.Kc == 0.0
    end

    # Float64: metric (FD + analytic), Weyl, zero sectors, chain residuals — a couple of δJ values.
    for δJ in (1.0, 0.5)
        run_metric_gate(Float64, δJ; atol_fd = 1e-10, atol_ana = 1e-12, use_fd = true)
        run_weyl_gate(Float64, δJ; atol = 1e-12)
        run_zero_gate(Float64, δJ; atol = 1e-12)
        run_chain_residuals(Float64, δJ; atol = 1e-11)
    end

    # G-F6B-4 {T}-generic: the whole set re-runs at BigFloat-256 (analytic anchors at 1e-40).
    setprecision(BigFloat, 256) do
        for δJ in (big"1.0", big"0.5")
            run_metric_gate(BigFloat, δJ; atol_fd = big(10.0)^(-9), atol_ana = big(10.0)^(-40), use_fd = false)
            run_weyl_gate(BigFloat, δJ; atol = big(10.0)^(-40))
            run_zero_gate(BigFloat, δJ; atol = big(10.0)^(-40))
            run_chain_residuals(BigFloat, δJ; atol = big(10.0)^(-38))
        end
    end
end
