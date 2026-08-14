# test_npghp.jl — Spacetime.NPGHP gate suite G-B1..G-B5 (NPGHP campaign Move B, PLAN.md binding).
#
# G-B1  HPHC NP_0 fixture match (rel ≤1e-13 on the baked TeukEvolution grid) + BigFloat
#        true-closed-form check for tau_0/pi_0 (sqrt(big(2)) — DEVIATION-Q4).
# G-B2  DECISIVE: BL table vs the first-principles oracle (test/npghp_oracle.jl) — sweep
#        a∈{0,.35,.7,.9,.99}×r∈{2.2,5,20}×θ∈{π/6,π/2,3π/4} isapprox(rtol=1e-8, atol=1e-9); the
#        atol is the oracle's absolute FD floor (max noise 2.99e-10, R_oracle.md §Step2). Plus the
#        primed table vs the LEG-SWAP oracle; a BigFloat-256 spot (rtol 1e-30); and ψ₂ = −1/ζ³.
# G-B3  a=0 equator reduces to the baked ModifiedGravity.schw_spincoeff literals (r∈{3,5,10}).
# G-B4  chart identity: psi2_0(R,cy) == ψ₂_BL(cl²/R,θ)·(r/cl²)³ on the fixture grid (cl=1 and √2).
# G-B5  ghp_weights vs the R_teuk.md §2 12-row (s,b,falloff,p,q) table + prime/bar weight maps.
# G-F1  HPHC↔BL GHP boost-map identities (F1_DERIVATION.md).
# G-F5  :In/:Out chart NP tables (F5_DERIVATION.md): chart-invariant BL table via the coords keyword;
#        DECISIVE first-principles oracle recomputed in :In/:Out with FD Christoffels of the analytic
#        chart metric (analytic-Γ shortcut is BL-only) + boost-law route + ψ₂/normalization + BigFloat.

using Test
using Spacetime
const NP = Spacetime.NPGHP

include("fixtures/npghp_np0_fixture.jl")
include("fixtures/npghp_mg_schw_literals.jl")
include("npghp_oracle.jl")
const F0 = NPGHPNp0Fixture
const MGS = NPGHPMGSchwLiterals
const OR = NPGHPOracle

_reldiff(x, y) = abs(x - y) / max(abs(x), abs(y), eps(Float64))
_clval(kind, ::Type{T}) where {T} = kind === :one ? one(T) : sqrt(T(2))
_θ(num, den, ::Type{T}) where {T} = (num * (T <: BigFloat ? big(π) : convert(T, π))) / den

# ---- G-F1 helpers (F1_DERIVATION.md; harness scratchpad/npghp_phase2/f1_leg_comparison.jl) --------
# The six NP_0 identity (LHS,RHS) pairs at (a,R,θ,cl): the module NP_0 field vs the boost/spin-rotation
# transformed BL field. A = hphc_boost = Δ/(2Σ); e^{iλ} = hphc_spin_phase = −ζ̄/ζ (ph below).
function _f1_six(a, R, θ, cl)
    cy = cos(θ); sy = sin(θ); r = cl^2 / R
    A  = NP.hphc_boost(a, R, cy; cl = cl)
    ph = NP.hphc_spin_phase(a, R, cy; cl = cl)
    H  = NP.np_spincoeff_hphc(R, cy, sy; a = a, cl = cl)
    B  = NP.np_spincoeff(a, r, θ)
    ((:rho,  H.rho_0,  (A / R) * B.ρ),          # ρ₀ = (A/R) ρ_BL          {1,1}   b=+1
     (:mu,   H.mu_0,   B.μ / (A * R)),          # μ₀ = μ_BL/(A·R)          {−1,−1} b=−1
     (:tau,  H.tau_0,  (ph / R^2) * B.τ),       # τ₀ = −(ζ̄/ζ) τ_BL/R²      {1,−1}  s=+1
     (:pi,   H.pi_0,   (conj(ph) / R^2) * B.π), # π₀ = −(ζ/ζ̄) π_BL/R²      {−1,1}  s=−1 (conj law)
     (:psi2, H.psi2_0, B.ψ₂ / R^3),             # ψ₂₀ = ψ₂_BL/R³           {0,0}   invariant
     (:eps,  H.eps_0,  B.γ / R^2))              # ε₀ = γ_BL/R² (boost connection lands on γ_BL)
end
# Chart Jacobian J = ∂(T,R,ϑ',φ')/∂(t,r,ϑ,φ_BL) [tex:3689/3695/3701/3924-3930] (M=1, Float64).
function _f1_chart_jacobian(a, r, cl)
    Δ = Spacetime.delta_kerr(a, r)
    J = zeros(Complex{Float64}, 4, 4)
    J[1, 1] = 1
    J[1, 2] = (r^2 + a^2) / Δ - 2 - 4 / r        # dT/dr = dr_*/dr − 2 − 4M/r
    J[2, 2] = -cl^2 / r^2                          # dR/dr
    J[3, 3] = 1
    J[4, 2] = a / Δ                                # dφ'/dr azimuthal shift
    J[4, 4] = 1
    J
end
# Reference HPHC legs, verbatim from eq:tetrad_IEF_HC [tex:3941-3968] (M=1, Float64).
function _f1_hphc_legs(a, R, θ, cl)
    cy = cos(θ); sy = sin(θ); L2 = cl^2; L4 = cl^4; rt2 = sqrt(2.0); CT = Complex{Float64}
    den_l = L4 + a^2 * R^2 * cy^2
    lH = (R^2 / den_l) .* CT[2 * (2 - (a^2 / L2) * R), -(L2 - 2R + (a^2 / L2) * R^2) / 2, 0, a]
    nH = CT[2 + 4R / L2, R^2 / L2, 0, 0]
    zm = L2 - im * a * R * cy
    mH = (R / (rt2 * zm)) .* CT[-im * a * sy, 0, -1, -im / sy]
    (lH, nH, mH)
end

# ---- G-F5 helpers (:In/:Out chart-map Jacobian; F5_DERIVATION.md) ---------------------------------
# BL→(ingoing/outgoing)-Kerr chart Jacobian J = ∂(x_chart)/∂(x_BL): dv = dt + sgn·(r²+a²)/Δ dr,
# dφ̃ = dφ + sgn·a/Δ dr (sgn=+1 :In, −1 :Out); r,θ unchanged. (M=1, Float64.)
function _f5_chart_jacobian(a, r; sgn::Int)
    Δ = Spacetime.delta_kerr(a, r)
    J = zeros(Complex{Float64}, 4, 4)
    J[1, 1] = 1; J[1, 2] = sgn * (r^2 + a^2) / Δ
    J[2, 2] = 1; J[3, 3] = 1
    J[4, 2] = sgn * a / Δ; J[4, 4] = 1
    J
end
# component-independent proportionality constant of leg P vs Q (over the nonzero components of Q)
function _f5_propconst(P, Q)
    k = argmax(abs.(Q)); c = P[k] / Q[k]
    (c, maximum(abs(P[j] - c * Q[j]) for j in 1:4))
end

@testset "NPGHP (Move B: G-B1..G-B5)" begin

    # ============================================================ G-B1 ==========================
    @testset "G-B1 HPHC NP_0 fixture match + BigFloat true-closed-form (DEVIATION-Q4)" begin
        np0keys = (:mu_0, :tau_0, :pi_0, :rho_0, :eps_0, :psi2_0)
        # -- Float64 grid: the module's actual output vs the baked upstream oracle -----------------
        worstF = 0.0; wF = ""
        for (ci, cfg) in enumerate(F0.CONFIGS), (ti, (num, den)) in enumerate(F0.THETAS),
                (ri, Rstr) in enumerate(F0.RLABELS)
            bhs = parse(Float64, cfg[1]); cl = _clval(cfg[2], Float64)
            R = parse(Float64, Rstr); θ = _θ(num, den, Float64); cy = cos(θ); sy = sin(θ)
            got = NP.np_spincoeff_hphc(R, cy, sy; a = bhs, cl = cl)
            ref = F0.F64[(ci, ri, ti)]
            for k in np0keys
                d = _reldiff(getfield(got, k), getfield(ref, k))
                d > worstF && (worstF = d; wF = "cfg$ci R=$Rstr θ=$num π/$den $k")
                @test isapprox(getfield(got, k), getfield(ref, k); rtol = 1e-13, atol = 1e-14)
            end
        end
        println("  [G-B1] Float64 vs baked fixture: worst rel = ", worstF, "  @ ", wF, "  (≤1e-13)")

        # -- BigFloat-256 grid: module (sqrt(big(2))) vs baked upstream (Float64 √2) ---------------
        setprecision(BigFloat, 256) do
            worstB = 0.0; wB = ""
            for (ci, cfg) in enumerate(F0.CONFIGS), (ti, (num, den)) in enumerate(F0.THETAS),
                    (ri, Rstr) in enumerate(F0.RLABELS)
                bhs = parse(BigFloat, cfg[1]); cl = _clval(cfg[2], BigFloat)
                R = parse(BigFloat, Rstr); θ = _θ(num, den, BigFloat); cy = cos(θ); sy = sin(θ)
                got = NP.np_spincoeff_hphc(R, cy, sy; a = bhs, cl = cl)
                refstr = F0.BF_STR[(ci, ri, ti)]
                for k in np0keys
                    re, ims = getfield(refstr, k)
                    refv = Complex{BigFloat}(parse(BigFloat, re), parse(BigFloat, ims))
                    d = Float64(_reldiff(getfield(got, k), refv))
                    d > worstB && (worstB = d; wB = "cfg$ci R=$Rstr θ=$num π/$den $k")
                    @test isapprox(getfield(got, k), refv; rtol = 1e-13, atol = big"1e-14")
                end
            end
            println("  [G-B1] BigFloat vs baked fixture: worst rel = ", worstB, "  @ ", wB,
                    "  (≤1e-13; tau_0/pi_0 carry the ~1e-16 Float64-√2 fixture delta)")

            # -- DEVIATION-Q4: at BigFloat the module tau_0/pi_0 ARE the true closed form ----------
            bhs = big"0.9"; cl = one(BigFloat); R = big"0.35"; θ = big(π) / 6
            cy = cos(θ); sy = sin(θ); bhm = one(BigFloat)
            got = NP.np_spincoeff_hphc(R, cy, sy; a = bhs, cl = cl)
            zm = cl^2 - im * bhs * R * cy; d2 = cl^4 + (bhs * R * cy)^2
            tau_true = (im * bhs * sy / sqrt(big(2))) / (zm^2)          # sqrt(big(2)) = TRUE
            pi_true = -(im * bhs * sy / sqrt(big(2))) / d2
            @test isapprox(got.tau_0, tau_true; rtol = big"1e-70")
            @test isapprox(got.pi_0, pi_true; rtol = big"1e-70")
            # and confirm the fixture (Float64 √2) differs from the true form at ~1e-16 (the artifact)
            fr = F0.BF_STR[(3, 2, 1)].tau_0
            tau_fix = Complex{BigFloat}(parse(BigFloat, fr[1]), parse(BigFloat, fr[2]))
            delta = Float64(_reldiff(tau_true, tau_fix))
            @test 1e-17 < delta < 1e-15
            println("  [G-B1] DEVIATION-Q4: module tau_0/pi_0 == true sqrt(big(2)) form to ≤1e-70; ",
                    "baked-fixture (Float64 √2) differs at rel = ", delta, "  (the honest artifact)")
        end
    end

    # ============================================================ G-B2 ==========================
    @testset "G-B2 BL table vs first-principles oracle (DECISIVE)" begin
        avals = [0.0, 0.35, 0.7, 0.9, 0.99]; rvals = [2.2, 5.0, 20.0]
        θvals = [π / 6, π / 2, 3π / 4]
        RTOL = 1e-8; ATOL = 1e-9
        zerokeys = (:κ, :σ, :ν, :λ, :ε)   # np structural type-D zeros

        # -- unprimed sweep -----------------------------------------------------------------------
        worst = 0.0; wpt = ""; npts = 0
        for a in avals, r in rvals, θ in θvals
            npts += 1
            my = NP.np_spincoeff(a, r, θ); orc = OR.np_spincoeffs(a, r, θ; h = 1e-6)
            ok = true
            for k in keys(orc)
                mv = getfield(my, k); ov = getfield(orc, k)
                isapprox(mv, ov; rtol = RTOL, atol = ATOL) || (ok = false)
                v = abs(ov) < 1e-6 ? abs(mv - ov) : _reldiff(mv, ov)   # metric for the printed worst
                v > worst && (worst = v; wpt = "(a=$a,r=$r,θ=$(round(θ,digits=3))) $k")
            end
            @test ok
            # ψ₂ closed form (exact) + oracle Riemann cross-check
            @test isapprox(my.ψ₂, -inv(Spacetime.zeta_kerr(a, r, θ)^3); rtol = 1e-13)
            @test isapprox(my.ψ₂, OR.psi2(a, r, θ; h = 1e-5); rtol = 1e-6, atol = 1e-8)
            # structural type-D exact zeros
            @test all(getfield(my, k) == 0 for k in zerokeys)
        end
        println("  [G-B2] unprimed sweep ($npts pts): all isapprox(rtol=$RTOL,atol=$ATOL); ",
                "worst rel/abs = ", worst)

        # -- primed sweep vs LEG-SWAP oracle (the [AXIOM-4] full-Kerr closure) ---------------------
        worstp = 0.0; wpp = ""
        for a in avals, r in rvals, θ in θvals
            my = NP.ghp_spincoeff(a, r, θ); orc = OR.np_spincoeffs(a, r, θ; h = 1e-6, primed = true)
            ok = true
            for k in keys(orc)
                mv = getfield(my, k); ov = getfield(orc, k)
                isapprox(mv, ov; rtol = RTOL, atol = ATOL) || (ok = false)
                v = abs(ov) < 1e-6 ? abs(mv - ov) : _reldiff(mv, ov)
                v > worstp && (worstp = v; wpp = "(a=$a,r=$r,θ=$(round(θ,digits=3))) $k")
            end
            @test ok
        end
        println("  [G-B2] primed vs leg-swap oracle: all isapprox; worst rel/abs = ", worstp,
                "  (ρ′=−μ, ε′=−γ, … validated at full Kerr, not the a=0 pin)")

        # -- BigFloat-256 spot (0.7,5.0,1.1) rtol 1e-30 -------------------------------------------
        setprecision(BigFloat, 256) do
            a = big"0.7"; r = big"5.0"; θ = big"1.1"
            my = NP.np_spincoeff(a, r, θ); orc = OR.np_spincoeffs(a, r, θ; h = big"1e-30")
            myp = NP.ghp_spincoeff(a, r, θ); orp = OR.np_spincoeffs(a, r, θ; h = big"1e-30", primed = true)
            w = 0.0
            for k in keys(orc)
                @test isapprox(getfield(my, k), getfield(orc, k); rtol = big"1e-30", atol = big"1e-30")
                @test isapprox(getfield(myp, k), getfield(orp, k); rtol = big"1e-30", atol = big"1e-30")
                w = max(w, Float64(_reldiff(getfield(my, k), getfield(orc, k))))
            end
            @test isapprox(my.ψ₂, -inv(Spacetime.zeta_kerr(a, r, θ)^3); rtol = big"1e-70")
            @test my.ρ isa Complex{BigFloat}
            println("  [G-B2] BigFloat-256 spot (0.7,5,1.1): np+primed isapprox(rtol=1e-30); ",
                    "nonzero worst rel ≈ ", w, " ; ψ₂ vs −1/ζ³ exact; eltype = Complex{BigFloat}")
        end
    end

    # ============================================================ G-B3 ==========================
    @testset "G-B3 a=0 equator == baked MG schw_spincoeff literals (r∈{3,5,10})" begin
        for (r, ref) in MGS.SCHW
            np = NP.np_spincoeff(0.0, r, π / 2)
            gp = NP.ghp_spincoeff(0.0, r, π / 2)   # primed: .ρ=ρ′, .ε=ε′
            @test isapprox(real(np.ρ), ref.ρ; rtol = 1e-13) && isapprox(imag(np.ρ), 0.0; atol = 1e-13)
            @test isapprox(real(-np.μ), ref.ρp; rtol = 1e-13)   # ρ′ = −μ  ([AXIOM-4])
            @test isapprox(real(-np.γ), ref.εp; rtol = 1e-13)   # ε′ = −γ
            @test isapprox(real(np.ψ₂), ref.ψ2; rtol = 1e-13)
            @test np.κ == 0 && np.σ == 0 && np.τ == 0 && np.ε == 0 && np.ν == 0 && np.λ == 0
            # primed view reads ρ′, ε′ directly
            @test isapprox(real(gp.ρ), ref.ρp; rtol = 1e-13)
            @test isapprox(real(gp.ε), ref.εp; rtol = 1e-13)
        end
        println("  [G-B3] a=0/equator np_spincoeff & primed view reduce to MG schw_spincoeff ",
                "(ρ,ρ′=−μ,ε′=−γ,ψ₂) at r∈{3,5,10}, rel ≤1e-13")
    end

    # ============================================================ G-B4 ==========================
    @testset "G-B4 chart identity psi2_0 == ψ₂_BL(cl²/R,θ)·(r/cl²)³" begin
        worst = 0.0; w = ""
        for cfg in F0.CONFIGS, (num, den) in F0.THETAS, Rstr in F0.RLABELS
            bhs = parse(Float64, cfg[1]); cl = _clval(cfg[2], Float64)
            R = parse(Float64, Rstr); θ = _θ(num, den, Float64); cy = cos(θ); sy = sin(θ)
            lhs = NP.np_spincoeff_hphc(R, cy, sy; a = bhs, cl = cl).psi2_0
            r = cl^2 / R
            rhs = NP.np_spincoeff(bhs, r, θ).ψ₂ * (r / cl^2)^3
            d = _reldiff(lhs, rhs); d > worst && (worst = d; w = "cl=$cl a=$bhs R=$Rstr θ=$num π/$den")
            @test isapprox(lhs, rhs; rtol = 1e-12)
        end
        println("  [G-B4] psi2_0 chart identity over the fixture grid (cl=1 and √2): worst rel = ",
                worst, "  @ ", w, "  (≤1e-12; ψ₂ is GHP {0,0} → boost-clean)")
    end

    # ============================================================ G-B5 ==========================
    @testset "G-B5 ghp_weights vs R_teuk §2 12-row (s,b,falloff,p,q) table + prime/bar maps" begin
        # (field, s, b, falloff, p, q) — R_teuk.md §2 (TeukEvolution GHP.jl:205,232,262-263)
        TABLE = [
            (:lin_f, -2, -2, 1, -4, 0), (:psi3, -1, -1, 2, -2, 0), (:psi2, 0, 0, 3, 0, 0),
            (:lam, -2, -1, 1, -3, 1), (:pi, -1, 0, 2, -1, 1), (:muhll, 0, 1, 3, 1, 1),
            (:hlmb, -1, 1, 2, 0, 2), (:hmbmb, -2, 0, 1, -2, 2), (:res_bianchi3, -2, -1, 2, -3, 1),
            (:res_bianchi2, -1, 0, 2, -1, 1), (:res_hll, 0, 2, 2, 2, 2), (:scd_f, -2, -2, 1, -4, 0),
        ]
        for (name, s, b, fall, p, q) in TABLE
            @test NP.ghp_weights(s, b) == (p, q)
            @test NP.ghp_prime_weights(p, q) == (-p, -q)     # prime: {p,q}→{−p,−q}
            @test NP.ghp_bar_weights(p, q) == (q, p)         # bar:   {p,q}→{q,p}
            # prime of weights == weights of primed (s,b)→(−s,−b)
            @test NP.ghp_weights(-s, -b) == NP.ghp_prime_weights(NP.ghp_weights(s, b)...)
        end
        println("  [G-B5] ghp_weights(s,b)=(s+b,−s+b) matches all 12 rows; prime/bar maps verified")
    end

    # ============================================================ G-F1 ==========================
    # F1: HPHC ↔ BL-Kinnersley GHP boost map (F1_DERIVATION.md). The Move-B/D cross-chart FIXTURE
    # gates are UPGRADED to cross-chart IDENTITIES: (a) the six NP_0 fields == the boost/spin-rotation
    # transformed BL fields; (b) the leg law l_H=A·(J l_K), n_H=A⁻¹·(J n_K), m_H=e^{iλ}·(J m_K);
    # (c) the D4-deferred MG þ′-ingredient cross-check, now unlocked. A=Δ/(2Σ), e^{iλ}=−ζ̄/ζ.
    @testset "G-F1 HPHC↔BL boost-map identities (F1_DERIVATION.md)" begin

        # ---- (a) the SIX NP_0 identity gates ----------------------------------------------------
        AS = (0.0, 0.35, 0.7, 0.95); RS = (0.08, 0.3, 0.55); TS = (π / 5, π / 2, 2π / 3)
        worstF = 0.0; wF = ""
        for a in AS, R in RS, θ in TS
            for (nm, lhs, rhs) in _f1_six(a, R, θ, 1.0)
                d = _reldiff(lhs, rhs); d > worstF && (worstF = d; wF = "a=$a R=$R θ=$(round(θ,digits=3)) $nm")
                @test isapprox(lhs, rhs; rtol = 1e-13, atol = 1e-14)
            end
        end
        # one cl = √2 triple (θ-sweep at a=0.7, R=0.3)
        for θ in TS
            for (nm, lhs, rhs) in _f1_six(0.7, 0.3, θ, sqrt(2.0))
                @test isapprox(lhs, rhs; rtol = 1e-13, atol = 1e-14)
            end
        end
        println("  [G-F1a] six NP_0 identities ρ₀=(A/R)ρ_BL, μ₀=μ_BL/(AR), τ₀=−(ζ̄/ζ)τ_BL/R², ",
                "π₀=conj-law, ψ₂₀=ψ₂_BL/R³, ε₀=γ_BL/R² (Float64, 108 samples + cl=√2 triple): ",
                "worst rel = ", worstF, "  @ ", wF, "  (≤1e-13)")

        # -- BigFloat-256 triple at rel ≤1e-40 + eltype checks ------------------------------------
        setprecision(BigFloat, 256) do
            a = big"0.95"; R = big"0.3"; cl = one(BigFloat)
            worstB = 0.0
            for θ in (big(π) / 5, big(π) / 2, 2 * big(π) / 3)
                for (nm, lhs, rhs) in _f1_six(a, R, θ, cl)
                    @test lhs isa Complex{BigFloat} && rhs isa Complex{BigFloat}
                    @test isapprox(lhs, rhs; rtol = big"1e-40", atol = big"1e-50")
                    worstB = max(worstB, Float64(_reldiff(lhs, rhs)))
                end
            end
            @test NP.hphc_boost(a, R, cos(big(π) / 5); cl = cl) isa BigFloat
            @test NP.hphc_spin_phase(a, R, cos(big(π) / 5); cl = cl) isa Complex{BigFloat}
            println("  [G-F1a] BigFloat-256 triple (a=0.95,R=0.3): worst rel = ", worstB,
                    " (≤1e-40); hphc_boost::BigFloat, hphc_spin_phase::Complex{BigFloat}")
        end

        # ---- (b) leg-level law: l_H==A·(J l_K), n_H==A⁻¹·(J n_K), m_H==e^{iλ}·(J m_K) -----------
        maxleg = 0.0; mln = 0.0; mm = 0.0
        for (a, R, θ) in ((0.9, 0.35, π / 2), (0.5, 0.1, 3π / 4), (0.7, 0.55, π / 5))
            cl = 1.0; r = cl^2 / R; cy = cos(θ)
            lH, nH, mH = _f1_hphc_legs(a, R, θ, cl)
            lK, nK, mK, _ = Spacetime.kinnersley_tetrad(a, r, θ; coords = :BL)
            J = _f1_chart_jacobian(a, r, cl)
            A = NP.hphc_boost(a, R, cy; cl = cl); ph = NP.hphc_spin_phase(a, R, cy; cl = cl)
            Jl = J * lK; Jn = J * nK; Jm = J * mK
            for μ in 1:4
                @test isapprox(lH[μ], A * Jl[μ]; rtol = 1e-12, atol = 1e-12)          # l_H = A·(J l_K)
                @test isapprox(nH[μ], (1 / A) * Jn[μ]; rtol = 1e-12, atol = 1e-12)    # n_H = A⁻¹·(J n_K)
                @test isapprox(mH[μ], ph * Jm[μ]; rtol = 1e-12, atol = 1e-12)          # m_H = e^{iλ}·(J m_K)
            end
            maxleg = max(maxleg, maximum(abs.(lH .- A .* Jl)),
                                 maximum(abs.(nH .- (1 / A) .* Jn)),
                                 maximum(abs.(mH .- ph .* Jm)))
            # const_l·const_n == 1 and |const_m| == 1 (measured from the largest-|P| component)
            kl = argmax(abs.(Jl)); kn = argmax(abs.(Jn)); km = argmax(abs.(Jm))
            cl_ = lH[kl] / Jl[kl]; cn_ = nH[kn] / Jn[kn]; cm_ = mH[km] / Jm[km]
            @test isapprox(cl_ * cn_, 1; rtol = 1e-12, atol = 1e-12)
            @test isapprox(abs(cm_), 1; rtol = 1e-12)
            mln = max(mln, abs(cl_ * cn_ - 1)); mm = max(mm, abs(abs(cm_) - 1))
        end
        println("  [G-F1b] leg law l_H=A·Jl_K, n_H=A⁻¹·Jn_K, m_H=e^{iλ}·Jm_K (3 samples): ",
                "max|leg diff| = ", maxleg, "  max|const_l·const_n−1| = ", mln,
                "  max||const_m|−1| = ", mm, "  (≤1e-12)")

        # ---- (c) MG þ′-ingredient cross-check — D4 G-D4f [FUTURE] PARTIALLY CLOSED --------------
        # HOME DECISION: ModifiedGravity does NOT depend on Spacetime (MG deps = Printf/NumericsCore/
        # Harmonics/Perturbation/Data — verified in ModifiedGravity/Project.toml), so per the task
        # rule this cross-chart check lives HERE in Spacetime's own suite, using the BAKED MG literals
        # (fixtures/npghp_mg_schw_literals.jl) rather than adding a live MG↔Spacetime dep (which would
        # invert the Tier ordering + trip the §7.7 hygiene gate). MEANING: at a=0 equator the HPHC
        # connection scalar ε₀, boost-mapped by ×R², reproduces MG's a=0 þ′ ingredient ε′_MG up to the
        # GHP prime sign — the COEFFICIENT-level cross-chart map. The FIXED derivation (F1_DERIVATION.md
        # line 33) is ε₀ = γ_BL/R² (POSITIVE), and the GHP prime gives ε′_MG = −γ_BL (G-B3), so the true
        # chain is ε₀·R² == γ_BL == −ε′_MG. (The task's inline shorthand "ε₀·R²==−γ_BL==ε′_MG" carries a
        # sign slip vs the FIXED law; asserted here to the derivation, cross-checked numerically.) The
        # full OPERATOR-level þ′ comparison remains [FUTURE] until F-later.
        for (r, ref) in MGS.SCHW
            R = 1.0 / r                                                   # cl=1 ⇒ R = 1/r
            eps0 = NP.np_spincoeff_hphc(R, 0.0, 1.0; a = 0.0, cl = 1.0).eps_0   # cy=cos(π/2)=0, sy=1
            γ_BL = NP.np_spincoeff(0.0, r, π / 2).γ
            @test isapprox(real(eps0) * R^2, real(γ_BL); rtol = 1e-13)   # ε₀·R² == γ_BL  (F1 law)
            @test isapprox(real(γ_BL), -ref.εp; rtol = 1e-13)            # γ_BL == −ε′_MG (GHP prime, G-B3)
            @test isapprox(real(eps0) * R^2, -ref.εp; rtol = 1e-13)      # ⟹ ε₀·R² == −ε′_MG
            @test isapprox(imag(eps0), 0.0; atol = 1e-14)                # a=0 ⇒ ε₀ real
        end
        println("  [G-F1c] MG þ′-ingredient bridge (home=Spacetime suite; MG⊅Spacetime): ",
                "ε₀·R² == γ_BL == −ε′_MG at r∈{3,5,10}, rel ≤1e-13 ",
                "(D4 G-D4f coefficient-level map CLOSED; operator-level [FUTURE])")
    end

    # ============================================================ G-F5 ==========================
    # F5: :In (ingoing-Kerr) / :Out (outgoing-Kerr) chart NP/GHP tables (F5_DERIVATION.md). The
    # geometry.jl :In/:Out Kinnersley tetrads are PURE coordinate push-forwards of the BL Kinnersley
    # legs (const_l=const_n=const_m=1 EXACTLY; no boost A=1, no phase e^{iλ}=1), so the NP spin
    # coefficients are chart-INVARIANT scalars and `np_spincoeff(:In)==np_spincoeff(:Out)==(:BL)`.
    # (a) DECISIVE oracle route: the first-principles oracle recomputed IN the :In/:Out charts from
    #     their own metric + tetrad + CENTRAL-FD Christoffels of the analytic chart metric (the
    #     BL analytic-Γ shortcut is unavailable — Spacetime.christoffel is BL-only) reproduces the
    #     chart-invariant table; h=hΓ ∈ {1e-6,1e-7} both pass, floor documented honestly.
    # (b) boost-law route: const_l/const_n/const_m measured from the geometry.jl legs pushed through
    #     the chart Jacobian == 1; the two independent routes AGREE.
    # (c) cross-chart consistency: table field-equal across BL/In/Out; ψ₂ chart-invariant (In/Out
    #     FD Riemann); tetrad normalization l·n=−1, m·m̄=+1 per chart; primed table invariant.
    # (d) BigFloat-256 spot. (e) type-D zero set κ=σ=ν=λ=ε preserved.
    @testset "G-F5 :In/:Out chart NP tables (F5_DERIVATION.md)" begin
        avals = [0.0, 0.35, 0.7, 0.9, 0.99]; rvals = [2.2, 5.0, 20.0]
        θvals = [π / 6, π / 2, 3π / 4]
        allkeys = (:κ, :ρ, :σ, :τ, :π, :ν, :μ, :λ, :ε, :γ, :β, :α)
        zerokeys = (:κ, :σ, :ν, :λ, :ε)
        RTOL = 1e-6; ATOL = 1e-7                       # measured In/Out FD-Γ floor (@1e-6: rel 3.6e-8,
                                                        # abs 1.3e-9), gate with ~30×/70× margin

        # ---- (a) DECISIVE: oracle in :In/:Out (FD Christoffels) == chart-invariant BL table --------
        for chart in (:In, :Out)
            sgn = chart === :In ? 1 : -1
            for h in (1e-6, 1e-7)
                worstR = 0.0; wR = ""; worstA = 0.0
                for a in avals, r in rvals, θ in θvals
                    bl = NP.np_spincoeff(a, r, θ)                 # chart-invariant reference
                    oc = OR.np_spincoeffs(a, r, θ; h = h, hΓ = h, coords = chart)
                    for k in allkeys
                        mv = getfield(bl, k); ov = getfield(oc, k)
                        @test isapprox(mv, ov; rtol = RTOL, atol = ATOL)
                        worstA = max(worstA, abs(mv - ov))
                        abs(mv) > 1e-3 && (v = abs(mv - ov) / max(abs(mv), abs(ov), eps());
                                           v > worstR && (worstR = v; wR = "(a=$a,r=$r,θ=$(round(θ,digits=3))) $k"))
                    end
                    # np_spincoeff(chart) IS the BL table (chart-invariant) — exact field equality
                    my = NP.np_spincoeff(a, r, θ; coords = chart)
                    @test all(getfield(my, k) === getfield(bl, k) for k in allkeys)
                    @test my.ψ₂ === bl.ψ₂
                end
                println("  [G-F5a] oracle($chart, FD-Γ) vs np_spincoeff, h=hΓ=$h (135 pts × 12 keys): ",
                        "worst REL(|ref|>1e-3) = ", worstR, " @ ", wR, " ; worst ABS = ", worstA,
                        "  (gate rtol=$RTOL, atol=$ATOL)")
            end
        end

        # ---- (b) boost-law route: const_l=const_n=const_m=1 (no boost/phase); routes AGREE --------
        for chart in (:In, :Out)
            sgn = chart === :In ? 1 : -1
            wcl = 0.0; wcn = 0.0; wcm = 0.0; wdev = 0.0; wln = 0.0; wmm = 0.0
            for a in avals, r in rvals, θ in θvals
                lB, nB, mB, _ = Spacetime.kinnersley_tetrad(a, r, θ; coords = :BL)
                lC, nC, mC, _ = Spacetime.kinnersley_tetrad(a, r, θ; coords = chart)
                J = _f5_chart_jacobian(a, r; sgn = sgn)
                (cl, dl) = _f5_propconst(lC, J * lB)
                (cn, dn) = _f5_propconst(nC, J * nB)
                (cm, dm) = _f5_propconst(mC, J * mB)
                @test isapprox(cl, 1; atol = 1e-12) && isapprox(cn, 1; atol = 1e-12) &&
                      isapprox(cm, 1; atol = 1e-12)
                wcl = max(wcl, abs(cl - 1)); wcn = max(wcn, abs(cn - 1)); wcm = max(wcm, abs(cm - 1))
                wdev = max(wdev, dl, dn, dm)
                wln = max(wln, abs(cl * cn - 1)); wmm = max(wmm, abs(abs(cm) - 1))
                @test isapprox(cl * cn, 1; atol = 1e-12) && isapprox(abs(cm), 1; atol = 1e-12)
            end
            println("  [G-F5b] $chart boost law: max|const_l−1|=", wcl, " max|const_n−1|=", wcn,
                    " max|const_m−1|=", wcm, " (no boost A=1, no phase e^{iλ}=1); leg maxdev=", wdev,
                    " ; |const_l·const_n−1|≤", wln, " ||const_m|−1|≤", wmm)
        end

        # ---- (c) cross-chart consistency ---------------------------------------------------------
        # table field-equal across the three charts (exact); primed table too.
        for a in avals, r in rvals, θ in θvals
            b = NP.np_spincoeff(a, r, θ); i = NP.np_spincoeff(a, r, θ; coords = :In)
            o = NP.np_spincoeff(a, r, θ; coords = :Out)
            @test all(getfield(i, k) === getfield(b, k) && getfield(o, k) === getfield(b, k)
                      for k in (allkeys..., :ψ₂))
            gb = NP.ghp_spincoeff(a, r, θ); gi = NP.ghp_spincoeff(a, r, θ; coords = :In)
            @test all(getfield(gi, k) === getfield(gb, k) for k in (allkeys..., :ψ₂))
        end
        # ψ₂ chart-invariance via the In/Out FD-of-FD Riemann (floor ~6e-7 @ h=1e-3) + normalization
        worstψ = 0.0; worstN = 0.0
        for chart in (:In, :Out)
            for a in avals, r in rvals, θ in θvals
                ref = -inv(Spacetime.zeta_kerr(a, r, θ)^3)
                got = OR.psi2(a, r, θ; h = 1e-3, hΓ = 1e-3, coords = chart)
                @test isapprox(got, ref; rtol = 1e-5, atol = 1e-6)
                worstψ = max(worstψ, abs(got - ref) / max(abs(ref), eps()))
                g = Spacetime.kerr_metric(a, r, θ; coords = chart)
                l, n, m, mb = Spacetime.kinnersley_tetrad(a, r, θ; coords = chart)
                ln = sum(g[p, q] * l[p] * n[q] for p in 1:4, q in 1:4)
                mmb = sum(g[p, q] * m[p] * mb[q] for p in 1:4, q in 1:4)
                @test isapprox(ln, -1; atol = 1e-12) && isapprox(mmb, 1; atol = 1e-12)
                worstN = max(worstN, abs(ln + 1), abs(mmb - 1))
            end
        end
        println("  [G-F5c] cross-chart: np/ghp tables field-equal BL≡In≡Out (exact ===); ",
                "ψ₂ chart-inv via In/Out FD-Riemann worst rel = ", worstψ, " (≤1e-5); ",
                "tetrad norm l·n=−1,m·m̄=1 worst dev = ", worstN)

        # ---- (d) BigFloat-256 spot: oracle(:In) reproduces the BL table to ~1e-40 -----------------
        setprecision(BigFloat, 256) do
            a = big"0.7"; r = big"5.0"; θ = big"1.1"
            bl = NP.np_spincoeff(a, r, θ)
            oc = OR.np_spincoeffs(a, r, θ; h = big"1e-20", hΓ = big"1e-20", coords = :In)
            w = 0.0
            for k in (:ρ, :τ, :π, :μ, :γ, :β, :α)
                @test isapprox(getfield(bl, k), getfield(oc, k); rtol = big"1e-30", atol = big"1e-30")
                w = max(w, Float64(abs(getfield(bl, k) - getfield(oc, k)) / max(abs(getfield(bl, k)), eps())))
            end
            @test NP.np_spincoeff(a, r, θ; coords = :In).ρ isa Complex{BigFloat}
            @test oc.ρ isa Complex{BigFloat}
            println("  [G-F5d] BigFloat-256 In spot (0.7,5,1.1) h=hΓ=1e-20: worst nonzero rel = ", w,
                    " (≤1e-30); eltype = Complex{BigFloat}")
        end

        # ---- (e) type-D zero set κ=σ=ν=λ=ε preserved (EXACT zeros; oracle ≤ FD floor) ------------
        for chart in (:In, :Out)
            for a in avals, r in rvals, θ in θvals
                my = NP.np_spincoeff(a, r, θ; coords = chart)
                @test all(getfield(my, k) == 0 for k in zerokeys)     # exact structural zeros
                oc = OR.np_spincoeffs(a, r, θ; h = 1e-6, hΓ = 1e-6, coords = chart)
                @test all(abs(getfield(oc, k)) < ATOL for k in zerokeys)   # oracle zeros ≤ FD floor
            end
        end
        @test_throws ArgumentError NP.np_spincoeff(0.5, 5.0, 1.0; coords = :Bogus)
        println("  [G-F5e] type-D zeros κ=σ=ν=λ=ε EXACT in :In/:Out; oracle zeros ≤ atol=$ATOL; ",
                "unknown coords throws ArgumentError")
    end

    # ---- grounding (not a numbered gate): HPHC GHP pre_* consistency + {T}-genericity -----------
    @testset "hphc GHP pre_* consistency + genericity" begin
        θ = π / 6; cy = cos(θ); sy = sin(θ)
        pc = NP.hphc_ghp_coeffs(0.7, cy, sy; a = 0.9, cl = 1)
        np0 = NP.np_spincoeff_hphc(0.7, cy, sy; a = 0.9, cl = 1)
        @test isapprox(pc.pre_edth, -0.7 * np0.tau_0; rtol = 1e-13)         # R_teuk [DERIVED]
        @test isapprox(pc.pre_thorn_prime_DR, complex(0.7); rtol = 1e-13)   # R/cl²
        @test isapprox(pc.pre_thorn_prime_DT, complex(4.8); rtol = 1e-13)   # 2+4·bhm·R/cl²
        # chart maps round-trip + horizon
        @test NP.hphc_r_of_R(NP.hphc_R_of_r(7.3; cl = 1.4); cl = 1.4) ≈ 7.3
        @test NP.hphc_horizon_R(0.0; cl = 1) ≈ 0.5                          # cl²/r₊ = 1/2 at a=0
        # {T}-genericity: BigFloat eltype flows through both HPHC entry points
        setprecision(BigFloat, 256) do
            g = NP.np_spincoeff_hphc(big"0.7", cos(big(π) / 6), sin(big(π) / 6); a = big"0.9", cl = one(BigFloat))
            p = NP.hphc_ghp_coeffs(big"0.7", cos(big(π) / 6), sin(big(π) / 6); a = big"0.9", cl = one(BigFloat))
            @test g.rho_0 isa Complex{BigFloat} && p.pre_edth isa Complex{BigFloat}
        end
        println("  [grounding] pre_edth=−R·tau_0, pre_thorn_prime_{DR,DT}, chart round-trip, ",
                "BigFloat eltype through both HPHC entry points — all pass")
    end
end
