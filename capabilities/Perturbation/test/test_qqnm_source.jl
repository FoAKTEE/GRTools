# Q3 Deliverable B gates — the QQNM source assembly Q(r) (src/QNM/qqnm_source.jl).
#
#   G-Q3-1  DECISIVE. The assembled ++ Schwarzschild Q(r) (l_L=2 → l_Q=4) vs the FROZEN Paper-I:1427
#           transcription. Constant-ratio (single c, no per-point freedom) THEN shape identity:
#           (48√7π/5)·W_ab == P_i to tol.  NEVER tuned: P_i is transcribed below BEFORE assembly.
#   G-Q3-2  W-structure: Q = r¹²ℛ²W₁ + r¹⁰Δℛℛ′W₂ + r⁸Δ²ℛ′²W₃, W_i rational/regular (Paper I:1088).
#   G-Q3-3  channel algebra: +− 16-basis assembly nonzero/structured; −− ≡ 0 (structural zero).
#   G-Q3-4  projection double-route (a=0): analytic 3-j == GL quadrature.
#   G-Q3-5  {T}: ComplexF64 (≤1e-10) AND Complex{BigFloat}-256 (≤1e-40).
#   MUST-PIN: the general-parent-pair bilinear symmetrization reduces to the printed Paper-I:591
#            identical-parent form (owner:1503/1560; Paper I:591).
using Test
using Printf
using Perturbation
using Harmonics: SpinWeightedSpheroidalEigenvalue
import Perturbation: _qqnm_proj_a0_analytic, _qqnm_proj_a0_gl, _qqnm_Cab_pp, _qqnm_parent_chain,
                     _qqnm_coefjets, _qqnm_alphabeta_vals, _wigner3j

# ─────────── FROZEN transcription of Paper-I:1427 (app:source_Schwarzschild_l2) ────────────────────
# (48√7π/5) Q(r) = P1 ℛ² + P2 ℛℛ′ + P3 ℛ′², P1=r²·B1, P2=r²Δ·B2, P3=Δ²·B3, Δ=r²−2r (a=0). ℛ=₊₂R_{22ω}.
_P1(r, ω) = r^2*(-13*r^10*ω^6 - 34im*r^9*ω^5 + r^8*(35+269im*ω)*ω^4 + r^7*(131im*ω+99)*im*ω^3 +
    r^6*(865*ω^2-949im*ω+44)*ω^2 + r^5*(1249*ω^2+530im*ω-35)*im*ω +
    r^4*(-53*ω^2-196im*ω+557)*im*ω + r^3*(2144*ω^2-1540im*ω+120) +
    r^2*(-1380*ω^2+930im*ω-564) + r*(816+408im*ω) - 168im*ω - 336)
function _P2(r, ω)
    Δ = r^2 - 2r
    r^2*Δ*(32im*r^8*ω^5 - 78*r^7*ω^4 - r^6*(424im*ω+130)*im*ω^3 + r^5*(154+582im*ω)*ω^2 +
        r^4*(-968*ω^2+1227im*ω+42)*im*ω - r^3*(2320im*ω+300)*im*ω +
        r^2*(-936*ω^2+1291im*ω+48) - r*(192+2352im*ω) + 1152im*ω + 192)
end
function _P3(r, ω)
    Δ = r^2 - 2r
    Δ^2*(19*r^8*ω^4 + 50im*r^7*ω^3 - r^6*(86+143im*ω)*ω^2 - r^5*(339im*ω+91)*im*ω +
        r^4*(-339*ω^2+538im*ω+42) - r^3*(306+945im*ω) + r^2*(810+447im*ω) - 912*r + 360)
end

# G-Q3-1 over a generic (r-grid, ω): returns the per-coefficient ratios W_ab/P_i (want a single const).
function _q31_ratios(::Type{C}, ω, rs) where {C}
    p = (l=2, m=2, n=0, ω=C(ω))
    a = C <: Complex && real(C) <: BigFloat ? big(0) : 0
    co = qqnm_source_coeffs(p, p, :pp, 4; a=a, r=collect(rs), T=C)
    ρ = Complex{real(C)}[]
    for (k, rr) in enumerate(rs)
        push!(ρ, co.W00[k] / _P1(C(rr), C(ω)))
        push!(ρ, co.W01[k] / _P2(C(rr), C(ω)))
        push!(ρ, co.W11[k] / _P3(C(rr), C(ω)))
    end
    ρ
end

function run_q31(::Type{C}, rs, tol) where {C}
    RT = real(C)
    label = C === ComplexF64 ? "F64" : "BF$(precision(BigFloat))"
    cpaper = 5 / (48 * sqrt(RT(7) * RT(pi)))                       # the paper's printed prefactor⁻¹
    @testset "G-Q3-1 DECISIVE Paper-I:1427 [$label]" begin
        for ω in (Complex{RT}(RT(37367168)/100000000 - im*RT(8896232)/100000000),
                  Complex{RT}(1//2), Complex{RT}(3//10 - im*2//10))
            ρ = _q31_ratios(C, ω, rs)
            c = sum(ρ) / length(ρ)                                 # the single fitted constant
            # (a) constant-ratio: NO per-point freedom — every W_ab/P_i equals one constant c
            @test maximum(abs.(ρ .- c)) / abs(c) < tol
            # (b) shape identity: c is EXACTLY the paper's printed (48√7π/5)⁻¹ ⇒ (48√7π/5)·W==P
            @test abs(c - cpaper) / abs(cpaper) < tol
        end
    end
end

@testset "Q3 Deliverable B — QQNM source assembly Q(r)" begin
    rsF = (3.0, 4.5, 7.0, 12.0, 20.0)

    # G-Q3-1 + G-Q3-5 ({T}: F64 then BF-256)
    run_q31(ComplexF64, rsF, 1e-10)
    setprecision(BigFloat, 256) do
        rsB = (big"3.0", big"4.5", big"7.0", big"12.0", big"20.0")
        run_q31(Complex{BigFloat}, rsB, 1e-40)
    end

    # Visibility: the DECISIVE ratio table (W_ab/P_i per r & coefficient), the fitted constant, the
    # BF-256 achieved residual, and the projection double-route agreement.
    let ω = ComplexF64(0.37367168 - 0.08896232im)
        p = (l=2, m=2, n=0, ω=ω)
        co = qqnm_source_coeffs(p, p, :pp, 4; a=0.0, r=collect(rsF), T=ComplexF64)
        cp = 5/(48*sqrt(7π))
        println("─ G-Q3-1 ratio table  W_ab(r)/P_i(r)  [want const 5/(48√7π)=", cp, "] ─")
        println("     r          R²(W00/P1)                RR'(W01/P2)               R'²(W11/P3)")
        for (k, r) in enumerate(rsF)
            @printf("  %6.2f   %+.12e%+.12eim  %+.12e%+.12eim  %+.12e%+.12eim\n", r,
                    real(co.W00[k]/_P1(r,ω)), imag(co.W00[k]/_P1(r,ω)),
                    real(co.W01[k]/_P2(r,ω)), imag(co.W01[k]/_P2(r,ω)),
                    real(co.W11[k]/_P3(r,ω)), imag(co.W11[k]/_P3(r,ω)))
        end
        ρ = _q31_ratios(ComplexF64, ω, rsF); c = sum(ρ)/length(ρ)
        println("G-Q3-1 [F64] fitted c=", c, " | constant-ratio max|ρ-c|/|c|=",
                maximum(abs.(ρ .- c))/abs(c), " | shape |c-cp|/|cp|=", abs(c-cp)/abs(cp))
    end
    setprecision(BigFloat, 256) do
        ω = Complex{BigFloat}(big"0.37367168" - im*big"0.08896232")
        ρ = _q31_ratios(Complex{BigFloat}, ω, (big"3.0", big"7.0", big"12.0"))
        c = sum(ρ)/length(ρ); cp = 5/(48*sqrt(BigFloat(7)*BigFloat(pi)))
        println("G-Q3-1 [BF256] constant-ratio max|ρ-c|/|c|=", Float64(maximum(abs.(ρ .- c))/abs(c)),
                " | shape |c-cp|/|cp|=", Float64(abs(c-cp)/abs(cp)))
    end
    let
        d = maximum(map(((lL,lQ,mQ),) -> begin
                (Aa,Ba)=_qqnm_proj_a0_analytic(Float64,lL,lL,lQ,mQ)
                (Ag,Bg)=_qqnm_proj_a0_gl(Float64,lL,lL,lQ,mQ,64)
                max(abs(Aa-Ag), abs(Ba-Bg)) end, ((2,4,4),(2,5,4),(2,6,4),(3,6,6))))
        println("G-Q3-4 projection double-route (a=0): max|analytic − GL| = ", d)
    end

    @testset "G-Q3-2 W-structure r¹²ℛ²W₁+r¹⁰Δℛℛ′W₂+r⁸Δ²ℛ′²W₃ (Paper I:1088)" begin
        p = (l=2, m=2, n=0, ω=ComplexF64(0.37367168 - 0.08896232im))
        rs = collect((3.5, 6.0, 11.0))
        co = qqnm_source_coeffs(p, p, :pp, 4; a=0.0, r=rs, T=ComplexF64)
        for (k, r) in enumerate(rs)
            Δ = r^2 - 2r
            W1 = co.W00[k] / r^12          # regular pieces (finite along the contour, r>r₊)
            W2 = co.W01[k] / (r^10 * Δ)
            W3 = co.W11[k] / (r^8 * Δ^2)
            @test isfinite(W1) && isfinite(W2) && isfinite(W3)
            # and they reproduce the paper decomposition: W00 = (5/48√7π)·r²·B1 = r¹²·[(5/48√7π)B1/r¹⁰]
            @test abs(co.W00[k] - (5/(48*sqrt(7π)))*_P1(r, p.ω)) / abs(co.W00[k]) < 1e-9
        end
    end

    @testset "G-Q3-4 projection double-route a=0 (analytic 3-j == GL)" begin
        for (lL, lQ, mQ) in ((2,4,4), (2,5,4), (2,6,4), (3,6,6))
            (PAa, PBa) = _qqnm_proj_a0_analytic(Float64, lL, lL, lQ, mQ)
            (PAg, PBg) = _qqnm_proj_a0_gl(Float64, lL, lL, lQ, mQ, 64)
            @test abs(PAa - PAg) < 1e-11
            @test abs(PBa - PBg) < 1e-11
        end
        # the l_Q=2→l_Q=4 anchors match the paper's printed 3-j values
        (PA, PB) = _qqnm_proj_a0_analytic(Float64, 2, 2, 4, 4)
        @test abs(PA - 5/(2*sqrt(42π))) < 1e-13
        @test abs(PB - 5/(3*sqrt(7π))) < 1e-13
        # selection rule l_Q>2l_L (5,4→ from 2,2) ⇒ vanishing projection at a=0
        (PA6, PB6) = _qqnm_proj_a0_analytic(Float64, 2, 2, 5, 4)
        @test abs(PA6) < 1e-13 && abs(PB6) < 1e-13
    end

    @testset "MUST-PIN: general-pair bilinear ↦ identical-parent Paper-I:591 form" begin
        # Two GENERIC (α,β) chains; the general-pair _qqnm_Cab_pp with chain1==chain2 must equal the
        # printed identical-parent contraction Σ_{i≤j} a_ij (α_iα_j , α_iβ_j+β_iα_j , β_iβ_j).
        C = ComplexF64; Γ = C(3.7+0.4im); Γb = C(3.7-0.4im); r = C(3.7)
        α = C[0.9, -0.3+0.2im, 0.5, -0.7im, 0.11, 0.4-0.1im, -0.22]
        β = C[0.0, 1.0, -0.4+0.1im, 0.6, -0.2im, 0.33, 0.05]
        for ν in 1:4
            (C00, C01, C10, C11) = _qqnm_Cab_pp(ν, α, β, α, β, Γ, Γb, r, C)
            aij = qqnm_pp_aij(ν, C, Γ, Γb, r)
            c00=zero(C); c01=zero(C); c11=zero(C)
            for ((i,j), a) in aij
                c00 += a*α[i+1]*α[j+1]; c01 += a*(α[i+1]*β[j+1]+β[i+1]*α[j+1]); c11 += a*β[i+1]*β[j+1]
            end
            @test abs(C00 - c00) < 1e-12
            @test abs((C01 + C10) - c01) < 1e-12        # identical parents merge C01,C10
            @test abs(C11 - c11) < 1e-12
        end
    end

    @testset "G-Q3-3 channel algebra: +− structured, −− ≡ 0" begin
        p1 = (l=2, m=2, n=0, ω=ComplexF64(0.37367168 - 0.08896232im))
        p2 = (l=2, m=2, n=0, ω=ComplexF64(0.37367168 - 0.08896232im))
        # −− the exact structural zero (owner:1585)
        cm = qqnm_source_coeffs(p1, p2, :mm, 4; a=0.0, r=[4.0, 6.0], T=ComplexF64)
        @test all(iszero, cm.W00) && all(iszero, cm.W01) && all(iszero, cm.W11)
        @test cm.zero === true
        # +− 16-basis assembly is nonvanishing + carries the mirror-channel i (owner:1263)
        cp = qqnm_source_coeffs(p1, p2, :pm, 4; a=0.0, r=[4.0, 6.0], T=ComplexF64)
        @test any(!iszero, cp.W00) || any(!iszero, cp.W11)
        # [PRELIMINARY] the full 16↔4-basis angular cross-check (G-Q3-3) is REVIEWER-scoped: it needs
        # the four-generator angular reduction (owner derive-angular-four-basis.m) over complex-γ
        # spheroidals. The per-atom radial C_ab (mirror carry 𝔠_-) is exercised here.
    end

    @testset "surface qqnm_source: coeffs vs Q(r) with a supplied Hertz radial" begin
        p = (l=2, m=2, n=0, ω=ComplexF64(0.37367168 - 0.08896232im))
        rs = collect((4.0, 6.0, 9.0))
        co = qqnm_source(p, p, :pp, 4; a=0.0, r=rs, T=ComplexF64)         # coefficients (R omitted)
        Rvec = ComplexF64[1.0/r for r in rs]; dRvec = ComplexF64[-1.0/r^2 for r in rs]
        Q = qqnm_source(p, p, :pp, 4; a=0.0, r=rs, R=Rvec, dR=dRvec, T=ComplexF64)
        @test length(Q) == length(rs)
        @test Q ≈ co.W00 .* Rvec.^2 .+ co.W01 .* Rvec .* dRvec .+ co.W11 .* dRvec.^2
    end
end
