# Grounding for the native GHP thorn/thorn-prime operators (ghp_operators.jl). NO Hertz here — the
# operators are grounded STANDALONE against Newman-Penrose/GHP algebra on equatorial Schwarzschild:
#
#  (1) schw_spincoeff literals — transcribed from Sourceincoordinates.nb Schwspincoeff (SC): ρ=−1/r,
#      ρ′=(r−2M)/2r², ε′=ε′†=−M/2r², ψ₂=−M/r³. Exact rationals at r=6, M=1.
#  (2) þ ≡ 𝒟† at a=0,m=0 — thorn is the Kinnersley l-leg derivative iωr²/Δ+∂_r, which is EXACTLY the
#      Chandrasekhar 𝒟†=∂_r+c (c=iK/Δ, K=ωr² at a=0,m=0) of Operators.jet_ops.opDd (in-GrTool cross-
#      check that thorn reuses the SAME radial jet calculus, including the ω time part).
#  (3) þ′ n-leg directional part — at ω=0, weight-0, þ′g|₀ = n^r g′ = −(Δ/2r²)g′(r0), and the extracted
#      derivative g′(r0) equals Operators.opD(jet)|₀ (opD=∂_r at ω=0).
#  (4) [þ,þ′]g = (p+q)(M/r³)g — the GHP commutator curvature identity on a type-{p,q} scalar in vacuum
#      Schwarzschild (τ=σ=0, curvature term −(p+q)ψ₂=(p+q)M/r³). Self-contained analytic oracle: the g′/g″
#      terms cancel operator-wise, leaving exactly (p+q)(M/r³)g. Checked for several (p,q) and radial g.
#  (5) þ′ρ = ρ̄ρ′ − ψ₂ — a Schwarzschild GHP Ricci identity (ρ real). At r=6: both sides = −1/216.
#
# All expected values are analytic (hand-derived / Operators cross-check); none taken from new code.

using ModifiedGravity, Test
using ModifiedGravity: schw_spincoeff, thorn, thorn_prime
using NumericsCore: jderiv, jmul, jscale, jadd   # Move C: independent jet-op reference (was MR.Lorenz.jet_ops)

# Independent closed-form tower of r²/Δ = r/(r−2M) (derivative VALUES [f,f′,f″,…] at r0) — the
# Chandrasekhar 𝒟† coefficient at a=0,m=0. Analytic (f^{(0)}=r0/(r0−2M); f^{(k)}=2M(−1)^k k!(r0−2M)^{−(k+1)},
# k≥1), so it grounds thorn's own jmul(r²,1/Δ) construction from OUTSIDE the module under test.
function _r2overΔ_tower(r0, N; M=1.0)
    CT = Complex{typeof(float(r0))}
    u = float(r0) - 2M
    out = zeros(CT, N+1)
    out[1] = CT(float(r0) / u)
    for k in 1:N
        out[k+1] = CT(2M * (-1)^k * factorial(k) / u^(k+1))
    end
    out
end

# radial Taylor tower of g(r)=r^s at r0, to order N (exact): g^{(k)} = s(s-1)…(s-k+1) r0^{s-k}.
function powjet(r0, s, N)
    CT = Complex{typeof(float(r0))}
    out = CT[]
    for k in 0:N
        ff = 1.0
        for i in 0:(k-1); ff *= (s - i); end
        push!(out, CT(ff * float(r0)^(s - k)))
    end
    out
end

@testset "ghp_operators" begin
    r0 = 6.0; N = 6

    @testset "schw_spincoeff literals (Schwspincoeff, r=6 M=1)" begin
        sc = schw_spincoeff(6.0)
        @test sc.ρ  ≈ -1/6
        @test sc.ρ′ ≈ (6 - 2)/(2*36)          # = 1/18
        @test sc.ε′ ≈ -1/(2*36)               # = −1/72
        @test sc.ε′dag ≈ -1/(2*36)
        @test sc.ψ₂ ≈ -1/216
        @test sc.κ == 0 && sc.σ == 0 && sc.τ == 0 && sc.ε == 0
    end

    @testset "þ ≡ 𝒟† at a=0,m=0 (incl. ω time part)" begin
        ω = 0.3
        jet = powjet(r0, 2, N) .+ powjet(r0, -1, N)   # arbitrary radial tower
        # 𝒟†g = ∂_r g + iω(r²/Δ)g on equatorial Schwarzschild a=0,m=0 (independent closed-form reference)
        Dd = jadd(jderiv(jet), jscale(im*ω, jmul(_r2overΔ_tower(r0, N), jet)))
        @test thorn(jet, 1, 1, r0, ω) ≈ Dd
    end

    @testset "þ′ n-leg directional part vs ∂_r (ω=0, weight 0)" begin
        jet = powjet(r0, 3, N) .+ powjet(r0, -2, N)
        f = 1 - 2/r0                            # Δ/r² = 1−2M/r
        # þ′g|₀ = −(Δ/2r²) g′(r0); at ω=0 the l-leg reduces to ∂_r so jderiv(jet)|₀ = g′(r0)
        @test thorn_prime(jet, 0, 0, r0, 0.0)[1] ≈ -(f/2) * jderiv(jet)[1]
        @test jderiv(jet)[1] ≈ jet[2]
    end

    @testset "[þ,þ′]g = (p+q)(M/r³)g  (GHP commutator curvature)" begin
        M = 1.0
        for (p, q, s) in [(1, 1, 2), (2, 1, -1), (0, 3, 3), (-1, -1, 4)]
            g = powjet(r0, s, N)
            # [þ,þ′]g = þ(þ′g) − þ′(þg); þ′g:{p,q}→{p−1,q−1}, þg:{p,q}→{p+1,q+1}
            t1 = thorn(thorn_prime(g, p, q, r0, 0.0; M=M), p-1, q-1, r0, 0.0; M=M)
            t2 = thorn_prime(thorn(g, p, q, r0, 0.0; M=M), p+1, q+1, r0, 0.0; M=M)
            comm0 = (t1 .- t2)[1]
            @test comm0 ≈ (p + q) * (M / r0^3) * g[1] atol=1e-12
        end
    end

    @testset "þ′ρ = ρ̄ρ′ − ψ₂  (Schwarzschild GHP Ricci identity)" begin
        sc = schw_spincoeff(r0)
        ρjet = powjet(r0, -1, N) .* (-1)        # ρ = −1/r ⇒ tower = −(r^{-1} tower)
        lhs = thorn_prime(ρjet, 1, 1, r0, 0.0)[1]     # ρ has GHP type {1,1}
        rhs = sc.ρ * sc.ρ′ - sc.ψ₂
        @test lhs ≈ rhs
        @test lhs ≈ -1/216
    end
end
