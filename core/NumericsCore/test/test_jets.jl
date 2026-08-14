# Characterization of the canonical 1D derivative-tower ("jet") algebra hoisted into
# NumericsCore.jets (Move A: collapse the duplicated Leibniz-binomial jet calculus).
#
# GROUNDING is two-pronged:
#   (G-A1) ANALYTIC — Leibniz product on jets of r², r³, exp(r); reciprocal identity jinv∘jmul==unit;
#          jderiv shift. Expected values are independent closed forms, never the code's own output.
#   (G-A2) CHARACTERIZATION — new == legacy BIT-EXACT (== on arrays; ±0 equal, no NaN/Inf here) vs the
#          two legacy implementations copied VERBATIM below as permanent regression oracles.
#
# LEGACY ORACLES (copied verbatim, cite source):
#   _legacy_lorenz_*   — MetricReconstruction/src/Lorenz/Operators.jl:96-115 (generic-length Vector).
#   _legacy_radiation_* — MetricReconstruction/src/Radiation/Operators.jl:20,55-92 (fixed NJ=5).
#
# CONVENTION VERDICT (verified against source + empirical probe): BOTH legacy homes use the SAME
# derivative-value tower [f,f′,…,f^(N)] with the Leibniz binomial product rule (NOT Taylor coeffs /
# plain Cauchy convolution). jmul/jadd/jscale/jderiv agree BIT-EXACT for Float64/ComplexF64/BigFloat/
# Complex{BigFloat}. The ONLY divergence is jinv's LEADING term for COMPLEX eltype: Lorenz uses
# `one(eltype)/f[1]`, Radiation uses `inv(u[1])`, and Julia's `inv(::Complex)` ≠ `one/z` in the last
# ulp (overflow-robust algorithm) — a floating-point library artifact, not a convention difference.
# Canonical follows LORENZ (the established, externally-consumed set), so it is bit-exact to the Lorenz
# oracle (real+complex) and agrees with the Radiation oracle to machine ε on complex (documented below).

using Test
using NumericsCore
using NumericsCore: constjet, jadd, jscale, jmul, jderiv, jinv

# ── legacy oracle #1: Lorenz/Operators.jl:96-115 (VERBATIM) ──────────────────────────────────────
_legacy_lorenz_binom(n, k) = factorial(n) ÷ (factorial(k)*factorial(n-k))
_legacy_lorenz_jmul(f, g) = (n=length(f)-1; [sum(_legacy_lorenz_binom(k,j)*f[j+1]*g[k-j+1] for j in 0:k) for k in 0:n])
_legacy_lorenz_jadd(f, g) = f .+ g
_legacy_lorenz_jscale(s, f) = s .* f
_legacy_lorenz_jderiv(f) = vcat(f[2:end], [zero(eltype(f))])
function _legacy_lorenz_jinv(f)
    n = length(f)-1
    g = zeros(eltype(f), n+1)
    g[1] = one(eltype(f))/f[1]
    for k in 1:n
        g[k+1] = -sum(_legacy_lorenz_binom(k,j)*f[j+1]*g[k-j+1] for j in 1:k)/f[1]
    end
    g
end
_legacy_lorenz_constjet(s, N::Integer, ::Type{CT}) where {CT} = vcat([CT(s)], zeros(CT, N))

# ── legacy oracle #2: Radiation/Operators.jl:20,55-92 (VERBATIM) ─────────────────────────────────
const _LEG_NJ = 5   # Radiation/Operators.jl:20  `const NJ = 5`
function _legacy_radiation_jet_mul(u::Vector{C}, v::Vector{C}) where {C}
    w = zeros(C, _LEG_NJ)
    @inbounds for k in 0:_LEG_NJ-1
        acc = zero(C)
        for j in 0:k
            acc += C(binomial(k, j)) * u[j+1] * v[k-j+1]
        end
        w[k+1] = acc
    end
    w
end
_legacy_radiation_jet_add(u::Vector{C}, v::Vector{C}) where {C} = u .+ v
function _legacy_radiation_jet_deriv(u::Vector{C}) where {C}
    d = zeros(C, _LEG_NJ)
    @inbounds for k in 1:_LEG_NJ-1
        d[k] = u[k+1]
    end
    d
end
function _legacy_radiation_jet_inv(u::Vector{C}) where {C}
    w = zeros(C, _LEG_NJ)
    w[1] = inv(u[1])
    @inbounds for k in 1:_LEG_NJ-1
        acc = zero(C)
        for j in 1:k
            acc += C(binomial(k, j)) * u[j+1] * w[k-j+1]
        end
        w[k+1] = -acc / u[1]
    end
    w
end

# Dependency-free deterministic PRNG (SplitMix64) — NumericsCore's test target is Test+TestKit+Data,
# so `Random` is intentionally unavailable; a fixed literal seed gives reproducibility with no dep.
mutable struct _SM; s::UInt64; end
function _next!(g::_SM)
    g.s += 0x9e3779b97f4a7c15
    z = g.s
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    z ⊻ (z >> 31)
end
# uniform real in (−3, 3): Float64 uses one 53-bit draw; BigFloat combines two for >53 bits of mantissa.
_randreal(g::_SM, ::Type{Float64}) = (Float64(_next!(g) >> 11) / 9007199254740992.0 * 2 - 1) * 3
function _randreal(g::_SM, ::Type{BigFloat})
    hi = BigFloat(_next!(g) >> 11)
    lo = BigFloat(_next!(g) >> 11)
    x = (hi + lo / (BigFloat(2)^53)) / (BigFloat(2)^53)     # in [0, ~1)
    (2x - 1) * 3
end
# deterministic random jet of eltype T and length n.
function _randjet(g::_SM, ::Type{T}, n) where {T}
    if T <: Complex
        R = real(T)
        return T[complex(_randreal(g, R), _randreal(g, R)) for _ in 1:n]
    else
        return T[_randreal(g, T) for _ in 1:n]
    end
end

@testset verbose = true "jets (derivative-tower Leibniz algebra)" begin

    # ───────────────────────── G-A1: analytic algebra identities ─────────────────────────
    @testset "G-A1 Leibniz product on r², r³, exp(r) vs closed form" begin
        C = ComplexF64
        r = 2.3
        u = C[r^3, 3r^2, 6r, 6, 0]      # jet of r³ (derivative values at r)
        v = C[r^2, 2r, 2, 0, 0]         # jet of r²
        # r³·r² = r⁵ with exact derivative values [r⁵,5r⁴,20r³,60r²,120r]
        @test maximum(abs.(jmul(u, v) .- C[r^5, 5r^4, 20r^3, 60r^2, 120r])) < 1e-10
        # length-generic: length taken from f
        @test length(jmul(u, v)) == length(u)
        # jderiv shift: d/dr of the r³ jet == the 3r² jet (top zero-padded)
        @test jderiv(u) == C[3r^2, 6r, 6, 0, 0]
        # jadd / jscale
        @test jadd(u, v) == u .+ v
        @test jscale(2.0, u) == 2.0 .* u

        # exp(r): every derivative value equals eʳ, so the jet is fill(eʳ, N+1). Leibniz gives
        # (eʳ·eʳ)^(k) = e^{2r} Σ_j C(k,j) = 2ᵏ e^{2r} — an analytic cross-check independent of code.
        N = 6
        er = exp(r)
        ej = fill(C(er), N+1)
        prod = jmul(ej, ej)
        for k in 0:N
            @test isapprox(prod[k+1], 2.0^k * exp(2r); rtol=1e-12)
        end
    end

    @testset "G-A1 reciprocal identity jinv∘jmul == unit jet" begin
        C = ComplexF64
        # 1/(r²−2r) via jinv of the Δ=r²−2r jet, then jmul back to the unit jet [1,0,0,…]
        r = 3.7
        Δ = C[r^2 - 2r, 2r - 2, 2, 0, 0, 0, 0]
        iΔ = jinv(Δ)
        @test isapprox(iΔ[1], 1/(r^2 - 2r); rtol=1e-13)
        unit = jmul(Δ, iΔ)
        @test isapprox(unit[1], 1; atol=1e-12)
        @test maximum(abs.(unit[2:end])) < 1e-10
        # complex leading value: 1/z reciprocal jet, leading == 1/z
        z = C(1.3, -0.7)
        fz = C[z, 1, 0, 0, 0]          # jet of (r + (z-r)) ≈ constant+linear; leading z
        @test isapprox(jinv(fz)[1], one(C)/z; rtol=1e-14)
    end

    # ───────────────────── G-A2: characterization — new == legacy bit-exact ─────────────────────
    @testset "G-A2 new == legacy Lorenz (bit-exact, ==) — Float64/BigFloat × Real/Complex × len 3/5/7" begin
        setprecision(BigFloat, 256) do
            for T in (Float64, ComplexF64, BigFloat, Complex{BigFloat})
                for n in (3, 5, 7)
                    g = _SM(0x9e3779b900000001)   # fixed literal seed
                    for _ in 1:200
                        f = _randjet(g, T, n)
                        gg = _randjet(g, T, n)
                        while abs(f[1]) < 0.1
                            f = _randjet(g, T, n)
                        end
                        s = _randjet(g, T, 1)[1]
                        @test jmul(f, gg)     == _legacy_lorenz_jmul(f, gg)
                        @test jadd(f, gg)     == _legacy_lorenz_jadd(f, gg)
                        @test jscale(s, f)    == _legacy_lorenz_jscale(s, f)
                        @test jderiv(f)       == _legacy_lorenz_jderiv(f)
                        @test jinv(f)         == _legacy_lorenz_jinv(f)
                        @test constjet(s, n-1, T) == _legacy_lorenz_constjet(s, n-1, T)
                    end
                end
            end
        end
    end

    @testset "G-A2 new == legacy Radiation (fixed NJ=5)" begin
        setprecision(BigFloat, 256) do
            for T in (Float64, ComplexF64, BigFloat, Complex{BigFloat})
                g = _SM(0x00000000deadbeef)   # fixed literal seed
                for _ in 1:400
                    u = _randjet(g, T, _LEG_NJ)
                    v = _randjet(g, T, _LEG_NJ)
                    while abs(u[1]) < 0.1
                        u = _randjet(g, T, _LEG_NJ)
                    end
                    # jmul / jadd / jderiv — bit-exact for every eltype
                    @test jmul(u, v)  == _legacy_radiation_jet_mul(u, v)
                    @test jadd(u, v)  == _legacy_radiation_jet_add(u, v)
                    @test jderiv(u)   == _legacy_radiation_jet_deriv(u)
                    # jinv — bit-exact for REAL eltype; for COMPLEX eltype the legacy Radiation leading
                    # term uses inv(u[1]) while canonical (Lorenz) uses one/u[1]; Julia inv(::Complex)
                    # ≠ one/z in the last ulp (documented FP-library artifact), so agreement is to ε.
                    if T <: Complex
                        @test isapprox(jinv(u), _legacy_radiation_jet_inv(u); rtol=1e-14, atol=1e-300)
                    else
                        @test jinv(u) == _legacy_radiation_jet_inv(u)
                    end
                end
            end
        end
    end

    # ───────────────────── BigFloat inheritance + 256-bit analytic value check ─────────────────────
    @testset "BigFloat inheritance: eltype preserved + ~1e-70 analytic agreement" begin
        setprecision(BigFloat, 256) do
            CB = Complex{BigFloat}
            r = big"3.7"; a = big"0.99"
            Δ = CB[r^2 - 2r + a^2, 2r - 2, big(2), 0, 0, 0, 0]
            iΔ = jinv(Δ)
            # eltype flows through untouched (no Float64 pin)
            @test eltype(iΔ) === CB
            @test eltype(jmul(Δ, iΔ)) === CB
            @test eltype(jderiv(Δ)) === CB
            @test eltype(jadd(Δ, iΔ)) === CB
            @test eltype(jscale(CB(2), Δ)) === CB
            @test eltype(constjet(1, 6, CB)) === CB
            # 256-bit value: reciprocal leading term vs analytic 1/(r²−2r+a²) to ~1e-70
            @test abs(iΔ[1] - 1/(r^2 - 2r + a^2)) < big"1e-70"
            # jinv∘jmul == unit jet to ~1e-70
            unit = jmul(Δ, iΔ)
            @test abs(unit[1] - 1) < big"1e-70"
            @test maximum(abs.(unit[2:end])) < big"1e-70"
            # exp(r) Leibniz identity in 256-bit: (eʳ·eʳ)^(k) = 2ᵏ e^{2r}, rel err ≤ 1e-70
            N = 6
            ej = fill(CB(exp(r)), N+1)
            prodj = jmul(ej, ej)
            for k in 0:N
                @test abs(prodj[k+1] - 2^k * exp(2r)) / abs(2^k * exp(2r)) < big"1e-70"
            end
        end
    end

    # ─────────── G-A2 (NTuple): allocation-free tuple methods == legacy length-4 j4* kernel ───────────
    # Third legacy home (Move F2): the hand-unrolled length-4 TUPLE kernel from MetricReconstruction
    # Lorenz/SchwarzschildComponents.jl:46-59, copied VERBATIM below as a permanent regression oracle.
    # The `NTuple` jmul/jinv methods + the untyped jadd/jscale acting on tuples (`f.+g`, `s.*f`) must
    # reproduce it BIT-EXACT and allocate ZERO heap. Comparison uses isequal (distinguishes ±0 and works
    # by value on the mutable BigFloat, where === is object identity), Float64/BigFloat-256 × Real/Complex.
    @inline _legacy_j4mul(f, g) = (f[1]*g[1],
                                   f[1]*g[2]+f[2]*g[1],
                                   f[1]*g[3]+2f[2]*g[2]+f[3]*g[1],
                                   f[1]*g[4]+3f[2]*g[3]+3f[3]*g[2]+f[4]*g[1])
    @inline _legacy_j4add(f, g)   = (f[1]+g[1], f[2]+g[2], f[3]+g[3], f[4]+g[4])
    @inline _legacy_j4scale(s, f) = (s*f[1], s*f[2], s*f[3], s*f[4])
    function _legacy_j4inv(f)
        g1 = inv(f[1])
        g2 = -(1*f[2]*g1)*g1
        g3 = -(2*f[2]*g2 + 1*f[3]*g1)*g1
        g4 = -(3*f[2]*g3 + 3*f[3]*g2 + 1*f[4]*g1)*g1
        (g1, g2, g3, g4)
    end
    _beq(x::Real, y::Real) = isequal(x, y)
    _beq(x::Complex, y::Complex) = isequal(real(x), real(y)) && isequal(imag(x), imag(y))
    _tupbeq(a, b) = length(a) == length(b) && all(_beq(x, y) for (x, y) in zip(a, b))
    _randtup4(g, ::Type{T}) where {T} = ntuple(_ -> _randjet(g, T, 1)[1], 4)

    @testset "G-A2 NTuple == legacy length-4 j4* kernel (bit-exact) — Float64/BigFloat × Real/Complex" begin
        setprecision(BigFloat, 256) do
            for T in (Float64, ComplexF64, BigFloat, Complex{BigFloat})
                g = _SM(0x0000000012345678)   # fixed literal seed
                for _ in 1:2000
                    f  = _randtup4(g, T)
                    gg = _randtup4(g, T)
                    while abs(f[1]) < 0.1
                        f = _randtup4(g, T)
                    end
                    s = _randjet(g, T, 1)[1]
                    @test _tupbeq(jmul(f, gg),  _legacy_j4mul(f, gg))
                    @test _tupbeq(jadd(f, gg),  _legacy_j4add(f, gg))
                    @test _tupbeq(jscale(s, f), _legacy_j4scale(s, f))
                    @test _tupbeq(jinv(f),      _legacy_j4inv(f))
                    # output stays a statically-sized NTuple (not promoted to a Vector)
                    @test jmul(f, gg) isa NTuple{4}
                    @test jinv(f)     isa NTuple{4}
                end
            end
        end
    end

    @testset "NTuple jmul/jinv/jadd/jscale allocate zero heap (ComplexF64, matched-eltype scalar)" begin
        f = (ComplexF64(1.7,-0.4), ComplexF64(-0.3,0.9), ComplexF64(0.11,0.02), ComplexF64(-0.02,0.05))
        g = (ComplexF64(0.9,0.2),  ComplexF64(0.4,-0.1), ComplexF64(-0.05,0.3), ComplexF64(0.007,-0.11))
        s = ComplexF64(0.0, 2.0)
        jmul(f, g); jinv(f); jadd(f, g); jscale(s, f)   # warm (compile)
        @test @allocated(jmul(f, g))   == 0
        @test @allocated(jinv(f))      == 0
        @test @allocated(jadd(f, g))   == 0
        @test @allocated(jscale(s, f)) == 0
    end

    @testset "NTuple methods leave the generic-length Vector jet path unchanged (dispatch guard)" begin
        @test jmul(Float64[1,2,3,4,5], Float64[1,1,1,1,1]) == Float64[1,3,8,20,48]
        @test jmul(Float64[1,2,3,4,5], Float64[1,1,1,1,1]) isa Vector{Float64}
        @test jinv(Float64[2,1,0,0,0]) isa Vector{Float64}
    end
end
