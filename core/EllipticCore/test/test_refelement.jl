# test_refelement.jl — grounded gates for src/refelement.jl (DAG node: mesh-basis-quadrature,
# elliptic-ssf GLOBAL_DAG row 8 + TOPOLOGICAL_SCHEDULE W2-A). Written BEFORE the implementation
# (TDD red run captured in evidence/mesh-basis-quadrature/pkg_test_red.log).
#
# Oracle gates (GLOBAL_DAG row 8 + node contract):
#   * polynomial interpolation exact through degree N−1: 1e-13 @ Float64, 1e-60 @ BigFloat(256)
#   * Gauss quadrature exact for degree ≤ 2N−1 (falsifier: degree 2N is NOT exact — Gauss, not magic)
#   * differentiation matrix D exact for degree ≤ N−1
#   * generalized-SBP identity (precise statement): with W = Diagonal(weights) and the face-lift
#     matrix B = face_hi*face_hi' − face_lo*face_lo', the EXACT matrix identity
#         W*D + D'*W == B
#     holds for Legendre-GAUSS collocation. Because Gauss nodes exclude ±1, the boundary operator
#     is carried by the face-interpolation vectors (generalized SBP); the pure antisymmetry
#     W*D == −(W*D)' is FALSE — its defect is exactly the rank-≤2 matrix B. Proof: nodal vectors
#     ↔ polynomials of degree ≤ N−1 bijectively; for such u,v the integrand u′v + uv′ has degree
#     ≤ 2N−3 ≤ 2N−1 so quadrature is exact and uᵀ(WD + DᵀW)v = ∫(u′v+uv′) = [uv]₋₁¹ = uᵀBv.
#   * face-interpolation exactness at ξ = ±1 through degree N−1
#   * Runge 1/(1+25x²): exponential error decay in N (log-err slope < −0.15; theory −log ρ ≈ −0.199,
#     ρ = 1/5 + √(26)/5), identical gate at BigFloat(256)
#   * cache keyed on (real scalar type, N, BigFloat precision): hit/miss counters, ===-identity on
#     hits, ComplexF64 normalizes to the Float64 entry, BigFloat 256 vs 512 are distinct entries
#   * tensor-product 2D index helpers: exhaustive round-trip + per-dimension line slicing grounded
#     against Julia's own column-major LinearIndices/reshape

using Test
using LinearAlgebra
using EllipticCore
using NumericsCore: gl_nodes_weights

# ---- tiny deterministic generator (keeps the stdlib Random OUT of the test env) ----------------
mutable struct LCG
    s::UInt64
end
# Knuth MMIX constants; top 53 bits → uniform in [-1, 1)
rand_unit!(r::LCG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f;
                      Int64(r.s >> 11) / 2.0^53 * 2 - 1)
randcoeffs!(r::LCG, deg::Int) = [rand_unit!(r) for _ in 0:deg]

# Horner evaluation of coeffs c (ascending powers) at x — promotes to the wider of the two types.
function polyval(c::AbstractVector, x)
    acc = c[end] * one(x)
    for k in length(c)-1:-1:1
        acc = acc * x + c[k]
    end
    return acc
end
polyderiv(c::AbstractVector) = [(k - 1) * c[k] for k in 2:length(c)]
# exact ∫₋₁¹ of the poly with coeffs c, in RT arithmetic (odd powers vanish)
polyint_m11(::Type{RT}, c::AbstractVector) where {RT} =
    sum(RT(2) * RT(c[k]) / RT(k) for k in 1:2:length(c))

# least-squares slope of y vs x
function lsq_slope(xs::Vector{Float64}, ys::Vector{Float64})
    n = length(xs)
    sx, sy = sum(xs), sum(ys)
    return (n * sum(xs .* ys) - sx * sy) / (n * sum(abs2, xs) - sx^2)
end

@testset "refelement (mesh-basis-quadrature)" begin

exactness_tol = Dict{DataType,Any}(Float64 => 1e-13, BigFloat => big"1e-60")

@testset "exactness gates [RT=$RT]" for RT in (Float64, BigFloat)
    setprecision(BigFloat, 256) do
        tol = exactness_tol[RT]
        rng = LCG(0x00c0ffee_c0ffee01)

        for N in 2:12
            el = refelement(RT, N)
            x, w = el.nodes, el.weights

            @testset "N=$N: nodes/weights delegate to gl_nodes_weights" begin
                xref, wref = gl_nodes_weights(RT, N)
                @test x == xref
                @test w == wref
                @test issorted(x) && x[1] > -1 && x[end] < 1   # interior (Gauss, not Lobatto)
                @test all(>(0), w)
                @test abs(sum(w) - 2) <= tol
            end

            @testset "N=$N: interpolation exact through degree N-1" begin
                # collision path: interpolating to the nodes themselves gives exact unit rows
                @test interpolation_matrix(el, x) == Matrix{RT}(I, N, N)
                pts = vcat(-one(RT), one(RT), RT[RT(rand_unit!(rng)) for _ in 1:15])
                P = interpolation_matrix(el, pts)
                @test size(P) == (length(pts), N)
                for k in 0:N-1                                   # monomial basis: complete gate
                    @test maximum(abs, P * (x .^ k) .- pts .^ k) <= tol
                end
                c = randcoeffs!(rng, N - 1)                      # random-poly property gate
                @test maximum(abs, P * polyval.(Ref(c), x) .- polyval.(Ref(c), pts)) <= tol
            end

            @testset "N=$N: quadrature exact through degree 2N-1" begin
                for k in 0:2N-1
                    exact = iseven(k) ? RT(2) / RT(k + 1) : zero(RT)
                    @test abs(dot(w, x .^ k) - exact) <= tol
                end
                c = randcoeffs!(rng, 2N - 1)
                @test abs(dot(w, polyval.(Ref(c), x)) - polyint_m11(RT, c)) <=
                      tol * (1 + sum(abs, c))
            end

            @testset "N=$N: differentiation exact through degree N-1" begin
                for k in 0:N-1
                    dk = k == 0 ? zero(x) : k .* x .^ (k - 1)
                    @test maximum(abs, el.D * (x .^ k) .- dk) <= tol * N^2
                end
                # coefficients must be RT BEFORE differentiating: (k-1)*c[k] rounds at Float64
                # (~1e-16 rel), which would corrupt the 1e-60 BigFloat oracle, not test D.
                c = RT.(randcoeffs!(rng, N - 1))
                dc = polyderiv(c)
                dref = isempty(dc) ? zero(x) : polyval.(Ref(dc), x)
                @test maximum(abs, el.D * polyval.(Ref(c), x) .- dref) <= tol * N^2
            end

            @testset "N=$N: face interpolation exact at xi=+-1" begin
                for k in 0:N-1
                    @test abs(dot(el.face_lo, x .^ k) - (-one(RT))^k) <= tol
                    @test abs(dot(el.face_hi, x .^ k) - one(RT)) <= tol
                end
            end

            @testset "N=$N: generalized SBP  W*D + D'*W == face-lift B" begin
                W = Diagonal(w)
                B = el.face_hi * el.face_hi' - el.face_lo * el.face_lo'
                @test maximum(abs, W * el.D + el.D' * W - B) <= tol * N^2
                # the PURE antisymmetry is false for Gauss: the defect is the rank-<=2 lift
                @test maximum(abs, W * el.D + (W * el.D)') > 1000 * tol
            end
        end
    end
end

@testset "quadrature falsifier: degree 2N is NOT exact (Gauss, no free lunch)" begin
    for N in 2:6
        el = refelement(Float64, N)
        err = abs(dot(el.weights, el.nodes .^ (2N)) - 2 / (2N + 1))
        @test err > 1e-8
    end
end

@testset "Runge 1/(1+25x^2): exponential decay [RT=$RT]" for RT in (Float64, BigFloat)
    setprecision(BigFloat, 256) do
        f = x -> 1 / (1 + 25 * x^2)
        xs = collect(range(-one(RT), one(RT), length=401))
        Ns = collect(8:4:32)
        errs = map(Ns) do N
            el = refelement(RT, N)
            maximum(abs, interpolation_matrix(el, xs) * f.(el.nodes) .- f.(xs))
        end
        @test all(diff(errs) .< 0)                        # strictly decreasing along the ladder
        slope = lsq_slope(Float64.(Ns), Float64.(log.(errs)))
        @test slope < -0.15                               # theory: -log(rho) ~ -0.199
    end
end

@testset "cache: keyed (RT, N, precision), hit/miss counters, ===-identity" begin
    EllipticCore.refelement_cache_reset!()
    @test EllipticCore.refelement_cache_stats() == (hits = 0, misses = 0)

    a = refelement(Float64, 7)
    @test EllipticCore.refelement_cache_stats() == (hits = 0, misses = 1)
    b = refelement(Float64, 7)
    @test b === a
    @test EllipticCore.refelement_cache_stats() == (hits = 1, misses = 1)

    c = refelement(ComplexF64, 7)                # geometry is real: normalizes to Float64 entry
    @test c === a
    @test EllipticCore.refelement_cache_stats() == (hits = 2, misses = 1)

    d = refelement(Float64, 8)
    @test d !== a
    @test EllipticCore.refelement_cache_stats() == (hits = 2, misses = 2)

    # BigFloat entries are precision-keyed: a 256-bit element must never serve a 512-bit caller
    e256 = setprecision(() -> refelement(BigFloat, 5), BigFloat, 256)
    e512 = setprecision(() -> refelement(BigFloat, 5), BigFloat, 512)
    @test e256 !== e512
    @test precision(e256.nodes[1]) == 256
    @test precision(e512.nodes[1]) == 512
    @test EllipticCore.refelement_cache_stats() == (hits = 2, misses = 4)
    e512b = setprecision(() -> refelement(Complex{BigFloat}, 5), BigFloat, 512)
    @test e512b === e512
    @test EllipticCore.refelement_cache_stats() == (hits = 3, misses = 4)
end

# G16 lock-regression gate (repro-audit G16 claimed this cache was racy; it is NOT — the lock has
# been there since the package was written — so this testset pins the property instead of "fixing"
# it). `refelement` is called from the THREADED DG matrix/apply layers, so a lost update would give
# two callers DIFFERENT element objects for one key and build the same key twice.
# MUTATION-TESTED: with `lock(_REFELEMENT_LOCK)` deleted from src/refelement.jl, the same 480 calls
# at -t4 gave misses = 32 (not 12) and the ===-identity check FALSE. Probabilistic — it detects a
# removed lock with high probability, it does not prove the absence of every interleaving.
@testset "cache: concurrent refelement() is race-free (G16 lock gate)" begin
    if Threads.nthreads() < 2
        @test_skip "concurrent cache gate needs -t2"   # single-threaded run: no interleaving to expose
    else
        EllipticCore.refelement_cache_reset!()
        Ns = collect(4:15)                                  # 12 distinct cache keys
        jobs = [Ns[((k - 1) % length(Ns)) + 1] for k in 1:(40 * length(Ns))]
        got = Vector{Any}(undef, length(jobs))
        Threads.@threads for k in eachindex(jobs)
            got[k] = refelement(Float64, jobs[k])
        end
        canonical = Dict{Int,Any}()
        for k in eachindex(jobs)
            canonical[jobs[k]] = get(canonical, jobs[k], got[k])
        end
        @test all(k -> got[k] === canonical[jobs[k]], eachindex(jobs))  # ONE object per key
        st = EllipticCore.refelement_cache_stats()
        @test st.misses == length(Ns)                       # every key built EXACTLY once
        @test st.hits + st.misses == length(jobs)           # every call went through the cache
        EllipticCore.refelement_cache_reset!()              # leave no state for later testsets
    end
end

@testset "tensor-product index helpers (exhaustive, grounded on column-major)" begin
    for Nx in 1:6, Ny in 1:6
        L = LinearIndices((Nx, Ny))
        A = reshape(1:Nx*Ny, Nx, Ny)
        for j in 1:Ny, i in 1:Nx
            k = tensor_lin(i, j, Nx)
            @test k == L[i, j]                        # column-major, i fastest
            @test tensor_ij(k, Nx) == (i, j)          # exact round trip
        end
        # per-dimension line slicing: ranges (viewable), matching reshape rows/columns
        for j in 1:Ny
            line = tensor_line(1, j, Nx, Ny)
            @test line isa AbstractRange
            @test collect(line) == A[:, j]
        end
        for i in 1:Nx
            line = tensor_line(2, i, Nx, Ny)
            @test line isa AbstractRange
            @test collect(line) == A[i, :]
        end
        # the lines of either dimension partition 1:Nx*Ny exactly
        @test sort(vcat([collect(tensor_line(1, j, Nx, Ny)) for j in 1:Ny]...)) == collect(1:Nx*Ny)
        @test sort(vcat([collect(tensor_line(2, i, Nx, Ny)) for i in 1:Nx]...)) == collect(1:Nx*Ny)
    end
    @test_throws ArgumentError tensor_line(3, 1, 2, 2)
end

end # testset refelement
