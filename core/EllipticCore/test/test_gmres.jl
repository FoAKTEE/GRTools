# test_gmres.jl — grounded gates for src/gmres.jl (DAG node: complex-gmres, elliptic-ssf
# GLOBAL_DAG row 16 + TOPOLOGICAL_SCHEDULE W3-B). Written BEFORE the implementation (TDD red
# run captured in evidence/complex-gmres/pkg_test_red.log).
#
# Target: parametric-T restarted GMRES (Arnoldi/MGS + complex Givens QR of the Hessenberg),
# ADAPTED from the SpecCore reference (specjl SpecCore/src/elliptic/solvers/gmres.jl), with
#   * complex Givens rotations in the standard complex-safe (LAPACK zlartg) form: c real,
#     s complex, rotation [c s; -conj(s) c];
#   * matrix-free A as any callable applyA!(out, v);
#   * right-preconditioning hook Minv!(z, v) (identity default) — the CRP's preconditioned
#     usage (main.tex:670-674); right preconditioning preserves the TRUE residual, so the
#     Givens estimate history tracks ‖b − A x‖ itself;
#   * the SpecCore SPD monotone-residual HARD ERROR REMOVED (audit-speccore risk 3), replaced
#     by diagnostics valid for non-Hermitian complex operators: happy-breakdown detection
#     (:happy), restart-cycle stagnation (:stagnation), residual-history reporting.
#
# Oracle gates (worker node contract; all reference numbers calibrated BEFORE this file ran,
# scratchpad calibrate_gmres{,_restart,_pc}.jl, independent textbook constructions only):
#   * complex Givens vs LinearAlgebra.givens on COMPLEX inputs: c/s/r parity ≤ 1e-13 (200
#     random pairs; calibrated max deviation 4.6e-16). For REAL inputs Julia's givens uses a
#     different SIGN convention (calibrated deviation 2.58 on sign flips), so real + BigFloat
#     inputs are gated on the DEFINING properties instead: c real ≥ 0, c²+|s|² = 1,
#     -conj(s)f + cg = 0, cf + sg = r — which is what GMRES correctness rests on.
#   * known small complex systems vs dense `\` at 1e-12: random well-conditioned
#     non-Hermitian (n=60, calibrated 14 iters), Hermitian PD (calibrated 22 iters), 1D
#     complex-shifted Helmholtz tridiagonal (n=64, k²=36, ε=10; calibrated 64 iters =
#     Krylov exhaustion — unpreconditioned Helmholtz is genuinely hard).
#   * reference-iterate parity vs an in-test textbook GMRES (dense Arnoldi + FRESH Hessenberg
#     least-squares via `\` each step — algorithmically independent of the incremental-Givens
#     implementation): iteration count ±1, per-iterate residual norms to 1e-6 rel + 1e-11·β.
#   * restart correctness: restart ≥ n ≡ unrestarted (identical arithmetic path); a genuinely
#     restarted run (restart=5) still converges (calibrated 14 total inner iters).
#   * complex-vs-2-real-fields gate (CRP main.tex:700-703 claim 1): complexified Helmholtz
#     converges in ≤ the iterations of the equivalent [Ar -Ai; Ai Ar] block system
#     (calibrated 64 vs 128 at n=64, 48 vs 96 at n=48; the real-block run also exercises the
#     real-T instantiation of the same parametric code).
#   * residual-history monotone nonincreasing (GMRES minimization property, which DOES hold,
#     unlike SPD-ness) with the identity right-preconditioner, and identity-hook invariance;
#     nontrivial right preconditioner (Laplacian-preconditioned Helmholtz, the classic
#     shifted-Laplacian setup): calibrated 12 iters vs 64 unpreconditioned, monotone.
#   * diagnostics: happy breakdown on a 3-distinct-eigenvalue operator at exactly 3
#     iterations (calibrated) and on αI at 1; stagnation on Saad's cyclic-shift example
#     (restart < n makes EXACTLY zero progress per cycle — calibrated ‖Δx‖ = 0) flagged
#     :stagnation instead of burning maxiter; plain maxiter truncation stays :none.
#   * tolerance semantics max(abstol, reltol·‖b‖) (yaml: AbsoluteResidual 1e-12,
#     RelativeResidual 0, MaxIterations 500 defaults); b = 0 and warm-start x0 edge cases.
#   * BigFloat(256) smoke: Complex{BigFloat} system solved to 1e-60 (Float64 arithmetic
#     cannot pass; dense `\` oracle residual calibrated 3.7e-77), precision preserved.
#
# Test helpers are deliberately UNtyped so this file loads before the implementation exists
# (each @testset then errors individually in the red run without aborting sibling files).

using Test
using LinearAlgebra
using EllipticCore

# ---- tiny deterministic generator (keeps stdlib Random OUT of the test env; same Knuth-MMIX
# ---- LCG convention as the sibling test files, distinct name/seed) --------------------------
mutable struct GMRNG
    s::UInt64
end
gmr_unit!(r::GMRNG) = (r.s = 0x5851f42d4c957f2d * r.s + 0x14057b7ef767814f;
                       Int64(r.s >> 11) / 2.0^53 * 2 - 1)
gmr_c!(r::GMRNG) = complex(gmr_unit!(r), gmr_unit!(r))

# matrix-free closure over a dense/structured matrix ("A as any callable")
applier(A) = (out, v) -> mul!(out, A, v)

# 1D complex-shifted Helmholtz tridiagonal: (1/h²)tridiag(-1,2,-1) − (k²+iε)I on n interior
# points of (0,1) — the CRP's Helmholtz-type model problem shape (main.tex:680 item (2)).
function helmholtz1d(n, k2, eps)
    h = 1.0 / (n + 1)
    A = zeros(ComplexF64, n, n)
    for i in 1:n
        A[i, i] = 2 / h^2 - (k2 + im * eps)
        i > 1 && (A[i, i - 1] = -1 / h^2)
        i < n && (A[i, i + 1] = -1 / h^2)
    end
    return A
end

# equivalent 2×-real block system of A z = b: [Ar -Ai; Ai Ar][zr; zi] = [br; bi]
realblock(A) = [real(A) -imag(A); imag(A) real(A)]

# the calibrated fixtures (single source in this file; seeds fixed)
function wellcond_fixture()
    rng = GMRNG(0x5eed_2026_0728_0003)
    n = 60
    R = [gmr_c!(rng) for _ in 1:n, _ in 1:n]
    A = Matrix{ComplexF64}(3I, n, n) + R / (2 * sqrt(n))
    B = [gmr_c!(rng) for _ in 1:n, _ in 1:n] / sqrt(n)   # same draw order as calibration
    b = [gmr_c!(rng) for _ in 1:n]
    b ./= norm(b)
    return A, B' * B + I, b        # (well-conditioned non-Hermitian, Hermitian PD, rhs)
end
function helmholtz_fixture(n)
    A = helmholtz1d(n, 36.0, 10.0)
    rngb = GMRNG(UInt64(n) * 7919 + 36)
    b = [gmr_c!(rngb) for _ in 1:n]
    b ./= norm(b)
    return A, b
end

# in-test textbook GMRES oracle: dense Arnoldi/MGS + a FRESH Hessenberg least-squares solve
# via `\` at every step (no incremental Givens — algorithmically independent reference).
# Returns (x, iters, resnorms) with resnorms[1] = ‖b‖, x0 = 0.
function textbook_gmres(A, b; tol, maxiter)
    T = eltype(b)
    n = length(b)
    β = norm(b)
    resnorms = [β]
    β ≤ tol && return (zeros(T, n), 0, resnorms)
    V = [b / β]
    H = zeros(T, maxiter + 1, maxiter)
    x = zeros(T, n)
    for j in 1:maxiter
        w = A * V[j]
        for i in 1:j
            H[i, j] = dot(V[i], w)
            w -= H[i, j] * V[i]
        end
        H[j + 1, j] = norm(w)
        lucky = abs(H[j + 1, j]) ≤ 1e-14 * β
        lucky || push!(V, w / H[j + 1, j])
        e1 = zeros(T, j + 1)
        e1[1] = β
        y = H[1:j + 1, 1:j] \ e1
        push!(resnorms, norm(e1 - H[1:j + 1, 1:j] * y))
        if resnorms[end] ≤ tol || lucky
            for i in 1:j
                x += y[i] * V[i]
            end
            return (x, j, resnorms)
        end
    end
    j = maxiter
    e1 = zeros(T, j + 1)
    e1[1] = β
    y = H[1:j + 1, 1:j] \ e1
    for i in 1:j
        x += y[i] * V[i]
    end
    return (x, j, resnorms)
end

# residual-history monotonicity (allow incremental-Givens/restart roundoff slack)
monotone(hist) = all(hist[k + 1] ≤ hist[k] * (1 + 1e-8) + 1e-15 * hist[1]
                     for k in 1:length(hist) - 1)

@testset "gmres (complex-gmres)" begin

@testset "complex Givens rotation: parity with LinearAlgebra.givens + defining properties" begin
    # COMPLEX inputs: full parity with Julia's givens (both follow the zlartg convention;
    # calibrated max deviation 4.6e-16 over 200 pairs)
    rng = GMRNG(0x0badc0de_2026_0728)
    for _ in 1:200
        f, g = gmr_c!(rng), gmr_c!(rng)
        c, s, r = EllipticCore.givens_rotation(f, g)
        G, rj = givens(f, g, 1, 2)
        @test c isa Real
        @test abs(G.c - c) ≤ 1e-13
        @test abs(G.s - s) ≤ 1e-13
        @test abs(rj - r) ≤ 1e-13
    end
    # REAL + BigFloat inputs: Julia's real givens uses a different sign convention
    # (calibrated deviation 2.58), so gate the DEFINING properties instead.
    function check_props(f, g; tol)
        c, s, r = EllipticCore.givens_rotation(f, g)
        @test c isa Real && c ≥ 0
        @test abs(c^2 + abs2(s) - 1) ≤ tol
        @test abs(-conj(s) * f + c * g) ≤ tol * max(abs(f), abs(g))   # annihilation
        @test abs(c * f + s * g - r) ≤ tol * max(abs(f), abs(g))      # rotated value
    end
    for _ in 1:100
        check_props(gmr_unit!(rng), gmr_unit!(rng); tol = 1e-14)
        check_props(gmr_c!(rng), gmr_c!(rng); tol = 1e-14)
    end
    setprecision(BigFloat, 256) do
        f = complex(big"1.3", big"-0.7")
        g = complex(big"0.4", big"2.1")
        check_props(f, g; tol = big"1e-70")
        c, _, _ = EllipticCore.givens_rotation(f, g)
        @test c isa BigFloat && precision(c) == 256
    end
    # edge cases (zlartg convention)
    let (c, s, r) = EllipticCore.givens_rotation(2.0 + 1.0im, 0.0 + 0.0im)
        @test c == 1 && s == 0 && r == 2.0 + 1.0im
    end
    let g0 = 2.0 - 1.0im, (c, s, r) = EllipticCore.givens_rotation(0.0 + 0.0im, g0)
        @test c == 0 && abs(s - conj(g0) / abs(g0)) ≤ 1e-15 && abs(r - abs(g0)) ≤ 1e-15
    end
    let (c, s, r) = EllipticCore.givens_rotation(0.0 + 0.0im, 0.0 + 0.0im)
        @test c == 1 && s == 0 && r == 0
    end
end

@testset "known systems vs dense backslash at 1e-12" begin
    Awc, Ah, b = wellcond_fixture()
    Ahm, bh = helmholtz_fixture(64)
    cases = (("well-conditioned non-Hermitian", Awc, b, 20),
             ("Hermitian positive definite", Ah, b, 30),
             ("complex-shifted Helmholtz n=64", Ahm, bh, 65))
    for (name, A, rhs, itcap) in cases
        res = gmres(applier(A), rhs)                       # yaml defaults: 1e-12 abs, 500 max
        @test res isa GMRESResult
        @test res.converged
        @test res.breakdown == :none || res.breakdown == :happy
        @test norm(rhs - A * res.x) ≤ 1e-12               # TRUE residual, recomputed here
        @test res.residual ≤ 1e-12                        # ... and as reported
        @test abs(res.residual - norm(rhs - A * res.x)) ≤ 1e-14
        @test norm(res.x - A \ rhs) ≤ 1e-10               # against the direct solution
        @test 0 < res.iterations ≤ itcap
        @test length(res.history) == res.iterations + 1
        @test abs(res.history[1] - norm(rhs)) ≤ 1e-14
        @test monotone(res.history)
    end
end

@testset "reference-iterate parity vs in-test textbook GMRES" begin
    Awc, _, b = wellcond_fixture()
    Ahm, bh = helmholtz_fixture(64)
    for (A, rhs) in ((Awc, b), (Ahm, bh))
        res = gmres(applier(A), rhs)
        _, itref, rnref = textbook_gmres(A, rhs; tol = 1e-12, maxiter = 500)
        @test abs(res.iterations - itref) ≤ 1              # iteration-count parity ±1
        β = norm(rhs)
        for k in 1:min(length(res.history), length(rnref))
            @test abs(res.history[k] - rnref[k]) ≤ 1e-6 * rnref[k] + 1e-11 * β
        end
    end
end

@testset "restart correctness" begin
    Awc, _, b = wellcond_fixture()
    n = length(b)
    full = gmres(applier(Awc), b)                          # restart defaults to maxiter ≥ n
    atn = gmres(applier(Awc), b; restart = n)              # restart = n ≥ Krylov need
    @test atn.iterations == full.iterations
    @test atn.history == full.history                      # identical arithmetic path
    @test norm(atn.x - full.x) ≤ 1e-13 * norm(full.x)
    # a genuinely restarted run still converges (calibrated: 14 total inner iters at m=5)
    small = gmres(applier(Awc), b; restart = 5)
    @test small.converged
    @test small.iterations ≥ full.iterations               # restarting can never need fewer
    @test small.iterations ≤ 30
    @test norm(b - Awc * small.x) ≤ 1e-12
    @test norm(small.x - Awc \ b) ≤ 1e-10
    @test_throws ArgumentError gmres(applier(Awc), b; restart = 0)
end

@testset "complex-vs-2-real-fields gate (CRP main.tex:700-703 claim 1)" begin
    for n in (48, 64)
        A, b = helmholtz_fixture(n)
        resc = gmres(applier(A), b)
        Ar = realblock(A)
        br = [real(b); imag(b)]
        resr = gmres(applier(Ar), br)                      # real-T instantiation, same code
        @test resc.converged && resr.converged
        @test eltype(resr.x) == Float64
        @test norm(br - Ar * resr.x) ≤ 1e-12
        # the CRP's stated motivation: complex arithmetic restricts the solution space —
        # the complexified solve needs no more iterations than the 2×-real block solve
        # (calibrated: 48 vs 96, 64 vs 128 — Krylov exhaustion at n vs 2n)
        @test resc.iterations ≤ resr.iterations
        # the block solution is the same field, split (checks the two runs solved the
        # same physics): zr + i·zi == z
        @test norm((resr.x[1:n] + im * resr.x[n + 1:2n]) - resc.x) ≤ 1e-10
    end
end

@testset "right preconditioning: identity invariance, monotone history, Laplacian PC" begin
    Ahm, bh = helmholtz_fixture(64)
    plain = gmres(applier(Ahm), bh)
    # explicit identity hook ≡ no hook (must not change the math)
    ident = gmres(applier(Ahm), bh; Minv! = (z, v) -> copyto!(z, v))
    @test ident.iterations == plain.iterations
    @test ident.history ≈ plain.history rtol = 1e-12
    @test norm(ident.x - plain.x) ≤ 1e-12 * norm(plain.x)
    @test monotone(ident.history)                          # GMRES minimization property
    # nontrivial right PC: shifted-Laplacian M = (1/h²)tridiag(-1,2,-1) (the classic
    # Helmholtz preconditioner; calibrated 12 iters vs 64 unpreconditioned, monotone)
    n = 64
    M = real(helmholtz1d(n, -0.0, 0.0))                    # pure discrete Laplacian
    F = lu(Matrix{ComplexF64}(M))
    pc = gmres(applier(Ahm), bh; Minv! = (z, v) -> ldiv!(z, F, v))
    @test pc.converged
    @test norm(bh - Ahm * pc.x) ≤ 1e-12                    # TRUE residual of the ORIGINAL system
    @test pc.iterations ≤ 20 < plain.iterations
    @test monotone(pc.history)                             # right PC preserves the residual —
    @test length(pc.history) == pc.iterations + 1          # history is the true-residual estimate
    @test norm(pc.x - Ahm \ bh) ≤ 1e-10
end

@testset "diagnostics: happy breakdown, stagnation, maxiter truncation" begin
    rng = GMRNG(0x00d1a6_2026_0728)
    # happy breakdown at exactly 3 iterations: diagonal operator with 3 distinct eigenvalues
    # (Krylov space closes at dimension 3; calibrated res 1.4e-15 at j=3)
    n = 30
    vals = ComplexF64[2 + im, 3 - 2im, 5 + 0.5im]
    d = ComplexF64[vals[mod1(i, 3)] for i in 1:n]
    bd = ComplexF64[gmr_c!(rng) for _ in 1:n]
    hb = gmres((out, v) -> (out .= d .* v), bd)
    @test hb.breakdown == :happy
    @test hb.iterations == 3
    @test hb.converged
    @test norm(bd - d .* hb.x) ≤ 1e-12
    # α·I: Krylov space closes at dimension 1
    α = 2.5 - 1.5im
    one_step = gmres((out, v) -> (out .= α .* v), bd)
    @test one_step.breakdown == :happy
    @test one_step.iterations == 1
    @test norm(one_step.x - bd / α) ≤ 1e-14 * norm(bd)
    # stagnation: Saad's cyclic-shift example — restart < n makes EXACTLY zero progress per
    # cycle (calibrated ‖Δx‖ = 0.0), so the restarted solver must FLAG it instead of
    # spinning to maxiter (the non-Hermitian replacement for SpecCore's SPD hard error)
    ns = 20
    P = zeros(ComplexF64, ns, ns)
    for i in 1:ns
        P[mod1(i + 1, ns), i] = 1
    end
    e1 = zeros(ComplexF64, ns)
    e1[1] = 1
    st = gmres(applier(P), e1; restart = 5, maxiter = 100)
    @test st.breakdown == :stagnation
    @test !st.converged
    @test st.iterations < 100                              # stopped early, did not burn budget
    @test abs(st.history[end] - 1.0) ≤ 1e-12               # residual truly never moved
    # ... while the UNRESTARTED run on the same operator converges exactly at j = n
    # (flat history then a plunge — monotone-nonincreasing throughout)
    unre = gmres(applier(P), e1)
    @test unre.converged
    @test unre.iterations == ns
    @test monotone(unre.history)
    # plain maxiter truncation (no stagnation misdiagnosis): Helmholtz cut off at 5
    Ahm, bh = helmholtz_fixture(64)
    cut = gmres(applier(Ahm), bh; maxiter = 5)
    @test !cut.converged
    @test cut.breakdown == :none
    @test cut.iterations == 5
    @test cut.residual ≥ 1e-12                             # honest: not converged
    @test abs(cut.residual - norm(bh - Ahm * cut.x)) ≤ 1e-14
end

@testset "tolerance semantics + edge cases" begin
    Awc, _, b = wellcond_fixture()
    # b = 0 → x = 0 in zero iterations
    z = gmres(applier(Awc), zeros(ComplexF64, length(b)))
    @test z.converged && z.iterations == 0 && norm(z.x) == 0 && z.residual == 0
    # warm start at the exact solution → zero iterations
    xex = Awc \ b
    warm = gmres(applier(Awc), b; x0 = xex)
    @test warm.converged && warm.iterations == 0
    @test norm(warm.x - xex) == 0
    # warm start from a perturbed point still converges to the right answer
    warm2 = gmres(applier(Awc), b; x0 = xex .+ 0.1)
    @test warm2.converged && norm(b - Awc * warm2.x) ≤ 1e-12 && warm2.iterations > 0
    # relative tolerance: ‖r‖ ≤ max(abstol, reltol·‖b‖) with a scaled rhs
    bigb = b * 1e6
    rel = gmres(applier(Awc), bigb; abstol = 0.0, reltol = 1e-10)
    @test rel.converged
    @test norm(bigb - Awc * rel.x) ≤ 1e-10 * norm(bigb)
    # a LOOSE abstol stops earlier than the tight default
    loose = gmres(applier(Awc), b; abstol = 1e-4)
    tight = gmres(applier(Awc), b)
    @test loose.converged && loose.iterations < tight.iterations
    @test norm(b - Awc * loose.x) ≤ 1e-4
    @test_throws ArgumentError gmres(applier(Awc), b; maxiter = -1)
    @test_throws ArgumentError gmres(applier(Awc), b; x0 = zeros(ComplexF64, 3))
end

@testset "BigFloat(256) smoke" begin
    setprecision(BigFloat, 256) do
        rng = GMRNG(0xb16_f10a7)
        n = 12
        gB!(r) = complex(BigFloat(gmr_unit!(r)), BigFloat(gmr_unit!(r)))
        A = Matrix{Complex{BigFloat}}(3I, n, n) + [gB!(rng) for _ in 1:n, _ in 1:n] / 8
        b = [gB!(rng) for _ in 1:n]
        res = gmres(applier(A), b; abstol = big"1e-60", maxiter = 50)
        @test res.converged
        @test eltype(res.x) == Complex{BigFloat}
        @test precision(real(res.x[1])) == 256
        @test res.residual isa BigFloat
        @test norm(b - A * res.x) ≤ big"1e-60"             # Float64 arithmetic cannot pass
        @test norm(res.x - A \ b) ≤ big"1e-55"             # dense-\ oracle (calibrated 3.7e-77)
        @test monotone(res.history)
    end
end

end # testset gmres
