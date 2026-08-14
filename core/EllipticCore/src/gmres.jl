# gmres.jl — parametric-T restarted GMRES: Arnoldi/MGS + incremental complex-Givens QR of
# the Hessenberg, matrix-free operator, right-preconditioning hook.
# DAG node: complex-gmres (elliptic-ssf GLOBAL_DAG row 16, TOPOLOGICAL_SCHEDULE W3-B).
#
# ORIGIN: adapted from SpecCore (the user's own specjl repo),
# specjl/packages/SpecCore/src/elliptic/solvers/gmres.jl (dependency-free full-memory GMRES,
# itself verified against SpECTRE Gmres.hpp + ResidualMonitorActions.hpp) and the
# right-preconditioned variant pgmres_solve in solvers/schwarz.jl:535-622.
# Deviations from the SpecCore reference, per the elliptic-ssf DESIGN.md gmres.jl row:
#   * parametric in the scalar type T <: Number — Complex{Float64} first-class,
#     Complex{BigFloat} capable (RT = real(T) tolerances/history); the reference was
#     Vector{Float64}-hardcoded (audit-speccore §2).
#   * real Givens (hypot form) → the standard complex-safe zlartg formulation via abs/conj
#     (cosine real, sine complex), written from the mathematics in `givens_rotation` and
#     gated against LinearAlgebra.givens on complex inputs in test_gmres.jl.
#   * the SPD monotone-residual HARD ERROR (gmres.jl:78-81) REMOVED — it encodes an SPD
#     assumption the complex-shifted non-Hermitian m-mode operator does not satisfy
#     (audit-speccore risk 3). Replaced by diagnostics valid for any operator: happy-breakdown
#     detection, restart-cycle stagnation detection, and residual-history reporting.
#   * restart cycles added (SpecCore was full-memory only); restart = maxiter (the default)
#     reproduces the unrestarted algorithm exactly.
# Defaults abstol = 1e-12, reltol = 0, maxiter = 500 match the CRP reference yaml
# (AbsoluteResidual 1e-12, RelativeResidual 0, MaxIterations 500 [SOLID per decomp-crp §2]).
#
# Complex discipline (EllipticCore module docstring): inner products conjugate on the left
# (LinearAlgebra.dot); the Givens cosine is real, the sine complex; the rotation is
# G = [c s; -conj(s) c], so row updates read (h_i, h_{i+1}) ← (c·h_i + s·h_{i+1},
# -conj(s)·h_i + c·h_{i+1}).

export gmres, GMRESResult

"""
    GMRESResult{T, RT}

Outcome of [`gmres`](@ref). Fields:

  - `x::Vector{T}` — the (approximate) solution.
  - `iterations::Int` — total inner (Arnoldi) iterations across all restart cycles.
  - `residual::RT` — the TRUE residual norm `‖b − A x‖`, recomputed explicitly from the
    operator (not the Givens estimate), guarding against roundoff drift.
  - `converged::Bool` — whether `residual ≤ max(abstol, reltol·‖b‖)`.
  - `breakdown::Symbol` — `:none`, `:happy` (Arnoldi lucky breakdown — the Krylov space
    closed, so the returned iterate is exact up to roundoff), or `:stagnation` (a restart
    cycle produced no meaningful residual reduction; the honest non-Hermitian replacement
    for the SPD monotone-residual hard error of the SpecCore reference).
  - `history::Vector{RT}` — residual-norm history: `history[1] = ‖b − A x₀‖`, then one
    Givens estimate per inner iteration (`length == iterations + 1`). With RIGHT
    preconditioning the estimate tracks the true residual of the original system, so the
    history is monotone nonincreasing up to roundoff (GMRES minimization property).
"""
struct GMRESResult{T<:Number,RT<:Real}
    x::Vector{T}
    iterations::Int
    residual::RT
    converged::Bool
    breakdown::Symbol
    history::Vector{RT}
end

"""
    givens_rotation(f, g) -> (c, s, r)

Complex-safe Givens rotation (LAPACK `zlartg` convention): returns a REAL cosine `c ≥ 0`,
a complex sine `s`, and `r` such that

    [ c        s ] [f]   [r]
    [ -conj(s) c ] [g] = [0],      c² + |s|² = 1.

For `f ≠ 0`: with `d = hypot(|f|, |g|)` and the phase `ph = f/|f|`,
`c = |f|/d`, `s = ph·conj(g)/d`, `r = ph·d`. Edge cases: `g = 0 ⇒ (1, 0, f)`;
`f = 0, g ≠ 0 ⇒ (0, conj(g)/|g|, |g|)`. Works for any `T <: Number` (real types reduce to
a real rotation with nonnegative cosine). Unexported; the GMRES kernel and tests use it.
Gated against `LinearAlgebra.givens` on complex inputs (identical convention); for REAL
inputs Julia's `givens` uses a different sign convention, so only the defining properties
above are shared.
"""
function givens_rotation(f::T, g::T) where {T<:Number}
    RT = real(T)
    if iszero(g)
        return one(RT), zero(T), f
    elseif iszero(f)
        ag = abs(g)
        return zero(RT), conj(g) / ag, T(ag)
    else
        af = abs(f)
        d = hypot(af, abs(g))
        ph = f / af
        return af / d, ph * conj(g) / d, ph * T(d)
    end
end

"""
    gmres(applyA!, b; abstol = 1e-12, reltol = 0, maxiter = 500, restart = maxiter,
          Minv! = nothing, x0 = nothing) -> GMRESResult

Solve `A x = b` for the matrix-free operator `applyA!(out, v)` (any callable writing `A·v`
into `out`) with restarted GMRES: Arnoldi with modified Gram-Schmidt and an incremental
complex-Givens QR of the Hessenberg. Parametric in the scalar type of `b` (`ComplexF64`
first-class; `Complex{BigFloat}` supported, tolerances in `RT = real(T)`).

Convergence when the recomputed true residual satisfies `‖b − A x‖ ≤ max(abstol,
reltol·‖b‖)` — defaults match the CRP reference yaml (absolute 1e-12, relative 0, 500
iterations). `restart` is the Arnoldi-basis cap per cycle; the default (`= maxiter`)
reproduces unrestarted full-memory GMRES exactly. `maxiter` caps TOTAL inner iterations
across cycles. `x0` warm-starts the iteration (default zero).

`Minv!(z, v)` is the RIGHT-preconditioning hook (identity when `nothing`): the Krylov
recursion runs on `A·M⁻¹` while the residual — and therefore the convergence test, the
reported `history`, and the minimization property — remains that of the ORIGINAL system
(the CRP's preconditioned-GMRES usage, main.tex:670-674; multigrid-V-cycle `Minv!` arrives
with the solver-stack node).

No SPD monotone-residual hard error (the SpecCore reference's falsifier assumes an SPD
operator the complex-shifted m-mode problem does not provide): non-convergence surfaces as
`converged = false` with `breakdown ∈ {:none, :stagnation}` plus the full residual
history — see [`GMRESResult`](@ref). Assumes `A` (and `M`) nonsingular.
"""
function gmres(applyA!, b::AbstractVector{T};
               abstol::Real = 1e-12, reltol::Real = 0,
               maxiter::Integer = 500, restart::Integer = maxiter,
               Minv! = nothing,
               x0::Union{Nothing,AbstractVector} = nothing) where {T<:Number}
    Tf = float(T)
    RT = real(Tf)
    n = length(b)
    maxiter ≥ 0 || throw(ArgumentError("maxiter must be ≥ 0, got $maxiter"))
    restart ≥ 1 || throw(ArgumentError("restart must be ≥ 1, got $restart"))
    x0 === nothing || length(x0) == n ||
        throw(ArgumentError("x0 has length $(length(x0)), expected $n"))

    x = x0 === nothing ? zeros(Tf, n) : Vector{Tf}(x0)
    bnorm = norm(b)
    tol = max(RT(abstol), RT(reltol) * bnorm)

    # initial residual (r0 = b for the default zero start)
    w = Vector{Tf}(undef, n)                    # A·z scratch
    r = Vector{Tf}(undef, n)
    if x0 === nothing
        copyto!(r, b)
    else
        applyA!(w, x)
        @. r = b - w
    end
    β = norm(r)
    history = [β]
    β ≤ tol && return GMRESResult{Tf,RT}(x, 0, β, true, :none, history)

    m_max = min(Int(restart), Int(maxiter))
    H = zeros(Tf, m_max + 1, m_max)             # Hessenberg (upper-triangularized in place)
    cs = zeros(RT, m_max)                       # Givens cosines (real by construction)
    sn = zeros(Tf, m_max)                       # Givens sines (complex)
    g = zeros(Tf, m_max + 1)                    # rotated RHS (‖residual‖ estimate = |g[k+1]|)
    stag_rtol = sqrt(eps(RT))                   # "no meaningful reduction" threshold
    iters = 0
    breakdown = :none
    converged = false
    res_true = β

    while true                                  # restart cycles
        cycle_start = res_true                  # == β == ‖b − A x‖ entering this cycle
        V = Vector{Vector{Tf}}()                # Arnoldi basis of K(A·M⁻¹, r)
        push!(V, r ./ β)
        Z = Vector{Vector{Tf}}()                # preconditioned directions z_j = M⁻¹ v_j
        fill!(g, zero(Tf))
        g[1] = β
        breaktol = 100 * eps(RT) * β            # happy-breakdown scale (1e-14·β at Float64)
        k = 0
        lucky = false
        m = min(m_max, Int(maxiter) - iters)
        for j in 1:m
            if Minv! === nothing
                applyA!(w, V[j])
            else
                z = Vector{Tf}(undef, n)
                Minv!(z, V[j])
                push!(Z, z)
                applyA!(w, z)                   # w = A M⁻¹ v_j
            end
            for i in 1:j                        # modified Gram-Schmidt (dot conjugates left)
                hij = dot(V[i], w)
                H[i, j] = hij
                @. w = w - hij * V[i]
            end
            hw = norm(w)
            H[j + 1, j] = hw
            lucky = hw ≤ breaktol
            lucky || push!(V, w ./ hw)
            for i in 1:(j - 1)                  # accumulated rotations on column j
                t = cs[i] * H[i, j] + sn[i] * H[i + 1, j]
                H[i + 1, j] = -conj(sn[i]) * H[i, j] + cs[i] * H[i + 1, j]
                H[i, j] = t
            end
            c, s, ρ = givens_rotation(H[j, j], H[j + 1, j])
            cs[j] = c
            sn[j] = s
            H[j, j] = ρ
            H[j + 1, j] = zero(Tf)
            g[j + 1] = -conj(s) * g[j]
            g[j] = c * g[j]
            k = j
            iters += 1
            push!(history, abs(g[j + 1]))
            (abs(g[j + 1]) ≤ tol || lucky) && break
        end

        # back-substitution R y = g (R = H[1:k,1:k] upper triangular), then x += Σ yᵢ zᵢ
        y = Vector{Tf}(undef, k)
        for i in k:-1:1
            acc = g[i]
            for l in (i + 1):k
                acc -= H[i, l] * y[l]
            end
            y[i] = acc / H[i, i]
        end
        dirs = Minv! === nothing ? V : Z
        for i in 1:k
            @. x = x + y[i] * dirs[i]
        end

        applyA!(w, x)                           # honest residual (recomputed, not estimated)
        @. r = b - w
        res_true = norm(r)
        converged = res_true ≤ tol
        lucky && (breakdown = :happy)
        (converged || lucky || iters ≥ maxiter) && break
        if res_true ≥ cycle_start * (1 - stag_rtol)
            breakdown = :stagnation             # a whole cycle bought nothing — stop honestly
            break
        end
        β = res_true                            # restart from the true residual
    end

    return GMRESResult{Tf,RT}(x, iters, res_true, converged, breakdown, history)
end
