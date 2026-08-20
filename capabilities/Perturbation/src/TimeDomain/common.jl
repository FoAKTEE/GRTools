# Shared TimeDomain helpers. In the reference repo these were duplicated per-module
# (Mode1D/Flux1D/KerrModes/KerrFlux each carried private copies); flattening the physics into
# ONE submodule collapses each to a single definition (Phase 7-D DRY discipline). Every body is
# the repo's, ported verbatim (deviations noted inline).

# ---- interpolation probe row -------------------------------------------------------------------
"""barycentric Lagrange row `L_j(ξ)` on the reference element (repo `_lagrange_row`, one copy for
Mode1D/Flux1D/KerrModes/KerrFlux — identical in all four). The kernel is the single
`NumericsCore.bary_row`; the probe `ξ` is a COMPUTED coordinate, so the node-collision radius is the
repo's `10 eps(T)` (`nodetol`) rather than exact equality. Weights are rebuilt per call, as in the
repo — `bary_weights` runs the same `w[j] /= (r[j]-r[k])` loop in the same order, so the row is
bit-identical."""
_lagrange_row(el::RefElement{T}, ξ) where {T} =
    bary_row(el.r, bary_weights(el.r), T(ξ); nodetol = eps(T) * 10)

# ---- single-mode 2×2 flux action ---------------------------------------------------------------
"`A·(π,φ)` for the single-mode flux matrix `A = 1/(1+H)[[H,1],[1,H]]` (repo `_AU`, one copy for
Mode1D/Flux1D)."
@inline _AU(H, π, φ) = (iH = 1 / (1 + H); (iH * (H * π + φ), iH * (π + H * φ)))

# ---- smooth source ramp ------------------------------------------------------------------------
# erf via tanh-based Abramowitz-Stegun approximation; |err|<2e-3, fine for a ramp. KEPT VERBATIM
# (not SpecialFunctions.erf) so the ramp — and therefore every transient — is bit-compatible with
# the reference repo, which deliberately avoided the SpecialFunctions dep.
@inline function erf_approx(x)
    t = 1 / (1 + 0.3275911 * abs(x))
    y = 1 - (((((1.061405429t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * exp(-x^2)
    return sign(x) * y
end

"smooth erf ramp (Paper 2 App A): 0 for t≤0, →1 over [0,τ] (repo `_ramp`, one copy for
Flux1D/KerrFlux — the two repo bodies are the same formula, one with erf inlined)."
@inline function _ramp(t, τ, δ)
    t <= 0 && return 0.0
    t >= τ && return 1.0
    return 0.5 * (erf_approx(sqrt(δ) * (t - τ / 2)) + 1)
end

# ---- coupled-sector block matrices (KerrModes/KerrFlux) -----------------------------------------
"symmetric-tridiagonal `E_ππ` block: `diag 1-f c^ℓ_{ℓm}`, `offdiag -f c^ℓ_{ℓ±2,m}` (Paper 1 §II.E;
repo `_Eππ`, one copy for KerrModes/KerrFlux — identical in both; type-generic in `f`)."
function _Eππ(ℓs, m, f, n)
    E = zeros(typeof(float(f)), n, n)
    for i in 1:n
        Cmm, C0, Cpp = mode_coupling(ℓs[i], m)
        E[i, i] = 1 - f * C0
        i < n && (E[i, i+1] = -f * Cpp)
        i > 1 && (E[i, i-1] = -f * mode_coupling(ℓs[i], m)[1])  # C_{--}(ℓ_i) = offdiag to ℓ_i-2
    end
    return E
end

"""
coupled `2n×2n` flux matrix `A=(1-H)[[H S⁻¹, S⁻¹],[I+H²S⁻¹, H S⁻¹]]`, `S=E_ππ-H²I`
(`A=½[[I,I],[I,I]]` at 𝓘⁺ where `|1-H|<1e-14`); returns `(A, S⁻¹)` with `S⁻¹=0` at scri
(the lower-order term vanishes there). Repo `KerrModes._Amatrix` / `KerrFlux._Amat` unified:
both compute the same matrix; this keeps `_Amat`'s scalar-first `((1-H)·H)` product order and
always returns the tuple (type-generic in `eltype(Eππ)`).
"""
function _Amat(Eππ::AbstractMatrix{T}, H, n) where {T}
    if abs(1 - H) < 1e-14            # 𝓘⁺: A = ½[[I,I],[I,I]]
        A = zeros(T, 2n, 2n)
        for i in 1:n
            A[i, i] = T(0.5); A[i, n+i] = T(0.5); A[n+i, i] = T(0.5); A[n+i, n+i] = T(0.5)
        end
        return A, zeros(T, n, n)
    end
    S = Eππ - H^2 * Matrix(one(T) * Id, n, n)
    Si = inv(S)
    A = zeros(T, 2n, 2n)
    A[1:n, 1:n] .= (1 - H) * H .* Si
    A[1:n, n+1:2n] .= (1 - H) .* Si
    A[n+1:2n, 1:n] .= (1 - H) .* (Matrix(one(T) * Id, n, n) + H^2 .* Si)
    A[n+1:2n, n+1:2n] .= (1 - H) * H .* Si
    return A, Si
end

"block matrix-vector flux `F = A·[π;φ]` into `(Fπ[n], Fφ[n])` (repo `_AUvec!`, one copy for
KerrModes/KerrFlux — element-type generic, identical accumulation order)."
@inline function _AUvec!(Fπ, Fφ, A, π, φ, n)
    @inbounds for r in 1:n
        s1 = zero(eltype(π)); s2 = zero(eltype(π))
        for c in 1:n
            s1 += A[r, c] * π[c] + A[r, n+c] * φ[c]
            s2 += A[n+r, c] * π[c] + A[n+r, n+c] * φ[c]
        end
        Fπ[r] = s1; Fφ[r] = s2
    end
    return nothing
end

# ---- late-time power-law tail fit --------------------------------------------------------------
"""
    fit_tail_window(τs, ψ, win) -> n

Local power-law exponent of `|ψ(τ)| ~ τⁿ` via least squares of `log|ψ|` vs `log τ` over the time
window `win = (τlo, τhi)` (late-time tail; Paper 1 §IV.B). Returns the fitted `n` (`NaN` when
fewer than 3 usable samples). Repo `KerrModes.fit_tail_window` ≡ `Mode1D.fit_tail` body.
"""
function fit_tail_window(τs, ψ, win)
    x = Float64[]; y = Float64[]
    for (t, v) in zip(τs, ψ)
        win[1] <= t <= win[2] && abs(v) > 0 && (push!(x, log(t)); push!(y, log(abs(v))))
    end
    length(x) < 3 && return NaN
    nn = length(x); sx = sum(x); sy = sum(y)
    return (nn * sum(x .* y) - sx * sy) / (nn * sum(abs2, x) - sx^2)
end

"""
    fit_tail(τs, ψ; window=(2000.0, 1e9)) -> n

Keyword-window form of [`fit_tail_window`](@ref) (repo `Mode1D.fit_tail`; the two repo bodies
are identical — collapsed to one core).
"""
fit_tail(τs, ψ; window=(2000.0, 1e9)) = fit_tail_window(τs, ψ, window)
