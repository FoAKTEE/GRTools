# trace_branch.jl — the spin-0 trace head: h_mm̄^{(1)} directly from T_nn^{(1)}.
#   A1  Δ_{[2,0,1,1]} μ^{(1)} = −Φ22^{(1)}      (eq:mu_transport_final, LY:229)
#   A2  Δ_{[-1,1,0,0]} h_mm̄^{(1)} = 2 μ^{(1)}   (eq:trace_transport,   LY:234)
# Shell specialization (eq:hmm_eq_shell, LY:783): ∂_r(r²∂_r h_mm̄) = −(8r²/f²)Φ22 = −(2M/f(R))δ(r−R),
# basis {1, 1/r}, [AXIOM-E] outside branch c0⁺=c1⁺=0, δ jump [r²∂_r h_mm̄]=−2M/f(R), continuity.

"""
    solve_hmm_shell(T, M, R, Φ22; c0p=0, c1p=0) -> ShellProfile

Solve the trace-branch Euler equation for `h_mm̄^{(1)}` on the static Schwarzschild shell. `Φ22` is
the shell source [`ShellProfile`](@ref) (pure δ). The homogeneous solution is `h = c0 − c1/r` on each
side; the exterior amplitudes `(c0p, c1p)` are the **[AXIOM-E]** outside-branch input (paper default
`c0⁺=c1⁺=0`, "outside stays Schwarzschild", LY:793). The δ-source sets the jump
`[r²∂_r h_mm̄] = k` with `k = −(8r²/f²)|_R · Φ22_δ = −2M/f(R)`, and continuity fixes the interior
amplitudes. Reproduces ANCHOR-8 (LY:799–802).
"""
function solve_hmm_shell(::Type{T}, M, R, Φ22::ShellProfile{T}; c0p = zero(T), c1p = zero(T)) where {T}
    Mt, Rt = T(M), T(R)
    fR = 1 - 2Mt / Rt
    # decoupled RHS δ-coefficient  k = −(8r²/f²)|_R · Φ22_δ   (eq:hmm_eq_shell)
    k = -(8Rt^2 / fR^2) * Φ22.dc
    # r²∂_r(c0 − c1/r) = c1 ; jump [r²∂_r h] = c1⁺ − c1⁻ = k  ⇒  c1⁻ = c1⁺ − k
    c1m = c1p - k
    # continuity h⁺(R)=h⁻(R): c0⁺ − c1⁺/R = c0⁻ − c1⁻/R  ⇒  c0⁻ = c0⁺ − c1⁺/R + c1⁻/R
    c0m = c0p - c1p / Rt + c1m / Rt
    inside = laurent(T, 0 => c0m, -1 => -c1m)
    outside = laurent(T, 0 => c0p, -1 => -c1p)
    shellprofile(inside, outside, Rt)
end

"""
    mu_from_trace(T, M, h_mm̄) -> ShellProfile

Reconstruct `μ^{(1)}` from the trace `h_mm̄^{(1)}` via `eq:reconstruct_mu_gamma_ORG` (LY:218):
`μ^{(1)} = ½ Δ_{[-1,1,0,0]} h_mm̄^{(1)}`. On the static shell `Δ_{[-1,1,0,0]} = −(f/2)∂_r` (the
background `−μ+μ̄=0`), so `μ^{(1)} = −(f/4) ∂_r h_mm̄^{(1)}` (LY:778). Reproduces ANCHOR-9 (LY:816).
"""
function mu_from_trace(::Type{T}, M, h_mm̄::ShellProfile{T}) where {T}
    f = schw_f(T, M)
    pscale(-one(T) / 4, psmooth(f, pderiv(h_mm̄)))
end

"""
    gamma_from_trace(T, M, h_mm̄) -> ShellProfile

Reconstruct `γ^{(1)}` from `h_mm̄^{(1)}` via `eq:reconstruct_gamma` (LY:220):
`γ^{(1)} = −¼(μ−μ̄) h_mm̄^{(1)}`. In the static real-μ Schwarzschild background `μ−μ̄=0`, so
`γ^{(1)}=0` (γ^{(1)} is imaginary, LY:228; the shell's real spin-0 trace carries none).
"""
function gamma_from_trace(::Type{T}, M, h_mm̄::ShellProfile{T}) where {T}
    ShellProfile{T}(lzero(T), lzero(T), zero(T), zero(T), h_mm̄.R)
end
