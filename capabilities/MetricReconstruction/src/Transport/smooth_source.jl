# smooth_source.jl — F4: the generic SMOOTH static-source transport path (the disc-program enabler).
#
# The thin-shell chain (radial_profile/trace_branch/transport_chain) is CLOSED-FORM over the
# distributional (δ,δ′,Θ) shell profiles. This file adds the [FUTURE]-seam sibling declared in
# transport_chain.jl: `transport_solve(::SmoothSource)`, integrating each Fig.-1 transport equation
# as a first-order radial ODE with a SMOOTH RHS via `NumericsCore.integrate_to_targets!` (vendored
# DP5, {T}-generic, BigFloat-clean). It is GATED by convergence to the thin-shell closed forms as the
# smoothing width w → 0 (test/transport_smooth_source.jl, G-F4a…G-F4e).
#
# ── The conserved smooth-shell family (the thin shell's smooth analogue) ──────────────────────────
# Static, spherically symmetric, Schwarzschild background, ZERO radial pressure. Energy density
#   σ(r) = A·exp(−(r−R)²/2w²)         (a Gaussian of width w centred at R),
# with A normalised so ∫_{rin}^{rout} σ dr = δM/(4πR²)  (the value that makes the transport sources
# reproduce the thin-shell δ-coefficients in the w→0 limit ⇒ exterior ADM mass shift = δM; see
# [AXIOM-F4] and G-F4b). The tangential pressure is FIXED BY CONSERVATION ∇_μ T^{μν}=0 on
# Schwarzschild (the only nontrivial component, ν=r):
#
#   ∇_μ T^{μr} = ½ f′ σ − (2f/r) p_t = 0   ⇒   p_t(r) = ¼ r f′(r)/f(r) · σ(r) = M/(2 r f(r)) · σ(r)
#
# (T^{μν} = σ u^μu^ν + p_t Π^{μν}, u^μ = f^{-1/2}∂_t, Π = angular projector; f=1−2M/r, f′=2M/r²).
# Its w→0 limit p_t/σ → M/(2Rf(R)) reproduces the thin-shell EOS ratio (ANCHOR-1). Verified
# numerically against `Spacetime.christoffel` in G-F4e.
#
# ── The transport ODE system (static monopole; state y=[μ, h_mm̄, Ψ2, N], N ≡ f² h_ll) ────────────
# Reduce the Fig.-1 chain (A1,A2,B6,B7; the spin-1/2 sectors vanish, ANCHOR-5) to first order along
# the ingoing n-leg (Δ_static = −(f/2)∂_r + weight·bg):
#   A1  μ′     = (2/f)[Φ22 − (r−3M)/r² · μ]                         (Δ_{[2,0,1,1]}μ = −Φ22, LY:229/774)
#   A2  h_mm̄′ = −4μ/f                                              (Δ_{[-1,1,0,0]}h_mm̄ = 2μ, LY:235/778)
#   B6  Ψ2′    = −(3/r)Ψ2 − 6Mμ/(r³f) − 2R_g/f                     (Δ_{[3,0,0,0]}Ψ2 = −3Ψ2μ+R_g, LY:284/805)
#   B7  N′     = −16μ − 8r(μ′ − Ψ2 − 2Λ)   [N=f²h_ll]             (LY:311/836; f²-integrating factor)
# with smooth NP-Ricci sources Φ22(r)=πσf, Φ11(r)=πσ(r−M)/(rf), Λ(r)=πσ(r−3M)/(3rf) (reusing
# [`np_ricci_from_T`](@ref)) and R_g = −Φ22′ − Φ22/r + f Λ′ + (f/r)Φ11 (the LY:505 reduction, all
# analytic from σ,σ′). [AXIOM-E] BCs at the outer boundary (outside stays Schwarzschild + mass shift):
#   μ(rout)=0, h_mm̄(rout)=0, Ψ2(rout)=−δM/rout³, N(rout)=8δM/rout ; integrate inward (dir=−1).

import NumericsCore: integrate_to_targets!, rk_tol, CubicSplineInterpolant
import Spacetime: christoffel
import QuadGK: quadgk

# ================================================================================================
# SmoothSource — replaces the [FUTURE] stub. Static Schwarzschild smooth-shell parameters.
# ================================================================================================
"""
    SmoothSource(; M, R, w, dM=M, rin=2M+0.05, rout=R+6, ngrid=257)

A generic SMOOTH (non-distributional) static, spherically symmetric O(ε) stress-energy source: a
Gaussian energy-density shell of width `w` centred at areal radius `R` on a Schwarzschild background
of mass `M`, with tangential pressure fixed by conservation and radial pressure zero. `dM` is the
mass shift the profile is normalised to (default `= M`, matching the thin shell's δM; `dM=0` gives
the vacuum source for G-F4a). The transport ODEs are integrated on `[rin, rout]` (DP5), with `ngrid`
spline nodes for the callable solution. This is the F4 realisation of the seam declared in
`transport_chain.jl`; the thin-shell (`ShellSource`) path stays the closed-form gate.
"""
struct SmoothSource
    M::Float64
    R::Float64
    w::Float64
    dM::Float64
    rin::Float64
    rout::Float64
    ngrid::Int
end
function SmoothSource(; M, R, w, dM = M, rin = 2 * M + 0.05, rout = R + 6.0, ngrid = 257)
    SmoothSource(Float64(M), Float64(R), Float64(w), Float64(dM),
                 Float64(rin), Float64(rout), Int(ngrid))
end

# ------------------------------------------------------------------------------------------------
# The conserved profile family (σ, σ′, p_t) — {T}-generic.
# ------------------------------------------------------------------------------------------------
"""
    smooth_norm(src, ::Type{T}) -> A

The Gaussian amplitude `A` such that `∫_{rin}^{rout} A·exp(−(r−R)²/2w²) dr = dM/(4πR²)` (the mass-shift
normalisation; the truncated integral is done by `QuadGK.quadgk` so the enclosed mass is exactly
`dM/(4πR²)` on the integration domain — the exterior is exactly the mass-`dM` solution for any `w`).
"""
function smooth_norm(src::SmoothSource, ::Type{T}) where {T}
    R, w = T(src.R), T(src.w)
    rin, rout = T(src.rin), T(src.rout)
    I, _ = quadgk(r -> exp(-(r - R)^2 / (2 * w^2)), rin, rout; rtol = T <: BigFloat ? T(10)^(-30) : 1e-13)
    (T(src.dM) / (4 * T(pi) * R^2)) / I
end

"σ(r) = A·exp(−(r−R)²/2w²), the smooth energy density ({T}-generic; `A` = [`smooth_norm`](@ref))."
smooth_sigma(src::SmoothSource, r, A) = A * exp(-(r - oftype(A, src.R))^2 / (2 * oftype(A, src.w)^2))
"σ′(r) = −(r−R)/w² · σ(r)."
smooth_sigma_p(src::SmoothSource, r, A) =
    -(r - oftype(A, src.R)) / oftype(A, src.w)^2 * smooth_sigma(src, r, A)

"""
    smooth_pt(src, r, A) -> p_t

The conservation-fixed tangential pressure `p_t(r) = M/(2 r f(r)) · σ(r)` (radial pressure zero).
"""
function smooth_pt(src::SmoothSource, r, A)
    M = oftype(A, src.M)
    f = 1 - 2M / r
    M / (2 * r * f) * smooth_sigma(src, r, A)
end

# ------------------------------------------------------------------------------------------------
# NP-Ricci sources (reuse np_ricci_from_T) — smooth pointwise values + analytic derivatives.
# ------------------------------------------------------------------------------------------------
"""
    smooth_ricci(src, r, A) -> (Φ00, Φ22, Φ11, Λ)

The smooth NP-Ricci scalars at `r`, built from the fluid NP stress components
`Tll=σ/f, Tnn=σf/4, Tln=σ/2, Tmm̄=p_t` through [`np_ricci_from_T`](@ref) (the SAME projection the
shell uses). Gives `Φ22=πσf`, `Φ11=πσ(r−M)/(rf)`, `Λ=πσ(r−3M)/(3rf)`.
"""
function smooth_ricci(src::SmoothSource, r, A)
    M = oftype(A, src.M)
    f = 1 - 2M / r
    σ = smooth_sigma(src, r, A)
    pt = smooth_pt(src, r, A)
    np_ricci_from_T(σ / f, σ * f / 4, σ / 2, pt)
end

# ------------------------------------------------------------------------------------------------
# The DP5 RHS builder (state y = [μ, h_mm̄, Ψ2, N]) and the BC / integration driver.
# ------------------------------------------------------------------------------------------------
"""
    smooth_rhs(src, A, ::Type{T}) -> f!

Build the in-place first-order RHS `f!(g, r, y)` of the reduced transport system (A1,A2,B6,B7). All
sources analytic in `σ,σ′` (no finite differences ⇒ BigFloat-clean). `y=[μ, h_mm̄, Ψ2, N=f²h_ll]`.
"""
function smooth_rhs(src::SmoothSource, A::T, ::Type{T}) where {T}
    M = T(src.M)
    πc = T(pi)
    function f!(g, r, y)
        μ, hmm, Ψ2, N = y
        f = 1 - 2M / r
        fp = 2M / r^2
        σ = smooth_sigma(src, r, A)
        σp = smooth_sigma_p(src, r, A)
        pt = M / (2 * r * f) * σ
        # NP-Ricci + analytic derivatives (rf = r−2M)
        Φ22 = πc * σ * f
        Φ22p = πc * (σp * f + σ * fp)
        Φ11 = πc * σ * (r - M) / (r * f)
        rm2 = r - 2M
        h1 = (r - 3M) / (3 * rm2)                 # Λ = πσ·h1
        h1p = M / (3 * rm2^2)
        Λ = πc * σ * h1
        Λp = πc * (σp * h1 + σ * h1p)
        Rg = -Φ22p - Φ22 / r + f * Λp + (f / r) * Φ11
        μp = (2 / f) * (Φ22 - (r - 3M) / r^2 * μ)
        g[1] = μp
        g[2] = -4μ / f
        g[3] = -(3 / r) * Ψ2 - 6M * μ / (r^3 * f) - 2 * Rg / f
        g[4] = -16μ - 8 * r * (μp - Ψ2 - 2Λ)
        return nothing
    end
end

"[AXIOM-E] outer-boundary state [μ,h_mm̄,Ψ2,N] at r=rout (outside stays Schwarzschild + mass shift dM)."
function smooth_bc(src::SmoothSource, ::Type{T}) where {T}
    dM = T(src.dM)
    rout = T(src.rout)
    T[zero(T), zero(T), -dM / rout^3, 8dM / rout]
end

"""
    smooth_solve_targets(src, radii, ::Type{T}=Float64) -> Matrix{T}

Integrate the transport system inward from `rout` (BCs [`smooth_bc`](@ref)) to the given `radii`,
returning a `length(radii)×4` matrix whose row `i` is `[μ, h_mm̄, Ψ2, N]` at `radii[i]` (rows aligned
to the INPUT order of `radii`). Direct DP5 landing on the exact radii — machine-accurate (used by the
gates that need precision; the callable solution uses a spline, see [`transport_solve`](@ref)).
"""
function smooth_solve_targets(src::SmoothSource, radii, ::Type{T} = Float64; tol = rk_tol(T)) where {T}
    A = smooth_norm(src, T)
    f! = smooth_rhs(src, A, T)
    y0 = smooth_bc(src, T)
    rout = T(src.rout)
    perm = sortperm(collect(radii); rev = true)      # descending for dir=−1
    tgt = T.(collect(radii)[perm])
    Y = integrate_to_targets!(f!, y0, rout, tgt, -1; tol = T(tol))
    out = Matrix{T}(undef, length(radii), 4)
    for (k, i) in enumerate(perm)
        out[i, :] .= Y[k, :]
    end
    out
end

# ================================================================================================
# The callable solution + the Fig.-1 DAG driver `transport_solve(::SmoothSource)`.
# ================================================================================================
"""
    SmoothSolution{T}

The DP5 smooth-source reconstruction, as CALLABLES over `[rin,rout]` (cubic-spline interpolants of
the DP5 grid) — the smooth analogue of the [`ShellSource`](@ref) solution's `ShellProfile` fields.
Fields are functions of `r`: `μ, h_mm̄, Ψ2, h_ll, γ, N` plus the sources `Φ22, Λ`. Also carries the
grid `rs` and the underlying `spline`.
"""
struct SmoothSolution{T,S}
    src::SmoothSource
    rs::Vector{T}
    spline::S           # CubicSplineInterpolant with rows [μ, h_mm̄, Ψ2, h_ll]
    A::T
end

# field-index in the spline
const _SMOOTH_ROW = (μ = 1, h_mm̄ = 2, Ψ2 = 3, h_ll = 4)

_smooth_eval(sol::SmoothSolution, row::Int, r) = sol.spline(r, row)

"μ^{(1)}(r)."                       smooth_mu(sol::SmoothSolution, r)   = _smooth_eval(sol, 1, r)
"h_mm̄^{(1)}(r) (= the trace)."      smooth_hmm(sol::SmoothSolution, r)  = _smooth_eval(sol, 2, r)
"Ψ2^{(1)}(r)."                      smooth_psi2(sol::SmoothSolution, r) = _smooth_eval(sol, 3, r)
"h_ll^{(1)}(r)."                    smooth_hll(sol::SmoothSolution, r)  = _smooth_eval(sol, 4, r)
"γ^{(1)}(r) = 0 (real spin-0 trace, LY:228 — as for the shell)."
smooth_gamma(::SmoothSolution{T}, r) where {T} = zero(T)

# traceful-ORG metric components (LY:893–895), for the G-F4c end-to-end gate. h_tt=(f²/4)h_ll,
# h_rr=h_ll/4, h_φφ/sin²θ = r² h_mm̄ (compare to exact_htt/exact_hrr/exact_hφφ_over_sin2).
"h_tt^{(1)}(r) = (f²/4) h_ll (LY:893)."
smooth_htt(sol::SmoothSolution, r) = (1 - 2 * oftype(smooth_hll(sol, r), sol.src.M) / r)^2 / 4 * smooth_hll(sol, r)
"h_rr^{(1)}(r) = h_ll/4 (LY:893)."
smooth_hrr(sol::SmoothSolution, r) = smooth_hll(sol, r) / 4
"h_φφ^{(1)}(r)/sin²θ = r² h_mm̄(r) (LY:895)."
smooth_hφφ_over_sin2(sol::SmoothSolution, r) = r^2 * smooth_hmm(sol, r)

"""
    transport_solve(src::SmoothSource, ::Type{T}=Float64) -> SmoothSolution

Run the Fig.-1 transport DAG for the SMOOTH static Schwarzschild shell via DP5 (F4). Integrates
`[μ,h_mm̄,Ψ2,N=f²h_ll]` inward on an `ngrid` grid, returns a [`SmoothSolution`](@ref) whose fields are
spline callables (mirrors the `ShellSource` NamedTuple shape, callables instead of `ShellProfile`s).
The spin-1/2 sectors are zero by spherical symmetry (ANCHOR-5); `γ^{(1)}=0` (real trace). Gated by
w→0 convergence to the thin-shell closed forms (test/transport_smooth_source.jl). `tol` is the DP5
step tolerance (default the precision-inherited `rk_tol(T)`; for `BigFloat` pass an explicit modest
value — the default `2^(20−bits)` is impractically tight for a smoke run).
"""
function transport_solve(src::SmoothSource, ::Type{T} = Float64; tol = rk_tol(T)) where {T}
    rin, rout = T(src.rin), T(src.rout)
    rs = collect(range(rin, rout; length = src.ngrid))     # ascending grid
    Y = smooth_solve_targets(src, rs, T; tol = tol)         # rows aligned to rs (ascending)
    M = T(src.M)
    # build the spline rows [μ, h_mm̄, Ψ2, h_ll]  (convert N→h_ll = N/f²)
    rows = Matrix{T}(undef, 4, length(rs))
    for j in eachindex(rs)
        f = 1 - 2M / rs[j]
        rows[1, j] = Y[j, 1]           # μ
        rows[2, j] = Y[j, 2]           # h_mm̄
        rows[3, j] = Y[j, 3]           # Ψ2
        rows[4, j] = Y[j, 4] / f^2     # h_ll = N/f²
    end
    spl = CubicSplineInterpolant(rs, rows)
    A = smooth_norm(src, T)
    SmoothSolution{T,typeof(spl)}(src, rs, spl, A)
end

# ================================================================================================
# Conservation residual (G-F4e) — ∇_μ T^{μν} via Spacetime.christoffel, analytic ∂T.
# ================================================================================================
"""
    conservation_residual(src, r, θ, ::Type{T}=Float64) -> value

`max_ν |∇_μ T^{μν}|` at `(r,θ)` for the conserved smooth-shell stress-energy, computed with
`Spacetime.christoffel(r,θ;M)` and ANALYTIC ∂_r,∂_θ of `T^{μν}` (no finite differences). With the
conservation-fixed `p_t = M/(2rf)σ` this is algebraically zero (only the ν=r component is nontrivial:
½f′σ − 2f p_t/r ≡ 0), so the residual is at roundoff. The G-F4e gate asserts `≤ 1e-12` over a grid.
"""
function conservation_residual(src::SmoothSource, r, θ, ::Type{T} = Float64) where {T}
    A = smooth_norm(src, T)
    M = T(src.M)
    rt, θt = T(r), T(θ)
    f = 1 - 2M / rt
    fp = 2M / rt^2
    s = sin(θt)
    σ = smooth_sigma(src, rt, A)
    σp = smooth_sigma_p(src, rt, A)
    pt = M / (2 * rt * f) * σ
    # p_t′  (rf = r−2M, (rf)′=1):  p_t = (M/2) σ/(r−2M)
    ptp = (M / 2) * (σp / (rt - 2M) - σ / (rt - 2M)^2)
    # T^{μν} (contravariant), diagonal
    Ttt = σ / f
    Tθθ = pt / rt^2
    Tφφ = pt / (rt^2 * s^2)
    Tup = zeros(T, 4, 4)
    Tup[1, 1] = Ttt; Tup[3, 3] = Tθθ; Tup[4, 4] = Tφφ
    # ∂_r T^{μν}
    dTtt = σp / f - σ * fp / f^2
    dTθθ = ptp / rt^2 - 2pt / rt^3
    dTφφ = (ptp / rt^2 - 2pt / rt^3) / s^2
    dTr = zeros(T, 4, 4)
    dTr[1, 1] = dTtt; dTr[3, 3] = dTθθ; dTr[4, 4] = dTφφ
    # ∂_θ T^{μν}: only T^{φφ} depends on θ (1/sin²θ)
    dTθ = zeros(T, 4, 4)
    dTθ[4, 4] = pt / rt^2 * (-2 * cos(θt) / s^3)
    Γ = christoffel(rt, θt; M = M)                # Γ[λ,μ,ν] = Γ^λ_{μν}
    res = zero(T)
    for ν in 1:4
        val = dTr[2, ν] + dTθ[3, ν]              # ∂_μ T^{μν} (μ=r,θ; T^{rν}=0)
        for μ in 1:4, λ in 1:4
            val += Γ[μ, μ, λ] * Tup[λ, ν]        # Γ^μ_{μλ} T^{λν}
        end
        for μ in 1:4, λ in 1:4
            val += Γ[ν, μ, λ] * Tup[μ, λ]        # Γ^ν_{μλ} T^{μλ}
        end
        res = max(res, abs(val))
    end
    res
end
