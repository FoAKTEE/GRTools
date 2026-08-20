# SpiralWaves/grid.jl — the adaptive radial grid of arXiv:2608.17003 Appendix C (appendix.tex:53-128):
# three blocks — a left WKB grid, a sinh-clustered corotation grid (C4)-(C5) and a right WKB grid —
# on a per-mode domain selected by the fiducial edges (C1) r_in = r_ISCO, r_out = 5a• refined by the
# WKB-trust heuristic (C2) k r > 10/h0 and the phase-extent criterion (C3) ∫k dr ≥ 30. The WKB
# blocks are generated from the trial spacing (C7)-(C8) Δr_trial = smoothMin(Δr_WKB, Δr_in),
# Δr_WKB = 1/(k̄ n_W), smoothMin (C9)-(C10), through a cubic-spline interpolant (NumericsCore) of
# Δr_trial on a geometric ancillary grid (A:78, ~10⁴ points), with n_W ∈ [n_L, n_L+1] fixed so that
# the block lands EXACTLY on its two edges (A:117). Type-generic {T}; BigFloat inherits.
#
# Implementation choices where the paper text is under-specified (decomposition derivation_part2.md
# App. C [OPEN] items) — each recorded here and in the docstrings, all [PRELIMINARY]:
#   (k)   "local WKB wavenumber": the Newtonian local dispersion relation ω_d² = κ² + c_s² k², i.e.
#         k² = (ω_d² − κ²)/c_s² with the geodesic κ² (`kappa2`); k ≡ 0 in the evanescent zone k² < 0.
#   (C6)  k̄ = max[k, k(r = 10 h0⁻¹ k⁻¹)] is read as the r-local floor k̄(r) = max(k(r), k_trust/(h0 r)):
#         k̄ = k wherever the WKB trust criterion (C2) holds and the trust-threshold wavenumber where
#         it does not (it equals k(r_WKB) at the trust boundary, continuous).
#   (Δr_in) the extrapolated central spacing is linear in the distance from the central edge:
#         Δr_in(r) = Δr_C,edge + (dΔr_C/dr)_edge · |r − r_edge|, with the continuous sinh-grid spacing
#         Δr_C(x) = σ s_c cosh(s_c x)/sinh(s_c) · Δx, Δx = 2/(n_C−1), so Δr_C,edge = σ s_c coth(s_c) Δx
#         and (dΔr_C/dr)_edge = s_c tanh(s_c) Δx.
#   (nodes) from a spacing function to nodes: nodes sit at INTEGER values of the continuous cell-count
#         coordinate Φ(r) = ∫ dr'/S(r') measured from the central edge (S = the cubic spline of
#         Δr_trial); N = ⌈Φ_total(n_L)⌉ cells and n_W* ∈ [n_L, n_L+1] solves Φ_total(n_W*) = N, so the
#         block lands exactly on both edges (the 1-D root-find of A:117).
#   (C3)  the domain edge is the linear-interpolated crossing of the trapezoid-accumulated phase on
#         the ancillary grid (not the nearest ancillary node).
#
# CONTROLLED OPTIONS (EGD::sw-grid-gr-refinement, both DEFAULTING to the literal App.-C behaviour so
# every grid built before that node is bitwise unchanged):
#   `wavenumber = :newtonian | :gr` — which "local WKB wavenumber" k of (C2), (C3), (C6), (C7) means:
#         the Newtonian dispersion relation of [`wkb_wavenumber`](@ref) (the literal appendix, default)
#         or the master system's own radial eigenvalue k_r = e^{μ}u^t√(ω_d²−κ²)/c_s of
#         [`wkb_wavenumber_gr`](@ref) (ERRATUM E-M5: the appendix under-resolves the inner wave zone by
#         e^{μ}u^t → √3 at the χ = 0 ISCO). NOT a free parameter of the physics — only of the mesh.
#   `extent = :adaptive | :fiducial | :inner_fiducial | :outer_fiducial` — which of the two (C2)-(C3)
#         cuts is replaced by its (C1) fiducial edge (r_ISCO resp. 5a•); `:inner_fiducial` releases the
#         inner one only, `:outer_fiducial` the outer one only. One variable at a time for the
#         domain-truncation controls.

using NumericsCore: CubicSplineInterpolant, brent, gl_nodes_weights

"""
    smoothmin(a, b) -> smoothMin(a, b)

The soft minimum of Appendix C (C9)-(C10):

    smoothMin(a, b) = [a e^{a/f} + b e^{b/f}] / [e^{a/f} + e^{b/f}],   f = −max(a, b)/8,

a softmin with sharpness 8 (exponents lie in [−8, 0], so it never overflows): symmetric,
smoothMin(a, a) = a, min(a,b) ≤ smoothMin(a,b) ≤ max(a,b), → min(a, b) as the arguments separate
(smoothMin(a, b) = a + b e^{−8} + O(a²/b) for b ≫ a).
"""
function smoothmin(a, b)
    f = -max(a, b) / 8
    wa = exp(a / f); wb = exp(b / f)
    (a * wa + b * wb) / (wa + wb)
end

"""
    corotation_radius(bg, m=1) -> r_c

Corotation radius (61): the root of ω_d(r) = m[Ω•(a•) − Ω(r)] (`omega_d`) on [r_ISCO, 5a•], found
with `NumericsCore.brent` to the precision of `T` (BigFloat inherits). Ω(r) is monotonically
decreasing on the exterior so the root is unique; it is INDEPENDENT of m ≠ 0 (m only scales ω_d):
r_c = a• exactly for 𝔭 = 𝔮 = 0 (geodesic rotation) and r_c < a• for a pressure-supported
(sub-Keplerian) disc. Throws if [r_ISCO, 5a•] does not bracket the root.
"""
function corotation_radius(bg::DiscBackground{T}, m::Integer = 1) where {T}
    m == 0 && throw(ArgumentError("corotation_radius: m must be nonzero"))
    f(r) = Omega_emri(bg) - Omega(bg, r)
    brent(f, r_isco(bg), 5 * bg.aEMRI)
end

"""
    wkb_wavenumber2(bg, m, r) -> k²

Signed square of the local WKB wavenumber from the Newtonian local dispersion relation
ω_d² = κ² + c_s² k² (the leading-order form; IMPLEMENTATION CHOICE for the undefined k(r) of
(C2), [PRELIMINARY]): k² = (ω_d² − κ²)/c_s² with ω_d = `omega_d`, κ² = `kappa2` (geodesic), c_s =
`csound`. Negative between the inner and outer Lindblad resonances (evanescent zone containing
corotation), zero at the Lindblad resonances, positive in the wave zones.
"""
wkb_wavenumber2(bg::DiscBackground, m::Integer, r) = (omega_d(bg, m, r)^2 - kappa2(bg, r)) / csound(bg, r)^2

"Local WKB wavenumber k = √max(k², 0) (zero in the evanescent zone); see [`wkb_wavenumber2`](@ref)."
function wkb_wavenumber(bg::DiscBackground, m::Integer, r)
    k2 = wkb_wavenumber2(bg, m, r)
    sqrt(max(k2, zero(k2)))
end

"""
    wkb_local(bg, m, r, wavenumber::Symbol) -> k

The local WKB wavenumber the App.-C grid recipe is evaluated with: `:newtonian` (the literal
appendix) → [`wkb_wavenumber`](@ref), `:gr` → [`wkb_wavenumber_gr`](@ref), the master system's own
radial eigenvalue (erratum E-M5). THE one dispatch point for the `wavenumber` grid option — every
other use in this file goes through it. (`wkb_wavenumber_gr` lives in `solver.jl`, which is included
after this file; the call is resolved at run time, in the same module.)
"""
function wkb_local(bg::DiscBackground, m::Integer, r, wavenumber::Symbol)
    wavenumber === :newtonian && return wkb_wavenumber(bg, m, r)
    wavenumber === :gr && return wkb_wavenumber_gr(bg, m, r)
    throw(ArgumentError("wavenumber must be :newtonian or :gr, got $(repr(wavenumber))"))
end

"""
    wkb_wavenumber_floored(bg, m, r; k_trust=10, wavenumber=:newtonian) -> k̄

The regularised wavenumber (C6), read as the r-local floor k̄(r) = max[k(r), k_trust/(h0 r)]: equal
to k wherever the WKB-trust criterion (C2) k r > k_trust/h0 holds and to the trust-threshold
wavenumber where it does not (continuous at the trust boundary, where it equals k(r_WKB)).
[PRELIMINARY] reading of A:107-109. `wavenumber` selects k through [`wkb_local`](@ref).
"""
function wkb_wavenumber_floored(bg::DiscBackground, m::Integer, r; k_trust = 10,
                                wavenumber::Symbol = :newtonian)
    k = wkb_local(bg, m, r, wavenumber)
    max(k, k_trust / (bg.h0 * r))
end

"""
    corotation_block(r_c, σ_wc, s_c, nC) -> Vector

The central corotation grid (C4): r_C,j = r_c + σ_wc sinh(s_c x_j)/sinh(s_c), x_j uniform on
[−1, 1] with `nC` points (nominal σ_wc = 10⁻³ r_c, s_c = 4, (C5)). Spacing σ_wc s_c Δx/sinh(s_c)
at r_c, σ_wc s_c coth(s_c) Δx at the edges (Δx = 2/(nC−1)); symmetric about r_c, which is a node
for odd `nC`.
"""
function corotation_block(r_c, σ_wc, s_c, nC::Integer)
    T = float(promote_type(typeof(r_c), typeof(σ_wc), typeof(s_c)))
    nC ≥ 3 || throw(ArgumentError("corotation_block: need nC ≥ 3"))
    rc, σ, s = T(r_c), T(σ_wc), T(s_c)
    sh = sinh(s)
    [rc + σ * sinh(s * (-1 + 2 * T(j - 1) / (nC - 1))) / sh for j in 1:nC]
end

# ------------------------------------------------------------------------------------------------
# WKB regions: ancillary data + cell-count coordinate
# ------------------------------------------------------------------------------------------------
"""
    WKBRegion{T}

Ancillary data of one WKB block (`side` = :L on [r_lo, r_hi] = [r_in, r_C,in], :R on
[r_C,out, r_out]): the geometric ancillary grid `r_anc` (ascending), the floored wavenumber `kbar`
(C6) and the extrapolated central spacing `Δr_in` at its nodes, the central edge `r_edge`
(= r_hi for :L, r_lo for :R) and the Gauss–Legendre rule used to integrate the cell density 1/S.
Built by [`wkb_region`](@ref); consumed by [`cell_count`](@ref) and the node placement.
"""
struct WKBRegion{T<:Real}
    side::Symbol
    r_anc::Vector{T}
    kbar::Vector{T}
    Δr_in::Vector{T}
    r_edge::T
    glx::Vector{T}
    glw::Vector{T}
end

function _wkb_region(bg::DiscBackground{T}, m::Integer, side::Symbol, r_lo, r_hi, r_edge,
                     Δr_edge, slope, n_anc::Integer, k_trust, wavenumber::Symbol) where {T}
    r_lo < r_hi || throw(ArgumentError("wkb_region: empty region [$r_lo, $r_hi] on side $side"))
    n_anc ≥ 4 || throw(ArgumentError("wkb_region: need n_anc ≥ 4"))
    lo, hi = T(r_lo), T(r_hi)
    r_anc = [lo * (hi / lo)^(T(i) / (n_anc - 1)) for i in 0:n_anc-1]
    r_anc[1] = lo; r_anc[end] = hi                      # exact edges
    kbar = [wkb_wavenumber_floored(bg, m, r; k_trust = k_trust, wavenumber = wavenumber) for r in r_anc]
    Δr_in = [T(Δr_edge) + T(slope) * abs(r - T(r_edge)) for r in r_anc]
    glx, glw = gl_nodes_weights(T, 4)
    WKBRegion{T}(side, r_anc, kbar, Δr_in, T(r_edge), glx, glw)
end

# trial spacing (C7)-(C8) at the ancillary nodes and its cubic spline S(r)
function _trial_spline(reg::WKBRegion{T}, nW) where {T}
    Δr = [smoothmin(1 / (k * T(nW)), d) for (k, d) in zip(reg.kbar, reg.Δr_in)]
    CubicSplineInterpolant(reg.r_anc, Δr)
end

# ∫_{a}^{b} dr / S(r) by the region's Gauss–Legendre rule (S > 0 required)
@inline function _gl_density_integral(S::CubicSplineInterpolant{T}, a::T, b::T, glx, glw) where {T}
    h = (b - a) / 2; c = (a + b) / 2
    acc = zero(T)
    @inbounds for q in eachindex(glx)
        s = S(c + h * glx[q], 1)
        s > 0 || throw(DomainError(c + h * glx[q], "grid: non-positive trial spacing spline"))
        acc += glw[q] / s
    end
    acc * h
end

# cumulative cell count Φ at the ancillary nodes measured from the central edge (Φ ≥ 0, increasing
# AWAY from the centre); returns (S, Φ) with Φ[i] = ∫ between r_anc[i] and the edge.
function _phi_nodes(reg::WKBRegion{T}, nW) where {T}
    S = _trial_spline(reg, nW)
    n = length(reg.r_anc)
    Φ = Vector{T}(undef, n)
    if reg.side === :L
        Φ[n] = zero(T)
        for i in n-1:-1:1
            Φ[i] = Φ[i+1] + _gl_density_integral(S, reg.r_anc[i], reg.r_anc[i+1], reg.glx, reg.glw)
        end
    else
        Φ[1] = zero(T)
        for i in 2:n
            Φ[i] = Φ[i-1] + _gl_density_integral(S, reg.r_anc[i-1], reg.r_anc[i], reg.glx, reg.glw)
        end
    end
    S, Φ
end

"""
    cell_count(reg::WKBRegion, nW) -> N(n_W)

Continuous cell count of a WKB block for `nW` cells per radian of WKB phase: N(n_W) = ∫ dr/S(r)
over the region, S the cubic-spline interpolant of the trial spacing (C8) Δr_trial =
smoothMin(1/(k̄ n_W), Δr_in). Increasing in n_W; the block realises the integer N = ⌈N(n_L)⌉ with
n_W* ∈ [n_L, n_L+1] solving N(n_W*) = N (C7, A:117).
"""
function cell_count(reg::WKBRegion{T}, nW) where {T}
    _, Φ = _phi_nodes(reg, nW)
    reg.side === :L ? Φ[1] : Φ[end]
end

# place N+1 nodes at Φ = 0, 1, …, N (exact edges); returns ascending vector
function _place_nodes(reg::WKBRegion{T}, nW, N::Integer) where {T}
    S, Φ = _phi_nodes(reg, nW)
    n = length(reg.r_anc)
    nodes = Vector{T}(undef, N + 1)
    # Φ is monotone increasing away from the edge; walk segments away from the edge
    order = reg.side === :L ? (n-1:-1:1) : (2:n)        # segment index i pairs (i, i+1) resp. (i-1, i)
    j = 1                                               # next integer level to place
    for i in order
        a, b = reg.side === :L ? (reg.r_anc[i], reg.r_anc[i+1]) : (reg.r_anc[i-1], reg.r_anc[i])
        Φa, Φb = reg.side === :L ? (Φ[i], Φ[i+1]) : (Φ[i-1], Φ[i])   # Φ at a, b (a < b)
        # on :L Φ decreases from a to b (Φa > Φb); on :R Φ increases (Φa < Φb)
        Φlo, Φhi = min(Φa, Φb), max(Φa, Φb)
        while j ≤ N - 1 && Φlo ≤ j ≤ Φhi
            nodes[j + 1] = _invert_phi(S, a, b, Φa, Φb, T(j), reg.glx, reg.glw, reg.side)
            j += 1
        end
        j > N - 1 && break
    end
    j ≤ N - 1 && throw(ErrorException("grid: placed $(j-1) of $(N-1) interior nodes (cell count mismatch)"))
    # exact edges
    if reg.side === :L
        nodes[1] = reg.r_edge; nodes[N + 1] = reg.r_anc[1]
        reverse!(nodes)
    else
        nodes[1] = reg.r_edge; nodes[N + 1] = reg.r_anc[end]
    end
    nodes
end

# solve Φ(r) = target on the ancillary segment [a, b] (Φ(a) = Φa, Φ(b) = Φb, monotone): with
# t = |target − Φa| the problem is I(r) ≡ ∫_a^r dr'/S = t, I increasing, I' = 1/S > 0 — a safeguarded
# Newton iteration with a bisection fallback inside the bracket [a, b]
function _invert_phi(S, a::T, b::T, Φa::T, Φb::T, target::T, glx, glw, side::Symbol) where {T}
    t = abs(target - Φa)
    lo, hi = a, b
    r = a + (b - a) * t / abs(Φb - Φa)                  # linear first guess
    ftol = 32 * eps(T) * max(one(T), t)
    xtol = 4 * eps(T) * max(abs(a), abs(b))
    for _ in 1:200
        g = _gl_density_integral(S, a, r, glx, glw) - t
        abs(g) ≤ ftol && return r
        if g > 0
            hi = r
        else
            lo = r
        end
        rn = r - g * S(r, 1)                              # Newton: g' = 1/S
        (lo < rn < hi) || (rn = (lo + hi) / 2)
        hi - lo ≤ xtol && return rn
        r = rn
    end
    r
end

# ------------------------------------------------------------------------------------------------
# domain selection (C1)-(C3)
# ------------------------------------------------------------------------------------------------
# returns (r_in, r_out, rWKB_L, rWKB_R) on the fiducial ancillary grid; NaN marks "no trusted point"
function _wkb_domain(bg::DiscBackground{T}, m::Integer, r_c; n_anc, k_trust, phase_min,
                    wavenumber::Symbol) where {T}
    lo, hi = r_isco(bg), 5 * bg.aEMRI
    r = [lo * (hi / lo)^(T(i) / (n_anc - 1)) for i in 0:n_anc-1]
    r[1] = lo; r[end] = hi
    k = [wkb_local(bg, m, x, wavenumber) for x in r]
    thr = k_trust / bg.h0
    nan = T(NaN)
    # left side: r < r_c
    iL = findlast(i -> r[i] < r_c && k[i] * r[i] > thr, eachindex(r))
    r_in = lo; rWKB_L = nan
    if iL !== nothing
        rWKB_L = r[iL]
        ph = zero(T); i = iL
        while i > 1
            seg = (k[i-1] + k[i]) / 2 * (r[i] - r[i-1])
            if ph + seg ≥ phase_min
                kav = (k[i-1] + k[i]) / 2
                r_in = r[i] - (phase_min - ph) / kav          # linear crossing inside the segment
                break
            end
            ph += seg; i -= 1
        end
        i == 1 && (r_in = lo)                                  # criterion not met ⇒ fiducial edge
    end
    # right side: r > r_c
    iR = findfirst(i -> r[i] > r_c && k[i] * r[i] > thr, eachindex(r))
    r_out = hi; rWKB_R = nan
    if iR !== nothing
        rWKB_R = r[iR]
        ph = zero(T); i = iR; n = length(r)
        while i < n
            seg = (k[i] + k[i+1]) / 2 * (r[i+1] - r[i])
            if ph + seg ≥ phase_min
                kav = (k[i] + k[i+1]) / 2
                r_out = r[i] + (phase_min - ph) / kav
                break
            end
            ph += seg; i += 1
        end
        i == n && (r_out = hi)
    end
    (r_in, r_out, rWKB_L, rWKB_R)
end

# ------------------------------------------------------------------------------------------------
# the grid
# ------------------------------------------------------------------------------------------------
"""
    SpiralGrid{T}

The adaptive three-block radial grid of Appendix C for one azimuthal mode `m`: `r` (ascending,
strictly increasing), the block index ranges `blocks = (L, C, R)` (the two shared edge nodes belong
to both neighbours, so `length(r) == NL + nC + NR`), the corotation radius `r_c`, the realised domain
`[r_in, r_out]`, the central edges `r_Cin`/`r_Cout`, the WKB-trust points `rWKB_L`/`rWKB_R` (NaN when
no ancillary point satisfies (C2) on that side), the targets `(nL, nR, nC)`, the realised WKB cell
counts `(NL, NR)` and cells-per-radian `(nW_L, nW_R) ∈ [n_L, n_L+1]`, and the construction
parameters `params` (σ_wc, s_c, k_trust, phase_min, n_anc, extent, wavenumber).
"""
struct SpiralGrid{T<:Real}
    r::Vector{T}
    m::Int
    r_c::T
    r_in::T
    r_out::T
    r_Cin::T
    r_Cout::T
    rWKB_L::T
    rWKB_R::T
    nL::Int
    nR::Int
    nC::Int
    NL::Int
    NR::Int
    nW_L::T
    nW_R::T
    blocks::NamedTuple{(:L, :C, :R), Tuple{UnitRange{Int}, UnitRange{Int}, UnitRange{Int}}}
    params::NamedTuple{(:σ_wc, :s_c, :k_trust, :phase_min, :n_anc, :extent, :wavenumber),
                      Tuple{T, T, T, T, Int, Symbol, Symbol}}
end

# extrapolated central spacing parameters (Δr_C,edge, dΔr_C/dr at the edge)
function _central_edge_spacing(σ_wc::T, s_c::T, nC::Integer) where {T}
    Δx = 2 / T(nC - 1)
    (σ_wc * s_c * coth(s_c) * Δx, s_c * tanh(s_c) * Δx)
end

"""
    wkb_region(bg, m, g::SpiralGrid, side) -> WKBRegion

Rebuild the ancillary WKB-region data (:L or :R) of a finished grid `g` from its stored domain and
parameters — for inspection and for the [`cell_count`](@ref) consistency checks.
"""
function wkb_region(bg::DiscBackground{T}, m::Integer, g::SpiralGrid{T}, side::Symbol) where {T}
    Δr_edge, slope = _central_edge_spacing(g.params.σ_wc, g.params.s_c, g.nC)
    if side === :L
        _wkb_region(bg, m, :L, g.r_in, g.r_Cin, g.r_Cin, Δr_edge, slope, g.params.n_anc,
                    g.params.k_trust, g.params.wavenumber)
    elseif side === :R
        _wkb_region(bg, m, :R, g.r_Cout, g.r_out, g.r_Cout, Δr_edge, slope, g.params.n_anc,
                    g.params.k_trust, g.params.wavenumber)
    else
        throw(ArgumentError("wkb_region: side must be :L or :R"))
    end
end

# N = ⌈N(n_L)⌉ and n_W* ∈ [n_L, n_L+1] with N(n_W*) = N (extends the bracket upward if the region
# gains < 1 cell per unit n_W, which only happens when Δr_in dominates everywhere)
function _solve_nW(reg::WKBRegion{T}, nL) where {T}
    N0 = cell_count(reg, T(nL))
    N = ceil(Int, N0)
    N0 == N && return (N, T(nL))
    hi = T(nL) + 1
    while cell_count(reg, hi) < N
        hi += 1
        hi > T(nL) + 64 && throw(ErrorException("grid: n_W bracket failed to reach N=$N cells"))
    end
    nW = brent(x -> cell_count(reg, x) - N, T(nL), hi)
    (N, nW)
end

"The four (C1)-(C3) domain selections `build_grid` accepts; see its `extent` keyword."
const GRID_EXTENTS = (:adaptive, :fiducial, :inner_fiducial, :outer_fiducial)

"""
    build_grid(bg, m; nL=240, nR=240, nC=24000, n_anc=10_000, extent=:adaptive,
               σ_wc_frac=1e-3, s_c=4, k_trust=10, phase_min=30, rc=nothing) -> SpiralGrid

Generate the Appendix-C adaptive radial grid for azimuthal mode `m` of the background `bg`:

1. `r_c` = [`corotation_radius`](@ref) (61) — or the explicit `rc` when given, for a model whose
   corotation radius is NOT the relativistic root of (61) (the Newtonian-limit pipeline of
   `newtonian.jl` passes its own r_c so that the sinh block is centred on ITS pole); central block
   (C4)-(C5) with σ_wc = σ_wc_frac·r_c, `s_c`, `nC` points on [r_c − σ_wc, r_c + σ_wc].
2. Domain: fiducial (C1) r_in = r_ISCO, r_out = 5a• (`extent = :fiducial`, used for the wake images);
   for `extent = :adaptive` (torques) refined by (C2)-(C3) on a geometric ancillary grid of `n_anc`
   points: r_WKB = the trusted point (k r > k_trust/h0) closest to corotation on each side, then
   the edge is moved to where ∫k dr from r_WKB reaches `phase_min` radians if that happens before
   the fiducial edge.
3. Left/right WKB blocks on [r_in, r_C,in] / [r_C,out, r_out]: ancillary grid of `n_anc` points,
   k̄ (C6) = [`wkb_wavenumber_floored`](@ref), trial step Δr_trial = smoothMin(1/(k̄ n_W), Δr_in)
   (C7)-(C10) with Δr_in the linearly extrapolated central spacing, cubic-spline interpolant S of
   Δr_trial (NumericsCore), nodes at integer values of the cell-count coordinate ∫dr/S from the
   central edge; N = ⌈N(n_L)⌉ cells with n_W ∈ [n_L, n_L+1] solved so the block lands exactly on
   both edges. (`nL`, `nR`) = target cells per radian of WKB phase; the paper's reference triple is
   (240, 240, 24000) (A:140), its convergence runs (30, 3000) … (120, 12000).

`wavenumber` selects the k of steps 2-3 through [`wkb_local`](@ref): `:newtonian` (default) is the
literal appendix, `:gr` is the master system's own radial eigenvalue e^{μ}u^t√(ω_d²−κ²)/c_s
(erratum E-M5) — a strictly finer mesh in the inner wave zone and a strictly SMALLER (C2)-(C3)
domain, since the larger k reaches the trust threshold and the (C3) phase budget sooner.

All arithmetic is in `T` (BigFloat inherits). See the file header for the recorded implementation
choices ([PRELIMINARY]) behind the under-specified steps of the appendix.
"""
function build_grid(bg::DiscBackground{T}, m::Integer; nL::Integer = 240, nR::Integer = 240,
                    nC::Integer = 24000, n_anc::Integer = 10_000, extent::Symbol = :adaptive,
                    wavenumber::Symbol = :newtonian, σ_wc_frac = 1e-3, s_c = 4, k_trust = 10,
                    phase_min = 30, rc = nothing) where {T}
    m ≥ 1 || throw(ArgumentError("build_grid: m ≥ 1 required"))
    extent in GRID_EXTENTS ||
        throw(ArgumentError("build_grid: extent must be one of $(GRID_EXTENTS), got $(repr(extent))"))
    wavenumber in (:newtonian, :gr) ||
        throw(ArgumentError("build_grid: wavenumber must be :newtonian or :gr, got $(repr(wavenumber))"))
    (nL ≥ 1 && nR ≥ 1 && nC ≥ 3) || throw(ArgumentError("build_grid: need nL, nR ≥ 1 and nC ≥ 3"))
    s_cT, k_trustT, phase_minT = T(s_c), T(k_trust), T(phase_min)
    r_c = rc === nothing ? corotation_radius(bg, m) : convert(T, rc)
    σ_wc = T(σ_wc_frac) * r_c
    blkC = corotation_block(r_c, σ_wc, s_cT, nC)
    r_Cin, r_Cout = blkC[1], blkC[end]
    if extent === :fiducial
        r_in, r_out, rWKB_L, rWKB_R = r_isco(bg), 5 * bg.aEMRI, T(NaN), T(NaN)
    else
        r_in, r_out, rWKB_L, rWKB_R = _wkb_domain(bg, m, r_c; n_anc = n_anc, k_trust = k_trustT,
                                                  phase_min = phase_minT, wavenumber = wavenumber)
        extent === :inner_fiducial && (r_in = r_isco(bg))       # release the (C2)-(C3) INNER cut only
        extent === :outer_fiducial && (r_out = 5 * bg.aEMRI)    # release the OUTER cut only
    end
    (r_in < r_Cin && r_Cout < r_out) ||
        throw(ErrorException("build_grid: domain [$r_in, $r_out] does not contain the central block [$r_Cin, $r_Cout]"))
    Δr_edge, slope = _central_edge_spacing(σ_wc, s_cT, nC)
    regL = _wkb_region(bg, m, :L, r_in, r_Cin, r_Cin, Δr_edge, slope, n_anc, k_trustT, wavenumber)
    regR = _wkb_region(bg, m, :R, r_Cout, r_out, r_Cout, Δr_edge, slope, n_anc, k_trustT, wavenumber)
    NL, nW_L = _solve_nW(regL, nL)
    NR, nW_R = _solve_nW(regR, nR)
    blkL = _place_nodes(regL, nW_L, NL)          # r_in … r_Cin   (NL+1 nodes)
    blkR = _place_nodes(regR, nW_R, NR)          # r_Cout … r_out (NR+1 nodes)
    r = vcat(blkL[1:end-1], blkC, blkR[2:end])
    all(diff(r) .> 0) || throw(ErrorException("build_grid: assembled grid is not strictly increasing"))
    blocks = (L = 1:(NL + 1), C = (NL + 1):(NL + nC), R = (NL + nC):(NL + nC + NR))
    params = (σ_wc = σ_wc, s_c = s_cT, k_trust = k_trustT, phase_min = phase_minT, n_anc = Int(n_anc),
              extent = extent, wavenumber = wavenumber)
    SpiralGrid{T}(r, Int(m), r_c, r_in, r_out, r_Cin, r_Cout, rWKB_L, rWKB_R, Int(nL), Int(nR), Int(nC),
                  NL, NR, nW_L, nW_R, blocks, params)
end
