# Solve/horizon_axis.jl — the horizon/axis regularity INSTRUMENT (ESSF node
# `fix-horizon-axis-layer`, GLOBAL_DAG row 58). It measures one thing: whether the solved
# m-mode `Ψ_m` is regular on the two DEGENERATE boundaries of the CRP domain — the horizon
# `r = r₊`, where `𝔄^rr = Δ/(r²+a²)` vanishes, and the polar axis `z² = 1`, where
# `𝔄^{z²} = 4z²(1−z²)/(r²+a²)` vanishes — the two faces the paper leaves without boundary
# data because "the polynomial trial space … inherently selects regular solutions"
# (main.tex:396-401 + 614-615).
#
# THE LAYER IS A CONJUNCTION of this radial branch and the polar `(1−z²)^{−|m|}` branch, both
# of which only the trial space rejects: suppressing EITHER alone removes it (measured both
# ways). No corner asymptotic analysis exists, so the file claims a conjunction, never a
# direction. `horizon_hold` DEFAULTS TO `false`; the rule is to gate the MEASURED
# `layer_ratio` and switch the hold on when it fails — never by spin, never by `ν`.
#
# WHAT REGULARITY MEANS FOR THIS VARIABLE. The mode ansatz factors `sin^{|m|}θ` OUT of the
# field (`ansatz.jl`: Φ = Σ_m Ψ_m sin^{|m|}θ e^{im(φ̄−Ωt)} / r), so `Ψ_m` is the variable that
# is ANALYTIC at the axis and approaches a finite, generically non-zero limit there — it is
# `Φ`'s m-mode that falls like `sin^{|m|}θ`, never `Ψ_m`. The healthy signature is therefore
# a BOUNDED profile along the innermost radial node line, with an axis:equator ratio of order
# one and a fitted `log|Ψ_m|` vs `log(1−z²)` slope near ZERO. (A published statement that the
# plane must fall like `(1−z²)^{m/2}` reads the ansatz the wrong way round; the reference's
# own Fig. 5 surfaces are flat in `z²`.)
#
# WHAT THIS FILE ADDS — two measurements and one refinement rule.
#   * [`horizon_branch_index`](@ref) — the dimensionless exponent `ν` of the SECOND radial
#     branch at the horizon. At `r → r₊` the operator degenerates to
#     `−c x u'' + (γ^r_m − c) u' + … = 0` with `x = r − r₊` and
#     `c = dA^rr/dr|_{r₊} = (r₊−r₋)/(r₊²+a²)`, whose indicial roots are `0` (the physical,
#     analytic, ingoing branch) and `γ^r_m(r₊)/c`. With the vtu boost `H = −1`,
#     `γ^r_m(r₊) = 2im(Ω − Ω_H)` is PURELY IMAGINARY (`Ω_H = a/(r₊²+a²)` is the horizon
#     angular velocity), so the second branch is `x^{p}` with the SIGNED root
#     `p = 2im(Ω − Ω_H)(r₊²+a²)/(r₊−r₋)` — sign = sign(Ω − Ω_H) — of unit modulus,
#     infinitely oscillatory at `r₊`, and excluded ONLY by smoothness. The function returns
#     `ν = |p|`, and it is valid at `H = −1` alone (it throws otherwise: off the v-domain
#     boost `β_m` diverges like `1/Δ` and enters the indicial polynomial).
#   * [`horizon_axis_profile`](@ref) — the measured profile of `|Ψ_m|` along the innermost
#     radial node line, its axis:equator ratio, its regularity exponent, and the amplitude of
#     the near-horizon layer relative to the rest of the plane.
#   * [`horizon_order_hold`](@ref) — the one refinement rule those two measurements forced:
#     the element layer that touches `r₊` must keep its polynomial order, because that order
#     IS the only thing rejecting the `x^{iν}` branch. It is the `hold` half of
#     [`amr_step!`](@ref)`(…; horizon_hold = true)`, and that keyword defaults to `false`.
#
# NO PHYSICS IS DEFINED HERE. `horizon_branch_index` reads the SEALED coefficient callables
# (`flux_coeff_r`, `gamma_r_m`) rather than re-typing them, and the profile reads the SEALED
# nodal layout (`element_matrices`, `tensor_lin`, `mesh.offsets`) exactly as the figure
# drivers do. No coefficient, no boundary condition and no solver setting is defined or
# changed here; `horizon_order_hold` only NAMES elements.

export horizon_branch_index, horizon_axis_profile, horizon_order_hold

"""
    horizon_branch_index(orbit::CircularOrbit, m; H = -1) -> Real

`ν = |p|`, the MODULUS of the exponent of the SECOND (non-analytic) radial branch of the
m-mode operator at the horizon. The branch itself is `Ψ_m ∼ (r − r₊)^{p}` with the SIGNED,
purely imaginary

    p = γ^r_m(r₊, H = −1) / c = 2im(Ω − Ω_H) (r₊²+a²)/(r₊−r₋) ,    ν = |p| ,

`c = (r₊−r₋)/(r₊²+a²) = dA^rr/dr|_{r₊}` and `Ω_H = a/(r₊²+a²)` the horizon angular velocity.
The sign is that of `Ω − Ω_H` and is physical (it is the sense in which the spurious mode
winds as `x → 0`): at `a = 0.998`, `r₀ = 2 M`, `m = 17` the root is `p = −118.94 i`, i.e.
`x^{−iν}` with `ν = 118.94`. This function returns `ν`; the sign is
`sign(Ω − Ω_H) = sign(imag(gamma_r_m(orbit, m, r₊, -1)))`.

DERIVATION, from the SEALED coefficients only: near `r₊` the flux coefficient is
`A^rr = Δ/(r²+a²) = c x + O(x²)` with `x = r − r₊`, so

    −∂_r(A^rr ∂_r Ψ) + γ^r_m ∂_r Ψ + … = −c x Ψ'' + (γ^r_m(r₊) − c) Ψ' + …

whose indicial equation `p(p − 1)c = p(γ^r_m(r₊) − c)` has roots `0` and `γ^r_m(r₊)/c`. At
`H = −1` the real part of `γ^r_m` vanishes at the horizon (`2a²Δ/(r(r²+a²)²) → 0`), leaving
`γ^r_m(r₊) = 2im(Ω − Ω_H)` — purely imaginary — so `|x^p| ≡ 1`: BOUNDED but infinitely
oscillatory as `x → 0`, and excluded from the solution only by SMOOTHNESS, never by
boundedness and never by the horizon face (where `𝔄^rr = 0` makes every DG boundary term
vanish identically, main.tex:614-615 and [`NoneGhost`](@ref)).

`ν` is a property of `(a, r₀, m)` alone — no mesh, no worldtube radius, no penalty — and it
grows linearly with `m` and diverges as `r₊ → r₋`. It orders the layer AT A FIXED
CONFIGURATION only; across configurations it does not (see [`horizon_order_hold`](@ref)), so
it is a diagnostic, never a switch.

`H` IS ACCEPTED ONLY AS `−1`, the vtu v-domain boost of the production slicing, and any other
value throws. The derivation above drops `β_m`, which is legitimate only there: at `H = −1`
the `1/Δ` terms of `β_m` cancel (main.tex:386, `beta_m`'s docstring), while at `H = 0` `β_m`
genuinely diverges like `1/Δ`, contributes at the SAME order as the two kept terms and changes
the indicial polynomial. The keyword is kept so the restriction is explicit rather than
implied. Type-generic (a `BigFloat` orbit gives a `BigFloat` index).
"""
function horizon_branch_index(orbit::CircularOrbit, m::Integer; H = -1)
    H == -1 || throw(ArgumentError(
        "horizon_branch_index: H = $H — the indicial roots 0 and γ^r_m(r₊)/c are derived " *
        "with β_m dropped, which holds ONLY for the vtu v-domain boost H = −1 where the " *
        "1/Δ terms of β_m cancel (main.tex:386). Off that slicing β_m ~ b/x enters at the " *
        "same order and the polynomial is a different one; this function does not model it"))
    rp, rm = orbit.rp, orbit.rm
    w = rp^2 + orbit.a^2
    c = (rp - rm) / w                       # dA^rr/dr at r₊ (flux_coeff_r factored form)
    return abs(gamma_r_m(orbit, m, rp, H)) / c
end

"""
    horizon_axis_profile(sol::ModeSolution; line = 1, margin = 0.02) -> NamedTuple
    horizon_axis_profile(dom, orbit, mesh, u; line = 1, margin = 0.02) -> NamedTuple

The axis-regularity instrument of GLOBAL_DAG row 58, measured on a solved m-mode: the
profile of `|Ψ_m|` along the `line`-th distinct radial collocation line counted OUTWARD from
the inner boundary (`line = 1` is the innermost, the one that samples the horizon), as a
function of the angular chart coordinate.

`sol` is a [`ModeSolution`](@ref); an hp-AMR caller passes its sweep's final one,
`horizon_axis_profile(hist.solution)`.

THE FOUR-ARGUMENT METHOD IS THE SEAM FOR POOLED RESULTS. A process fan-out
([`solve_modes_parallel`](@ref)) returns per-mode records that carry a mesh and a nodal vector
but NO live `ModeSolution` — the operator, the stack and the caches do not survive the worker.
Such a caller measures the layer directly from what it has,
`horizon_axis_profile(dom, orbit, res.levels[end].mesh, res.u)`, with the SAME arithmetic: the
one-argument method delegates to it, so the two agree bitwise (gated). `dom` and `orbit` supply
only the angular chart, the lattice check and `r₊`; nothing is re-solved. Returned fields:

  - `r` — the radius of the sampled node line, and `x2` / `zeta` — the angular collocation
    coordinates of that line in the mesh's own chart and as `z²` (they coincide in the
    equatorial `:z2` chart; in the full `:z` chart `zeta = z²` and the profile is taken over
    `|z|`), sorted by `zeta`;
  - `absu` — `|Ψ_m|` at those nodes (inside the worldtube this is `|Ψ^ℛ_m|`, exactly as
    `sol.u` holds it);
  - `ratio` — `absu[end] / absu[1]`, the axis:equator ratio. **Healthy is O(1)**; the mode
    variable carries `sin^{|m|}θ` factored out, so it must approach a finite limit at the
    axis (file header);
  - `exponent` — the least-squares slope of `log(absu)` against `log(1 − zeta)` over the
    nodes with `zeta > 0.5`, i.e. the fitted `q` in `|Ψ_m| ∼ (1−z²)^q` near the axis.
    **Healthy is ≈ 0**; a strongly NEGATIVE `q` is the `(1−z²)^{−|m|}` polar branch leaking
    in;
  - `layer_max` — `max|Ψ_m|` over the nodes within `margin` of the inner boundary
    (`r ≤ r₊ + margin`), `bulk_max` — the maximum over all the others, and `layer_ratio` =
    `layer_max / bulk_max`. **Healthy is ≤ 1**: the near-horizon strip must never carry the
    largest value of the plane.

The instrument makes no claim about WHY a profile is unhealthy; the mechanism and its
controlled variations live in the node's evidence, and the control parameter is
[`horizon_branch_index`](@ref).
"""
function horizon_axis_profile(dom::BlockDomain, orbit::CircularOrbit, mesh::Mesh2D,
                              u::AbstractVector; line::Integer = 1, margin::Real = 0.02)
    return _hax_profile(mesh, u, dom, orbit, Int(line), margin)
end

horizon_axis_profile(sol::ModeSolution; line::Integer = 1, margin::Real = 0.02) =
    horizon_axis_profile(sol.dom, sol.orbit, sol.mesh, sol.u; line = line, margin = margin)

function _hax_nodal(mesh::Mesh2D{RT}, u::AbstractVector, chart::Symbol) where {RT}
    n = mesh.ndofs
    length(u) == n || throw(ArgumentError(
        "horizon_axis_profile: the nodal vector has length $(length(u)) but the mesh has " *
        "$n dofs — it belongs to a different level"))
    rs = RT[]; xs = RT[]; vs = float(real(eltype(u)))[]
    for e in eachindex(mesh.elements)
        el = mesh.elements[e]
        blk = mesh.blocks[mesh.blockindex[el.block]]
        em = element_matrices(RT, blk, el, mesh.orders[e])
        nx, ny = mesh.orders[e]
        off = mesh.offsets[e]
        for i in 1:nx, j in 1:ny
            push!(rs, em.x[1][i]); push!(xs, em.x[2][j])
            push!(vs, abs(u[off + tensor_lin(i, j, nx)]))
        end
    end
    chart === :z2 || chart === :z || throw(ArgumentError(
        "horizon_axis_profile: unknown angular chart $(repr(chart)) (expected :z2 or :z)"))
    zeta = chart === :z2 ? copy(xs) : xs .^ 2           # ζ = z² in both charts
    return rs, xs, zeta, vs
end

function _hax_profile(mesh::Mesh2D{RT}, u::AbstractVector, dom::BlockDomain,
                      orbit::CircularOrbit, line::Int, margin::Real) where {RT}
    line >= 1 || throw(ArgumentError("horizon_axis_profile: line must be ≥ 1, got $line"))
    rs, xs, zeta, vs = _hax_nodal(mesh, u, dom.chart)
    radii = sort(unique(rs))
    line <= length(radii) || throw(ArgumentError(
        "horizon_axis_profile: line = $line but the mesh has only $(length(radii)) distinct " *
        "radial collocation lines"))
    r_line = radii[line]
    idx = findall(==(r_line), rs)
    p = sortperm(zeta[idx])
    sel = idx[p]
    z = zeta[sel]; a = vs[sel]
    ratio = a[end] / a[1]
    hi = findall(>(0.5), z)
    # the polar half of the line, but never fewer than three nodes: with the horizon layer
    # held at P = 4 the polar block puts only two Gauss nodes above z² = 0.5
    length(hi) < 3 && (hi = collect(max(1, length(z) - 2):length(z)))
    q = _hax_slope(z[hi], a[hi])
    strip = findall(<=(orbit.rp + margin), rs)
    rest = findall(>(orbit.rp + margin), rs)
    lmax = isempty(strip) ? zero(eltype(vs)) : maximum(vs[strip])
    bmax = isempty(rest) ? zero(eltype(vs)) : maximum(vs[rest])
    return (r = r_line, x2 = xs[sel], zeta = z, absu = a, ratio = ratio, exponent = q,
            layer_max = lmax, bulk_max = bmax, layer_ratio = lmax / bmax)
end

# least-squares slope of log|u| against log(1 − ζ); NaN when fewer than three usable nodes
function _hax_slope(z::AbstractVector, a::AbstractVector)
    length(z) >= 3 || return NaN
    x = [log(1 - zi) for zi in z]
    y = [log(ai) for ai in a]
    keep = findall(k -> isfinite(x[k]) && isfinite(y[k]), eachindex(x))
    length(keep) >= 3 || return NaN
    xm = sum(x[keep]) / length(keep); ym = sum(y[keep]) / length(keep)
    sxx = sum((x[k] - xm)^2 for k in keep)
    sxx == 0 && return NaN
    return sum((x[k] - xm) * (y[k] - ym) for k in keep) / sxx
end

"""
    horizon_order_hold(dom::BlockDomain, mesh::Mesh2D) -> Dict{ElementId,NTuple{2,Bool}}

The elements of `mesh` whose lower radial face IS the horizon `r = r₊`, each mapped to
`(true, true)` — the elements on which the hp-AMR loop must NOT raise the polynomial order.
It is the `hold` half of [`amr_step!`](@ref)`(…; horizon_hold = true)`.

WHY. At `r₊` the operator has TWO bounded local solutions and the boundary imposes nothing on
either: the indicial roots are `0` (analytic, retarded) and `p = γ^r_m(r₊)/c`, purely
imaginary under the vtu boost, so `|x^p| ≡ 1` ([`horizon_branch_index`](@ref)). The face
cannot reject `x^p`, because every DG boundary term there carries `n_i𝔄^{ij}` and
`𝔄^rr(r₊) = 0` exactly (main.tex:614-615), and the discrete `n·F` of a polynomial is
identically zero whichever branch it approximates. The only rejection is the paper's own
"they expand fields in a polynomial series that inherently selects regular solutions"
(main.tex:399-401) — a statement about the SIZE of the trial space. Raising the order of the
element that touches `r₊` enlarges exactly the space in which `x^p` is representable.

WHY MESH REFINEMENT CANNOT DO IT AND ORDER CAN. On an element `(0, h)` the substitution
`x = h s` gives `x^p = h^p s^p` with `|h^p| = 1`: the branch's SHAPE on the horizon element is
scale-invariant, so h-refinement reproduces it unchanged and cannot reject it, while the
trial-space size is not scale-invariant.

BOTH DEGENERATE FACES ARE NECESSARY, NEITHER SUFFICIENT — an evidence-backed CONJUNCTION, not
a direction: no corner asymptotic analysis exists, and suppressing EITHER the radial `x^p`
branch or the polar `(1−z²)^{−|m|}` branch alone removes the layer (measured both ways). The
hold is on both dimensions for the same reason: holding only the radial order leaves it.

MEASURED, in one line (a = 0.998, r₀ = r_ISCO, m = 17, uniform mesh p-refined from P = 4): the
axis:equator ratio on the innermost radial line DIVERGES 9.5 → 1.5e5 over P = 4…12 without the
hold and CONVERGES 9.5 → 5.5 with it, while the maximum over the rest of the plane is
unchanged to every printed digit. Every table — that sweep, the uniform-h and hp-AMR sweeps,
the six-spin `ν` list, the coefficient variations, the blast radius — is in the node's
evidence, `evidence/fix-horizon-axis-layer/ab_table.md`.

THE SWITCH RULE IS MEASURED, NEVER INFERRED. `horizon_hold` DEFAULTS TO `false` and stays
there: the reference has no counterpart to this hold, the paper's own observables are
insensitive to the layer (`F_r` moves 2.6e-5 relative when it is switched on), the flip would
move sealed AMR bookkeeping, and the hold's own truncation cost on the held layer is not
separately bounded. `ν` does not order the layer across configurations either. So a node that
publishes a PLANE calls [`horizon_axis_profile`](@ref), gates its `layer_ratio ≤ 1.05`, and
enables `horizon_hold` when that fails — never by spin, never by `ν`. This function is the
documented cure for axis-contaminated planes, not a default the solver takes on its own.

The hold costs the elements it names their p-convergence, which is why it is confined to the
single element layer that touches `r₊`, where the retarded field is smooth and featureless in
the horizon-penetrating chart.

PRECONDITIONS AND SHAPE. The named elements are the horizon layer only because
[`block_domain`](@ref) REFUSES a lattice whose inner radial bound is not the outer horizon
(`blocks.jl`, 8 ulp relative); `crp_paper_config` takes `inner_radius` freely, so the
`BlockDomain` is what pins it. `mesh` must be a mesh OF THAT DOMAIN's lattice — the block-map
endpoints are checked against `dom`'s own bounds, so a same-shaped mesh from a different
`BlockDomain` is refused rather than silently held (the same guard
[`boundary_conditions`](@ref) uses). It may be any such mesh: a coarsened multigrid level or
an AMR-refined one. Under refinement only the INNER radial child keeps the horizon face, and
only it is named; the two clauses `block_lattice_index(dom, block)[1] == 0` (inner lattice
column) and `segments[1].index == 0` (inner radial child) are independent. Angular splits and
2:1-balance-induced splits produce several held children, each inheriting its parent's order.
The result is exactly the shape [`p_refine`](@ref)'s `hold` and [`h_refine`](@ref)'s `flags`
consume, and is EMPTY for a mesh with no element on that edge.
"""
function horizon_order_hold(dom::BlockDomain, mesh::Mesh2D)
    _bc_require_same_lattice(dom, mesh)     # block-map endpoints, not merely id range
    hold = Dict{ElementId,NTuple{2,Bool}}()
    for el in mesh.elements
        i, _ = block_lattice_index(dom, el.block)      # throws off-lattice
        i == 0 && el.segments[1].index == 0 && (hold[el] = (true, true))
    end
    return hold
end

# ─────────────────────────────────────────────────────────────────────────────────────────────
# THE SHARED HOLD-SELECTION RULE.
#
# Two figure nodes publish an object whose regularity this file's `layer_ratio` gates, and they
# had grown two different policies for turning a measurement into a decision. `fig_modes.jl`
# holds iff the default sweep's maximum over `m` exceeds 1.05 and then plots the held sweep;
# that is correct on its three panels, where the hold always measured smaller — but that is an
# EMPIRICAL fact about those panels, not a property of the hold, and a figure that spans a
# compressed near-extremal ISCO cannot assume it. `select_hold_sweep` is the one rule both use:
# the same trigger, plus the assertion that the held sweep actually reduces the ratio.

"""
    HoldSelection{S}

What [`select_hold_sweep`](@ref) decided, and every number it decided it on.

| field | meaning |
|---|---|
| `sweep::S` | the sweep to publish — the very object passed in, never a copy |
| `hold::Bool` | is that the HELD sweep? |
| `default_ratio::Float64` | the default sweep's `layer_ratio`, maximised over `m` |
| `held_ratio::Float64` | the held sweep's, or `NaN` when no held sweep was consulted |
| `gate::Float64` | the threshold the trigger used |
| `cleared::Bool` | does the PUBLISHED sweep's ratio sit at or below the gate? |
| `hole::Bool` | could the rule not do its job? (no held sweep to fall back on, a hold that did not reduce the ratio, or one that reduced it without clearing the gate) |
| `reason::String` | the sentence to put in the caller's log and on its figure |

`hole` and `!cleared` are different questions: a hold that halves a ratio from 10 to 2 has done
its job (`hole` is still true because 2 exceeds the gate), while a hold that raises it has not
(`hole` true, and the DEFAULT is published).
"""
struct HoldSelection{S}
    sweep::S
    hold::Bool
    default_ratio::Float64
    held_ratio::Float64
    gate::Float64
    cleared::Bool
    hole::Bool
    reason::String
end

function Base.show(io::IO, ::MIME"text/plain", s::HoldSelection)
    println(io, "HoldSelection  publish the ", s.hold ? "HELD" : "DEFAULT", " sweep")
    println(io, "  layer_ratio (max over m): default ", s.default_ratio,
            isnan(s.held_ratio) ? "" : string("   held ", s.held_ratio))
    println(io, "  gate ", s.gate, "  cleared ", s.cleared, "  hole ", s.hole)
    print(io,   "  ", s.reason)
end

Base.show(io::IO, s::HoldSelection) =
    print(io, "HoldSelection(", s.hold ? "held" : "default", ", ", s.default_ratio,
          isnan(s.held_ratio) ? "" : string(" -> ", s.held_ratio), ")")

"maximum of `v` over its finite entries; a mode that failed carries `NaN` and is skipped."
function _hax_ratio_max(v, what::AbstractString)
    isempty(v) && throw(ArgumentError(
        "select_hold_sweep: the $(what) sweep carries no layer_ratio at all; there is nothing " *
        "to maximise over and the rule cannot decide"))
    m = -Inf
    for x in v
        y = Float64(x)
        isnan(y) && continue
        y > m && (m = y)
    end
    isfinite(m) || throw(ArgumentError(
        "select_hold_sweep: every layer_ratio of the $(what) sweep is NaN (no mode of it was " *
        "measured); the rule cannot decide"))
    return m
end

"""
    select_hold_sweep(default::Pair, held = nothing; gate = 1.05) -> HoldSelection

Decide which of two already-solved sweeps a figure should publish, from their measured
[`horizon_axis_profile`](@ref) `layer_ratio`s. Each argument is a `sweep => layer_ratios` pair,
where `layer_ratios` is that sweep's ratio for EVERY mode; `held` may be `nothing` when no held
sweep was run.

THE RULE, in one sentence: publish the held sweep exactly when the default's ratio crosses the
`gate` AND the held sweep's ratio is strictly smaller; otherwise publish the default, and say
so.

  * **The maximum over `m` is the input, never the top mode.** `layer_ratio` is not monotone in
    `m` — measured at `a = +0.998`, `r₀ = 2 M` it is 0.58 at `m = 0`, crosses one at `m = 6`,
    peaks at 3.889 at `m = 12` and falls back to 3.6e-3 by `m = 20`, so a top-mode probe reads
    a contaminated panel as healthy.
  * **The default is not consulted against the held sweep when it is already inside the gate.**
    A caller that only runs the held sweep on demand (both figure drivers do) then never pays
    for it, and `held_ratio` comes back `NaN`.
  * **A hold that does not reduce the ratio is refused, not silently substituted.** The hold
    costs the horizon layer its p-convergence; taking it without a measured improvement would
    spend that for nothing. `hole` is set and `reason` says which case it was.
  * `NaN` ratios (a mode that failed) are skipped; a sweep whose ratios are all `NaN`, or empty,
    is refused rather than defaulted.

The gate's landed value is the 1.05 of this file's own switch rule. Nothing here solves, reads a
mesh or touches the field: it is bookkeeping over two measurements, so a caller can — and both
callers do — apply it to whatever object it is holding.

```jldoctest
julia> s = select_hold_sweep(:default => [0.5, 23.94, 3.2], :held => [0.83, 0.03, 2.0e-4]);

julia> s.sweep, s.hold, s.default_ratio, s.held_ratio, s.cleared, s.hole
(:held, true, 23.94, 0.83, true, false)
```
"""
function select_hold_sweep(default::Pair, held = nothing; gate::Real = 1.05)
    g = Float64(gate)
    (isfinite(g) && g > 0) || throw(ArgumentError(
        "select_hold_sweep: the gate must be a positive finite ratio, got $(gate)"))
    dsweep, dratios = first(default), last(default)
    dmax = _hax_ratio_max(dratios, "default")
    if dmax <= g
        return HoldSelection(dsweep, false, dmax, NaN, g, true, false,
                             "the default sweep's layer_ratio $(dmax) is inside the gate " *
                             "$(g); the held sweep was not consulted")
    end
    if held === nothing
        return HoldSelection(dsweep, false, dmax, NaN, g, false, true,
                             "the default sweep's layer_ratio $(dmax) crosses the gate $(g) " *
                             "and no held sweep was measured [HOLE]")
    end
    hsweep, hratios = first(held), last(held)
    hmax = _hax_ratio_max(hratios, "held")
    if hmax < dmax
        cleared = hmax <= g
        return HoldSelection(hsweep, true, dmax, hmax, g, cleared, !cleared,
                             cleared ?
                             "the hold takes the layer_ratio from $(dmax) to $(hmax), inside " *
                             "the gate $(g)" :
                             "the hold takes the layer_ratio from $(dmax) to $(hmax), which " *
                             "still crosses the gate $(g) [HOLE]")
    end
    return HoldSelection(dsweep, false, dmax, hmax, g, false, true,
                         "the hold did not reduce the layer_ratio ($(dmax) -> $(hmax)), so the " *
                         "DEFAULT sweep is published and the gate $(g) stays crossed [HOLE]")
end
