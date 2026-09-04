# Solve/horizon_axis.jl — the horizon/axis regularity INSTRUMENT (ESSF node
# `fix-horizon-axis-layer`, GLOBAL_DAG row 58). It measures one thing: whether the solved
# m-mode `Ψ_m` is regular on the two DEGENERATE boundaries of the CRP domain — the horizon
# `r = r₊`, where `𝔄^rr = Δ/(r²+a²)` vanishes, and the polar axis `z² = 1`, where
# `𝔄^{z²} = 4z²(1−z²)/(r²+a²)` vanishes — the two faces the paper leaves without boundary
# data because "the polynomial trial space … inherently selects regular solutions"
# (main.tex:396-401 + 614-615).
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
# THE TWO NUMBERS THIS FILE ADDS.
#   * [`horizon_branch_index`](@ref) — the dimensionless exponent `ν` of the SECOND radial
#     branch at the horizon. At `r → r₊` the operator degenerates to
#     `−c x u'' + (γ^r_m − c) u' + … = 0` with `x = r − r₊` and
#     `c = dA^rr/dr|_{r₊} = (r₊−r₋)/(r₊²+a²)`, whose indicial roots are `0` (the physical,
#     analytic, ingoing branch) and `γ^r_m(r₊)/c`. With the vtu boost `H = −1`,
#     `γ^r_m(r₊) = 2im(Ω − Ω_H)` is PURELY IMAGINARY (`Ω_H = a/(r₊²+a²)` is the horizon
#     angular velocity), so the second branch is `x^{iν}` with
#     `ν = 2m|Ω − Ω_H|(r₊²+a²)/(r₊−r₋)`: unit modulus, infinitely oscillatory at `r₊`, and
#     excluded ONLY by smoothness. It is the quantity the layer tracks.
#   * [`horizon_axis_profile`](@ref) — the measured profile of `|Ψ_m|` along the innermost
#     radial node line, its axis:equator ratio, its regularity exponent, and the amplitude of
#     the near-horizon layer relative to the rest of the plane.
#
# NO PHYSICS IS DEFINED HERE. `horizon_branch_index` reads the SEALED coefficient callables
# (`flux_coeff_r`, `gamma_r_m`) rather than re-typing them, and the profile reads the SEALED
# nodal layout (`element_matrices`, `tensor_lin`, `mesh.offsets`) exactly as the figure
# drivers do. This file only measures.

export horizon_branch_index, horizon_axis_profile

"""
    horizon_branch_index(orbit::CircularOrbit, m; H = -1) -> Real

The dimensionless exponent `ν` of the SECOND (non-analytic) radial branch of the m-mode
operator at the horizon, `Ψ_m ∼ (r − r₊)^{iν}`:

    ν = |γ^r_m(r₊, H)| · (r₊² + a²) / (r₊ − r₋) = 2|m| |Ω − Ω_H| (r₊²+a²)/(r₊−r₋) ,

with `Ω_H = a/(r₊²+a²)` the horizon angular velocity (`H = −1`, the vtu v-domain boost of
the production slicing). DERIVATION, from the SEALED coefficients only: near `r₊` the flux
coefficient is `A^rr = Δ/(r²+a²) = c x + O(x²)`, `x = r − r₊`,
`c = (r₊−r₋)/(r₊²+a²) = dA^rr/dr|_{r₊}`, so

    −∂_r(A^rr ∂_r Ψ) + γ^r_m ∂_r Ψ + … = −c x Ψ'' + (γ^r_m(r₊) − c) Ψ' + …

whose indicial equation `p(p − 1)c = p(γ^r_m(r₊) − c)` has roots `p = 0` and
`p = γ^r_m(r₊)/c`. At `H = −1` the real part of `γ^r_m` vanishes at the horizon
(`2a²Δ/(r(r²+a²)²) → 0`), leaving `γ^r_m(r₊) = 2im(Ω − Ω_H)` — purely imaginary — so the
second root is `iν`: `|x^{iν}| ≡ 1`, BOUNDED but infinitely oscillatory as `x → 0`, and
excluded from the solution only by SMOOTHNESS, never by boundedness and never by the
horizon face (where `𝔄^rr = 0` makes every DG boundary term vanish identically,
main.tex:614-615 and [`NoneGhost`](@ref)).

`ν` is the control parameter of the spurious horizon–axis layer this node measures: it is a
property of `(a, r₀, m)` alone — no mesh, no worldtube radius and no penalty enters — and it
grows linearly with `m` and diverges as `r₊ → r₋` (extremal spin). MEASURED at the CRP's own
`m = 17` (a = 0.5 ISCO 1.86, a = 0 ISCO 4.63, a = 0.9 ISCO 9.85, a = −0.5 ISCO 13.4,
a = +0.998 ISCO 27.5, a = −0.998 ISCO 290); the layer is worst in the middle of that range —
see [`horizon_axis_profile`](@ref).

Type-generic (a `BigFloat` orbit gives a `BigFloat` index). `H` is the caller's boost slot,
exactly as in [`gamma_r_m`](@ref); the production v domain is `H = −1`.
"""
function horizon_branch_index(orbit::CircularOrbit, m::Integer; H = -1)
    rp, rm = orbit.rp, orbit.rm
    w = rp^2 + orbit.a^2
    c = (rp - rm) / w                       # dA^rr/dr at r₊ (flux_coeff_r factored form)
    return abs(gamma_r_m(orbit, m, rp, H)) / c
end

"""
    horizon_axis_profile(sol; line = 1, margin = 0.02) -> NamedTuple

The axis-regularity instrument of GLOBAL_DAG row 58, measured on a solved m-mode: the
profile of `|Ψ_m|` along the `line`-th distinct radial collocation line counted OUTWARD from
the inner boundary (`line = 1` is the innermost, the one that samples the horizon), as a
function of the angular chart coordinate.

`sol` is a [`ModeSolution`](@ref); an hp-AMR caller passes its sweep's final one,
`horizon_axis_profile(hist.solution)`. Returned fields:

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
function horizon_axis_profile(sol::ModeSolution; line::Integer = 1, margin::Real = 0.02)
    return _hax_profile(sol.mesh, sol.u, sol.dom, sol.orbit, Int(line), margin)
end

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

The elements of `mesh` whose lower radial face IS the horizon `r = r₊` — the elements on which
the hp-AMR loop must NOT raise the polynomial order (`(true, true)`, both dimensions). It is the
`hold` half of [`amr_step!`](@ref)`(…; horizon_hold = true)`; that keyword DEFAULTS TO `false`,
which is the pre-row-58 sweep bitwise, so nothing sealed moves until a caller opts in (see
`amr_step!` for what flipping the default costs).

WHY. The horizon face is the one place in the CRP domain where the elliptic operator has TWO
bounded local solutions and the boundary imposes nothing on either. Writing `x = r − r₊` and
`c = dA^rr/dr|_{r₊}`, the radial part degenerates to `−c x Ψ'' + (γ^r_m(r₊) − c) Ψ' + … = 0`
with indicial roots `0` and `γ^r_m(r₊)/c = iν` ([`horizon_branch_index`](@ref)): the second
root is PURELY IMAGINARY under the vtu boost `H = −1`, so `x^{iν}` has unit modulus — bounded,
infinitely oscillatory, and NOT the retarded solution. The horizon face cannot reject it,
because every DG boundary term there carries `n_i𝔄^{ij}` and `𝔄^rr(r₊) = 0` exactly
(main.tex:614-615): the face is inert for ANY choice of `Ψ_ext` and `(n·F)_ext`, and the
discrete `n·F` computed from a polynomial is identically zero whichever branch the polynomial
approximates. The ONLY thing that rejects `x^{iν}` is the paper's own argument — "they expand
fields in a polynomial series that inherently selects regular solutions" (main.tex:399-401) —
and that argument is a statement about the SIZE of the trial space. Raising the order of the
element that touches `r₊` enlarges exactly the space in which `x^{iν}` becomes representable,
and the discrete solution then slides onto it.

MEASURED (a = 0.998, r₀ = r_ISCO, m = 17, `ν` = 27.5, uniform mesh, no AMR, `L = 0`), the
axis:equator ratio on the innermost radial line and the near-horizon layer amplitude relative
to the rest of the plane, WITHOUT this hold and WITH it:

| P | without: ratio / layer | with: ratio / layer |
|---|---|---|
| 4 | 9.5 / 1.00 | 9.5 / 1.00 |
| 6 | 3.9e1 / 1.00 | 8.4 / 1.00 |
| 8 | 2.8e3 / 1.62 | 7.7 / 1.00 |
| 10 | 1.9e4 / 118.6 | 7.4 / 1.00 |
| 12 | 1.5e5 / 45.0 | 5.5 / 1.00 |

Without the hold the plane's maximum migrates to the horizon–axis corner and grows without
bound; with it the ratio converges to its physical `O(1)` value and the maximum over the rest
of the plane is UNCHANGED to every printed digit (6.072e-3 at P = 10, 8.848e-3 at P = 12).
Holding only the radial order is not enough (ratio 5.6e2 at P = 10): the corner mode is
two-dimensional — the `x^{iν}` radial branch drives the polar `(1−z²)^{−|m|}` branch, whose
representation needs the angular order too.

The hold costs the elements it names their p-convergence, which is why it is confined to the
single element layer that touches `r₊`: there the retarded field is smooth and featureless
(the whole point of the horizon-penetrating chart), and `P = 4` resolves it. Elements one
layer out keep every order the sweep gives them.

The returned dictionary is in exactly the shape [`p_refine`](@ref)'s `hold` and
[`h_refine`](@ref)'s `flags` consume; it is EMPTY for a mesh with no element on the inner
lattice column's inner edge. `mesh` may be any mesh of `dom`'s block lattice (a coarsened
multigrid level or an AMR-refined one); a block outside the lattice raises `ArgumentError`
through [`block_lattice_index`](@ref).
"""
function horizon_order_hold(dom::BlockDomain, mesh::Mesh2D)
    hold = Dict{ElementId,NTuple{2,Bool}}()
    for el in mesh.elements
        i, _ = block_lattice_index(dom, el.block)      # throws off-lattice
        i == 0 && el.segments[1].index == 0 && (hold[el] = (true, true))
    end
    return hold
end
