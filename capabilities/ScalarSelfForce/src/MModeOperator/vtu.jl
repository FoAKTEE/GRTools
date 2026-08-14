# MModeOperator/vtu.jl — vtu piecewise slicing + junction DATA (ESSF node
# `vtu-slicing-jumps`, W5-D; arXiv:2606.04998 Sec. II.B.2, main.tex:315-347). The paper's
# production slicing is the piecewise-CONSTANT boost ("vtu"): H = −1 in the v domain
# (r_* ≤ r^v_*, ingoing null), 0 in the t domain (r^v_* ≤ r_* ≤ r^u_*, the source region,
# s = t), +1 in the u domain (r_* ≥ r^u_*, outgoing null), with the height function made
# continuous by the domain constants of integration (main.tex:316-330). The jump of H at
# each transition radius puts a jump into ∂_{r_*}Ψ_m — the two junction conditions
# [∂_{r_*}Ψ_m] = −imΩ Ψ_m (eq:jump_vt main.tex:337-339, eq:jump_tu main.tex:341-343) —
# which this node exposes as pure interface DATA ([`VtuJunction`](@ref)).
#
# SCOPE BOUNDARIES (predecessors own the rest; this file only adds the vtu instance):
#   - hyperboloidal.jl (node hyperboloidal-transformation) owns the HeightFunction carrier,
#     rstar, and t_slicing(); the vtu height function is built AS a HeightFunction;
#   - ansatz.jl (node mode-ansatz-decomposition) owns the phase map Ψ_m = Ψ̃_m e^{−imΩh};
#   - tsliced.jl (node elliptic-operator-tsliced) owns the (H, H′) coefficient slots —
#     inside each open domain the vtu instance feeds them the constants (∓1|0, 0);
#   - config.jl (node runtime-config-fixture) owns the yaml-twin RunConfig whose
#     HyperboloidalSlicingTransitions carry the production transition radii.
# CONSUMER (a later node, NOT here): dg-jump-conditions turns the VtuJunction data into an
# exterior-state modification under the dg_operator.jl header contract N6 (a PURE FUNCTION
# FaceState → FaceState applied to the received neighbor state) — inherited verbatim.

"""
    VtuSlicing(r_v, r_u, rstar_v, rstar_u)

The vtu slicing geometry: the two transition radii in BOTH radial charts — Boyer–Lindquist
`r_v`/`r_u` (the values the yaml carries; they coincide with radial block bounds, so the
transitions are element faces sampled from both sides, main.tex:323) and tortoise
`rstar_v = r_*(r_v)`/`rstar_u = r_*(r_u)` (the `r^v_*`/`r^u_*` of the printed height and
boost tables, main.tex:316-330). INVARIANT (enforced by the [`vtu_slicing`](@ref)
builders, not the raw constructor): `r₊ < r_v < r_0 < r_u` — the tortoise chart is defined
at both transitions and the t domain contains the worldline (`r^v_* < r_{*0} < r^u_*`,
main.tex:323). Type-generic (BigFloat inherited).
"""
struct VtuSlicing{T<:AbstractFloat}
    r_v::T
    r_u::T
    rstar_v::T
    rstar_u::T
end

"""
    vtu_slicing(orbit, r_v, r_u)          -> VtuSlicing
    vtu_slicing(orbit, config::RunConfig) -> VtuSlicing
    vtu_slicing(orbit, oc::OrbitConfig)   -> VtuSlicing

Build the [`VtuSlicing`](@ref) geometry from the Boyer–Lindquist transition radii,
converting to the tortoise chart via [`rstar`](@ref) once. Throws `ArgumentError` unless
`r₊ < r_v < r_0 < r_u` (t domain contains the worldline, main.tex:323).

The config forms read the yaml-twin `HyperboloidalSlicingTransitions` 4-tuple
(`v→t, v→t, t→u, t→u` — vtu is the DEGENERATE sharp limit of a transition band; a
non-degenerate tuple is refused) and require the config's `(M, a, r_0)` to match `orbit`
(the tortoise conversion and the junction frequency Ω both belong to that orbit).
"""
function vtu_slicing(orbit::CircularOrbit, r_v::Real, r_u::Real)
    T = float(promote_type(typeof(orbit.r0), typeof(r_v), typeof(r_u)))
    rv, ru = T(r_v), T(r_u)
    rv < ru || throw(ArgumentError("transition radii must satisfy r_v < r_u; got ($rv, $ru)"))
    rv > orbit.rp || throw(ArgumentError(
        "r_v = $rv must lie outside the horizon r₊ = $(orbit.rp) (tortoise chart undefined)"))
    rv < orbit.r0 < ru || throw(ArgumentError(
        "the t domain must contain the worldline, r_v < r₀ < r_u (main.tex:323); " *
        "got r_v = $rv, r₀ = $(orbit.r0), r_u = $ru"))
    VtuSlicing{T}(rv, ru, rstar(orbit, rv), rstar(orbit, ru))
end

vtu_slicing(orbit::CircularOrbit, config::RunConfig) = vtu_slicing(orbit, config.orbit)

function vtu_slicing(orbit::CircularOrbit, oc::OrbitConfig)
    t = oc.hyperboloidal_slicing_transitions
    (t[1] == t[2] && t[3] == t[4]) || throw(ArgumentError(
        "HyperboloidalSlicingTransitions $t has a smooth transition band; " *
        "vtu slicing is the degenerate sharp form (v→t, v→t, t→u, t→u)"))
    (oc.black_hole_mass == orbit.M && oc.black_hole_spin == orbit.a &&
     oc.orbital_radius == orbit.r0) || throw(ArgumentError(
        "orbit/config mismatch: config (M, a, r₀) = ($(oc.black_hole_mass), " *
        "$(oc.black_hole_spin), $(oc.orbital_radius)) vs orbit " *
        "($(orbit.M), $(orbit.a), $(orbit.r0))"))
    vtu_slicing(orbit, t[1], t[3])
end

"""
    vtu_height_function(v::VtuSlicing) -> HeightFunction

The vtu instance of the [`HeightFunction`](@ref) carrier (main.tex:316-330):

    H(r_*) = −1 (v domain) | 0 (t domain) | +1 (u domain),
    h(r_*) = r^v_* − r_*   | 0            | r_* − r^u_* ,

with `h` continuous across both interfaces (the domain integration constants of
main.tex:324-330) and `Hprime = Base.zero` — H′ vanishes IDENTICALLY pointwise
(main.tex:347): the distributional H′ concentrated at the transition radii is carried by
the junction data ([`vtu_junctions`](@ref)), never by this function. The asymptotic
hyperboloidal condition H(±∞) = ±1 (main.tex:278) holds by construction.

Pointwise convention at the (double-valued) interfaces: `H` takes the t-domain value on
the CLOSED `[r^v_*, r^u_*]`. The paper samples interfaces from both sides (main.tex:323);
per-side values belong to the per-element assembly, which classifies each element by
[`vtu_domain`](@ref) and never evaluates H on a face (Legendre-Gauss volume nodes are
strictly interior). `h` is continuous, so its interface values are convention-free.
Evaluations are type-generic (BigFloat argument ⇒ BigFloat result).
"""
function vtu_height_function(v::VtuSlicing)
    HeightFunction(
        rs -> rs <= v.rstar_v ? v.rstar_v - rs :
              rs >= v.rstar_u ? rs - v.rstar_u : zero(rs - v.rstar_v),
        rs -> rs < v.rstar_v ? -one(rs - v.rstar_v) :
              rs > v.rstar_u ? one(rs - v.rstar_u) : zero(rs - v.rstar_v),
        zero,
    )
end

"""
    vtu_domain(v::VtuSlicing, r) -> Symbol

Classify a Boyer–Lindquist radius `r` into the slicing domain: `:v` for `r < r_v`, `:t` on
the closed `[r_v, r_u]`, `:u` for `r > r_u` (equivalent to the printed r_*-table,
main.tex:316-322, since r_*(r) is strictly increasing). Matches the closed-t interface
convention of [`vtu_height_function`](@ref); elements are classified by any interior point
(the yaml `NullSlicingBlocks` — first radial column v, last column u).
"""
function vtu_domain(v::VtuSlicing, r)
    r < v.r_v && return :v
    r > v.r_u && return :u
    :t
end

"""
    VtuJunction(r, rstar, coeff)

One vtu junction condition as pure interface DATA (eq:jump_vt/eq:jump_tu,
main.tex:337-343): at the interface located at Boyer–Lindquist `r` (a radial block bound —
the mesh-side key) = tortoise `rstar`, the mode function is continuous while its tortoise
derivative jumps LINEARLY IN Ψ,

    Ψ_m|⁺ − Ψ_m|⁻ = 0 ,
    ∂_{r_*}Ψ_m|⁺ − ∂_{r_*}Ψ_m|⁻ = coeff · Ψ_m|_iface ,      coeff = −imΩ,

with ⁺/⁻ the limits from the greater/lesser-r_* side (the printed LHS order). Because
`F^r = A^rr ∂_rΨ_m` IS `∂_{r_*}Ψ_m` (A^rr = dr/dr_*, main.tex:354-358), the same data
reads `[F^r] = coeff·Ψ_m` — directly consumable by the DG jump hooks in `FaceState`
(φ = n·F) language under the dg_operator.jl N6 exterior-modification contract
(node dg-jump-conditions; not built here).
"""
struct VtuJunction{T<:AbstractFloat}
    r::T
    rstar::T
    coeff::Complex{T}
end

"""
    vtu_junctions(v::VtuSlicing, orbit, m) -> (vt::VtuJunction, tu::VtuJunction)

The two junction conditions of the vtu slicing as data: `vt` at `(r_v, r^v_*)`
(eq:jump_vt) and `tu` at `(r_u, r^u_*)` (eq:jump_tu), both with the SAME linear-in-Ψ jump
coefficient `−imΩ` (Ω = `orbit.Omega`, the frequency of the mode phase `e^{im(φ̄−Ωs)}`).
At `m = 0` both coefficients vanish — the slicing is degenerate (h drops out of the phase
map) and there are no jumps. `orbit` must be the orbit the slicing was built from.
"""
function vtu_junctions(v::VtuSlicing{T}, orbit::CircularOrbit, m::Integer) where {T}
    c = -im * (m * orbit.Omega)
    P = promote_type(T, real(typeof(c)))
    (vt = VtuJunction{P}(v.r_v, v.rstar_v, c),
     tu = VtuJunction{P}(v.r_u, v.rstar_u, c))
end
