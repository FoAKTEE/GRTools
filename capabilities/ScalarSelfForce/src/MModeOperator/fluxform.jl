# MModeOperator/fluxform.jl — horizon-penetrating compactification: the chart-converted
# coefficient bundle (ESSF node `horizon-penetrating-compactification`, W6-D;
# arXiv:2606.04998 Sec. II.B.3, main.tex:354-360). The hyperboloidal slicing eliminates
# asymptotic oscillations, which enables compactifying the radial chart; because
# Boyer–Lindquist r is inherently compact near the horizon, the paper discretizes IN r and
# expresses the differential operator in r derivatives through
#
#     ∂/∂r_* = Δ/(r² + a²) · ∂/∂r        (main.tex:356-358),
#
# placing the horizon at the finite radius r = r₊ instead of r_* → −∞. Combined with the
# vtu boost this yields ingoing Eddington–Finkelstein coordinates near the horizon (v),
# Boyer–Lindquist near the puncture (t), and outgoing EF in the wavezone (u), where a ~1/r
# map (EllipticCore `InverseRadialMap1D`, eq:compaction) partially compactifies large radius.
#
# SCOPE BOUNDARIES (predecessors own the substance; this file is the CLOSURE):
#   - tsliced.jl (node elliptic-operator-tsliced) owns the r-chart coefficient functions —
#     the printed A^rr = Δ/(r²+a²), A^zz, β_m, γ^r_m, γ^z_m (main.tex:366-383) ARE the
#     chart-converted forms, with pointwise (H, H′) slots;
#   - hyperboloidal.jl (node hyperboloidal-transformation) owns r_*(r) and the
#     HeightFunction carrier + slicing_kinematics (and sealed dr_*/dr = 1/A^rr);
#   - vtu.jl (node vtu-slicing-jumps) owns the piecewise production boost H ∈ {−1, 0, +1}.
# This node adds the one composition the DISCRETIZATION consumes — the pointwise bundle
# binding slicing kinematics into the coefficient slots — plus the horizon-penetration
# closure gates (test_fluxform.jl): r-chart ≡ r_*-chart under per-chart numerical
# differentiation, all-coefficient finiteness at r = r₊(1 + 1e-12) with H = −1 and down the
# collocation-reachable approach, u-domain smoothness across the 1/r-compaction map with
# H = +1, a = 0 degeneracy, BigFloat(256) inheritance.
#
# Part 2 (node `first-order-flux-form`, W6-E, GLOBAL_DAG row 30): `fluxform_callables` —
# the ASSEMBLED per-region form, coefficient callables over (r, z) with the constant vtu
# slab boost, in exactly the shape EllipticCore's DG residual layer consumes.
#
# CONSUMERS (later nodes): effective-source-direct (operator↔source identity),
# domain-blocks-map / solver-stack-ssf (per-element coefficient evaluation),
# equatorial-symmetry-z2 (the z² variant bundle).

"""
    fluxform_coefficients(orbit, m, hf::HeightFunction, r, z)
        -> (Arr, Azz, beta, gamma_r, gamma_z)

The m-mode elliptic operator coefficient bundle IN THE DISCRETIZATION CHART q = (r, z)
(main.tex:354-383): everything a nodal DG discretization — which differentiates in
Boyer–Lindquist r, never in r_* — needs at one quadrature point, with the slicing boost
already composed in. The first-order flux form reads

    △_m Ψ_m = −∂_r(Arr ∂_r Ψ_m) − ∂_z(Azz ∂_z Ψ_m) + beta Ψ_m
              + gamma_r ∂_r Ψ_m + gamma_z ∂_z Ψ_m ,

where the horizon-penetrating conversion ∂_{r_*} = Δ/(r²+a²) ∂_r (main.tex:356-358) is
fully absorbed: `Arr` = Δ/(r²+a²) = dr/dr_* is both the radial flux coefficient and the
chart Jacobian (so `Arr ∂_r Ψ` IS `∂_{r_*}Ψ`, and the horizon sits at finite r = r₊ where
`Arr` vanishes exactly — no explicit horizon boundary condition, main.tex:385). One
`slicing_kinematics` evaluation (one r_*(r)) feeds the boost `(H, H′)` of `hf` into the
tsliced.jl slots; with the vtu v-domain boost H = −1 every returned coefficient is finite
down the entire collocation-reachable horizon approach (the 1/Δ pieces of `beta` cancel,
main.tex:386 — gated at r = r₊(1 + 1e-12) in the node tests), and with H = +1 the bundle
stays smooth across the wavezone 1/r-compaction map (eq:compaction).

Pure and pointwise (evaluate at each quadrature node); type-generic — BigFloat orbits,
slicings, and coordinates propagate precision through every field. Field types follow
tsliced.jl: `Arr`, `Azz`, `gamma_z` real; `beta`, `gamma_r` complex. The z² equatorial
variant and the junction data are NOT bundled here (nodes equatorial-symmetry-z2 and
vtu-slicing-jumps / dg-jump-conditions).
"""
function fluxform_coefficients(orbit::CircularOrbit, m::Integer, hf::HeightFunction, r, z)
    kin = slicing_kinematics(hf, orbit, r)
    (Arr     = flux_coeff_r(orbit, r),
     Azz     = flux_coeff_z(orbit, r, z),
     beta    = beta_m(orbit, m, r, z, kin.H, kin.Hprime),
     gamma_r = gamma_r_m(orbit, m, r, kin.H),
     gamma_z = gamma_z_m(orbit, m, r, z))
end

"""
    fluxform_callables(orbit, m, v::VtuSlicing, region::Symbol)
        -> (A = (Arr, Azz), beta, gamma = (gamma_r, gamma_z))

The ASSEMBLED first-order flux form of one vtu region (node `first-order-flux-form`,
GLOBAL_DAG row 30): the coefficient CALLABLES `(r, z) -> value` of the discretization
chart, in exactly the shape EllipticCore's DG residual layer consumes —
`element_coefficients(em, blk, id, fc.A, fc.beta, fc.gamma)` (eq:dgres_strong coefficients,
main.tex:565-575). The operator they assemble is eq:flux/eq:first-order (main.tex:365-383):

    △_m Ψ_m = −∂_r(A[1] ∂_r Ψ_m) − ∂_z(A[2] ∂_z Ψ_m) + beta Ψ_m
              + gamma[1] ∂_r Ψ_m + gamma[2] ∂_z Ψ_m .

`region` ∈ `:v` | `:t` | `:u` (the [`vtu_domain`](@ref) vocabulary of `v`) selects the
CONSTANT slab boost H = −1 | 0 | +1 with H′ = 0 (main.tex:316-330 + 347); any other symbol
throws `ArgumentError`. This is the per-ELEMENT classification contract of the vtu assembly
(vtu.jl docstring): an element's coefficients carry its region's H everywhere INCLUDING its
interface faces — where `element_coefficients` samples 𝔄 (convention N2) and where the
pointwise piecewise boost would instead take the closed-t interface value. Inside each open
region the callables agree bitwise with [`fluxform_coefficients`](@ref) over
[`vtu_height_function`](@ref); because the slab boost is constant, no tortoise evaluation
is needed — the callables are pure algebra in (r, z), cheap at every quadrature node.

Flux-vanishing closure (main.tex:385, gated in the node tests): `A[1]` is EXACTLY zero at
the cached horizon r₊ (Δ factored, tsliced.jl) and `A[2]` exactly zero at z = ±1, so
F^r(r₊) = 0 and F^z(±1) = 0 for any finite ∂Ψ — no explicit horizon/pole boundary
conditions. The H = −1 (`:v`) coefficients are horizon-finite (main.tex:386; sealed
tsliced/W6-D gates). H is carried as an exact `Int`, so every value inherits the precision
of `(orbit, r, z)` (BigFloat orbits ⇒ BigFloat/Complex{BigFloat} values). The junction data
at the region interfaces is NOT here ([`vtu_junctions`](@ref) → node dg-jump-conditions);
the z² equatorial variant is node equatorial-symmetry-z2.
"""
function fluxform_callables(orbit::CircularOrbit, m::Integer, v::VtuSlicing, region::Symbol)
    H = region === :v ? -1 :
        region === :t ? 0 :
        region === :u ? 1 :
        throw(ArgumentError("unknown vtu region $region (expected :v, :t, or :u)"))
    (A     = ((r, z) -> flux_coeff_r(orbit, r),
              (r, z) -> flux_coeff_z(orbit, r, z)),
     beta  = (r, z) -> beta_m(orbit, m, r, z, H, 0),
     gamma = ((r, z) -> gamma_r_m(orbit, m, r, H),
              (r, z) -> gamma_z_m(orbit, m, r, z)))
end
