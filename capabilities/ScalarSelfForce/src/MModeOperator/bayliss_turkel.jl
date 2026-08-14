# MModeOperator/bayliss_turkel.jl — the second-order Bayliss–Turkel outer radiation
# condition (ESSF node `bayliss-turkel-bc`, W7-D, GLOBAL_DAG row 32; arXiv:2606.04998
# Sec. II.B.5 main.tex:396-412 + Sec. III.B.3 main.tex:612-618, following Osburn:2022bby).
# The outer radial boundary sits at a large but finite R (production R = 1e4, eq:compaction)
# instead of null infinity; there the paper imposes, in t-sliced form (main.tex:403-407),
#
#     ∂²_{r_*}Ψ̃_m − 2imΩ ∂_{r_*}Ψ̃_m − m²Ω²Ψ̃_m = 0   ⟺   (∂_{r_*} − imΩ)² Ψ̃_m = 0 ,
#
# the annihilator of the two leading outgoing multipoles (kernel (c₁ + c₂r_*)e^{imΩr_*};
# reflection error O(1/R²) relative — gate (b) of the node tests). U-SLICED EQUIVALENCE
# (main.tex:408-411, eq:bayliss_turkel_bc — documented API knowledge of this node): pushing
# the t-form through the hyperboloidal replacement rules ∂_{r_*} → ∂_{r_*} + imΩH,
# ∂²_{r_*} → ∂²_{r_*} + 2imΩH∂_{r_*} + imΩH′ − m²Ω²H² (main.tex:306-308, sealed
# hyperboloidal-transformation node) with the u-slab boost H = +1, H′ = 0 collapses it
# IDENTICALLY to
#
#     ∂²_{r_*}Ψ_m = 0 ,
#
# because the u-sliced phase map Ψ = Ψ̃e^{−imΩh} with dh/dr_* = 1 conjugates
# (∂_{r_*} − imΩ) into plain ∂_{r_*} (algebraically gated in test (c); the H = −1 v-domain
# counterpart is (∂_{r_*} − 2imΩ)², the horizon-side ingoing annihilator).
#
# DG REALIZATION (main.tex:614-618, as detailed in Vu:2024): ∂²_{r_*}Ψ = 0 is equivalent to
# ∂_r F^r = 0 (F^r = ∂_{r_*}Ψ = A^rr∂_rΨ, A^rr = Δ/(r²+a²) — sealed fluxform layer), and
# substituting that into the z-reduced interior equation −∂_iF^i + β_mΨ + γ^i_m∂_iΨ = 0
# solves for the boundary flux:
#
#     (n·F)_N = −(Δ/(r²+a²))·(β_m/γ^r_m)·Ψ ,                       (main.tex:616)
#
# imposed through the numerical flux by the ghost data (n·F)_ext = 2(n·F)_N − (n·F)_int and
# Ψ_ext = Ψ_int (main.tex:618) — a Neumann-type condition whose data is LINEAR-HOMOGENEOUS
# in the trial field, so its linearized application equals the full one and the BT term
# lives inside the GMRES/Schwarz/multigrid matrix (spectre Sommerfeld.cpp:71-81,
# audit-spectre-ssf §1.4).
#
# SCOPE BOUNDARIES (predecessors own the substance; later nodes own the wiring):
#   - tsliced.jl / fluxform.jl (rows 26/29/30): β_m, γ^r_m, A^rr and the assembled
#     per-region callables — the coefficient SOURCE of the Robin data;
#   - hyperboloidal.jl / vtu.jl (rows 27/28): replacement rules, boost slabs, r_*(r);
#   - EllipticCore boundary_hooks.jl (row 14): the ghost seam ghost_state /
#     ghost_state_linearized + conventions N1-N7 — this file only ADDS the GhostBoundary
#     subtype the seam anticipated ("Bayliss-Turkel via flux ... app-side consumers");
#   - Domain/config.jl (row 7): SommerfeldBC and the m = 0 ⇒ no-outer-BC yaml branch;
#   - dg-boundary-conditions (row 41, LATER): Solve/bc_glue.jl wires the ghost into the
#     Mesh2D sweep per boundary face and runs the reduced-ODE contamination oracle.
# Oracle: spectre Sommerfeld (ref-code/spectre .../Scalar/BoundaryConditions/Sommerfeld.cpp
# :34-69, hard regression values Test_Sommerfeld.cpp:91-93 + 123-124 — gate (d)).

export bt_tsliced_residual, bt_robin_coefficient, BaylissTurkelGhost, outer_boundary_ghost

"""
    bt_tsliced_residual(orbit, m, psi, dpsi, d2psi) -> Complex

The t-sliced second-order Bayliss–Turkel operator (main.tex:403-407) applied to boundary
values of the t-sliced mode function: with `(psi, dpsi, d2psi)` = (Ψ̃_m, ∂_{r_*}Ψ̃_m,
∂²_{r_*}Ψ̃_m) at the evaluation point and k = mΩ from the orbit,

    B[Ψ̃] = ∂²_{r_*}Ψ̃ − 2ik ∂_{r_*}Ψ̃ − k²Ψ̃  =  (∂_{r_*} − ik)² Ψ̃ ,

which vanishes on the outgoing kernel (c₁ + c₂r_*)e^{ikr_*} exactly and on the outgoing
asymptotics e^{ikr_*}/rⁿ to O(1/r^{n+2}) (gate (b)); the ingoing mode e^{−ikr_*} is NOT
annihilated (residual −4k²Ψ̃ — the falsifier). U-sliced equivalence (eq:bayliss_turkel_bc,
main.tex:408-411): under the H = +1 replacement rules this operator IS ∂²_{r_*} (file
header; gate (c) pushes the probed coefficient triple through the substitution). Linear in
`(psi, dpsi, d2psi)` — probing with unit vectors recovers the printed coefficient triple
(1, −2imΩ, −m²Ω²). Type-generic (BigFloat inherited through `orbit.Omega`).
"""
function bt_tsliced_residual(orbit::CircularOrbit, m::Integer, psi, dpsi, d2psi)
    k = m * orbit.Omega
    d2psi - 2im * k * dpsi - k^2 * psi
end

"""
    bt_robin_coefficient(bc::SommerfeldBC, alpha, beta, gamma_r) -> Complex

Pointwise Robin coefficient of the Bayliss–Turkel/Sommerfeld outer boundary condition: the
imposed normal flux at one outer-boundary face point is `(n·F)_N = robin · Ψ`, with
`(alpha, beta, gamma_r)` = (A^rr = Δ/(r²+a²), β_m, γ^r_m) evaluated AT that face point
(the operator's own coefficients — sealed tsliced/fluxform layers). Faithful port of the
spectre Sommerfeld oracle (Sommerfeld.cpp:34-69; hard regression values gated in test (d)),
all four `(hyperboloidal_slicing, order)` branches:

  - hyperboloidal, order 1: `robin = 0` — plain outflow; u-slicing already propagates
    purely outgoing radiation (Sommerfeld.cpp:43).
  - hyperboloidal, order 2 (THE production branch, yaml:118-124): `robin = −αβ_m/γ^r_m` —
    the u-sliced BT condition ∂²_{r_*}Ψ = 0 ⟺ ∂_rF^r = 0 substituted into the z-reduced
    interior equation and solved for the flux (main.tex:616; Sommerfeld.cpp:45-46).
  - no slicing (t-sliced), order 1: `robin = ik` with k = mω, ω = 1/(a + √(r0³/M))
    recomputed from the BC's own orbit fields exactly as spectre does
    (Sommerfeld.cpp:54-60; Osburn:2022bby Eqs. 4.10-4.11).
  - no slicing, order 2: `robin = (k² − β_m)/(γ^r_m − 2ik)` — the t-form substituted the
    same way (Sommerfeld.cpp:61-64, transcribed AS-BUILT: spectre leaves β un-multiplied
    by α in this branch, an α → 1 large-R simplification consistent with the BC's own
    O(1/R²) accuracy; unused in production, kept verbatim for oracle fidelity).

Any other order throws `ArgumentError` (spectre ERROR). The coefficient is trial-field
INDEPENDENT — the BC data `robin·Ψ` is linear-homogeneous in Ψ, which is why the ghost
built from it is applied identically in full and linearized operator applications.
Type-generic; m = 0 production runs never call this (the yaml selects no outer BC —
[`outer_boundary_ghost`](@ref)).
"""
function bt_robin_coefficient(bc::SommerfeldBC, alpha, beta, gamma_r)
    if bc.hyperboloidal_slicing
        bc.order == 1 && return zero(im * alpha * beta)
        bc.order == 2 && return -alpha * beta / gamma_r          # main.tex:616
        throw(ArgumentError("Order $(bc.order) not implemented for Sommerfeld boundary " *
                            "condition with hyperboloidal slicing (spectre parity)"))
    end
    a = bc.black_hole_spin * bc.black_hole_mass
    omega = 1 / (a + sqrt(bc.orbital_radius^3 / bc.black_hole_mass))  # Sommerfeld.cpp:57
    k = bc.m_mode_number * omega
    bc.order == 1 && return im * k + zero(im * alpha * beta)
    bc.order == 2 && return (k^2 - beta) / (gamma_r - 2im * k)   # Sommerfeld.cpp:61-64
    throw(ArgumentError("Order $(bc.order) not implemented for Sommerfeld boundary " *
                        "condition (spectre parity)"))
end

"""
    BaylissTurkelGhost(robin::AbstractVector) <: EllipticCore.GhostBoundary

The Bayliss–Turkel outer-boundary ghost policy on ONE outer radial face — the app-side
GhostBoundary instance the EllipticCore seam anticipated (boundary_hooks.jl): `robin` holds
the pointwise [`bt_robin_coefficient`](@ref) values at the face points, and the ghost
exterior implements the DG imposition of main.tex:614-618,

    (n·F)_N   = robin ⊙ ψ_int ,
    Ψ_ext     = ψ_int ,                       (main.tex:618)
    (n·F)_ext = 2(n·F)_N − (n·F)_int          [receiver's outward normal]

which in the SENDER convention (N3: the ghost's own outward normal is −n_int) reads

    ghost_state(bc, int) = FaceState(ψ_int, φ_int − 2·robin ⊙ ψ_int) .

Through the numerical flux (any penalty): Ψ* = ψ_int (zero jump — the penalty term is
inert) and (nF)* = robin ⊙ ψ_int, the classical weak Neumann imposition with FIELD-linear
data. Because the data is linear-homogeneous in the trial field, `ghost_state_linearized`
IS `ghost_state` (unlike Dirichlet/Neumann fixed data there is nothing to zero): the BT
term belongs to the GMRES/Schwarz/multigrid matrix, exactly spectre's
`apply_linearized = apply` (Sommerfeld.cpp:71-81). Face identity and physics are captured
at CONSTRUCTION (the spectre argument_tags pattern, W5-B ledger convention) — EllipticCore
stays zero-physics.
"""
struct BaylissTurkelGhost{T<:Number} <: GhostBoundary
    robin::Vector{T}
end

EllipticCore.ghost_state(bc::BaylissTurkelGhost, int::FaceState) =
    FaceState(int.psi, int.nflux .- 2 .* (bc.robin .* int.psi))
EllipticCore.ghost_state_linearized(bc::BaylissTurkelGhost, int::FaceState) =
    ghost_state(bc, int)

"""
    outer_boundary_ghost(cfg::RunConfig, alpha, beta, gamma_r)
        -> Union{BaylissTurkelGhost, Nothing}

Build the outer radial boundary ghost of one run from its `RunConfig` and the operator
coefficient values `(A^rr, β_m, γ^r_m)` at the outer-face points — honoring the m = 0
no-outer-BC flag semantics of the production input (yaml:115-125, sealed W4-E node):

  - `m = 0` ⇒ `cfg.domain.outer_radial_bc === nothing` ⇒ returns `nothing` — the m = 0
    mode is static (k = mΩ = 0: nothing radiates) and γ^r_0 can vanish (at a = 0 it is
    identically zero for ANY boost — every γ^r term carries a or m — so the Robin data
    −αβ₀/γ^r_0 is ill-defined); NO outer condition is imposed, and the outer face gets
    whatever no-op treatment the solve layer applies (row-41 node).
  - `m ≠ 0` ⇒ the config's order-2 hyperboloidal [`SommerfeldBC`](@ref) drives
    [`bt_robin_coefficient`](@ref) pointwise → [`BaylissTurkelGhost`](@ref).

Horizon and angular boundaries need no counterpart at all: their perpendicular flux
vanishes identically (A^rr(r₊) = 0, A^zz(z² = 1) = 0 exactly — sealed row-30/31 gates;
main.tex:385 + 399-401 + 614), so both boundary terms of the numerical flux drop out.
"""
outer_boundary_ghost(cfg::RunConfig, alpha::AbstractVector, beta::AbstractVector,
                     gamma_r::AbstractVector) =
    cfg.domain.outer_radial_bc === nothing ? nothing :
    BaylissTurkelGhost(bt_robin_coefficient.(Ref(cfg.domain.outer_radial_bc),
                                             alpha, beta, gamma_r))
