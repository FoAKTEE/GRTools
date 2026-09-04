"""
    ScalarSelfForce

ScalarSelfForce (GrToolsScalarSelfForce) — the m-mode ELLIPTIC Kerr scalar self-force
capability: the GrTools reproduction of arXiv:2606.04998 (Vu et al., "Self-force
calculations with numerical relativity methods"). It owns the PHYSICS of the elliptic
formulation and consumes the physics-free hp-DG solver stack from `core/EllipticCore`.

Dependency budget (binding, enforced by `test/hygiene.jl` §7.7 gates): EXACTLY
`{LinearAlgebra, Printf, Distributed, NumericsCore, BigSpecial, Spacetime, EllipticCore,
Data}` — `Distributed` is the stdlib IPC of the row-47 m-mode process fan-out, admitted by
the wave-12 mission decision (stdlib only, mirroring the SparseArrays decision).
Special functions come ONLY through BigSpecial (no direct SpecialFunctions /
HypergeometricFunctions); the EffectiveSource C oracles are baked test fixtures, NEVER
live deps (no GSL/ccall channel); GPU stacks are banned (CUDA/KernelAbstractions).

Conventions (locked package-wide):
  - Geometric units `G = c = 1`, background mass `M = 1`; Kerr spin `a ∈ (−1, 1)`,
    retrograde orbits encoded via `a < 0` with `m ≥ 0` and the ×2 real-part mode-sum
    convention.
  - Mode ansatz `Φ = Σ_m Ψ_m(r, z) · sin^{|m|}θ / r · e^{i m (φ̄ − Ω t)}` with `z = cosθ`
    and the horizon-regular azimuthal chart `φ̄ = φ + a/(r₊−r₋)·ln((r−r₊)/(r−r₋))`
    (main.tex:223-236).
  - Radial chart: horizon-penetrating Boyer–Lindquist `r` via
    `∂_{r_*} = Δ/(r²+a²) ∂_r` (main.tex:354-358); slicing is vtu (v→t→u null/Boyer
    slabs) with piecewise boost `H ∈ {−1, 0, +1}` and the two junction conditions
    `[∂_{r_*}Ψ] = −imΩΨ` (main.tex:337-343).
  - Equatorial symmetry via the `z²` variant on `z ∈ [0, 1]` (main.tex:389-394).

Submodules (each = one elliptic-ssf GLOBAL_DAG node, TDD, added as it lands):
  - `Orbits/`          — `CircularOrbit{T}` thin adapters over Spacetime (E, L, Ω, u^t,
                         u^φ, r_ISCO, r_*, φ̄ offset).
  - `MModeOperator/`   — ansatz → t-sliced △̃_m → hyperboloidal → vtu jumps →
                         horizon-penetrating first-order flux form (F^i, β_m, γ^i_m) →
                         z² variant → Bayliss–Turkel outer BC.
  - `Puncture/`        — 4th-order Detweiler–Whiting Φ^P (A_ijk), ₃P/2F1 m-mode
                         integrals (Bourg Eq. 115 over `BigSpecial.legendre_P3`),
                         analytic derivatives.
  - `EffectiveSource/` — S^eff_m direct evaluation + near-worldline series switch.
  - `Worldtube/`       — regularized-region partition, Ψ^P mortar jump data,
                         fixed-source assembly (operator-on-zero trick).
  - `Domain/`          — `SSFConfig` (yaml-twin runtime config) + EllipticCore `Mesh2D`
                         block instantiation.
  - `Solve/`           — BC glue, single-mode solver-stack driver, hp-AMR loop,
                         process-parallel m-mode fan-out.
  - `Extraction/`      — F_r^m extraction, Heffernan regularization, mode-sum + tail.

Type policy: type-generic `{T}` throughout so BigFloat precision is inherited from
NumericsCore/BigSpecial (never hard-Float64'd); the production path is ComplexF64.
"""
module ScalarSelfForce

using LinearAlgebra
using Printf
using Distributed          # stdlib; the row-47 process fan-out (see Solve/parallel.jl)
using NumericsCore
using BigSpecial
using Spacetime
using EllipticCore
using Data

# Submodule files are included as each elliptic-ssf GLOBAL_DAG node is TDD-built and
# grounded against its oracle (baked EffectiveSource-C fixtures, spectre-verified
# formulas, Warburton–Barack tables). Each include + its exports land with the node.

include("Orbits/orbits.jl")                  # orbits-adapters
export CircularOrbit

include("MModeOperator/ansatz.jl")           # mode-ansatz-decomposition
# Locked conventions (see the mode_term docstring for the normative statement):
# z = cosθ chart, sin^{|m|}θ/r factoring, e^{+im(φ̄−Ωt)} phase, e^{−imΩh} hyperboloidal
# factor, and Ψ_{−m} = conj(Ψ_m) for real Φ (guards the spectre doxygen sign trap).
export sintheta_of_z, sinm_theta, sinm_theta_dz, sinm_theta_d2z
export kerr_varphi, mode_term, mode_amplitude, reconstruct_field
export psi_from_tilde, tilde_from_psi

include("MModeOperator/tsliced.jl")          # elliptic-operator-tsliced
# First-order flux-form coefficient layer (eq:flux + eq:first-order, main.tex:365-383; z²
# equatorial variant main.tex:389-394): A^rr, A^zz, A^{z²} + β_m, γ^r_m, γ^z_m, γ^{z²}_m with
# caller-supplied boost (H, H′ = dH/dr_*) — t-slicing = the (0, 0) defaults; the vtu piecewise
# assembly + junction conditions are a later node. |m| reading of the angular-origin terms
# (coefficient_{−m} = conj(coefficient_m)); Δ factored off cached horizons (exact zero at r₊).
export flux_coeff_r, flux_coeff_z, flux_coeff_z2
export beta_m, gamma_r_m, gamma_z_m, gamma_z2_m

include("MModeOperator/hyperboloidal.jl")    # hyperboloidal-transformation
# Slicing kinematics (main.tex:255-262 + 267-314): Kerr tortoise r_*(r) (wraps
# Spacetime.tortoise_kerr; dr_*/dr = 1/A^rr) and the HeightFunction carrier
# (h, H = dh/dr_*, H′ = dH/dr_*) of r_*, evaluated at r via slicing_kinematics. Composes
# with the ansatz.jl phase map and the tsliced.jl (H, H′) coefficient slots; closure
# oracle = the conjugation identity △_m[Ψ̃e^{−imΩh}] = e^{−imΩh}△̃_m[Ψ̃]. t_slicing() is
# the exact-zero t-domain instance; the vtu piecewise instance is a later node.
export rstar, HeightFunction, t_slicing, slicing_kinematics

include("Puncture/field4d.jl")               # puncture-field-4d
# 4D Detweiler-Whiting puncture: A_ijk transcribed from kerr-circular.c
# effsource_set_particle (generated field4d_coeffs.jl); Φ^P = R⁻⁷ Σ A_ijk Δr^i Δθ^j Q^k.
export Puncture, puncture_phi

include("Puncture/mmode_integrals.jl")       # mmode-puncture-integrals
# m-mode puncture via the CRP ₃P route (main.tex:459-486; Bourg Eq. 115 over
# BigSpecial.legendre_P3/pochhammer): mode_integral_Rn = ∫R^n e^{−imΔφ}dΔφ closed form,
# puncture_phi_m = the BARE azimuthal integral of Φ^P (the C effsource_PhiS_m up to its
# e^{−imφ_p} rotation), puncture_psi_m = THE locked-convention mode function Ψ^P_m(r, z)
# (sin^{|m|}θ/r factoring included via mode_amplitude; Ψ^P_{−m} = conj(Ψ^P_m)).
export mode_integral_Rn, puncture_phi_m, puncture_psi_m

include("Puncture/derivatives.jl")           # puncture-derivatives
# Analytic first derivatives of the m-mode puncture (main.tex:481-487):
# mode_integral_Rn_drho2 = exact power reduction ∂J_n/∂ϱ² = (n/2)J_{n−2} (extends n to −9),
# puncture_dphi_m = bare (I_m, ∂I_m/∂Δr, ∂I_m/∂Δθ) bundle (C effsource_calc_m dPhiS slots
# up to the phase rotation), puncture_dpsi_m = Ψ-level bundle with the trivial −imΩ/+im
# t/φ̄ factors on the same API.
export mode_integral_Rn_drho2, puncture_dphi_m, puncture_dpsi_m

include("Puncture/second_derivatives.jl")    # effective-source-direct (additive puncture layer)
# Analytic SECOND derivatives of the m-mode puncture (main.tex:481-490), the effective
# source's single puncture call: mode_integral_Rn_d2rho2 = the twice-applied power reduction
# ∂²J_n/∂(ϱ²)² = (n/2)((n−2)/2)J_{n−4} (extends mode_integral_Rn to n = −11), and
# puncture_d2phi_m = the FUSED bundle (I_m, ∂I_m/∂Δr, ∂I_m/∂Δθ, ∂²I_m/∂Δr², ∂²I_m/∂Δθ²) over
# ONE J ladder. No MIXED second derivative exists: the first-order flux form is diagonal in
# (r, z) — the same five slots the C oracle records as NaN. The c_j ladder identity
# ∂c_j/∂ϱ² = −(j+1)c_{j+1} keeps `_puncture_cs` the only home of the binomial algebra.
export mode_integral_Rn_d2rho2, puncture_d2phi_m

include("Domain/config.jl")                  # runtime-config-fixture
# Typed twin of the CRP production input anc/ScalarSelfForce.yaml (baked + sha256-pinned at
# test/fixtures/runtime_config/): crp_paper_config(; spin, orbital_radius, m_mode, …) = the
# yaml's Jinja2 parameters (defaults = the header production example), reproducing every
# derived number (v→t = r0 − 1.5wr, t→u = r0 + 3wr, block bounds, L+1 source pre-refinement,
# m = 0 no-outer-BC branch). CRP_PAPER_MMODES = 0:20 is the paper's per-mode run sweep.
export RunConfig, crp_paper_config, SommerfeldBC, BlockRefinement, CRP_PAPER_MMODES

include("MModeOperator/vtu.jl")              # vtu-slicing-jumps
# vtu piecewise slicing (main.tex:315-347; after config.jl — reads the RunConfig/OrbitConfig
# transition radii): VtuSlicing (r^v/r^u in both BL and tortoise charts), vtu_height_function
# = the piecewise HeightFunction instance (boost H ∈ {−1, 0, +1}, continuous h, H′ ≡ 0
# pointwise — the interface deltas live in the junction data), vtu_domain classification, and
# the two junction conditions [∂_{r_*}Ψ] = −imΩΨ (eq:jump_vt/eq:jump_tu) as pure VtuJunction
# data for the DG jump hooks (node dg-jump-conditions; dg_operator.jl N6 contract).
export VtuSlicing, vtu_slicing, vtu_height_function, vtu_domain, VtuJunction, vtu_junctions

include("MModeOperator/fluxform.jl")         # horizon-penetrating-compactification + first-order-flux-form
# Chart-converted coefficient bundle in the DISCRETIZATION chart q = (r, z) (main.tex:354-360):
# ∂_{r_*} = Δ/(r²+a²) ∂_r fully absorbed — fluxform_coefficients composes one slicing_kinematics
# evaluation into the tsliced.jl slots, the pointwise (Arr, Azz, beta, gamma_r, gamma_z) the DG
# quadrature consumes. Horizon-penetration closure gates (H = −1 finiteness at r₊(1+1e-12),
# u-domain 1/r-compaction smoothness at H = +1) live in test_fluxform.jl. Part 2
# (first-order-flux-form, W6-E): fluxform_callables — the assembled per-region form, constant
# vtu slab boost H ∈ {−1, 0, +1} baked into (A, beta, gamma) callables over (r, z), the exact
# shape EllipticCore's element_coefficients consumes; closure gates (F^r(r₊) = 0 / F^z(±1) = 0
# exact, spectral −∂F + βΨ + γ∂Ψ == △_mΨ per region) live in test_first_order_flux_form.jl.
export fluxform_coefficients
export fluxform_callables

include("MModeOperator/z2.jl")               # equatorial-symmetry-z2
# z² equatorial-symmetry chart (main.tex:389-394): fluxform_callables_z2 — the assembled
# per-region bundle in the half-domain chart q = (r, w = z²) with the printed chain-rule
# factors 𝓕^{z²} = 4z²𝓕^z, γ^{z²} = 2zγ^z + 2(1−z²)/(r²+a²) (coefficient functions sealed in
# tsliced.jl) — plus the RunConfig chart-selector method of fluxform_callables driven by
# impose_equatorial_symmetry (yaml:36). Closure gates ((r, w)-operator == (r, z)-operator on
# even fields per-chart spectrally, A^{z²} exact zeros at w ∈ {0, 1}, equator regularity)
# live in test_equatorial_symmetry_z2.jl.
export fluxform_callables_z2

include("MModeOperator/bayliss_turkel.jl")   # bayliss-turkel-bc
# Second-order Bayliss–Turkel outer radiation condition (main.tex:403-411 + DG realization
# main.tex:614-618): bt_tsliced_residual = the printed t-form (∂_{r_*} − imΩ)², u-sliced
# equivalence ∂²_{r_*}Ψ = 0 documented + gated through the replacement rules;
# bt_robin_coefficient = the pointwise (n·F)_N = robin·Ψ data (all four spectre Sommerfeld
# branches, hard-regression-pinned); BaylissTurkelGhost = the EllipticCore GhostBoundary
# policy (linear-homogeneous: linearized == full, in the solver matrix);
# outer_boundary_ghost = the RunConfig entry honoring the m = 0 no-outer-BC yaml branch.
export bt_tsliced_residual, bt_robin_coefficient, BaylissTurkelGhost, outer_boundary_ghost

include("Extraction/regularization.jl")      # heffernan-regularization
# Heffernan m-mode regularization parameters of the RADIAL self-force for circular equatorial
# Kerr orbits (arXiv:1211.6446 KerrRegularization.tex:1675-1731; CRP main.tex:520-522 + 723):
# reg_param_Fr(order, orbit) = the m-independent circular-Kerr coefficients F_{r[2]}/F_{r[4]}
# (purely algebraic — the m-mode parameters need NO elliptic integrals, unlike the ℓ-mode ones
# of the same paper's §IV, which share the symbol names); mmode_reg_factor(_over_pi) = the
# printed m-dependence for orders 2, 4, 6, 8, exact in Rational{BigInt}; mmode_reg_Fr = the
# assembled per-m SUBTRACTION term F^m_{r[order]}; mmode_reg_factor_tail(_over_pi) /
# mmode_reg_tail_Fr = the closed-form analytic TAIL Σ_{m≥mstart}, the completion node
# mode-sum-total-force adds to its truncated numerical sum. FOLDED m ≥ 0 ×2 real-part
# convention (see the mmode_reg_factor_over_pi docstring, which cross-references the normative
# mode_term statement) — a < 0 is retrograde, and the closed form at m = 0 carries the ×2
# unconditionally. Extraction of F_r^m and the sum itself are LATER nodes.
export HEFFERNAN_MMODE_ORDERS, HEFFERNAN_CIRCULAR_FR_ORDERS
export reg_param_Fr, mmode_reg_Fr, mmode_reg_tail_Fr, mmode_reg_term_Fr
export mmode_reg_factor, mmode_reg_factor_over_pi
export mmode_reg_factor_tail, mmode_reg_factor_tail_over_pi

include("Domain/blocks.jl")                  # domain-blocks-map
# The CRP block lattice as an EllipticCore Mesh2D (main.tex:541-554): block_domain(orbit,
# config) instantiates the 5 radial × 2 angular (r, z²) production lattice — or 5 × 3 in the
# full (r, z) chart — from the sealed RunConfig, with LinearMap1D everywhere except the outer
# wavezone column's 1/r InverseRadialMap1D (eq:compaction), the yaml block numbering
# (i + nr·j, radial index fastest — pinned by NullSlicingBlocks), the RefinedLevels initial
# h refinement (source column L+1), and one BlockConnection per interior lattice face. Carries
# the metadata the later nodes bind to: block_region (the sealed vtu_domain classification),
# worldtube_blocks (the regularized region), boundary_face_classes (:horizon/:outer/:equator/
# :pole — node dg-boundary-conditions' ghost table), and the per-block / whole-domain
# coefficient bundles over the sealed fluxform_callables. ZERO duplicated coefficient algebra.
export BlockDomain, block_domain, block_lattice_bounds, block_lattice_index, block_region
export worldtube_blocks, boundary_face_classes, boundary_class, boundary_class_rule
export block_coefficients, domain_coefficients

include("EffectiveSource/direct.jl")         # effective-source-direct
# The DIRECT pointwise m-mode effective source (main.tex:414-460 + 481-490):
# puncture_psi_chart = the CONVENTION TRANSFORMATION CHAIN — the chart-complete Ψ^P_m and its
# (r, z) first/second derivatives, including the FIELD-POINT horizon-regular rotation
# e^{−im(φ_p + χ(r))} that a constant `varphi_p` cannot carry (χ′ = a/Δ, spectre
# CircularOrbit.cpp:322-341); effective_source_m = S^eff_m = −△_mΨ^P_m over the SEALED
# t-region (H = 0) flux-form coefficients — the DG fixed source of △_mΨ^ℛ_m = S^eff_m;
# effective_source_bundle additionally exposes the worldtube's (Ψ^P, ∂Ψ^P) data and the
# near-worldline roundoff-loss DIAGNOSTIC (term_scale, cancellation ≈ the λ^{-2} amplification
# of main.tex:488-490) that node `effective-source-series` will switch on — data, not a switch.
export puncture_psi_chart, effective_source_m, effective_source_bundle

include("EffectiveSource/series_coeffs.jl")  # effective-source-series (GENERATED, see the
include("EffectiveSource/series.jl")         # generator header) + effective-source-series
# The NEAR-WORLDLINE SERIES and the route SWITCH (main.tex:488-504):
# effective_source_series_coefficients = the paper's eleven coefficients S^[0], S^[1]_1,
# S^[2]_{2,0,2L,0L,4,6,8,10,12} (mechanically transcribed from kerr-circular-comoving.c
# `effsource_calc_m_series`, with the half-integer harmonic number H_{−m−1/2} completing
# S^[2]_0 / S^[2]_2 exactly); effective_source_series_m = S^eff_m from that series in the
# SAME convention as effective_source_m, accurate through O(λ²) with an O(λ³ log ϱ)
# remainder and FINITE on the worldline circle; effective_source_switched_bundle /
# effective_source_m_switched = the combined evaluator that trades the series truncation
# against the direct route's eps·cancellation roundoff — THE evaluator the DG fixed-source
# assembly calls. The criterion (and its precision-aware eps^{1/5} crossover) is stated once,
# in the effective_source_switched_bundle docstring.
export effective_source_series_coefficients, effective_source_series_m
export effective_source_switched_bundle, effective_source_m_switched

include("Worldtube/worldtube.jl")            # worldtube-formulation
# The worldtube formulation (main.tex:505-523 + 596-608) as PURE DATA for the DG layer:
# worldtube_data(dom, orbit, config) assembles, in one pass over the sealed
# effective_source_bundle, (i) the regularized-region partition of a BlockDomain
# (inside_worldtube / worldtube_elements / worldtube_boundary_faces — the interior faces
# between an inside and an outside element, with their mortar geometry from
# worldtube_face_points), (ii) the eq:jump_worldtube data ±Ψ^P / ∓n·𝓕(Ψ^P) as one FaceState
# per mortar (worldtube_jump_data; worldtube_jump_signs is the ONE normative statement of the
# sender-normal sign table, spectre Equations.cpp:80-91), packable as EllipticCore
# FixedDataHooks (worldtube_hooks — trial-independent, so they reach the RHS through
# operator_on_zero), and (iii) the fixed source +effective_source_m inside / exactly zero
# outside in the nodal and the massive M·S form (worldtube_source /
# worldtube_massive_source). worldtube_jump_values is the pointwise kernel, carrying the
# printed z² chain rule 𝓕^{z²} = A^{z²}∂_wΨ, ∂_w = ∂_z/(2z). NO solve, NO operator, NO vtu
# junction wiring — node dg-jump-conditions binds this into the hooks table.
export WorldtubeData, worldtube_data, inside_worldtube, worldtube_elements
export worldtube_boundary_faces, worldtube_face_points, worldtube_jump_signs
export worldtube_jump_values, worldtube_jump_data, worldtube_hooks
export worldtube_source, worldtube_massive_source

include("Solve/bc_glue.jl")                  # dg-boundary-conditions
# The physical boundary conditions of the CRP domain bound to the EllipticCore ghost seam
# (main.tex:396-401 + 612-618): boundary_conditions(dom, orbit, config; mesh) is the
# MeshOperator `bcs` table — the sealed order-2 Bayliss–Turkel ghost on the `:outer` faces
# for m ≠ 0 (Robin data evaluated at each face's own quadrature points) and NoneGhost
# everywhere else. NoneGhost is the port of spectre's `None` boundary condition (empty
# apply ⇒ the convention-N7 transmissive exterior FaceState(ψ_int, −φ_int)): the numerical
# flux reduces to the interior flux bitwise, so the face imposes nothing — which is what
# the vanishing perpendicular flux coefficient at the horizon/equator/poles demands, and
# what the yaml's m = 0 outer branch selects. Callable per LEVEL mesh (the multigrid
# build_op contract): the face class is recomputed from the lattice perimeter.
export NoneGhost, boundary_conditions

include("Solve/jump_glue.jl")                # dg-jump-conditions
# The two INTERIOR jump conditions of main.tex:596-608 bound onto the EllipticCore
# ModifyBoundaryData seam: jump_hooks(dom, wt, orbit, config, m; mesh) is the MeshOperator
# `hooks` table — one hook PER MORTAR (the receiver's mortar order), which is what the
# production lattice's 2:1 worldtube faces need. The sealed worldtube_jump_data become
# FixedDataHooks (trial-independent ⇒ they leave the matrix untouched and reach the RHS
# through operator_on_zero), while vtu_junction_faces tags every interior face that crosses
# the vtu slicing (spectre's NullSlicingBlocks rule over the sealed block_region, cross-checked
# against the transition radii) with its sealed VtuJunction, and vtu_junction_hook turns
# `[∂_{r_*}Ψ] = −imΩΨ` into a LinearizedHook — the term that IS in the GMRES matrix. A face
# carrying both composes as ComposedHook(linearized, fixed) in the sealed W5-B order. NO
# chart conversion: A^rr = dr/dr_* already makes n·𝓕 the tortoise derivative. Callable per
# LEVEL mesh (multigrid/AMR); the mesh-bound worldtube data must be rebuilt by the caller.
export vtu_junction_faces, vtu_junction_hook, jump_hooks

include("Solve/assemble.jl")                 # fig-matrix (the row-42 assembly contract)
# The ONE composition of the sealed CRP pieces into the production DG operator:
# assemble_operator(T, dom, wt, orbit, config, m; mesh) = MeshOperator over
# domain_coefficients + jump_hooks + boundary_conditions + config.dg.penalty_parameter,
# BITWISE equal to the hand-written spelling and callable on ANY mesh of the block lattice.
# operator_builder is the MultigridHierarchy `build_op` factory over a caller-supplied
# per-level WorldtubeData factory (rebuilding that bundle stays worldtube_data's job), and
# coarsest_operator walks coarsen_mesh to the multigrid COARSEST level — by default the yaml
# twin of InitialCoarseLevels: L (L + 1 levels in total), which at the production L = 0 is the
# initial 20-element / 320-dof lattice arXiv:2606.04998 Fig. 2 plots (main.tex:622-626).
# linearized_rhs folds the worldtube fixed source with the operator's own fixed-data lift,
# b = worldtube_massive_source(wt) − operator_on_zero(op), extending EllipticCore's general
# method. NO solve, NO multigrid cycle, NO AMR — rows 42/43 own those.
export assemble_operator, operator_builder, coarsest_operator
export linearized_rhs

include("Solve/driver.jl")                   # solver-stack-ssf (the row-42 JOIN)
# THE FIRST REAL m-MODE SOLVE: solve_mode(orbit, config, m) builds the block domain, the
# mesh-bound worldtube bundle, the operator and the SolverStack (GMRES o multigrid V-cycle o
# 2 additive-Schwarz passes o cached bottom LU) over operator_builder's per-level factory,
# assembles the PHYSICAL right-hand side and solves - returning a ModeSolution holding
# (dom, mesh, wt, op, stack, b, u, monitor, record, result) with u = Psi^R inside the
# worldtube / Psi outside. solver_options is the field-by-field SolverStackOptions twin of
# the sealed RunConfig (the yaml solver block) with the two multigrid LEVEL-POLICY knobs
# (coarse_levels, coarse_orders) exposed for the AMR driver. mode_fixed_source is the VOLUME
# source: the SWITCHED evaluator inside the tube, exactly zero outside, massive M.S form -
# the worldtube MORTAR data stay on the direct bundle (the one normative statement of that
# rule lives in its docstring). Read-out: mode_residual (the DG residual against the FULL
# table), reconstruct_full_field (Psi everywhere, Psi^R + Psi^P inside) and field_at
# (Lagrange interpolation on the containing element, one-sided via `element`).
# rebuild_for_mesh! is the ONE-STEP per-level rebuild an hp-AMR sweep calls; THE LOOP,
# its criterion and its stopping rule are node amr-driver (row 43). NO extraction, NO figure.
export solver_options, ModeSolution, solve_mode, solve_mode!, rebuild_for_mesh!
export mode_fixed_source, mode_residual, reconstruct_full_field, field_at

include("Solve/amr.jl")                      # amr-driver (the row-43 production loop)
# THE hp-AMR SWEEP: solve_mode_amr(orbit, config, m) runs the paper's Amr.Iterations levels of
# "h-refine the elements containing the puncture, p-refine everywhere else" over ONE
# ModeSolution (one stack, one epoch, one pair of caches, each level warm-started from the
# transferred previous solution) and returns an AmrHistory - per level the mesh, the puncture
# elements, the SolveRecord (iterations, residual staircase, build_time, cache counters) and
# the particle diagnostics Psi^R(x_p), d_r Psi^R(x_p) - plus the final ModeSolution.
# amr_flags is the CRITERION: spectre's RefineAtPuncture (block-logical snap + CLOSED element
# containment, so a puncture on a corner flags every element that touches it) merged with
# IncreaseResolution, under the yaml Limits; the one normative statement of the criterion and
# its snapping rule lives in its docstring. amr_level_options is the multigrid LEVEL POLICY
# (InitialCoarseLevels at level 0 - the direct bottom solve - then coarsening to the fixed
# block structure so the bottom LU stays a cache hit, with shrink = 0.5 so a graded hp cascade
# reaches it in a handful of levels). ScalarSelfForce.amr_step! is one level (not exported).
# The particle diagnostics are read through the SEALED row-42 field_at / field_grad_at - this
# node adds no read-out of its own - and resolution_error / relative_resolution_error are THE
# one definition of the paper's numerical-resolution error (main.tex:730) that rows 44/50/53
# share: the difference between the highest and second-highest AMR level.
# NO self-force extraction (row 44), NO process fan-out (row 47), NO figure.
export puncture_chart_position, amr_flags, amr_level_options
export solve_mode_amr, AmrLevel, AmrHistory
export resolution_error, relative_resolution_error

include("Extraction/selfforce.jl")           # mode-self-force-extraction (row 44)
# F_r^m AT THE PARTICLE, the number the paper is about (Osburn-Nishimura arXiv:2206.07031
# Eq. 3.10, main.tex:178-184; CRP main.tex:428-430 + 701-703): mode_self_force assembles
# F_alpha = q d_alpha Phi^R from Psi^R_m and its chart derivatives at x_p in the locked ansatz
# conventions - UNFOLDED at m = 0, and for m > 0 the x2 real part of e^{+im(varphi_p - Omega t)}
# times [d_r Psi/r - Psi/r^2 + i m (a/Delta) Psi/r], where the last term is the DEPHASING
# derivative that d_r picks up from the horizon-regular azimuth (spectre adds it by hand at
# ObserveSelfForce.cpp:70-72). F_t = -Omega F_phi per mode; F_theta is a BITWISE zero in the
# production z^2 chart. Ingredients come from the SEALED field_at / field_grad_at or the
# sealed AmrLevel diagnostics - no stencil is re-derived. mode_self_force_history and
# mode_force_convergence are the per-AMR-level sequence and the resolution floor
# (main.tex:719-724); regularized_mode_force subtracts the Heffernan F^m_{r[order]} with the
# sealed m = 0 HALVING rule and is the per-m term the mode sum adds up. THE SUM ITSELF is node
# mode-sum-total-force (row 46); NO figure, NO fan-out.
export mode_self_force, mode_self_force_history, mode_force_convergence
export regularized_mode_force

include("Solve/parallel.jl")                 # mmode-parallel-driver (the row-47 fan-out)
# THE m-MODE FAN-OUT: solve_modes_parallel(orbit, config, ms = CRP_PAPER_MMODES) runs the
# paper's independent modes (main.tex:716) as ONE PROCESS PER MODE - each worker rebuilds its
# own mode_config, runs the sealed solve_mode_amr with its own stack and its own two caches,
# pins BLAS explicitly, and ships back a ModeResult (every AmrLevel + the final sol.u; the
# live ModeSolution is process-local by row 42's rule and deliberately does not travel).
# solve_modes_serial is the same loop in-process and is the BITWISE ORACLE of the fan-out.
# THE COMPUTE POLICY of the compute-usage skill is stated once, normatively, in the
# solve_modes_parallel docstring and enforced mechanically: mode_pool_size caps
# nworkers x max(threads_per_worker, blas_threads) at PARALLEL_CPU_BUDGET = 24 (a bigger
# request is clamped AND logged; a CPU spent inside OpenBLAS is a CPU), system_load reads the
# LIVE whole-system usage so the pre-launch guard can REFUSE with the numbers when
# live + planned would reach PARALLEL_LOAD_GUARD = 0.90, and a monitor task samples the machine
# into the run record (and a log file) for the whole run. THE GUARD GOVERNS LAUNCH ONLY:
# run.max_load is a record of what the shared box then did, warned but never enforced.
# mode_pool is the REUSABLE pool the figure drivers need: a fresh process pays ~200 s of JIT
# before its first solve, so several sweeps must share one pool; isopen asks the cluster
# whether the workers are alive, a dead worker is reported by id and pruned, and timeout_s
# bounds a hung sweep by removing its workers. mode_config is the per-m RunConfig (only
# MModeNumber and the m = 0 no-outer-BC yaml branch move) and mode_stamp is the per-mode record
# of the knobs used (the shift is 0 for every m: main.tex:680 / yaml:194,203).
# NO extraction, NO figure, NO GPU.
export PARALLEL_CPU_BUDGET, PARALLEL_MAX_THREADS_PER_WORKER, PARALLEL_LOAD_GUARD
export PARALLEL_LOAD_PROBES
export LoadSample, system_load, mode_pool_size, mode_config
export ModeStamp, mode_stamp
export ModePool, mode_pool
export ModeResult, ModeFailure, ModeRunRecord, ModeSet
export solve_modes_serial, solve_modes_parallel

end # module
