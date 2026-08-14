# MModeOperator/z2.jl — equatorial-symmetry z² chart of the assembled first-order flux form
# (ESSF node `equatorial-symmetry-z2`, W6-F; arXiv:2606.04998 main.tex:389-394). For
# equatorial orbits the field is even in z = cosθ, so the paper optionally formulates the
# equations in the coordinate w = z² on the HALF-domain w ∈ [0, 1] — the production lattice
# is (r, z²) (`DomainCreator...BlockBounds[z²] = [0, dz², 1]`, config.jl), selected by
# `Background.CircularOrbit.ImposeEquatorialSymmetry` (yaml:36). The printed modifications
#
#     𝓕^{z²}_m = 4z² 𝓕^z_m ,      γ^{z²}_m = 2z γ^z_m + 2(1−z²)/(r²+a²)
#
# are the chain rule ∂_z = 2z ∂_w applied to the (r, z)-chart operator on even fields:
# −∂_z(A^zz ∂_zΨ) + γ^z ∂_zΨ = −∂_w(4z²A^zz ∂_wΨ) + (2zγ^z + 2A^zz) ∂_wΨ with 2A^zz =
# 2(1−z²)/(r²+a²) — no w^{−1/2} term survives, and the angular flux coefficient vanishes
# EXACTLY at both ends (pole w = 1 as before, equator w = 0 in the even sector), so the
# half-domain needs no explicit angular boundary data either (main.tex:385 + 399-400).
#
# SCOPE BOUNDARIES (predecessors own the substance; this file is the chart CLOSURE):
#   - tsliced.jl (node elliptic-operator-tsliced) owns the z² coefficient FUNCTIONS
#     flux_coeff_z2 / gamma_z2_m (factors read from main.tex:389-394, gated pointwise vs the
#     z form at 4.3e-42);
#   - fluxform.jl part 2 (node first-order-flux-form) owns the assembled z-chart
#     `fluxform_callables` — the region → slab-H map, its validation, and the row-30
#     operator-identity/flux-vanishing closure;
#   - config.jl (node runtime-config-fixture) owns RunConfig/OrbitConfig incl. the
#     impose_equatorial_symmetry flag; vtu.jl owns the RunConfig → VtuSlicing chain.
# This node adds the two entries the domain/solve nodes consume when the flag is on: the
# (r, w)-chart bundle `fluxform_callables_z2` and the RunConfig-driven chart selector method
# of `fluxform_callables` — plus the row-31 closure gates (test_equatorial_symmetry_z2.jl):
# (r, w)-chart operator == (r, z)-chart operator on even fields under per-chart spectral
# differentiation, F^{z²}(w = 1) = 0 exact, equator regularity at w = 0, BigFloat(256).
#
# CONSUMERS (later nodes): domain-blocks-map (the (r, z²) Mesh2D pairs with this bundle),
# solver-stack-ssf / effective-source-direct (flag-consistent operator evaluation).

"""
    fluxform_callables_z2(orbit, m, v::VtuSlicing, region::Symbol)
        -> (A = (Arr, Az2), beta, gamma = (gamma_r, gamma_z2))

The equatorial-symmetry variant of [`fluxform_callables`](@ref) (main.tex:389-394): the
same assembled per-region first-order flux form, but in the HALF-domain chart
q = (r, w = z²), w ∈ [0, 1] — coefficient callables `(r, w) -> value` in the exact bundle
shape EllipticCore's `element_coefficients(em, blk, id, fc.A, fc.beta, fc.gamma)` consumes.
The operator they assemble is

    △_m Ψ_m = −∂_r(A[1] ∂_r Ψ_m) − ∂_w(A[2] ∂_w Ψ_m) + beta Ψ_m
              + gamma[1] ∂_r Ψ_m + gamma[2] ∂_w Ψ_m ,

which equals the (r, z)-chart operator on even-z fields via ∂_z = 2z ∂_w (the row-31 gate).
The radial sector never touches the angular chart: `A[1]` and `gamma[1]` are the SAME
closures as the z-chart bundle's, and `beta`'s z-dependence enters only through
1 − z² = 1 − w. The angular sector carries the printed chain-rule factors (tsliced.jl,
read from main.tex:390-393):

    A[2]     = A^{z²} = 4z² A^zz          = 4w(1−w)/(r²+a²) ,
    gamma[2] = γ^{z²} = 2zγ^z + 2A^zz     = (4|m|w + 2(1−w))/(r²+a²) ,

polynomials in w — no w^{−1/2} survives in the even sector, so w = 0 is an interior-regular
equator boundary. `A[2]` vanishes EXACTLY at BOTH chart ends (pole w = 1 and equator w = 0):
neither angular boundary of the half-domain needs explicit data (`BoundaryConditions[z²] =
None`, config.jl), mirroring F^z(±1) = 0 of the full chart (main.tex:385). `region` and the
slab boost H ∈ {−1, 0, +1} follow the sealed z-chart constructor (which also rejects unknown
regions); precision is inherited from `(orbit, r, w)` exactly as there (`sqrt` preserves
BigFloat precision).
"""
function fluxform_callables_z2(orbit::CircularOrbit, m::Integer, v::VtuSlicing,
                               region::Symbol)
    fc = fluxform_callables(orbit, m, v, region)   # sealed z-chart bundle; validates region
    (A     = (fc.A[1],                             # radial sector: chart-untouched
              (r, w) -> flux_coeff_z2(orbit, r, sqrt(w))),
     beta  = (r, w) -> fc.beta(r, sqrt(w)),        # z-dependence only via 1 − z² = 1 − w
     gamma = (fc.gamma[1],
              (r, w) -> gamma_z2_m(orbit, m, r, sqrt(w))))
end

"""
    fluxform_callables(orbit, config::RunConfig, region::Symbol)

Chart selector for the production configuration: builds the vtu slicing from the config's
transition radii ([`vtu_slicing`](@ref)) and returns, for `config.orbit.m_mode_number`, the
(r, z²)-chart bundle [`fluxform_callables_z2`](@ref) when
`config.orbit.impose_equatorial_symmetry` is set (`ImposeEquatorialSymmetry`, yaml:36 — the
production default, matching the (r, z²) block lattice of the yaml domain) and the
(r, z)-chart [`fluxform_callables`](@ref) otherwise. Both bundles have the same shape; the
domain/solve nodes pair the selected bundle with the matching mesh chart.
"""
function fluxform_callables(orbit::CircularOrbit, config::RunConfig, region::Symbol)
    v = vtu_slicing(orbit, config)
    m = config.orbit.m_mode_number
    config.orbit.impose_equatorial_symmetry ?
        fluxform_callables_z2(orbit, m, v, region) :
        fluxform_callables(orbit, m, v, region)
end
