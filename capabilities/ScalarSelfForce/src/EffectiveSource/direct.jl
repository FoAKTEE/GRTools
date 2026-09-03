# EffectiveSource/direct.jl — the DIRECT pointwise m-mode effective source (ESSF node
# `effective-source-direct`, W7-C; arXiv:2606.04998 Sec. III, main.tex:414-460 + 481-490).
# This file owns TWO things and nothing else:
#
#   1. the m-mode CONVENTION TRANSFORMATION CHAIN — the map from the bare azimuthal integral
#      I_m of the 4D puncture (Puncture/, in Boyer–Lindquist φ) to the mode function Ψ^P_m in
#      the LOCKED package chart (ansatz.jl: horizon-regular φ̄, sin^{|m|}θ/r factoring,
#      z = cosθ, m ≥ 0 with Ψ_{−m} = conj(Ψ_m)) — see the [`puncture_psi_chart`](@ref)
#      docstring, the single normative statement, cross-referencing `mode_term`;
#   2. the assembly S^eff_m = −△_m Ψ^P_m over the SEALED first-order flux-form coefficients
#      (tsliced.jl / fluxform.jl) with the worldtube slab boost H = 0.
#
# The near-worldline SERIES switch (main.tex:488-504) is a LATER node
# (`effective-source-series`); this file only EXPOSES the roundoff-loss diagnostic that node
# will switch on, as data ([`effective_source_bundle`](@ref) `cancellation`).
#
# WHY THE DIRECT ROUTE CANCELS (main.tex:488-490, measured in the node evidence): each of the
# seven assembled terms diverges like 1/λ² towards the worldline while S^eff_m stays finite —
# four orders of singularity cancel by construction — so the roundoff amplification grows like
# λ^{-2}: measured 1e3…2e6 on a near-puncture box (Δr ∈ [0.2, 1.2]) and 1e5…1.7e10 at
# λ ∈ [1e-2, 1e-1). Float64 therefore loses ~eps·cancellation of relative accuracy, which is
# exactly why the paper switches to the series near the worldline.

"""
    puncture_psi_chart(p::Puncture, m::Integer, r, z; phi_p = 0)
        -> (psi, dr, dz, drr, dzz)

The CHART-COMPLETE m-mode puncture bundle in the discretization chart `q = (r, z)`: the mode
function `Ψ^P_m` the elliptic solver's operator acts on, together with its first and second
`(r, z)` derivatives (no mixed derivative — the first-order flux form is diagonal).

NORMATIVE CONVENTION CHAIN (the single statement for this package; the ansatz side is locked
in `MModeOperator/ansatz.jl`, see [`mode_term`](@ref)). The 4D puncture is a function of the
BOYER–LINDQUIST offsets `Φ^P(Δr, Δθ, φ − φ_p(t))` with `φ_p(t) = Ω t + φ_p`, while the mode
integral (eq:modeIntegral, main.tex:233-236) analyzes in the HORIZON-REGULAR azimuth
`φ̄ = φ + χ(r)`, `χ = kerr_varphi_offset(a, r)`. At fixed `(t, r, θ)` the two differ by a
constant, so

    Ψ^P_m(r, z) = r e^{+imΩt}/(2π sin^{|m|}θ) ∫ Φ^P e^{−imφ̄} dφ̄
                = r e^{−im(φ_p + χ(r))} / (2π sin^{|m|}θ) · I_m(Δr, Δθ) ,

with `I_m` the bare integral of [`puncture_phi_m`](@ref) (the co-rotating puncture's
`e^{−imΩt}` cancels eq:modeIntegral's `e^{+imΩt}`, so `Ψ^P_m` is time-independent). The
rotation is evaluated at the FIELD POINT `r`, not at the particle: this is the piece a
constant `varphi_p` cannot carry, and it is what makes `∂_r` in the φ̄ chart differ from `∂_r`
in the φ chart by the exact term `−im χ′(r) Ψ^P_m`, `χ′ = a/Δ`
(`Spacetime.kerr_varphi_offset_dr`) — the same term spectre adds by hand
(CircularOrbit.cpp:338-341, `− i a m/Δ · singular_field`). Passing
`varphi_p = kerr_varphi(a, r, phi_p)` to the sealed [`puncture_psi_m`](@ref) reproduces `psi`
exactly; the sealed [`puncture_dpsi_m`](@ref) `dr` slot differs from this one by precisely
that chain term (both gated in the node tests).

Slots (all `Complex`; `Ψ^P_{−m} = conj(Ψ^P_m)` slot by slot):
  * `psi`  — `Ψ^P_m(r, z)`;
  * `dr`, `drr` — `∂_r`, `∂²_r` at fixed `z`, including the `e^{−imχ(r)}` chain
    (`χ″ = −a(2r−2M)/Δ²`, the exact derivative of the sealed `χ′ = a/Δ`);
  * `dz`, `dzz` — `∂_z`, `∂²_z` at fixed `r`, through `Δθ = −asin z` (`∂Δθ/∂z = −1/sinθ`,
    `∂²Δθ/∂z² = −z/sin³θ`) and the ansatz factor `sin^{|m|}θ` whose z-derivatives come from
    the sealed [`sinm_theta_dz`](@ref) / [`sinm_theta_d2z`](@ref) — the single ansatz home.

SLICING CONTRACT (binding on every caller): the puncture lives in the worldtube, which the
production domain places strictly inside the vtu `t` slab (`r^v < r₀ − w_r`, `r₀ + w_r < r^u`;
main.tex:316-330 + 545-548), where `H = 0` and `h = 0` — so `Ψ^P_m ≡ Ψ̃^P_m`, no `e^{−imΩh}`
factor applies, and the coefficients used downstream are exactly
`fluxform_callables(orbit, m, v, :t)`. There is deliberately no `region` argument: an
effective source evaluated with a `v` or `u` boost would not be the source of any equation
this package solves, and the worldtube never leaves the `t` slab. A caller that ever needs a
boosted region must raise it as a formulation change, not pass a flag.

DOMAIN (checked, not silently NaN): requires `ϱ > 0` (off the worldline circle — the mode
function diverges logarithmically towards it, main.tex:424) and `|z| < 1` (the z-chart
Jacobian `∂Δθ/∂z = −1/sinθ` diverges at the poles for every `m`, and for `m ≠ 0` so does the
ansatz factor). Both raise `DomainError`; node `amr-driver`'s puncture-snap can otherwise put
a collocation point exactly on the worldline and propagate a silent NaN into the fixed
source. Type-generic: BigFloat orbits/coordinates give `Complex{BigFloat}` slots.
"""
function puncture_psi_chart(p::Puncture, m::Integer, r::Real, z::Real; phi_p::Real = 0)
    o = p.orbit
    # DOMAIN guards (see the docstring): the worldline circle and the poles are genuine
    # singularities of this chart, and a silent NaN would propagate into the DG fixed source.
    abs(z) < 1 || throw(DomainError(z,
        "puncture_psi_chart: |z| < 1 required (z = cosθ; the ∂Δθ/∂z = −1/sinθ chart Jacobian, \
         and for m ≠ 0 the sin^{|m|}θ ansatz factor, diverge at the poles)"))
    dr = r - o.r0
    dtheta = -asin(z)                                    # Δθ = θ − π/2 (locked branch)
    rho2 = o.rho2_dr2 * dr^2 + o.rho2_dtheta2 * dtheta^2
    rho2 > 0 || throw(DomainError((r, z),
        "puncture_psi_chart: ϱ² > 0 required — (r, z) is ON the worldline circle (r₀, 0), \
         where the m-mode puncture diverges logarithmically (main.tex:424)"))
    b = puncture_d2phi_m(p, m, dr, dtheta)
    tp = 2 * oftype(b.phi, pi)

    # ansatz factors (ansatz.jl is the single home): N = r / (2π sin^{|m|}θ)
    s = sinm_theta(m, z)
    s1 = sinm_theta_dz(m, z)
    s2 = sinm_theta_d2z(m, z)
    sig = sintheta_of_z(z)
    N = r / (tp * s)
    N_r = one(r) / (tp * s)
    N_z = -r * s1 / (tp * s^2)
    N_zz = r * (2 * s1^2 - s * s2) / (tp * s^3)

    # z-chart of the bare integral: Δθ = −asin z
    th_z = -one(z) / sig
    th_zz = -z / sig^3
    I_z = b.dtheta * th_z
    I_zz = b.dthetatheta * th_z^2 + b.dtheta * th_zz

    # the REAL φ-chart mode function f = N·I_m and its (r, z) derivatives (N_rr ≡ 0)
    f = N * b.phi
    f_r = N_r * b.phi + N * b.dr
    f_rr = 2 * N_r * b.dr + N * b.drr
    f_z = N_z * b.phi + N * I_z
    f_zz = N_zz * b.phi + 2 * N_z * I_z + N * I_zz

    # horizon-regular chart rotation e^{−im φ̄(r; φ_p)} and its exact r-derivatives
    chi1 = kerr_varphi_offset_dr(o.a, r)                 # χ′ = a/Δ (sealed closed form)
    chi2 = -o.a * (2 * r - 2 * o.M) / delta_kerr(o.a, r)^2   # χ″, exact derivative of χ′
    E = cis(-m * kerr_varphi(o.a, r, phi_p))

    # psi through the SEALED ansatz stripper, so it is the landed puncture_psi_m at
    # varphi_p = kerr_varphi(a, r, phi_p) (single assembly home)
    psi = mode_amplitude(m, b.phi * E / tp; Omega = o.Omega, ts = zero(b.phi), r = r, z = z)

    return (psi = psi,
            dr = E * (f_r - im * m * chi1 * f),
            dz = E * f_z,
            drr = E * (f_rr - 2 * im * m * chi1 * f_r +
                       (-im * m * chi2 - m^2 * chi1^2) * f),
            dzz = E * f_zz)
end

# d/dr of the SEALED radial flux coefficient A^rr = Δ/(r²+a²) (tsliced.jl `flux_coeff_r`) —
# the one coefficient derivative the divergence expansion ∂_r(A^rr ∂_rΨ) needs and the
# operator layer does not own. Δ from the orbit's cached horizons, exactly as there.
function _flux_coeff_r_dr(orbit::CircularOrbit, r)
    w = r^2 + orbit.a^2
    D = (r - orbit.rp) * (r - orbit.rm)
    return ((2 * r - 2 * orbit.M) * w - 2 * r * D) / w^2
end

# d/dz of the SEALED angular flux coefficient A^zz = (1−z²)/(r²+a²) (`flux_coeff_z`).
_flux_coeff_z_dz(orbit::CircularOrbit, r, z) = -2 * z / (r^2 + orbit.a^2)

"""
    effective_source_bundle(p::Puncture, m::Integer, r, z; phi_p = 0)
        -> (source, psi, dpsi_r, dpsi_z, term_scale, cancellation)

The m-mode effective source with its roundoff diagnostic and the puncture data the worldtube
assembly needs, from ONE puncture evaluation.

  * `source` — `S^eff_m(r, z)`, see [`effective_source_m`](@ref);
  * `psi`, `dpsi_r`, `dpsi_z` — the [`puncture_psi_chart`](@ref) slots (the worldtube mortar
    data of node `worldtube-formulation` is `±Ψ^P_m` and `∓n·𝓕^i(Ψ^P_m)` with
    `𝓕^r = flux_coeff_r(orbit, r)·dpsi_r`, `𝓕^z = flux_coeff_z(orbit, r, z)·dpsi_z`);
  * `term_scale` — `max_i |t_i|` over the seven assembled terms `t_i` of `−△_mΨ^P_m`;
  * `cancellation` — `term_scale / |source|`, the ROUNDOFF-LOSS DIAGNOSTIC of
    main.tex:488-490: four orders of singularity cancel in the assembly, so the achievable
    relative accuracy of `source` in a type `T` is `≈ eps(T)·cancellation` (measured
    prediction ratio ≤ 2e2 against BigFloat(256) truth). It grows like `λ^{-2}` towards the
    worldline (measured ~1e4 per two decades of λ). This is DATA, not a switch: the node
    `effective-source-series` owns the criterion that trades this against the near-worldline
    power series. `Inf` where `source` vanishes exactly.

Allocation-free on the Float64 path.
"""
function effective_source_bundle(p::Puncture, m::Integer, r::Real, z::Real;
                                 phi_p::Real = 0)
    o = p.orbit
    b = puncture_psi_chart(p, m, r, z; phi_p = phi_p)
    # −△_mΨ = ∂_r(A^rr ∂_rΨ) + ∂_z(A^zz ∂_zΨ) − βΨ − γ^r ∂_rΨ − γ^z ∂_zΨ, coefficients from
    # the sealed t-region (H = 0) flux form — identical to fluxform_callables(orbit, m, ·, :t)
    t1 = flux_coeff_r(o, r) * b.drr
    t2 = _flux_coeff_r_dr(o, r) * b.dr
    t3 = flux_coeff_z(o, r, z) * b.dzz
    t4 = _flux_coeff_z_dz(o, r, z) * b.dz
    t5 = -beta_m(o, m, r, z) * b.psi
    t6 = -gamma_r_m(o, m, r) * b.dr
    t7 = -gamma_z_m(o, m, r, z) * b.dz
    src = t1 + t2 + t3 + t4 + t5 + t6 + t7
    scale = max(abs(t1), abs(t2), abs(t3), abs(t4), abs(t5), abs(t6), abs(t7))
    return (source = src, psi = b.psi, dpsi_r = b.dr, dpsi_z = b.dz,
            term_scale = scale, cancellation = scale / abs(src))
end

"""
    effective_source_m(p::Puncture, m::Integer, r, z; phi_p = 0) -> Complex

The m-mode EFFECTIVE SOURCE `S^eff_m(r, z)` — the fixed source of the residual-field elliptic
problem the solver actually solves (main.tex:414-460):

    □Φ^{ℛ} = S^eff = −4πρ − □Φ^{𝒫}   ⟹   △_m Ψ^{ℛ}_m = S^eff_m ,
    S^eff_m := −(Σ/(r²+a²))·P_m[S^eff] = −△_m Ψ^P_m       (off the worldline),

with `P_m` the eq:modeIntegral projector and `−Σ/(r²+a²)` the mode-projection scale fixed by
eq:elliptic_eqn (main.tex:247-262). Equivalently: `S^eff_m` is the m-mode of `−□Φ^𝒫` in the
normalization that makes the SEALED first-order flux form the left-hand side, so

    −∂_i𝓕^i(Ψ^P_m) + β_m Ψ^P_m + γ^i_m ∂_i Ψ^P_m  ==  −S^eff_m

pointwise (the spectre `Test_CircularOrbit` identity, gated on a 20×20 near-puncture grid in
the node tests). Downstream consumers therefore feed `+S^eff_m` as the DG fixed source of
`△_m` inside the worldtube and zero outside; this is the same sign as spectre's
`Tags::FixedSource<Tags::MMode>` (CircularOrbit.cpp:322-332).

`phi_p` is the particle's BOYER–LINDQUIST azimuth (the baked fixtures use 0.31; spectre uses
0); the horizon-regular chart rotation at the field point is applied internally — see the
[`puncture_psi_chart`](@ref) convention chain. Evaluated with the vtu `t`-slab boost `H = 0`
(the puncture lives in the worldtube). Type-generic; `Complex{BigFloat}` under BigFloat input.

ACCURACY (main.tex:488-490): the direct assembly cancels four orders of singularity, so its
relative accuracy degrades like `eps·λ^{-2}` towards the worldline —
[`effective_source_bundle`](@ref) exposes that amplification as `cancellation`. Use the
near-worldline series (node `effective-source-series`) where it exceeds the caller's budget.
"""
effective_source_m(p::Puncture, m::Integer, r::Real, z::Real; phi_p::Real = 0) =
    effective_source_bundle(p, m, r, z; phi_p = phi_p).source
