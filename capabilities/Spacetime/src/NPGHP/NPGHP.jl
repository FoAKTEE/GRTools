"""
    Spacetime.NPGHP

NP/GHP coefficient kernel for the Kerr background. This module owns **coefficients**
(pointwise closed-form geometry, chart-parameterized); operator **application** (jets,
SepField, grid matrices, mode iω) stays with each consumer — so the þ/þ′ operators are
deliberately absent. Namespaced `Spacetime.NPGHP.*` (symbols like `ρ`, `μ` would collide
flat).

Convention pins ([AXIOM], gated by the first-principles oracle test/npghp_oracle.jl):
ρ = −1/ζ, μ = −Δ/(2Σζ), τ = −ia sinθ/(√2 Σ), ψ₂ = −1/ζ³; primed view ρ′ = −μ, ε′ = −γ
(full-Kerr leg-swap gated). Spin coefficients are frame scalars, so the table is
chart-invariant over `coords ∈ (:BL, :In, :Out)`.

HPHC chart (R = cl²/r, scri⁺ at R=0): the six peeling-stripped NP₀ scalars, the eleven
GHP `pre_*` coefficient fields, and the HPHC↔BL boost map `A = Δ/(2Σ)`,
`e^{iλ} = −ζ̄/ζ` (cross-chart identity gates G-F1).

Deliberate deviation from the TeukEvolution source: `sqrt(T(2))` replaces its Float64
`sqrt(2)` literals so BigFloat reaches full precision.

Derivations and history: progress/emri-disk/npghp_reformulation/ (F1/F5 notes);
user docs: docs/capabilities/spacetime.md.
"""
module NPGHP

include("weights.jl")
export ghp_weights, ghp_prime_weights, ghp_bar_weights

include("spincoeff_bl.jl")
export np_spincoeff, ghp_spincoeff

include("hphc.jl")
export hphc_R_of_r, hphc_r_of_R, hphc_horizon_R, np_spincoeff_hphc, hphc_ghp_coeffs
export hphc_boost, hphc_spin_phase

end # module NPGHP
