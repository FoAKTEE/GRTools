"""
    Harmonics

Angular-eigenfunction FOUNDATION for GrTools — the second spine root (Tier-2, design §4.2).
It owns the full angular harmonic surface of the Kerr perturbation problem: spin-weighted
spheroidal eigenvalues λ(s,l,m,γ), the spheroidal harmonic S_lm(θ) carrier and its
derivatives, the spin-weighted spherical-harmonic Y seed, and the UNIFIED m<0 parity surface
that collapses the per-package FEW/KerrBackwards/real-Y forks into one tested flag.

Heavy-dep-FREE by construction: its closure is EXACTLY {LinearAlgebra, NumericsCore} (the
§3.2 binding-DAG spine-root pin). It carries NO local Kummer 1F1 /
continued-fraction copy — those promoted light-tier primitives are consumed from NumericsCore
(design §15.1 H-CORE-STRAT). Wigner is deliberately NOT folded in (would leak coupling onto the
in-degree spine, §4.2).

Submodules (per design §4.2 + §6.2 crosswalk rows dest=Harmonics):
  - `SpheroidalEigenvalues` — λ(s,l,m,γ) + ∂_γλ behind ONE dispatcher over three backends:
    pentadiagonal-LAPACK spectral Float64 fast path, type-generic Leaver continued-fraction
    BigFloat path, small-γ Taylor + large-γ Laurent asymptotics. The BigFloat eigenvalue home
    (H-CORE-STRAT); fold candidate for KerrLorenzCirc `swsh_eigenvalue_big` once its signature
    is confirmed.
  - `SpheroidalHarmonics` — the callable S_lm(θ) carrier (spherical-expansion + Leaver reps)
    with mixed ∂_θ/∂_φ/∂_γ; the Float64-VALUE tier.
  - `SphericalHarmonics` — the γ=0 canonical spin-weighted Y seed (NumericsCore.legendre),
    half-angle BigFloat high-l path, exact θ/θθ derivatives, and ONE
    `SpinWeightedSphericalHarmonicY(s,l,m,θ,φ; parity)` whose `parity` flag (`:none` default,
    `:few` applies the (-1)^l FEW/KerrBackwards flip) unifies the parity surface — one tested
    flag, never a silent per-call fork.
  - `Precision` — `angular_precision_floor()` surfacing the Float64 value-ceiling vs the
    BigFloat eigenvalue-floor.

Every primitive is type-generic `{T}` so BigFloat precision (256/512/1024/2048 bits) is
inherited; the spectral LAPACK eigensolve is wrapped in the BLAS-guard (BigFloat suffers
GC+BLAS contention under threads).

Grounding: test-only characterization against the LIGHT `SpinWeightedSpheroidalHarmonics`
oracle (deps = LinearAlgebra only) at its tolerances, PLUS analytic identities — Y
orthonormality, l=|m| closed forms, the γ→0 spheroidal→spherical limit, the eigenvalue
l(l+1)-s(s+1) at γ=0, and the Condon-Shortley/Goldberg m↔-m relations.
"""
module Harmonics

using LinearAlgebra
using NumericsCore

# Modules are built TDD (test-first, grounded vs the SWSH oracle + analytic identities) in the
# capability BUILD rounds. Each include + its exports land with its grounded test.

include("SphericalHarmonics.jl")          # γ=0 spin-weighted Y seed front door + unified parity
export SpinWeightedSphericalHarmonicY, spin_weighted_spherical_Y, angular_parity
export ylm, ylm_few, ylms_few
# DiscModes real-Y path: SphericalHarmonics is the public home (§6 crosswalk) — re-export the
# NumericsCore seed primitives so callers depend on Harmonics, not NumericsCore directly.
export Ylm_real, Ylm_eq

include("SpheroidalEigenvalues.jl")       # λ(s,l,m,γ) + ∂_γλ: spectral (Float64) + Leaver (BigFloat)
export SpinWeightedSpheroidalEigenvalue, SpinWeightedSpheroidalEigenvalueDeriv

include("SpheroidalEigenvalueSeries.jl")  # small-γ Taylor + large-γ Laurent eigenvalue coeffs
export SpinWeightedSpheroidalEigenvalueSeries, SpinWeightedSpheroidalEigenvalueSeriesLargeGamma

include("SpheroidalHarmonics.jl")         # callable S_lm(θ,φ): SphericalExpansion + Leaver reps
export SWSHSFunction, SWSHSLeaverFunction, SpinWeightedSpheroidalHarmonicS
export SpinWeightedSpheroidalHarmonicSSeries, expansion_coefficients
export swsh_derivative, swsh_gamma_derivative

include("SphereCollocation.jl")           # GL-collocation swal + ð/ð′ ladder / Laplacian / filter / killtop matrices (TeukEvolution Sphere.jl port, D2)
export num_l, swal, Y_vals, cos_vals, sin_vals, inner_product
export swal_vals, swal_laplacian_matrix, swal_raising_matrix, swal_lowering_matrix
export swal_filter_matrix, swal_killtop_matrix, angular_matrix_mult!

include("precision.jl")                   # angular precision policy: Float64 value-ceiling vs BigFloat eigenvalue-floor
export angular_precision_floor

end # module Harmonics
