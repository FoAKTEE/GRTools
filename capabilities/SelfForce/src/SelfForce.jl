"""
    SelfForce

Tier-5 1st+2nd-order gravitational self-force capability (design §4.5): the GSF corpus that fuses
the legacy `H1Lorenz`, `Punctures`, and `SecondOrderRicci` packages over the shared h1 corpus
(intrinsic HDF5, mediated by `Data.H1Schema`). This is HDF5 HOME #1 of the capability layer. The
public surface is held in three disjoint submodule namespaces:

  - [`LorenzField`](@ref) — the type-generic Lorenz-gauge BLS first-order field solver (legacy
    `H1Lorenz`). It CONSUMES `MetricReconstruction.Lorenz.LorenzRadial`'s radial field-equation RHS
    (`even_ode!`/`odd_ode!`/`monopole_ode!`/`odd_dipole_ode!`) — it does not own them — assembles
    the R1..R10 sources, the in/out boundary-condition series, integrates with the NumericsCore
    vendored DP5+GBS, builds the gauge field, and WRITES the h1 corpus via `Data.write_mode_h5`.

  - [`SecondOrder`](@ref) — the second-order Ricci source (legacy `SecondOrderRicci`). It reads the
    h1 corpus via `Data`, evaluates the 10 quadratic combinations `r2_1..r2_10` on the
    `Kernels.r2` GPU host, couples 3 harmonics through `Wigner`, supplies the analytic
    `h1_mass`/`h1_spin` and `d2hdr2`, and assembles `second_order_ricci` / `ricci_source` then
    `write_src`.

  - [`Puncture`](@ref) — the puncture / effective-source machinery (legacy `Punctures`): static h0
    AGM elliptic integrals, auto-generated coefficients, terminating-`2F1` α-integration via
    `BigSpecial`, Wigner-D β=π/2 mode mixing, digamma regularisation sums, the singular/retarded
    field subtraction, and the worldline self-field.

Grounding: test-only characterization against CAPTURED fixtures (the legacy oracles are HEAVY and
are never live test-deps) — the compiled g++/GSL `r2_oracle.cpp`+`h1.cc` oracles, wolframscript /
Mathematica for the puncture, analytic identities, and a round-trip h1-corpus write/read. Every
primitive is type-generic `{T}` so BigFloat precision (256..2048-bit) is inherited.

Deps: NumericsCore, BigSpecial, Wigner, Kernels, Data, Harmonics, MetricReconstruction; external
HDF5 (the h1 corpus).
"""
module SelfForce

# ── Submodule: Lorenz-gauge first-order field solver + h1-corpus writer (legacy H1Lorenz) ────────
include("LorenzField.jl")

# ── Submodule: second-order Ricci source over the h1 corpus (legacy SecondOrderRicci) ────────────
include("SecondOrder.jl")

# ── Submodule: puncture / effective-source machinery (legacy Punctures) ──────────────────────────
include("Puncture.jl")

using .LorenzField
using .SecondOrder
using .Puncture

# Public API is re-exported through the submodule namespaces (LorenzField.*/SecondOrder.*/
# Puncture.*). Concrete exports land with each TDD module.

end # module SelfForce
