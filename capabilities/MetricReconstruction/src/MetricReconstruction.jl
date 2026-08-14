"""
    MetricReconstruction

Tier-4 metric-perturbation capability (design §4.4). Four API-separated submodule namespaces let
concept-name overlaps (`weyl_scalar`, `to_BLS`) coexist without renaming:

  - [`Radiation`](@ref) — radiation/ingoing-gauge CCK (Chrzanowski–Cohen–Kegeles) reconstruction:
    the Hertz potential from the Teukolsky ψ4 amplitudes, the regenerated ψ0/ψ4 Weyl scalars, the
    10 reconstructed metric components, and the `reconstruct_kerr` / `to_BLS_kerr` / `closure_kerr`
    assembly. Absorbs the legacy `MetricReconstruction` package.

  - [`Lorenz`](@ref) — Dolan-2024 Lorenz-gauge reconstruction for a Kerr (a≠0) circular orbit, and
    `Lorenz.LorenzRadial`, the CANONICAL home of the a=0 first-order Lorenz-gauge radial
    field-equation RHS (`even_ode!`/`odd_ode!`/`monopole_ode!`/`odd_dipole_ode!`) that
    `SelfForce.LorenzField` and `AccretionDisc.Force` CONSUME but do not own. Absorbs the legacy
    `KerrLorenzCirc` package.

  - [`Transport`](@ref) — Li--Yunes traceful outgoing-radiation-gauge transport for the bounded
    static-source specializations implemented by this package.

  - [`GHZ`](@ref) — sourced Green--Hollands--Zimmerman carriers and bounded correctness kernels in
    outgoing Schwarzschild slicing (`a=0`). The current window includes hardened local residuals;
    the finite-ray corrector head/middle and `l>=1` paired tail; a separate head-only `l=m=0`
    monopole tail; an exact-grid `RayMomentCache`; the finite-ray eta equation for a supplied,
    already-adjusted source; and positive-quarter Hertz inversion followed by arbitrary-frequency
    fixed-outgoing-`u` tetrad-metric injection. The mode, solution, and carrier APIs retain explicit
    provenance and free-IVP policies; the pure residual helpers are local scalar oracles. They do
    not construct `J-E*x`, solve `O†Φ=η`, select physical/retarded data,
    interpolate cached rays, reconstruct `a≠0` or a completed coordinate metric, or provide
    a full second-order self-force.

Grounding includes test-only characterization against the LIGHT live oracles
`MetricReconstruction` and `KerrLorenzCirc` (verified runtime-HDF5-free) at their frozen tolerances,
plus analytic identities and the bounded executable tests owned by each namespace.

Deps: NumericsCore, Wigner (clebsch_gordan, promoting KLC's hand-rolled cg), Spacetime, Harmonics,
Perturbation; external QuadGK.
"""
module MetricReconstruction

# ── Submodule: radiation/ingoing-gauge CCK reconstruction (legacy MetricReconstruction) ──────────
include("Radiation.jl")

# ── Submodule: Dolan-2024 Lorenz-gauge reconstruction + LorenzRadial (legacy KerrLorenzCirc) ─────
include("Lorenz.jl")

# ── Submodule: Li–Yunes traceful radiation-gauge (ORG) reconstruction for generic type-D sources ─
include("Transport.jl")

# ── Submodule: sourced GHZ Hertz-plus-corrector reconstruction ─────────────────────────────────
include("GHZ.jl")

using .Radiation
using .Lorenz
using .Transport
using .GHZ

# Public API is re-exported through the submodule namespaces (Radiation.*/Lorenz.*/Transport.*) to
# preserve the disambiguated weyl_scalar/to_BLS concept names. Concrete exports land with each TDD
# module.
export Radiation, Lorenz, Transport, GHZ

end # module MetricReconstruction
