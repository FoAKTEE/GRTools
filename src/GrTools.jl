"""
    GrTools

Umbrella module for the GrTools suite — Julia packages for black-hole perturbation
theory: Kerr geometry, spin-weighted harmonics, the Teukolsky equation, metric
reconstruction, gravitational and scalar self-force, post-Newtonian expansions,
EMRI waveforms, and relativistic accretion-disc modes.

`using GrTools` loads every package in the suite; each is exposed as a qualified
submodule binding (`GrTools.Spacetime`, `GrTools.Perturbation`, …). No symbols are
re-exported at the top level, so subpackage APIs never collide — bring the ones you
want into scope yourself:

```julia
using GrTools
using GrTools.Spacetime, GrTools.Perturbation
```

or depend on the individual packages directly (each is an independent Julia package
under `core/` or `capabilities/` with its own `Project.toml` and test suite).

# Suite layout

Core tier (physics-free numerics):

  - [`NumericsCore`](@ref) — special functions, jets, root-finding, interpolation,
    collocation, ODE, DG primitives (dep: LinearAlgebra only)
  - [`BigSpecial`](@ref) — Complex{BigFloat}-correct Γ/ψ/₂F₁/₁F₁
  - [`Data`](@ref) — data root, Pkg artifacts, format readers
  - [`TestKit`](@ref) — test-only tolerance oracle and CI hygiene gates
  - [`Wigner`](@ref) — 3j/6j/Clebsch–Gordan/Wigner-d, exact rational core
  - [`Kernels`](@ref) — second-order Ricci quadratic kernels
  - [`EllipticCore`](@ref) — complex-parametric 2D hp-DG elliptic solver stack

Capability tier (physics over core):

  - [`Spacetime`](@ref) — Kerr geometry, charts, tetrads, NP/GHP fields
  - [`Harmonics`](@ref) — spin-weighted spherical/spheroidal harmonics
  - [`Perturbation`](@ref) — the Teukolsky engine: radial solutions, fluxes, QNM/QQNM
  - [`MetricReconstruction`](@ref) — h_μν from Weyl scalars (CCK, Lorenz gauge)
  - [`SelfForce`](@ref) — gravitational self-force (Lorenz gauge, punctures, 2nd order)
  - [`ScalarSelfForce`](@ref) — m-mode elliptic Kerr scalar self-force
  - [`PostNewtonian`](@ref) — PN expansions and flux catalogues
  - [`ModifiedGravity`](@ref) — modified-Teukolsky ring-forced sources
  - [`EMRIWaveforms`](@ref) — adiabatic EMRI trajectories and waveforms
  - [`Hybrid`](@ref) — hybrid PN/Teukolsky inspiral models
  - [`AccretionDisc`](@ref) — relativistic disc modes and forcing

# Conventions

Geometric units `G = c = 1`, black-hole mass `M = 1` unless stated. Boyer–Lindquist
coordinates; Kinnersley tetrad; fields `∝ e^{i(mφ − ωt)}`. Every numerical path is
type-generic: BigFloat inputs inherit the active `setprecision(BigFloat, bits)`
end-to-end.
"""
module GrTools

# Core tier
import NumericsCore
import BigSpecial
import Data
import TestKit
import Wigner
import Kernels
import EllipticCore

# Capability tier
import Spacetime
import Harmonics
import Perturbation
import MetricReconstruction
import SelfForce
import ScalarSelfForce
import PostNewtonian
import ModifiedGravity
import EMRIWaveforms
import Hybrid
import AccretionDisc

# The subpackages, exposed as qualified bindings (GrTools.Spacetime, …).
export NumericsCore, BigSpecial, Data, TestKit, Wigner, Kernels, EllipticCore,
       Spacetime, Harmonics, Perturbation, MetricReconstruction, SelfForce,
       ScalarSelfForce, PostNewtonian, ModifiedGravity, EMRIWaveforms, Hybrid,
       AccretionDisc

"""
    GrTools.PACKAGES

NamedTuple of the suite's packages grouped by tier — `PACKAGES.core` and
`PACKAGES.capabilities`, each a tuple of the loaded subpackage modules.
"""
const PACKAGES = (
    core = (NumericsCore, BigSpecial, Data, TestKit, Wigner, Kernels, EllipticCore),
    capabilities = (Spacetime, Harmonics, Perturbation, MetricReconstruction,
                    SelfForce, ScalarSelfForce, PostNewtonian, ModifiedGravity,
                    EMRIWaveforms, Hybrid, AccretionDisc),
)

end # module GrTools
