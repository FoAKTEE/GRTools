"""
    ModifiedGravity

Beyond-GR (modified Teukolsky) source machinery (Tier-4, design §4.7) — the THIN, independently
versioned research surface absorbing legacy `ModifiedTeukolsky` (Polcar & Witzany 2025,
arXiv:2507.15720, "Towards relativistic inspirals into black holes surrounded by matter").

Physics: a thin equatorial matter ring back-reacts on the Kerr background; the reconstructed (0,1)
metric builds an EFFECTIVE SOURCE so the perturbed Weyl scalar ψ4 obeys a Teukolsky equation with
that source. This capability owns the source-table I/O, the per-mode radial Hertz front-end, the
closed-form NP/GHP source kernels, the assembled source amplitudes, the perturbed (1,1) ring
dynamics, and the closed-form ψ4←Hertz (CCK) reconstruction.

The CCK psi4←Hertz identity is hard-baked as a closed-form kernel (`_psi4_hertz_kernel`). The Hertz
front-end REUSES Perturbation's Teukolsky radial engine; the finite-r h44 GHP kernel REUSES the
canonical NumericsCore.jets jet calculus and Perturbation.radial_jet (the s=+2 Teukolsky ODE
tower). HDF5-free + Symbolics-free.

Submodules (design §4.7; built TDD across the BUILD rounds — each include + its grounded test
lands together):
  - `SourceTables`    — SourceMode POD + 12-col `.dat` ring-forcing table reader/writer (via Data /
    NumericsCore.dataschema; stdlib parse path preserved).
  - `HertzPotential`  — radial Hertz per (l,m,n,k) from the -2 ψ4 amplitude × homogeneous Up, with
    the load-bearing `kloc = (x==1 ? 0 : k)` convention; reuses Perturbation Teukolsky radial.
  - `SourceKernels`   — the closed-form NP/GHP `Sourceace`/`Sourcenorm` kernels (rational functions
    at equatorial Schwarzschild θ=π/2), type-generic {T}.
  - `SourceAmplitudes`— per-mode Sl/Sn/Sm source amplitudes + ω + Snorm assembly, and the `.dat`
    generator (`rhs_amplitudes`).
  - `RingDynamics`    — perturbed (1,1) ring geodesic + density/continuity response, NP→BL tetrad
    map, `ring_urho_amplitudes`, `ring_functions`.
  - `GHPOperators`    — native GHP `thorn`/`thorn_prime` (þ,þ′≡PDNP[{2,−NP}]) directional derivatives +
    `schw_spincoeff`, reusing the canonical NumericsCore.jets jet calculus (the h44 leaf).
  - `WeylFromHertz`   — `psi4_from_hertz` via the closed-form CCK `_psi4_hertz_kernel`, and
    `h44_from_hertz` — the finite-r reconstructed (4,4)=m̄m̄ metric component assembled NATIVELY from
    the recovered GHP operator (SC l.368), no xAct.

Grounding: frozen `ModifiedTeukolsky` characterization literals, the 3 shipped `.dat` source-table
oracles (round-trip via Data), the analytic CCK ψ4 reconstruction identity (reldiff ~4e-9), the
standalone GHP algebra (commutator/Ricci identities) for þ/þ′, and — for h44 — the necessary
`h44 → 2·ψ4/ω²` asymptotic plus the sufficient finite-r match to the independent wolframscript
operator oracle (`oracle/h44_from_hertz_oracle.wl`).

Public API: stubs only at scaffold stage; real entry points land in the per-module BUILD rounds.
"""
module ModifiedGravity

using Printf
using NumericsCore
using Harmonics
using Perturbation
using Data

# Submodule includes land in the BUILD rounds (TDD, each grounded), in dependency order:
# SourceTables → SourceKernels → HertzPotential → SourceAmplitudes → (RingDynamics, WeylFromHertz).
include("ghp_operators.jl")       # thorn/thorn_prime GHP derivatives + schw_spincoeff (h44 leaf)
include("source_tables.jl")       # SourceMode, read/write_source_table, source_mode_line
include("source_kernels.jl")      # Sourceace*/Sourcenorm* NP/GHP kernels (yv-callable, no Harmonics dep)
include("hertz_potential.jl")     # hertz_potential_mode (RadHZ/dRadHZ)
include("source_amplitudes.jl")   # source_amplitudes (Sl/Sn/Sm/Smdag + ω + Snorm)
include("rhs_amplitudes.jl")      # rhs_amplitudes (.dat source-mode generator: Σ_l at ringrad)
include("ring_dynamics.jl")       # ring_urho_amplitudes, ring_functions (perturbed (1,1) ring)
include("weyl_from_hertz.jl")     # psi4_from_hertz + h44_from_hertz (native GHP thorn_prime kernel, oracle-grounded)

end # module ModifiedGravity
