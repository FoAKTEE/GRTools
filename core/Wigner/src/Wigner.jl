"""
    Wigner (GrToolsWigner)

Hand-rolled angular-momentum coupling — the Tier-1 core stratum (C2) that ELIMINATES the
`WignerSymbols` package from every GrTools runtime path (design §5.3 / §7.4 / L4). Its only
dependency is `NumericsCore`; it sits BESIDE `Harmonics`, never on the in-degree light-numeric
spine.

Owns (per design §5.3 + §6 crosswalk rows dest=Wigner):
  - `wigner3j` / `wigner6j` — Wigner 3j + 6j symbols (Racah closed form, BigInt-exact)
  - `clebsch_gordan` — Racah Clebsch-Gordan (promotes KerrLorenzCirc's hand-rolled `cg`)
  - `coupling` — 3-harmonic spin-weighted angular coupling (the canonical home of the
    SecondOrderRicci `coupling`, with its WignerSymbols `[deps]` dropped)
  - `wignerD` — reduced Wigner-d (general β) + the β=π/2 special case (promoted from
    Punctures `_wignerd_halfpi` / `_wignerD_fixed`)

Every primitive is type-generic `{T}` so BigFloat precision is inherited; the exact
combinatorics run in `BigInt`/rational arithmetic and downcast only at the boundary.

Grounding: analytic 3j/6j identities (orthogonality, Regge symmetries, closed forms,
triangle/selection rules, 3j↔CG) PLUS test-only characterization against the
`WignerSymbols.jl` oracle. The runtime package has NO `WignerSymbols` dependency.
"""
module Wigner

using NumericsCore

# Modules built TDD (test-first, grounded vs analytic identities + WignerSymbols oracle) in the
# capability BUILD rounds. Stubs now; each include + its exports land with its grounded test.
include("wigner3j.jl")        # wigner3j, wigner6j (Racah)
export wigner3j, wigner6j
include("clebsch_gordan.jl")  # clebsch_gordan (Racah CG)
export clebsch_gordan
include("coupling.jl")        # coupling (3-harmonic angular coupling)
export coupling
include("wignerD.jl")         # wignerd, wignerd_halfpi, wignerD_fixed (reduced Wigner-d)
export wignerd, wignerd_halfpi, wignerD_fixed

end # module Wigner
