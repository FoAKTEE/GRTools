"""
    R3Derivation

Route **R3** — the Julia-native symbolic derivation of the relativistic master system (47)–(49), the
eliminated velocities and the diagnostics integrands of arXiv:2608.17003, starting from the paper's
equations only:  line element (4)/(A1), fluid equations (3), linearisation (15)–(22), Fourier ansatz
(45)–(46) and the reduction scheme of M:604–608.

This is a STANDALONE project (`capabilities/AccretionDisc/derivation/`), **not** a dependency of
`AccretionDisc`: it produces the generated sources, it never runs inside the package.  Route R2 (our
independent Mathematica derivation) is the ORACLE that gates this route; no expression of it is read
in here.  The notebook-derived oracle of route R1 was removed from this repository on 2026-08-22
(owner decision, EGD::jd-notebook-purge) and is archived in the mission repository; route R3 never
read it either, so nothing in this project changed.

Layout — `symbols.jl` (alphabet + the radial derivation `ddr`), `geometry.jl` (the one generic
Christoffel helper), `jets.jl` (first-order λ-jets carrying the Fourier phase), `fluid.jl` (Eq. (3),
written once and used for both the background and its linearisation), `background.jl` (Sec. II.A),
`certify.jl` (the symbolic and exact-rational admission modalities), `stage1.jl` (the background
certificate suite), `perturbation.jl` + `reduce.jl` (stage 2), `diagnostics.jl` (stage 3).
"""
module R3Derivation

using Symbolics
using LinearAlgebra
using Random
using Printf

include("symbols.jl")
include("geometry.jl")
include("jets.jl")
include("fluid.jl")
include("background.jl")
include("certify.jl")
include("series.jl")
include("stage1.jl")
include("perturbation.jl")
include("reduce.jl")
include("newtonian.jl")
include("stage2.jl")
include("diagnostics.jl")
include("serialize.jl")
include("codegen.jl")
include("locally_isothermal.jl")
include("thermal_balance.jl")

end # module
