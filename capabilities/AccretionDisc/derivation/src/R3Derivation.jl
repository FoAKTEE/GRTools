"""
    R3Derivation

Route **R3** — the Julia-native symbolic derivation of the relativistic master system (47)–(49), the
eliminated velocities and the diagnostics integrands of arXiv:2608.17003, starting from the paper's
equations only:  line element (4)/(A1), fluid equations (3), linearisation (15)–(22), Fourier ansatz
(45)–(46) and the reduction scheme of M:604–608.

This is a STANDALONE project (`capabilities/AccretionDisc/derivation/`), **not** a dependency of
`AccretionDisc`: it produces the generated sources, it never runs inside the package.  The authors'
supplementary notebook (route R1) and the independent Mathematica derivation (route R2) are ORACLES
that gate this route; no expression of theirs is read in here.

Layout — `symbols.jl` (alphabet + the radial derivation `ddr`), `geometry.jl` (the one generic
Christoffel helper), `jets.jl` (first-order λ-jets carrying the Fourier phase), `fluid.jl` (Eq. (3),
written once and used for both the background and its linearisation), `background.jl` (Sec. II.A),
`certify.jl` (the symbolic and exact-rational admission modalities), `series.jl` (truncated Laurent
series over ℚ) and `stage1.jl` (the background certificate suite).  Stages 2 and 3 land in the
following nodes of the series.
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

end # module
