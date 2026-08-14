"""
    Hybrid.HybridC

The WaSABI-C SF-PN-hybrid forcing engine (design §4.9): the 9-var closed-form forcing channels
(Fω/Fν/Fm/Fχ1/Fχ2 + the δm/δν/δχ shift combinations), the δx waveform-frequency correction and
the δx-coupled trajectory evolution (`Vern9` + `ContinuousCallback` stop).

The δx / dδxdχ / forcing closed forms are a **deterministic machine-generated artifact**
(`HybridC/wasabic_codegen.jl`, the WL InputForm translation `codegen_wasabic.wl`) over the
`WaSABICModel` interpolant bundle (`HybridC/wasabic_deps.jl`). Routes the complex digamma to the
canonical `BigSpecial.PolyGamma` (AHA — the legacy hand-rolled `_digamma` fork is collapsed into
the one shared primitive), interpolation to `Hybrid.Interpolation`, the Kerr 1SF fluxes to
`Hybrid.Models.forcing_0PAKerr`, and the `.m` SF tables to `Data`. Extends the `Drivers`
`int_inspiral`/`duration`/`omega_at`/`phi_at` and `Models.omega_crit` generics for the
`WaSABICModel`.
"""
module HybridC

using Spacetime: kerr_isco
using Data: read_wl_table, read_wl_table_complex, read_wl_data, read_wl_grid2d, read_wl_grid3d
using BigSpecial: PolyGamma
using ..Interpolation: Interp1, Interp2Hermite, Interp3Lagrange,
    interp1_order8, interp2_hermite, interp3_lagrange, partial, derivative
using ..Models: forcing_0PAKerr
import ..Models: omega_crit
import ..Drivers: int_inspiral, duration, omega_at, phi_at

using OrdinaryDiffEqVerner: Vern9
using SciMLBase: ODEProblem, ContinuousCallback, terminate!, solve

include("HybridC/wasabic_deps.jl")
include("HybridC/wasabic_codegen.jl")
include("HybridC/wasabic_inspiral.jl")

export WaSABICModel, wasabic_model, PolyGamma,
    wc_x, wc_s, wc_sigma, wc_δx, wc_dδxdχ,
    wc_Fφ, wc_forcing, omega_crit, WaSABICSolution,
    int_inspiral, duration, omega_at, phi_at, wc_state_at

end # module HybridC
