"""
    Hybrid.Drivers

N-var inspiral ODE evolution engine (design §4.9): the single-trajectory `Vern9` integrator, the
threaded batch driver, the `*_at` trajectory evaluators and the `r0`-stop / `ω`-crit termination
callbacks, over the 2-var (Ω, φ) Schwarzschild, the 4-arg (ω, φ; χ1) Kerr and the 7-var
(Ω, φ, ν, M, χt1, χt2, δm) spin/δm states. This is where the modern `OrdinaryDiffEqVerner` /
`SciMLBase` precompile floor is confined (swapped off the legacy `OrdinaryDiffEq` monolith).
"""
module Drivers

using OrdinaryDiffEqVerner: Vern9
using SciMLBase: ODEProblem, ContinuousCallback, terminate!, solve
using NumericsCore: blas_single_thread

using ..Models: Forcing0PASchwarz, Forcing1PAT1, Forcing0PAKerr, Model1PAT1ea, Model1PAT1R,
    inspiral_rhs!, kerr_inspiral_rhs!, _spin_eom!, stop_r0, omega_crit

include("Drivers/inspiral_driver.jl")

export InspiralSolution, SpinInspiralSolution, int_inspiral, int_inspiral_batch,
    duration, omega_at, r0_at, phi_at, chit1_at, dm_at

end # module Drivers
