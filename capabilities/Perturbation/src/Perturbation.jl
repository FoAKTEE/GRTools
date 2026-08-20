"""
    Perturbation

Teukolsky perturbations of Kerr: the homogeneous radial solutions (In/Up) over four backends
(MST, confluent-Heun, Sasaki–Nakamura, GSN), the renormalized angular momentum ν, point-particle
amplitudes and fluxes, the MST scattering phase, spinning-secondary linear-in-spin fluxes, the
QNM/QQNM surface, and the time-domain stack (`Perturbation.TimeDomain`).

LIGHT-DEP FIREWALL: zero HDF5/JLD2/Symbolics/OrdinaryDiffEq/CUDA/KernelAbstractions edges.
"external-ODE-free" rests on the VENDORED DP5/GBS integrators — never refactor them onto
OrdinaryDiffEq. Heavy math reaches a user only as the HypergeometricFunctions/SpecialFunctions
closure carried through BigSpecial; a closure-load gate asserts this.

Every primitive is type-generic `{T}`: BigFloat precision (256–2048 bits) is inherited; the
BigFloat MST/CF tier stays on CPU; heavy flux sweeps scale via processes under the BLAS guard.
"""
module Perturbation

using LinearAlgebra
using Printf
using NumericsCore
using BigSpecial
using Spacetime
using Harmonics

# --- Radial: the MST homogeneous Teukolsky radial engine (dependency-ordered includes) ---------
include("Radial/renorm_angmom.jl")   # RenormalizedAngularMomentum ν-solve + MST/monodromy machinery
include("Radial/mst.jl")             # MST minimal-solution coefficients f_n (cf_convergent tier)
include("Radial/hypergeom.jl")       # exact 2F1 (In) / U (Up) MST terms + recurrence arrays
include("Radial/radial.jl")          # raw MST homogeneous solutions R^in/R^up + r-derivatives
include("Radial/carrier.jl")         # TeukolskyRadialFunction carrier (callable R(r); :Static/:MST/:HeunC)
include("Radial/heunc.jl")           # confluent-Heun radial backend (:HeunC value + analytic deriv)
include("Radial/teukolsky_radial_static.jl")  # static ω=0 closed-form 2F1 path internals

# --- AmplitudesFluxes: Z/C amplitudes + energy/L fluxes + point-particle source convolution -----
include("AmplitudesFluxes/amplitudes.jl")      # MST asymptotic amplitudes (In/Up trans/inc/ref)
include("Radial/teukolsky_radial.jl")          # public TeukolskyRadial + teukolsky_radial_deriv (needs amplitudes)
include("Radial/radial_sn.jl")                 # :SasakiNakamura backend (s=-2 Up; vendored DP5, in-firewall)
# --- :GSN backend: Generalized Sasaki-Nakamura (arXiv:2306.16469; vendored ref-code v0.8.5
# commit 982b087). s∈{0,±1,±2}×{In,Up} at real ω, {T}-generic, vendored DP5/GBS.
include("Radial/gsn_coordinates.jl")           # Kerr primitives + {T}-typed tortoise map (a42/Newton)
include("Radial/gsn_transform.jl")             # SN transformation α, β, η (c₀..c₋₄)
include("Radial/gsn_potentials.jl")            # GSN sF/sU, Teukolsky V_T, SN V_SN
include("Radial/gsn_conversions.jl")           # Teukolsky↔SN amplitude factors + TS constants
include("Radial/gsn_asymptotics.jl")           # in-house Taylor arithmetic + boundary-series coefficients
include("Radial/gsn_boundary.jl")              # asymptotic-ansatz boundary values (X, X′)
include("Radial/gsn_solve.jl")                 # unified real-ω solve + amplitude extraction + gates
include("Radial/gsn_backend.jl")               # SN→Teukolsky conversion + TeukolskyRadial carrier glue
include("AmplitudesFluxes/fluxes.jl")          # energy / angular-momentum fluxes from Z
include("AmplitudesFluxes/pointparticle.jl")   # GHP source convolution → Z + TeukolskyPointParticleMode*
include("Radial/radial_jet.jl")                # per-mode radial derivative tower (shared home; MR re-exports)

# --- Scattering: MST scattering-phase pipeline (REUSES the Radial ν/MST/amplitude tier) ----------
include("Scattering/scattering.jl")

# --- SpinFluxes: spinning-secondary linear-in-spin GW fluxes -------------------------------------
include("SpinFluxes/circular.jl")    # circular-equatorial s=-2 flux + shared FabEq helper basis
include("SpinFluxes/eccentric.jl")   # eccentric-equatorial s=-2 flux (wr Mino sum, 2 velocity branches)
# spinning-secondary kinematics: the numeric Correction engine + orbit-correction assemblies live in
# Spacetime.SpinningSecondary (single owner); only the closed-form inclined branch and the spherical
# engine are still local (scoped out of Spacetime so far)
include("SpinFluxes/spinning_secondary_local.jl")        # KerrSpinOrbit (equatorial+inclined) + pex1CTilde + δx3
include("SpinFluxes/spinning_secondary_bindings.jl")     # KerrSpinOrbitCorrection/Num + _oc_* ← Spacetime.SpinningSecondary
include("SpinFluxes/spinning_orbit_spherical_local.jl")  # spherical engine + KerrSpinOrbitCorrectionSpherical + r/x FD
# the 10 TeukolskySpinMode* entrypoints
include("SpinFluxes/inclined.jl")                # generic angular ops + Fab + _inclined_spin_mode + inclined_spin_flux
include("SpinFluxes/fab_basis.jl")               # fab^(i) basis + _angular_fromcorr + dJz*/_KCorrection
include("SpinFluxes/from_correction.jl")         # TeukolskySpinModeFromCorrection core (wr×wθ, numeric correction)
include("SpinFluxes/spin_mode_num.jl")           # TeukolskySpinMode core (full finite-s numeric orbit)
include("SpinFluxes/spin_mode_correction.jl")    # TeukolskySpinModeCorrection core (monopole + dE¹/dLz¹)
include("SpinFluxes/spherical.jl")               # TeukolskySpinModeSpherical core (θ-only wz sum)
include("SpinFluxes/spherical_correction_num.jl")        # …SphericalCorrectionNum core (δω-FD ω-derivative)
include("SpinFluxes/spherical_correction.jl")            # …SphericalCorrection core (+ r/x FluxesDerivatives)
include("SpinFluxes/spherical_correction_analytical.jl") # …SphericalCorrectionAnalytical core (+ CarterConstantK)
include("SpinFluxes/darwin.jl")                  # …EquatorialCorrectionAnalyticalDarwin core (χr grid)
include("SpinFluxes/entrypoints.jl")             # the 10 public TeukolskySpinMode* wrappers

# --- QNM: the Kerr quasinormal-mode frequency solver ---------------------------------------------
include("QNM/qnm_frequency.jl")      # :leaver reference (radial CF ⟂ complex-γ spheroidal λ)
include("QNM/qnm_fast.jl")           # :fast — shift-invert angular + spin-sequence continuation
include("QNM/qnm_hyperboloidal.jl")  # :hyperboloidal — Chebyshev collocation QEP (F64→BF RQI)

# --- QQNM: the second-order source, solve, and amplitude-ratio surface ---------------------------
include("QNM/qqnm_source_tables.jl") # the baked certified a_ij/b_ij tables (sha256-provenanced)
include("QNM/qqnm_source.jl")        # the quotient fold + angular structures + child projection
include("QNM/qqnm_angular.jl")       # general-a angular reduction (after qqnm_source: reuses its chains)
export qqnm_pp_aij, qqnm_pp_aij_gammaform, qqnm_pm_bij, qqnm_mm_is_zero
export qqnm_lambda_plus, qqnm_source_coeffs, qqnm_source
export qqnm_angular_structures, qqnm_angular_projections, qqnm_source_coeffs_gl
include("QNM/qqnm_solve.jl")         # the sourced hyperboloidal solve → the scri amplitude E
include("QNM/qqnm_ratio.jl")         # the R = A_Q/(A₁A₂) assembly
export qqnm_solve, qqnm_hphc_cond, qqnm_channel_omega, qqnm_ratio, qqnm_ratio_np

# --- TimeDomain: the time-domain DG/spectral Teukolsky stack (own namespace) ---------------------
include("TimeDomain.jl")

export RenormalizedAngularMomentum, TeukolskyRadial, TeukolskyRadialFunction, teukolsky_radial_deriv
export radial_jet, clear_radial_jet_caches!
export energy_flux, angular_momentum_flux
export TeukolskyPointParticleModeCircular, TeukolskyPointParticleModeEccentric,
       TeukolskyPointParticleModeSpherical, TeukolskyPointParticleModeGeneric
export starobinsky_constant, scattering_phase, scattering_phase_numerical, scattering_phase_analytic,
       scattering_phase_table, mst_coefficient, mst_convergence
export qnm_frequency, qnm_cf_residual, qnm_sequence_fast
export hphc_radial_operator, qnm_radial_hphc, hphc_mirror_carry
export TeukolskySpinModeCircularCorrectionAnalytical,
       TeukolskySpinModeEquatorialCorrectionAnalyticalPhases
export TeukolskySpinMode, TeukolskySpinModeFromCorrection, TeukolskySpinModeAnalytical,
       TeukolskySpinModeCorrection, TeukolskySpinModeCorrectionAnalytical,
       TeukolskySpinModeSpherical, TeukolskySpinModeSphericalCorrectionNum,
       TeukolskySpinModeSphericalCorrection, TeukolskySpinModeSphericalCorrectionAnalytical,
       TeukolskySpinModeEquatorialCorrectionAnalyticalDarwin
# the s=0 time-domain DG stack lives in its own namespace (use Perturbation.TimeDomain.<name>)
export TimeDomain

end # module Perturbation
