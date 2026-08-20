# Migration guide: from the predecessor codes

GrTools regroups the functions of the earlier research codes by physics
capability. This page is the lookup table (old package → new home): for every
predecessor package / file / function group, its new GrTools home. Look up the
old package you were using, find the capability (and submodule) that now owns
the function.

**★ = crosses a predecessor-package boundary** (an actionable reorg move — the
function no longer lives with its historical siblings). Rows without ★ moved
wholesale with their package into the named capability.

The predecessor codes are not part of this repository; their outputs live on as
baked characterization fixtures under each package's `test/fixtures/`.

## → Spacetime

| Legacy pkg / file / fn | New home | ★ |
|---|---|---|
| KerrGeodesics / kerrgeo_constants.jl / KerrGeoEnergy,AngularMomentum,CarterConstant,ConstantsOfMotion | Spacetime.Geodesics | |
| KerrGeodesics / kerrgeo_frequencies.jl / KerrGeoRadialRoots,KerrGeoFrequencies,ProperFrequencyFactor,_kgMinoFreq* | Spacetime.Geodesics | |
| KerrGeodesics / kerrgeo_trajectory.jl / KerrGeoTrajectory(+Equatorial,Radial,Polar),KerrGeoProperTime | Spacetime.Geodesics | |
| KerrBackwards / utility.jl / constants_of_motion | Spacetime.Geodesics (delegating wrapper) | ★ |
| KerrBackwards / utility.jl / separatrix,_sep_poly_*,_sep_root | Spacetime.Separatrix.separatrix (CANONICAL) | ★ |
| FastEMRIWaveforms / separatrix.jl / get_separatrix,_separatrix_polynomial_*,_brentq | Spacetime.Separatrix (alias of separatrix) | ★ |
| WaSABI / kerr_isco.jl + DiscModes _risco + DiscEMRIForce _risco_disc + FEW kerrflux isco | Spacetime.Separatrix.kerr_isco | ★ |
| MetricReconstruction / background.jl / kerr_metric_BL/In/Out(+inv),kinnersley_* ; DiscModes/kerr_geometry.jl ; KerrLorenzCirc/operators.jl | Spacetime.Geometry (coords-method) | ★ |
| MetricReconstruction / {einstein,weyl_scalars}.jl / christoffel_BL,_dg,_dginv ; DiscModes ; Punctures/worldline.jl | Spacetime.Geometry | ★ |
| MetricReconstruction / background.jl / Δ/Σ/ζ/K/Q ; KerrLorenzCirc/operators.jl / geom,Kfn,Qfn,Dcoef | Spacetime.Geometry.bg_scalars | ★ |
| H1Lorenz / orbit.jl / rs_of_r,r_of_rs,_f,_df,_Vl,orbit_energy/ang_mom/phi_freq ; DiscModes/lorenz_metric.jl | Spacetime.Geometry | ★ |
| MetricReconstruction / weyl_amplitudes.jl / _horizons,_Omega_plus | Spacetime.Geometry | ★ |
| SpinningOrbit / spin_orbit.jl, linear_parts.jl, spin_deviations.jl, numeric_correction*.jl, orbit_*.jl | Spacetime.SpinningSecondary | |
| WaSABI / pn_kinematics.jl / kepler_Ω,r0_of_Ω,pn_x (+aliases) | Spacetime.PNKinematics | ★ |

## → Harmonics

| Legacy pkg / file / fn | New home | ★ |
|---|---|---|
| SpinWeightedSpheroidalHarmonics / spectral.jl / SpinWeightedSpheroidalEigenvalue(+Deriv),_swsh_eigenvalue_spectral | Harmonics.SpheroidalEigenvalues | |
| SWSH / leaver.jl / _swsh_eigenvalue_leaver,_swsh_leaver_coeffs | Harmonics.SpheroidalEigenvalues | |
| SWSH / series_coeffs.jl / EigenvalueSeries(+LargeGamma),_swsh_alpha | Harmonics.SpheroidalEigenvalues | |
| SWSH / harmonic.jl, leaver.jl / SpinWeightedSpheroidalHarmonicS,SWSHSFunction,SWSHSLeaverFunction,swsh_derivative | Harmonics.SpheroidalHarmonics | |
| SWSH / spherical_Y.jl / SpinWeightedSphericalHarmonicY,_y_terms,_y_dtheta | Harmonics.SphericalHarmonics | |
| FastEMRIWaveforms / harmonics.jl / ylm | Harmonics.SphericalHarmonics (convention=:few shim) | ★ |
| KerrBackwards / ylm.jl / ylm_few,ylms_few | Harmonics.SphericalHarmonics (convention=:few shim) | ★ |
| DiscModes / metric_recon.jl / Ylm_real,Ylm_eq | Harmonics.SphericalHarmonics (real-Y path) | ★ |
| KerrLorenzCirc / precision.jl / swsh_eigenvalue_big | Harmonics.SpheroidalEigenvalues (fold candidate) | ★ |
| NEW | Harmonics: spin_weighted_spherical_Y, angular_parity, angular_precision_floor | |

## → Perturbation

| Legacy pkg / file / fn | New home | ★ |
|---|---|---|
| Teukolsky / {continued_fraction,renorm_angmom,mst,radial,hypergeom,heunc,teukolsky_radial(_static),radial_numint,radial_sn}.jl | Perturbation.Radial | |
| Teukolsky / {amplitudes,fluxes,source,point_particle}.jl | Perturbation.AmplitudesFluxes | |
| ScatteringPhaseMST / *.jl (scattering_phase,starobinsky_constant,table,brefl_over_binc,mst_coefficient,mst_convergence) | Perturbation.Scattering | ★ |
| TeukolskySpinFluxes / *.jl (13 TeukolskySpinMode*, circular/eccentric/inclined flux + corrections) | Perturbation.SpinFluxes | ★ |
| Teukolsky / hypergeom.jl,amplitudes.jl / _f21_maclaurin,_f21_zinf,_safe_2f1,_gamma_big,_pochhammer | → BigSpecial (delegated, gate) | ★ |

## → MetricReconstruction

| Legacy pkg / file / fn | New home | ★ |
|---|---|---|
| MetricReconstruction / api.jl,operators.jl,modes.jl,tetrad_components*.jl,weyl_amplitudes.jl,projection*.jl,einstein.jl,weyl_scalars.jl | MetricReconstruction.Radiation | |
| KerrLorenzCirc / assembly.jl,operators.jl,angular.jl,background.jl,jumps.jl,closure.jl,radial.jl,radial_kerr.jl,components*.jl | MetricReconstruction.Lorenz | ★ |
| KerrLorenzCirc / angular.jl / cg | → Wigner.clebsch_gordan | ★ |
| MetricReconstruction / background.jl,einstein.jl / kerr_metric,kinnersley,christoffel_BL,_dg | → Spacetime.Geometry (delegate) | ★ |

## → SelfForce

| Legacy pkg / file / fn | New home | ★ |
|---|---|---|
| H1Lorenz / {orbit,sources,boundary_conditions,integrator,scaling,static,gauge_fields,coupled_sets,assemble,driver,io_hdf5}.jl | SelfForce.LorenzField | |
| H1Lorenz / ode_systems.jl / even_ode!,odd_ode!,monopole_ode!(g,rs,y,l,ω),odd_dipole_ode! | → MetricReconstruction.Lorenz.LorenzRadial (CANONICAL; LorenzField CONSUMES) | ★ |
| Punctures / {static,h0_coefficients,alpha_integration,puncture_modes*,mode_mixing,regsum,tensor_basis,subtraction*,singular_field,worldline_self_field}.jl | SelfForce.Puncture | |
| SecondOrderRicci / r2_1..r2_10.jl,r2terms.jl,h1_mass.jl,h1_spin.jl,d2hdr2.jl,modes.jl,assembly.jl,read_h1.jl | SelfForce.SecondOrder | |
| H1Lorenz / math_funcs.jl / associated_legendre_P,_fact | → NumericsCore.legendre + BigSpecial.biggamma | ★ |
| H1Lorenz / driver.jl / gamma_to_lm,lm_to_gamma,mode_list ; SecondOrderRicci/modes.jl/a_il | → NumericsCore.modeindex | ★ |
| SecondOrderRicci / r2_1..r2_10 | → Kernels.r2 (GPU host) | ★ |
| SecondOrderRicci / coupling.jl / coupling | → Wigner.coupling | ★ |
| Punctures / mode_mixing.jl / _wignerd_halfpi,_wignerD_fixed | → Wigner.wignerD | ★ |
| Punctures / alpha_integration.jl / _poch,_f21_term ; regsum.jl/_regsum_tail | → BigSpecial.{pochhammer,hyp2f1,digamma} | ★ |
| Punctures / worldline_self_field.jl / _christoffel(r,θ,M)(+dr,dtheta) | → Spacetime.Geometry (a=0/M-only METHOD of canonical christoffel) | ★ |

## → PostNewtonian

| Legacy pkg / file / fn | New home | ★ |
|---|---|---|
| TeukolskyPN / pnseries.jl / PNSeries,Normal,arithmetic,change_series_parameter,power_counting,expand_log/gamma | PostNewtonian.Series | |
| TeukolskyPN / mst_coefficients.jl / MSTCoefficients,_solve_nu,_eps_to_*,imagunit_sym | PostNewtonian.MST | |
| TeukolskyPN / radial_pn.jl,mode_pn*.jl,mode_pn_assemble.jl / TeukolskyRadialPN,TeukolskyPointParticleModePN | PostNewtonian.RadialMode | |
| TeukolskyPN / mode_pn_source_data.jl / inline oracle Dicts | → Data.artifacts | ★ |
| PostNewtonianSelfForce / {postnewtonian_data,coefficients,resummed_series}.jl | PostNewtonian.SelfForce | ★ |
| PostNewtonianSelfForce / catalog_selector.jl,series_data.jl | PostNewtonian.Catalog | ★ |
| TeukolskyPN/mode_pn_source.jl _eval_mma_ast + PNSF/series_data.jl _eval_pn_ast | PostNewtonian.MMAInputForm.read_mma_inputform_symbolic (unified, gate) | ★ |

## → ModifiedGravity

| Legacy pkg / file / fn | New home | ★ |
|---|---|---|
| ModifiedTeukolsky / source_table_io.jl,hertz_potential.jl,source_kernels.jl,source_amplitudes.jl,rhs_amplitudes.jl,ring_geodesic.jl,ring_functions.jl,weyl_from_hertz.jl | ModifiedGravity.{SourceTables,HertzPotential,SourceKernels,SourceAmplitudes,RingDynamics,WeylFromHertz} | |
| ModifiedTeukolsky / source_table_io.jl / _mma_sci | → Data.dlm (only if 2nd user) | ★ |

## → EMRIWaveforms

| Legacy pkg / file / fn | New home | ★ |
|---|---|---|
| KerrBackwards / trajectory.jl,pn5fluxes.jl ; FEW/trajectory/{emri_inspiral,dopr853,ode_base,pn5_rhs,pn5}.jl | EMRIWaveforms.Trajectory (KB.* / FEW.* forks) | |
| FEW/trajectory/{flux,kerrflux}.jl ; KerrBackwards/{schwarzeccflux,fastschwarz_waveform,slowschwarz_waveform}.jl | EMRIWaveforms.FluxModels (forked) | |
| FEW/amplitude/{ampinterp2d,romannet}.jl ; KerrBackwards/{romanamplitude,interp2damplitude}.jl | EMRIWaveforms.AmplitudeModels (forked) | |
| FEW/modeselector.jl ; KerrBackwards/modeselector.jl | EMRIWaveforms.ModeSelect (forked) | |
| FEW/summation/{directmodesum,interpolatedmodesum,fdinterpolatedmodesum}.jl ; KerrBackwards/{summation,fdsummation}.jl | EMRIWaveforms.ModeSum (forked; kernels→Kernels) | |
| FEW/summation/aaksummation.jl ; KerrBackwards/{summation,bessel_ak,parametermap_aak}.jl | EMRIWaveforms.AAK (forked; →Kernels.aak) | |
| FEW/mappings/{kerrecceq,pn,jacobian,schwarzecc}.jl,geodesic_elq.jl ; KerrBackwards/inclination.jl | EMRIWaveforms.Mappings | |
| FEW/waveform/generate.jl ; KerrBackwards/waveform.jl | EMRIWaveforms.Generate (forked frame fns) | |
| FEW/geodesic.jl,separatrix.jl ; KerrBackwards/utility.jl geometry | → Spacetime (migrate-out) | ★ |
| both: interpolant/gslbicubic/modeindex/constants/elliptic | → NumericsCore (hoist, single copy) | ★ |
| FEW: import Teukolsky/TeukolskySpinFluxes/SpinningOrbit/TeukolskyPN | DELETE (dead, 0 call sites) | ★ |

## → Hybrid

| Legacy pkg / file / fn | New home | ★ |
|---|---|---|
| WaSABI / inspiral_{0PASchwarz,1PAT1,1PAT1ea,0PAKerr}.jl,amplitude_*.jl,interpolation.jl(interp_logexp,kerr_grid_coords) | Hybrid.Models | |
| WaSABI / inspiral_WaSABI-C{,_deps,_codegen}.jl | Hybrid.HybridC | |
| WaSABI / binding_energy.jl | Hybrid.BindingEnergy | |
| WaSABI / inspiral_driver.jl,inspiral.jl | Hybrid.Drivers | |
| WaSABI / registry.jl | Hybrid.Registry | |
| WaSABI / binary_inspiral.jl,waveform.jl | Hybrid.BinaryWaveform | |
| WaSABI / utilities.jl | Hybrid.Introspection | |
| WaSABI / kerr_isco.jl,pn_kinematics.jl | → Spacetime | ★ |
| WaSABI / interpolation.jl Interp1/2/3 ; inspiral_WaSABI-C_deps.jl _digamma,PolyGamma | → NumericsCore.interp + BigSpecial.digamma | ★ |
| WaSABI / wl_parser.jl read_wl_* | → Data.wl | ★ |

## → AccretionDisc

| Legacy pkg / file / fn | New home | ★ |
|---|---|---|
| DiscModes / newtonian_disc.jl,master_equation_cool_newt.jl,solver_newt.jl,heightaverage.jl,relativistic_disc.jl,off_equator.jl,relativistic_master(PDE).jl,lindblad_torque.jl,metric_recon.jl(bl_angular_eq) | AccretionDisc.Modes | |
| DiscEMRIForce / loader.jl,assemble.jl,kerr_source.jl,source.jl,fluid.jl,pipeline.jl | AccretionDisc.Force | |
| DiscModes / {kerr_geometry,metric_recon,relativistic_disc}.jl geometry+_risco | → Spacetime.Geometry | ★ |
| DiscModes / metric_recon.jl / Ylm_real,Ylm_eq,assoc_legendre | → Harmonics | ★ |
| DiscModes / lorenz_metric.jl / even_ode!,odd_ode!,monopole_ode!(g,rs,y,ω),odd_dipole_ode! | → MetricReconstruction.Lorenz.LorenzRadial (CANONICAL) | ★ |
| DiscModes / lorenz_metric.jl / r_of_rs | → Spacetime.Geometry | ★ |
| DiscModes / integrator.jl / DP5 ; DiscEMRIForce/numerics.jl / fornberg/cheb/bary | → NumericsCore (ODE / interp) | ★ |
| DiscEMRIForce / loader.jl / read_h1_mode path | → Data.h1schema (datadir injected) | ★ |
