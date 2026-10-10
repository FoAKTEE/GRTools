"""
    NumericsCore (GrToolsBase)

Dimensionless math + infrastructure spine for GrTools — the unique DAG sink whose
ONLY dependency is `LinearAlgebra`. It carries no physics geometry; that
lives in `Spacetime`. Everything here is type-generic `{T}` so BigFloat precision is
inherited, never hard-Float64'd on a generic path.

Owns (per design §5.1):
  - Carlson `RF`/`RC`/`RD`/`RJ` + complete/incomplete Legendre `K`/`E`/`Π` + Jacobi am/sn/cn/dn
  - the interpolation library (cubic/bicubic/tricubic, GSL-bicubic, barycentric, Chebyshev,
    Fornberg-FD, cubic-Hermite, WL-faithful Interp1/2/3)
  - normalized associated-Legendre `P` seed
  - `(l,m,n,k)` `ModeIndexTable` + ladder factors + `a_il`
  - BigFloat precision policy + `bc_tol`/`rk_tol`/`fd_step` selectors
  - split-named units + convention helpers (`p_to_y`, ε↔ω, q↔a, imagunit)
  - `blas_guard`, `backend_iface`
  - vendored adaptive-RK integrator (DP5/GBS/DOPR853), continued-fraction primitive,
    self-contained Kummer 1F1, root-finders (secant/bisect/brent + A42), dep-free DataSchema structs
  - 1D nodal-DG spectral-element core (LGL `RefElement` with self-contained Newton-on-P'_N
    Gauss-Lobatto rule, multi-domain `Mesh`, `rk4!`/`cfl_dt`) — DG-Teulolsky port, FGQ-free

Legacy packages absorbed (as dedup'd primitives): KerrGeodesics, SpinWeightedSpheroidalHarmonics,
Teukolsky, H1Lorenz, KerrBackwards, FastEMRIWaveforms, DiscModes, KerrLorenzCirc, SecondOrderRicci,
DiscEMRIForce, Punctures, WaSABI, SpinningOrbit. Those stay green as the characterization oracle.
"""
module NumericsCore

using LinearAlgebra

# Module files are added incrementally as each primitive is TDD-ported (carlson, precision,
# modeindex, units, conventions, blas_guard, interp, legendre, ode, rootfind, ...).

include("specialfun_carlson.jl")
export RC, RF, RD, RJ
export ellipK, ellipE, ellipPi, ellipFinc, ellipEinc, ellipPiinc
export jacobiAM, jacobiSN, jacobiCN, jacobiDN

include("precision.jl")
export bc_tol, rk_tol, fd_step, bigprec_for_mag, with_precision

include("modeindex.jl")
export gamma_to_lm, lm_to_gamma, mode_list
export ModeIndexTable, mode_index_table
export a_il, lambdahat, Lambdahat

include("units.jl")
export MTSUN_SI, YRSID_SI, MRSUN_SI, AUsec, Gpc_in_m, Gpc_in_s, GPCINSEC, PI, PI_2, PI_4

include("conventions.jl")
export p_to_y, eps_to_omega, omega_to_eps, q_to_a, a_to_q, imagunit

include("blas_guard.jl")
export with_blas_threads, blas_single_thread

include("rootfind.jl")
export secant, bisect, brent, a42

include("jets.jl")
export constjet, jadd, jscale, jmul, jderiv, jinv

include("continued_fraction.jl")
export cf_finite, cf_convergent, cf_lentz

include("dataschema.jl")
export H1Schema, AmpGridSchema, FluxGridSchema, validate_schema, schema_version, SCHEMA_VERSION

include("interp.jl")
export CubicSplineInterpolant, BicubicNotAKnot, TricubicNotAKnot, GslBicubic
export natural_spline_node_derivs, gsl_bsearch, cubic_hermite
export fornberg_weights, highorder_deriv, bary_weights, bary_eval, bary_row, bary_interp, local_interp
export cheb_nodes, cheb_diffmatrix, cheb_deriv, cheb_smooth_deriv, fd1c, fd2c
export Interp1WL, interp1_movingwindow, Interp2Hermite, interp2_hermite, partial
export Interp3Lagrange, interp3_lagrange
export LogChebInterp, log_cheb_interp, BigLogRInterp, biglogr_interp

include("chebyshev.jl")   # after interp.jl: CGL nodes + D1/D2 (Weideman–Reddy) + dense solve for
export cheb_lobatto, cheb_collocation, cheb_solve   # the Q2-H hyperboloidal collocation operator

include("fft.jl")
export fft_radix2!, fft_radix2

include("ode.jl")
export dp5_step!, dp5_step_buffered!, integrate_grid, integrate_to_targets!,
       gragg_midpoint!, gbs_kmax, gbs_integrate_to_targets!
export DOPR853, take_step, integrate, dense_eval

include("legendre.jl")
export factorial_int, assoc_legendre_P, legendre_P_thetaderiv, assoc_legendre
export Ylm_real, Ylm_eq, pl_deriv, spin_weighted_Y

include("kummer.jl")
export kummer_1F1, kummer_U

include("dg.jl")
export RefElement, Mesh, build_mesh, element_of, lagrange_interp, node_coords, rk4!, cfl_dt

include("fd_uniform.jl")   # after dg.jl: reuses _legendre_PdP-free; uniform-grid FD sweep + KO filter
export set_d1!, set_d2!, ko_filter!, KO6_KERNEL

include("gauss_legendre.jl")   # after dg.jl: reuses _legendre_PdP for the Newton-on-P_n root solve
export gl_nodes_weights

include("backend_iface.jl")
export AbstractBackend, CPUBackend, select_backend, backend_array, backend_collect

end # module NumericsCore
