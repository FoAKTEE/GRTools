"""
    AccretionDisc.Force

The circular-equatorial Kerr thin-disc torque front-end (from legacy `DiscEMRIForce`), the full
design-§4.10 pipeline: (i) the a=0 Lorenz BLS corpus loader (`H1Mode` over `Data.H1Schema`) and the
a≠0 front-end CONSUMING `MetricReconstruction.Lorenz.{closure_kerr,smatrices_kerr,reconstruct_kerr,
to_BLS_kerr}` (`kerr_mode`/`kerr_modes` — NEVER through SelfForce, §7.3); (ii) the parity-aware
tetrad→BLS `h_μν` assembly; (iii) the fluid-frame effective force → 𝒬-projected source →
radial forcing 𝔉(r) (+ the wing-precision `source_F_grid` array pass); (iv) fluid reconstruction
(δU,δp,δε) + the advected angular-momentum flux; (v) the forced master-ODE solve with the F8
torque-balance / flux cross-checks and the F5 source-structure driver. HDF5 is pulled directly ONLY
here, activating the `Data` HDF5 extension (design §5.4 Rev-2 extension-activation rule).
Type-generic `{T}`; BigFloat-inherited. Threaded per-(l)/per-grid-point sweeps are race-free
per-slot + blas_single_thread-guarded, bit-identical to serial (§8).
"""
module Force

import HDF5   # activates Data's HDF5 weakdep ext (read_h1/write_mode_h5 bodies) for this session

include("Loader.jl")
export H1Mode, load_h1_mode, load_modes, available_modes

include("Assemble.jl")
export kerr_parity, tetmodes_from_bls, hmunu, hmunu_eq

include("Source.jl")
export disc_a0, disc, force_mode, Vvec, source_F, source_F_grid

include("Fluid.jl")
export FluidVars, reconstruct, Jadv

include("Pipeline.jl")
export ForcingInterp, build_forcing, forced_solve
export f8_torque_balance, f8_torque_balance_ilr, f8_flux_crosscheck, f5_source

include("KerrFrontend.jl")
export kerr_mode, kerr_modes

end # module Force
